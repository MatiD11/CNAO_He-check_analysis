import logging
from contextlib import redirect_stdout
import os
from pathlib import Path
from shutil import make_archive
from signal import SIGINT
import threading
from time import sleep

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image, ImageSequence
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit
from sklearn.cluster import DBSCAN
import sympy as sp
import pandas as pd

from flask import Flask, render_template, request, redirect, url_for


#refractive index
n_2 = 1.58
#scintillator and mirror length [mm]
l_s = 200.
l_m = 200.
#total number of pixels
n_pz = 2560
n_py = 2160
#lens-camera distance
L = 650
#mirror-scintillator distance
S = 22.
#mm per pixels at the beginning and end of the scintillator
k_o=0.162
k_s = 0.21

#define pixel corresponding to S   
def pix_S(k_s):
    return int(S/k_s)

#transform pixels of the total image into number of pixels corresponding to apparent positions
def pix_local_scint(pix_z):
    return pix_z - n_pz/2

def pix_local_scinty(pix_y):
    return n_py/2 - pix_y

def pix_local_mirr(pix_x):
    return n_pz - pix_x-n_pz/2-pix_S(k_s)

#geometrical optics formulas
def t2(t1):
        return sp.asin(sp.sin(t1 / n_2))
def th1(za,L_):
    return sp.atan(za / (L_ + l_s))
def ph1(xa, zt,L_):
    return sp.atan((xa + S) / (L_ + l_s + zt + S))
def ztt(za, xt,L_):
    return (L_ * za) / (L_ + l_s - xt) + (l_s - xt) * sp.tan(t2(th1(za,L_)))
def xtt(xa, zt,L_):
    return xa - zt * (sp.tan(ph1(xa, zt,L_)) - sp.tan(t2(ph1(xa, zt,L_))))
def xat(zt,n_x):
    return ((l_s+S+zt)*(k_s-k_o)/l_s + k_o)*n_x
def zat(xt,n_z):
    return ((l_s-xt)*(k_s-k_o)/l_s + k_o)*n_z

#numerical solution of the gemetrical optics formulas
def find_position(nx, ny, nz, L_):
    n_x = pix_local_mirr(nx)
    n_y = pix_local_scinty(ny)
    n_z = pix_local_scint(nz)

    ztf, xtf, xaf, zaf = sp.symbols('ztf xtf xaf zaf')

    eq1 = ztf - ztt(zaf, xtf, L_)
    eq2 = xtf - xtt(xaf, ztf, L_)
    eq3 = xaf - xat(ztf, n_x)
    eq4 = zaf - zat(xtf, n_z)
    
   
    sol = sp.nsolve((eq1, eq2, eq3, eq4), (xaf, xtf, zaf, ztf), (1, 1, 1, 1))
   
    sol_arr = np.array(sol)
    positions = sol_arr.flatten()
    positions[0] = l_s - sol[0]
    positions[1] = l_s - sol[1]
    positions[2] = l_s - sol[2]
    positions[3] = l_s - sol[3]
    positions = np.append(positions,zat(sol[1],n_y)) 
    positions= np.append(positions, ztt(positions[4],sol[1],L))
    
    return positions

#define a gaussian fit
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def find_biggest_clust(half_points, half_labels):
    unique_clusters = np.unique(half_labels)
    unique_clusters = unique_clusters[unique_clusters != -1]

    largest_cluster_size = 0
    largest_cluster_points = None

    for cluster_label in unique_clusters:
        cluster_points = half_points[half_labels == cluster_label]
        cluster_size = len(cluster_points)
        if cluster_size > largest_cluster_size:
            largest_cluster_size = cluster_size
            largest_cluster_points = cluster_points
    return largest_cluster_points

def tif2array(page):
    imarray = np.array(page)
    
    #take bkg images
    bkg_images = [f'bkg/bkg_img{i}.tif' for i in range(1, 6)]
    bkg_arrays = [np.array(Image.open(bkg)) for bkg in bkg_images]
    bkg = np.mean(bkg_arrays, axis=0)
    
    #bkg subtraction
    img_int32 = imarray.astype(np.int32)
    bkg_int32 = bkg.astype(np.int32)
    img = img_int32 - bkg_int32
    img[img < 0] = 0
    img = img.astype(np.uint16)

    #median filter to clean hot spots
    img = median_filter(img, size = 4).astype("uint16")
    
    #remove noisy borders
    img[:5, :] = 0
    img[-5:, :]=0
    
    return img

def find_bragg_peak(image):
    scint_row = np.sum(image[2:2158, 1280:], axis=0)
    scint_row = scint_row[::-1]
    max_index = np.argmax(scint_row)
    max = scint_row[max_index]
    post_max_array = scint_row[max_index + 1:]
    closest_index_post_max = (np.abs(post_max_array - 0.9*max)).argmin()
    closest_index = max_index + 1 + closest_index_post_max  
    bragg_peak_x = n_pz - closest_index
    return bragg_peak_x

def find_y_range(image):
    y_scint_column = np.sum(image[2:2158, 1280:],axis=1)
    threshold_fraction = 0.5
    threshold_value = threshold_fraction * y_scint_column[5:2154].max()
    bright_indices = np.where(y_scint_column[5:2154] > threshold_value)[0]

    if len(bright_indices) == 0:
        raise ValueError("Interval not found")
    start_bright = bright_indices[0]
    end_bright = bright_indices[-1]
    y_range = (start_bright + end_bright) // 2
    return y_range

def find_x_mirror(image):
    mirror_row = np.sum(image[2:2158, :1280], axis=0)
    spot_mirror = np.argmax(mirror_row)
    half_max_value = 0.5 * np.max(mirror_row)    	
    left_index = np.where(mirror_row[:spot_mirror] < half_max_value)[0]
    right_index = np.where(mirror_row[spot_mirror:] < half_max_value)[0]
    if len(left_index) == 0:
        start_fit = 0
    else:
        start_fit = left_index[-1]
    if len(right_index) == 0:
        end_fit = len(mirror_row) - 1
    else:
        end_fit = spot_mirror + right_index[0] 
    x_fit = np.arange(start_fit, end_fit)
    y_fit = mirror_row[start_fit:end_fit]
    par_fit, _ = curve_fit(gaussian, x_fit, y_fit, p0=[y_fit.max(), spot_mirror, 1])
    x_mirror = int(np.round(par_fit[1]))
    return x_mirror

def cluster_img(img):
    # Perform some additional filtering to help clusterization
    img[img < (img.max()*0.1)] = 0

    # Clusterization
    x = np.array(np.nonzero(img)).astype(int).T
    db = DBSCAN(eps=20, min_samples=1000).fit(x)
    labels = db.labels_

    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    print(f'Found {n_clusters_} clusters and {n_noise_} noise points.')

    left_half_points = x[x[:, 1] < img.shape[1] // 2]
    left_half_labels = labels[x[:, 1] < img.shape[1] // 2]
    right_half_points = x[x[:, 1] > img.shape[1] // 2]
    right_half_labels = labels[x[:, 1] > img.shape[1] // 2]
    largest_clust_points_spot = find_biggest_clust(left_half_points, left_half_labels)
    largest_clust_points_range = find_biggest_clust(right_half_points, right_half_labels)

    if largest_clust_points_spot is not None:
        centroid_x = np.median(largest_clust_points_spot[:, 1])
    else:
        print("There are no clusters in the left half of the image.")
        
    if largest_clust_points_range is not None:
        centroid_y = np.median(largest_clust_points_range[:, 0])
        cluster_image = np.zeros_like(img, dtype=np.uint16)
        for point in largest_clust_points_range:
            cluster_image[point[0], point[1]] = img[point[0], point[1]]
        clust_range_position = find_bragg_peak(cluster_image)
        cluster_position = find_position(centroid_x, centroid_y, clust_range_position, L)
        return cluster_position

    else:
        print("There are no clusters in the right half of the image.")


def process_img(proton_range):

    out_directory = Path(f'static/proton_{proton_range}')
    out_directory.mkdir(parents=True, exist_ok=True)
    tif = Image.open(f'images/immaginep_{proton_range}.tif')
    positions = []
    cluster_positions = []    
    img_rows = []
    with open(out_directory / f'results_{proton_range}.txt', 'w') as outfile:
        with open("status.txt", "a") as status_file:
            with redirect_stdout(status_file):
                print(f"Selected {proton_range} mm proton range data. Image processing has started.\n")
                status_file.flush()
                outfile.write('X\tY\tRange_position\tX_cluster\tY_cluster\tRange_position_cluster\n')
                for i, page in enumerate(ImageSequence.Iterator(tif)):
                    img = tif2array(page)
                    print(f"Image no. {i+1}:")
                    status_file.flush()

                    x_mirror = find_x_mirror(img)
                    bragg_peak = find_bragg_peak(img)
                    y_range = find_y_range(img)
                  
                    position = find_position(x_mirror, y_range, bragg_peak, L)
                    positions.append([position[1], position[5], position[3]])

                    cluster_position = cluster_img(img)    
                    cluster_positions.append([cluster_position[1], cluster_position[5], cluster_position[3]])

                    print("Pixel-position mapping completed.\n")
                    status_file.flush() 
                    
                    line = f"{position[1]}\t{position[5]}\t{position[3]}\t{cluster_position[1]}\t{cluster_position[5]}\t{cluster_position[3]}\n"            
                    outfile.write(line)       

                    '''img_rgb = cv2.cvtColor(img, cv2.COLOR_GRAY2RGB)
                    img_dot = cv2.circle(img_rgb, (bragg_peak, y_range), radius=15, color=(255, 0, 0), thickness=-1)
                    img_dot = cv2.circle(img_dot, (x_mirror, y_range), radius=15, color=(255, 0, 0), thickness=-1)'''
                    plt.figure(figsize=(12, 6))
                    plt.subplot(1, 2, 1)
                    
                    plt.imshow(np.array(page), cmap = 'gray',vmin=0, vmax=255)
                    plt.title("Original image")
                    plt.subplot(1, 2, 2)
                    plt.imshow(img, cmap = 'gray',vmin=0, vmax=255)
                    plt.title("Processed image")
                    
                    plt.tight_layout()
                    plt.savefig(out_directory / f'images_{i}.png')
                    plt.close() 

                    img_rows.append(np.sum(img[2:2158,:], axis=0))

                print(f"Analysis completed. Redirecting to result page.\n")
                status_file.flush()

    positions_x = [pos[0] for pos in positions]
    positions_y = [pos[1] for pos in positions]

    cluster_positions_x = [cpos[0] for cpos in cluster_positions]
    cluster_positions_y = [cpos[1] for cpos in cluster_positions]

    fig, ax = plt.subplots(figsize=(8, 8)) 

    ax.set_facecolor('lightblue')
    ax.scatter(positions_x, positions_y, c='blue', label='Integration method')
    ax.scatter(cluster_positions_x, cluster_positions_y, c='red', label='Clusterization method')

    ax.set_xlim(0, 200)
    ax.set_ylim(-100, 100)
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_title('Hit positions of proton beam in the scintillator cube')
    ax.legend()
    plt.savefig(out_directory / f'hit_map.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    #bragg peaks for different x positions   
    plt.figure()
    plt.plot(np.arange(1600, len(img_rows[0])), img_rows[0][1600:])
    plt.plot(np.arange(1600, len(img_rows[5])), img_rows[5][1600:])        
    plt.plot(np.arange(1600, len(img_rows[6])), img_rows[6][1600:])
    plt.legend(['x = 130, y = -30','x = 130, y = 0','x = 130, y = 30'])
    plt.title("Bragg peak")
    plt.xlabel("Pixel")
    plt.ylabel("Intensity")
    plt.savefig(out_directory / f'x_130.png')
    plt.close()

    plt.figure()
    plt.plot(np.arange(1600, len(img_rows[1])), img_rows[1][1600:])
    plt.plot(np.arange(1600, len(img_rows[4])), img_rows[4][1600:])        
    plt.plot(np.arange(1600, len(img_rows[7])), img_rows[7][1600:])
    plt.legend(['x = 100, y = -30','x = 100, y = 0','x = 100, y = 30'])
    plt.title("Bragg peak")
    plt.xlabel("Pixel")
    plt.ylabel("Intensity")
    plt.savefig(out_directory / f'x_100.png')
    plt.close()

    plt.figure()
    plt.plot(np.arange(1600, len(img_rows[2])), img_rows[2][1600:])
    plt.plot(np.arange(1600, len(img_rows[3])), img_rows[3][1600:])        
    plt.plot(np.arange(1600, len(img_rows[8])), img_rows[8][1600:])
    plt.legend(['x = 70, y = -30','x = 70, y = 0','x = 70, y = 30'])
    plt.title("Bragg peak")
    plt.xlabel("Pixel")
    plt.ylabel("Intensity")
    plt.savefig(out_directory / f'x_70.png')
    plt.close()


app = Flask('Image_analysis')

status_file = "status.txt"
Path(status_file).write_text("") 
valid_ranges = [30, 60, 101]

shutdown_requested = False

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == "POST":
        proton_range = request.form.get("proton_range")
        return redirect(url_for('processing', proton_range=proton_range))
    return render_template('index.html', valid_ranges=valid_ranges)

@app.route('/processing')
def processing():
    proton_range = request.args.get('proton_range')
    if not hasattr(app, 'processing_thread') or not app.processing_thread.is_alive():
        app.processing_thread = threading.Thread(target=process_img, args=(proton_range,))
        app.processing_thread.start()    
    return render_template('processing.html', proton_range=proton_range)

@app.route('/status')
def status():
    with open(status_file, "r") as f:
        status = f.read()
    if "Analysis completed" in status:
        sleep(2)
        return redirect(url_for('result', proton_range=request.args.get('proton_range')))
    return status, 200, {'Content-Type': 'text/plain'}

@app.route('/result')
def result():
    proton_range = request.args.get('proton_range')

    Path(status_file).write_text("") 

    out_directory = Path(f'static/proton_{proton_range}')
    out_directory.mkdir(parents=True, exist_ok=True)

    # Zip the results directory after processing
    static_dir = os.path.join(app.root_path, out_directory)
    zip_path = os.path.join(app.root_path, f'{out_directory}.zip')
    make_archive(zip_path.replace('.zip', ''), 'zip', static_dir)

    # Read the results file
    results_file = Path(out_directory / f'results_{proton_range}.txt')
    df = pd.read_csv(results_file, sep='\t', skiprows=0)

    # Extract the required columns
    ranges_integration = df.iloc[:, 2].values
    ranges_clusterization = df.iloc[:, 5].values

    # Calculate means and standard deviations
    mean_integration = np.mean(ranges_integration) if ranges_integration.size > 0 else None
    std_integration = np.std(ranges_integration) if ranges_integration.size > 0 else None
    mean_clusterization = np.mean(ranges_clusterization) if ranges_clusterization.size > 0 else None
    std_clusterization = np.std(ranges_clusterization) if ranges_clusterization.size > 0 else None

    if not hasattr(app, 'shutdown_thread'):
        print("\nAnalysis completed. Click button on results page to shut down.")
        app.shutdown_thread = threading.Thread(target=wait_for_shutdown)
        app.shutdown_thread.start()

    return render_template(
        'result.html',
        proton_range=proton_range,
        ranges_integration=ranges_integration,
        ranges_clusterization=ranges_clusterization,
        mean_integration=mean_integration,
        std_integration=std_integration,
        mean_clusterization=mean_clusterization,
        std_clusterization=std_clusterization,
        zip=zip
    )

@app.route('/shutdown_server', methods=['POST'])
def shutdown_server():
    global shutdown_requested
    shutdown_requested = True
    return '', 204  # No Content

@app.route('/shutdown_page', methods=['GET'])
def shutdown_page():
    return render_template('shutdown.html')

def wait_for_shutdown():
    global shutdown_requested
    while not shutdown_requested:
        sleep(1)
    print("Shutting down server...")
    os.kill(os.getpid(), SIGINT)

if __name__ == '__main__':
    host = os.getenv('FLASK_RUN_HOST', '127.0.0.1')
    port = int(os.getenv('FLASK_RUN_PORT', 5000))
    print(f"Starting Image Analysis interface. Access it at http://{host}:{port}")
    
    logging.getLogger('werkzeug').setLevel(logging.ERROR)  # Reduce verbosity
    
    app.run(host=host, port=port, debug=False)




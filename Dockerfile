FROM python:3.10-slim
WORKDIR /usr/app
COPY . .
RUN pip install --no-cache-dir -r requirements.txt 
EXPOSE 80
CMD [ "python3", "app.py" ]
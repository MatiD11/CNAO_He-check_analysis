#!/bin/bash

export FLASK_RUN_HOST=0.0.0.0
export FLASK_RUN_PORT=80

docker build -t image_analysis .
docker run -p 80:80 --env FLASK_RUN_HOST --env FLASK_RUN_PORT --rm image_analysis

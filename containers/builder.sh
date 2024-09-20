#!/bin/bash

# Base directory to search from (replace `.` with your specific path if needed)
BASE_DIR="."

# Find all Dockerfiles and build images
find "$BASE_DIR" -name 'Dockerfile' | while read dockerfile; do
    # Get the directory of the Dockerfile
    docker_dir=$(dirname "$dockerfile")
    
    # Extract the directory name (used as the image name)
    image_name=$(basename "$docker_dir")
    
    # Build the Docker image
    echo "Building image: $image_name from $dockerfile"
    docker build -t "$image_name" "$docker_dir"
done
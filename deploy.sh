#!/bin/bash

echo "Starting deployment"

# Optional: fail if any command fails
set -e

# Pull the latest changes from the repo
echo "Pulling latest code"
git pull

# Stop and remove old containers
echo "Shutting down old containers"
docker-compose -f docker-compose-prod.yml down

# Rebuild and restart the containers
echo "Building and launching new containers"
docker-compose -f docker-compose-prod.yml up --build -d

echo "Deployment complete!"
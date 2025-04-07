#!/bin/bash

# static file collection not needed in development with whitenoise

python manage.py makemigrations
# Apply database migrations
echo "Apply database migrations"
python manage.py migrate

# Start server
echo "Starting server"
python manage.py runserver 0.0.0.0:8000
# gunicorn cfchem.wsgi:application --bind 0.0.0.0:8000
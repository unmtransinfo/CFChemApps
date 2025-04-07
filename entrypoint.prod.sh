#!/bin/bash

echo "Applying database migrations"
python manage.py migrate --noinput

echo "Collecting static files"
python manage.py collectstatic --noinput

# Start server
echo "Starting server"
gunicorn cfchem.wsgi:application --bind 0.0.0.0:8000
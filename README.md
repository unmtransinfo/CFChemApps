## Installation without Docker

    conda create --name <env_name> python=3.8
    pip install Flask
    conda install -c rdkit rdkit

To run the application, use 

    python manage.py makemigrations
    python manage.py migrate
    python manage.py runserver 

## Installation with Docker
If you have docker installed on your system, you can run the app directly using docker
Use the following command to do so

    docker-compose up -d --build 
You can leave out `-d` flag if you do not wish to run in detached mode

> On successfully running the server, go to http://127.0.0.1:5000/depict/ to check out the portal.

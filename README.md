## Installation without Docker

    conda create --name env_name python=3.8
    pip install Flask
    conda install -c rdkit rdkit

To run the application, use 

    flask --app app run --debug
You can skip the `--debug` flag if you do not want to run the app in debug mode

## Installation with Docker
If you have docker installed on your system, you can run the app directly using docker
Use the following command to do so

    docker-compose up -d
You can leave out `-d` flag if you do not wish to run in detached mode
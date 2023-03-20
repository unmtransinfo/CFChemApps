# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster

# Set environment variables
ENV PIP_DISABLE_PIP_VERSION_CHECK 1
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
# PIP_DISABLE_PIP_VERSION_CHECK disables an automatic check for pip updates each time
# PYTHONDONTWRITEBYTECODE means Python will not try to write .pyc files
# PYTHONUNBUFFERED ensures our console output is not buffered by Docker

WORKDIR /app

RUN apt update && apt install -y libsm6 libxext6
RUN apt-get install -y libxrender-dev

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .

CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0"]
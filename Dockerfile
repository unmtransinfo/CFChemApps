# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster

# Set environment variables
ENV PIP_DISABLE_PIP_VERSION_CHECK 1
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
# PIP_DISABLE_PIP_VERSION_CHECK disables an automatic check for pip updates each time
# PYTHONDONTWRITEBYTECODE means Python will not try to write .pyc files
# PYTHONUNBUFFERED ensures our console output is not buffered by Docker

RUN apt update && apt install -y libsm6 libxext6
RUN apt-get install -y libxrender-dev

WORKDIR /cfchem

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

# copy entrypoint.sh
COPY ./entrypoint.sh /cfchem/entrypoint.sh
RUN sed -i 's/\r$//g' /cfchem/entrypoint.sh
RUN chmod +x /cfchem/entrypoint.sh

COPY . .

# ENTRYPOINT ["/cfchem/entrypoint.dev.sh"]
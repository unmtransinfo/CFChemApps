# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster

# Set environment variables
ENV PIP_DISABLE_PIP_VERSION_CHECK 1
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
# PIP_DISABLE_PIP_VERSION_CHECK disables an automatic check for pip updates each time
# PYTHONDONTWRITEBYTECODE means Python will not try to write .pyc files
# PYTHONUNBUFFERED ensures our console output is not buffered by Docker

RUN apt update && apt install -y libsm6 libxext6 libcairo2 libpango1.0-0 libpangocairo-1.0-0 libgdk-pixbuf2.0-0 libffi-dev
RUN apt-get install -y gcc libxrender-dev

WORKDIR /cfchemapps

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

# copy entrypoint.sh
COPY ./entrypoint.dev.sh /entrypoint.dev.sh
RUN sed -i 's/\r$//g' /entrypoint.dev.sh
RUN chmod +x /entrypoint.dev.sh

COPY . .

ENTRYPOINT ["/entrypoint.dev.sh"]
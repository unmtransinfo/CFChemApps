# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster

WORKDIR /cfchem

COPY requirements.txt requirements.txt
RUN pip3 install --upgrade pip
# RUN pip3 install wheel setuptools manimce
RUN pip3 install -r requirements.txt

COPY . /cfchem

CMD [ "python3", "-m", "flask", "run", "--host=0.0.0.0:3000" ]
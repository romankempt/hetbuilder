FROM python:3.9

WORKDIR /hetbuilder-app

COPY requirements.txt .

ADD bin/hetbuilder

RUN pip install -r requirements.txt
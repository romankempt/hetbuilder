#-----------------------------------------------------------------------------
FROM ubuntu:focal
LABEL Author, Roman Kempt

#  $ docker build . -t continuumio/miniconda3:latest -t continuumio/miniconda3:4.5.11
#  $ docker run --rm -it continuumio/miniconda3:latest /bin/bash
#  $ docker push continuumio/miniconda3:latest
#  $ docker push continuumio/miniconda3:4.5.11

ENV APP_HOME /app
WORKDIR $APP_HOME
COPY . $APP_HOME

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate" >> ~/.bashrc  

#---------------- Prepare the environnment
RUN conda instal --file environment.yaml

RUN pip install -vv .

#ENTRYPOINT ["conda", "run", "--name", "hetbuilder", "python", "-m", "flask", "run", "app/app.py"]
WORKDIR $HOME
#ENTRYPOINT [ "hetbuilder --help" ]
#COPY requirements.txt requirements.txt

#ENTRYPOINT [ "hetbuilder"]
#CMD [ "--help" ]


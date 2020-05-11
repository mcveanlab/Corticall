FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/sv/bin/:${PATH}

# install packages
RUN apt update
RUN apt-get install -y --no-install-recommends curl \
    && rm -rf /var/lib/apt/lists/*

FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/mccortex/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install packages
RUN apt update
RUN apt-get install -y --no-install-recommends curl make gcc g++ autoconf libc6-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# install cd-hit
RUN git clone https://github.com/weizhongli/cdhit.git \
    && cd cdhit \
    && make MAX_SEQ=10000000 \
    && make install

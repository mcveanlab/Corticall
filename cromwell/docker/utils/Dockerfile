FROM ubuntu:latest

# install base packages
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends git ssh ca-certificates libssl-dev curl python python3.6 time wget gcc make libc6-dev zlib1g-dev && \
	rm -rf /var/lib/apt/lists/*

# install Flash paired-end merger
RUN wget http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11-Linux-x86_64.tar.gz && \
    gunzip -c FLASH-1.2.11-Linux-x86_64.tar.gz | tar xvf - && \
    cd FLASH-1.2.11-Linux-x86_64 && \
    mv flash /usr/local/bin/

# install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh \
	&& bash /miniconda.sh -b -p /miniconda
ENV PATH=/miniconda/bin/:/miniconda/envs/utils/bin/:/root/google-cloud-sdk/bin/:$PATH

# install conda packages
COPY environment.yml /
RUN conda env create -f /environment.yml \
	&& conda clean -a

# install gsutil
RUN curl https://sdk.cloud.google.com | bash

# install readsim
RUN git clone https://github.com/noporpoise/readsim.git \
    && cd readsim \
    && make \
    && mv readsim /usr/local/bin

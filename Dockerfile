FROM nvidia/cuda:11.6.2-base-ubuntu20.04

# ENV TZ=America/Chicago
# RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# https://forums.developer.nvidia.com/t/18-04-cuda-docker-image-is-broken/212892/8
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/3bf863cc.pub
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/machine-learning/repos/ubuntu2004/x86_64/7fa2af80.pub

RUN apt-get update \
    && apt-get install -y wget git libnvidia-compute-510-server uuid \
        build-essential bzip2 zlib1g-dev libbz2-dev libcurl4-gnutls-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev \
        # python3 python3-pip lsb-release apt-transport-https \
    && apt-get clean \
    && apt purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# ENV PLATFORM=$(lsb_release -cs)
# RUN bash -l -c 'echo export PLATFORM="$(lsb_release -cs)" >> /etc/bash.bashrc'
# ENV PLATFORM=focal

# RUN wget -O- https://europe.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add -
# RUN echo "deb http://europe.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list
# RUN apt-get update

# RUN apt update
# RUN apt install -y ont-guppy
# RUN pip install megalodon
# RUN git clone https://github.com/nanoporetech/rerio
# RUN rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases_5mC_CpG_v001 rerio/basecall_models/res_dna_r941_prom_modbases_5mC_CpG_v001
# RUN chmod a+r /rerio/basecall_models/*


RUN wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy_6.1.2_linux64.tar.gz
RUN tar -xf ont-guppy_6.1.2_linux64.tar.gz
RUN rm ont-guppy_6.1.2_linux64.tar.gz
ENV PATH=/ont-guppy/bin:$PATH

# RUN apt install -y cuda-toolkit-11-6
# RUN apt-get clean
# RUN apt install -y nvidia-cuda-toolkit 

WORKDIR /opt
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar -xvf htslib-1.13.tar.bz2 && \
    rm -r /opt/htslib-1.13.tar.bz2 && \
    cd htslib-1.13 && \
    ./configure && \
    make install
WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar xvf samtools-1.13.tar.bz2 && \
    rm -r /opt/samtools-1.13.tar.bz2 && \
    cd samtools-1.13/ && \
    ./configure && \
    make install

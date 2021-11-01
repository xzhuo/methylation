FROM nvidia/cuda:11.4.2-devel-ubuntu20.04

ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update \
    && apt-get install -y git python3 python3-pip \
        wget lsb-release \
    && apt-get clean

# ENV PLATFORM=$(lsb_release -cs)
# RUN bash -l -c 'echo export PLATFORM="$(lsb_release -cs)" >> /etc/bash.bashrc'
ENV PLATFORM=focal

RUN wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add -
RUN echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list
RUN apt-get update

RUN apt update
RUN apt install -y ont-guppy
RUN pip install megalodon
RUN git clone https://github.com/nanoporetech/rerio
RUN rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases_5mC_CpG_v001 rerio/basecall_models/res_dna_r941_prom_modbases_5mC_CpG_v001

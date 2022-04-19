FROM nvidia/cuda:11.4.2-devel-ubuntu20.04

ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update \
    && apt-get install -y git python3 python3-pip \
        wget lsb-release libnvidia-compute-460-server apt-transport-https \
    && apt-get clean

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

RUN apt install -y nvidia-cuda-toolkit cuda-toolkit-11-2

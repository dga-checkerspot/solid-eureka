FROM ubuntu:18.04
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

RUN conda install -c anaconda svn

RUN conda install -c conda-forge unzip

RUN wget https://sourceforge.net/projects/popoolation2/files/latest/download/popoolation2_1201.zip --no-check-certificate

RUN unzip popoolation2_1201.zip

RUN chmod 777 popoolation2_1201/mpileup2sync.pl

RUN chmod 777 popoolation2_1201/fst-sliding.pl

RUN conda install -c anaconda perl

RUN conda install -c bioconda perl-pod-usage

RUN PATH=$PATH:/root/miniconda3/lib/5.26.2/

RUN conda install -c bioconda samtools

RUN conda install -c bioconda minimap2

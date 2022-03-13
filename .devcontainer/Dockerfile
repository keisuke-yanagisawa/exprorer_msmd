FROM gromacs/gromacs

SHELL ["/bin/bash", "-c"]

# install conda
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

# install packages
RUN conda install -c conda-forge mamba
RUN mamba create -n exprorer
RUN mamba install -n exprorer -c conda-forge python=3.7 rdkit ambertools=20.15 \
    pandas numpy biopython jinja2 griddataformats parmed openbabel==2.4.1 git pyyaml

RUN conda init \
    && echo "conda activate exprorer" >> ~/.bashrc


RUN mamba install -y tqdm scikit-learn
RUN mamba install -y -c conda-forge gputil

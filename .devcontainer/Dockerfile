FROM gromacs/gromacs:2021.5

# install conda
ENV PATH="/miniconda3/bin:${PATH}"
ARG PATH="/miniconda3/bin:${PATH}"

# for update of cuda
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh \
    && bash Miniconda3-py39_4.11.0-Linux-x86_64.sh -b \
         -p /miniconda3 \
    && rm -f Miniconda3-py39_4.11.0-Linux-x86_64.sh 
RUN conda --version

# install packages
# libarchive is required for mamba
RUN conda install -y -c conda-forge libarchive
RUN conda install -y -c conda-forge mamba
RUN mamba update conda
RUN mamba install -y -c conda-forge \
    ambertools numpy biopython \
    jinja2 griddataformats \
    parmed git pyyaml \
    tqdm scikit-learn gputil joblib
RUN pip install coverage pytest

RUN mamba install -y -c conda-forge flake8 autopep8

# ====== settings for git / ssh =======
RUN echo "mkdir -p /root/.ssh" >> /etc/bash.bashrc
RUN echo "cp /root/.ssh_host/* /root/.ssh/" >> /etc/bash.bashrc
RUN echo "chmod 600 /root/.ssh/*" >> /etc/bash.bashrc

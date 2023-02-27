FROM gromacs/gromacs:2021.5

# install conda
ENV PATH="/miniconda3/bin:${PATH}"
ARG PATH="/miniconda3/bin:${PATH}"

# for update of cuda
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
RUN apt-get update
RUN apt-get install -y make
# RUN apt-get upgrade

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh \
    && bash Miniconda3-py39_4.11.0-Linux-x86_64.sh -b \
         -p /miniconda3 \
    && rm -f Miniconda3-py39_4.11.0-Linux-x86_64.sh 
RUN conda --version

# install packages
RUN conda install -y -c conda-forge mamba
RUN mamba update conda
RUN mamba install -y -c conda-forge \
    ambertools=20.15 numpy=1.22.3 biopython=1.79 \
    jinja2=3.0.3 griddataformats=0.7.0 \
    parmed=3.4.3 git=2.35.1 pyyaml=6.0 \
    tqdm=4.62.3 scikit-learn=1.0.2 gputil=1.4.0 joblib=1.1.0 

# ====== settings for git / ssh =======
# change permission to read ssh settings if exists
RUN chown root:$USER /root/.ssh/config || true
RUN chmod 600        /root/.ssh/config || true
FROM kyanagis/exprorer_msmd:2025.01.25

RUN mamba install -y pytest
RUN pip install mypy pytest isort ipykernel
RUN apt install -y openssh-client

# ====== settings for git / ssh =======
RUN echo "mkdir -p /root/.ssh" >> /etc/bash.bashrc
RUN echo "cp /root/.ssh_host/* /root/.ssh/" >> /etc/bash.bashrc
RUN echo "chmod 600 /root/.ssh/*" >> /etc/bash.bashrc

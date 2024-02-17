FROM kyanagis/exprorer_msmd:2024.02.12

RUN mamba install -y pytest

# ====== settings for git / ssh =======
RUN echo "mkdir -p /root/.ssh" >> /etc/bash.bashrc
RUN echo "cp /root/.ssh_host/* /root/.ssh/" >> /etc/bash.bashrc
RUN echo "chmod 600 /root/.ssh/*" >> /etc/bash.bashrc

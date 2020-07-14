FROM ubuntu:18.04
MAINTAINER Yu Huang <polyactis@gmail.com>
ADD bashrc /root/.bashrc
ADD pyflow /root/pyflow_source
ADD samtools /usr/local/bin/
ADD strelka-2.9.10.centos6_x86_64 /usr/local/strelka
RUN apt-get update
RUN apt-get install -y apt-utils iputils-ping less libboost-program-options1.65.1 libboost-iostreams1.65.1 libbz2-1.0 libhdf5-dev liblzma5 libssl1.0.0 lynx man python-h5py python-matplotlib python-numpy python-pandas telnet tmux vim
# install pyflow, use absolute path. Default path is /.
# do not split it in two commands (cd and then run this. as every command starts from /)
RUN cd /root/pyflow_source/ && python setup.py build install
RUN rm -rf /root/pyflow_source
ENV HOME /root
CMD source /root/.bashrc

# Development docker image
FROM nvidia/cuda:11.1-devel

RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt -y upgrade
RUN DEBIAN_FRONTEND=noninteractive apt -y install intel-mkl cmake libopenblas-dev gdb libgtest-dev build-essential extra-cmake-modules cmake-extras

ENV MKL_INCLUDE_PATH /usr/include/mkl
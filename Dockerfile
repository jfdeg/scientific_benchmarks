# Development docker image
FROM nvidia/cuda:11.1-devel

RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt -y install intel-mkl cmake libopenblas-dev


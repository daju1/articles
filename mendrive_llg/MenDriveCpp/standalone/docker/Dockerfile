# Образ для сборки и запуска MenDrive standalone (Ubuntu 22.04)
FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive LANG=C.UTF-8 LC_ALL=C.UTF-8 MPLBACKEND=Agg
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential g++ make libomp-dev time procps \
        python3 python3-pip python3-numpy python3-matplotlib python3-pandas python3-scipy \
        python3-notebook jupyter-nbconvert ca-certificates \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /work/standalone
# исходники монтируются томом (см. run.sh), поэтому в образ не копируются
CMD ["bash"]

FROM ubuntu:18.04

ARG USER_NAME
ARG GROUP_NAME
ARG USER_ID
ARG GROUP_ID

ARG DEBIAN_FRONTEND=noninteractive

# Configure a 'docker' sudo user without password
RUN apt-get update && apt-get -y install sudo
RUN addgroup --gid ${GROUP_ID} ${GROUP_NAME}
RUN adduser --disabled-password --gecos '' --uid ${USER_ID} --gid ${GROUP_ID} ${USER_NAME}
RUN usermod -aG sudo ${USER_NAME}
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# Install build-essential package (for cross compilation)
RUN apt-get update
RUN apt-get -y install build-essential

# Install vim, nano
RUN apt-get -y install vim nano

# Install a bunch of packages for sage project
# https://doc.sagemath.org/html/en/installation/source.html

RUN apt -y install software-properties-common apt-utils 
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update

# Unable to locate package libbraiding-dev


RUN apt-get -y install bc binutils bzip2 ca-certificates cliquer curl eclib-tools fflas-ffpack flintqs g++ g++ gcc gcc gfan gfortran git glpk-utils gmp-ecm lcalc libatomic-ops-dev libboost-dev libbz2-dev libcdd-dev libcdd-tools libcliquer-dev libcurl4-openssl-dev libec-dev libecm-dev libffi-dev libflint-arb-dev libflint-dev libfreetype6-dev libgd-dev libgf2x-dev libgivaro-dev libglpk-dev libgmp-dev libgsl-dev libiml-dev liblfunction-dev liblrcalc-dev liblzma-dev libm4rie-dev libmpc-dev libmpfi-dev libmpfr-dev libncurses5-dev libntl-dev libopenblas-dev libpari-dev libpcre3-dev libplanarity-dev libppl-dev libpython3.7-dev libreadline-dev librw-dev libsqlite3-dev libsymmetrica2-dev libz-dev libzmq3-dev m4 make nauty pari-doc pari-elldata pari-galdata pari-galpol pari-gp2c pari-seadata patch perl pkg-config planarity ppl-dev python3 python3-distutils python3.7 r-base-dev r-cran-lattice sqlite3 tachyon tar xz-utils yasm flex

RUN apt-get -y install cmake coinor-cbc coinor-libcbc-dev libboost-dev libfile-slurp-perl libisl-dev libjson-perl libmongodb-perl libperl-dev libsvg-perl libterm-readkey-perl libterm-readline-gnu-perl libterm-readline-gnu-perl libxml-libxslt-perl libxml-writer-perl libxml2-dev ninja-build pandoc pari-gp2c

RUN apt-get -y install openssh-server openssh-client

RUN apt-get -y install tk tk-dev

RUN apt-get -y install texlive       # debian latex

RUN apt-get -y install texlive-generic-extra #(to generate pdf documentation)
RUN apt-get -y install texlive-latex-extra #(to generate pdf documentation)
RUN apt-get -y install texlive-xetex #(to convert Jupyter notebooks to pdf)
RUN apt-get -y install latexmk #(to generate pdf documentation)
RUN apt-get -y install pandoc #(to convert Jupyter notebooks to pdf)
RUN apt-get -y install dvipng #(to render text with LaTeX in Matplotlib)
RUN apt-get -y install default-jdk #(to run the Jmol 3D viewer from the console and generate images for 3D plots in the documentation)
RUN apt-get -y install ffmpeg #(to produce animations)
RUN apt-get -y install libavdevice-dev #(to produce animations)
RUN apt-get -y install texlive-lang-cyrillic
# RUN apt-get -y install texlive-full

RUN apt-get -y install texlive-publishers #(revtex, etc)

RUN apt-get -y install net-tools

USER ${USER_NAME}
# WORKDIR ${PROJECT_PATH}

CMD ["/bin/bash"]











FROM ubuntu:16.04

MAINTAINER Sergey Cherkashyn <sergey.cherkashyn@iwedia.com>

ARG USER_NAME
ARG GROUP_NAME
ARG USER_ID
ARG GROUP_ID

# Configure a 'docker' sudoer user without password
RUN apt-get update && apt-get -y install sudo
RUN addgroup --gid ${GROUP_ID} ${GROUP_NAME}
RUN adduser --disabled-password --gecos '' --uid ${USER_ID} --gid ${GROUP_ID} ${USER_NAME}
RUN usermod -aG sudo ${USER_NAME}
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# Install build-essential package (for cross compilation)
RUN apt-get update
#RUN apt-get -y install build-essential aptitude

# Install vim, nano
#RUN apt-get -y install vim nano

# Install a bunch of packages for teatro30-ui5 project
RUN apt-get -y install gcc g++ make
#bison flex gettext texinfo patch wget mercurial stow telnet xterm libncurses5-dev expect ddd dos2unix cpio rpm squashfs-tools attr libcap-dev libcap2-bin xinetd tftpd tftp perl xdotool python3-minimal u-boot-tools default-jre realpath nfs-kernel-server gperf gawk git python python-lxml pkg-config bc unzip zip make subversion

# Inteprete sh as bash
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

USER ${USER_NAME}










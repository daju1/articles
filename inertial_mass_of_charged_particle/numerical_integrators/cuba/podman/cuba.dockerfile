FROM ubuntu:16.04

ARG USER_NAME
ARG GROUP_NAME
ARG USER_ID
ARG GROUP_ID

ARG DEBIAN_FRONTEND=noninteractive

# Configure a 'sudo' user without password (so the in-container user matches the host user)
RUN apt-get update && apt-get -y install sudo
RUN if ! getent group ${GROUP_ID} > /dev/null; then \
        addgroup --force-badname --gid ${GROUP_ID} ${GROUP_NAME}; \
    fi
RUN adduser --force-badname --disabled-password --gecos '' --uid ${USER_ID} --gid ${GROUP_ID} ${USER_NAME}
RUN usermod -aG sudo ${USER_NAME}
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# Toolchain for building the Cuba library and compiling the integrator code.
#   gcc / g++ / make      -> compile Cuba 4.2 and the calc_* targets (see cuba.mk)
#   wget / ca-certificates -> build_cuba.sh downloads Cuba-4.2.tar.gz over https
#   tar                   -> unpack the Cuba tarball
RUN apt-get update && apt-get -y install \
    gcc \
    g++ \
    make \
    wget \
    ca-certificates \
    tar

# Update package lists and install Python 3 + Pip
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# (Optional) Verify installation by checking the Python version
RUN python3 --version

# Interpret sh as bash (same as the docker/ environment next to this folder)
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

USER ${USER_NAME}

CMD ["/bin/bash"]

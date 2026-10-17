We use CONTAINER_ENGINE (default podman)
optionally we can use docker

This folder is the podman counterpart of ../docker. It builds an Ubuntu image
with a C/C++ toolchain, then you build the Cuba 4.2 library and the calc_*
integrator targets *inside* the container, so the host stays clean.

# Install podman
```
brew install podman
```

# init and start podman engine
```
podman machine init
podman machine start
```

# check podman installation
```
podman --version
podman info
```

# (optional, corporate network) copy certificate inside the podman machine VM
```
podman machine ssh
```

```
sudo cp /Users/your.name/Downloads/ZscalerRootCertificate-2048-SHA256.crt /etc/pki/ca-trust/source/anchors/
sudo update-ca-trust
```

```
podman machine stop
podman machine start
```

# build the image and start the container
```
./build.sh
./run.sh
```

# If you have docker installed in the system (and a Docker license)
```
CONTAINER_ENGINE=docker ./build.sh
CONTAINER_ENGINE=docker ./run.sh
CONTAINER_ENGINE=docker ./exec.sh
```

# inside the container: build Cuba and run an integrator
`run.sh` drops you into the project root (the "cuba" folder, one level above
this one) with your host user, and the project is mounted at the same path.

```
# 1. build the Cuba 4.2 library (downloads Cuba-4.2.tar.gz, configures, make install -> ./local)
./build_cuba.sh

# 2. build and run one of the integrator targets
make -f cuba.mk clean && make -f Makefile-calc_proton_mass all && ./calc_proton_mass
```

# open an extra shell in the running container
```
./exec.sh
```

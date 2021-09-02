FROM rootproject/root:6.24.02-conda
# FROM uofscphysics/root

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y cmake

ENV PATH /usr/local/physics/bin:$PATH

WORKDIR /usr/local/src

COPY . /usr/local/src/physics_code
RUN mkdir compile \
    && cd compile \
    && cmake /usr/local/src/physics_code -DCMAKE_INSTALL_PREFIX=/usr/local/physics \
    && make -j2 \
    && make install

WORKDIR /data
ENTRYPOINT /bin/bash

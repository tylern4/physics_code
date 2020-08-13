FROM uofscphysics/root

ENV PATH /usr/local/physics/bin:$PATH

WORKDIR /usr/local/src

COPY . /usr/local/src/physics_code
RUN mkdir compile \
    && cd compile \
    && cmake /usr/local/src/physics_code -DCMAKE_INSTALL_PREFIX=/usr/local/physics \
    && make -j$(nproc) \
    && make install

WORKDIR /data
ENTRYPOINT /bin/bash

FROM tylern4/rootv6:latest

ENV PATH /usr/local/physics/bin:$PATH
ENV LD_LIBRARY_PATH /usr/local/physics/lib:/usr/local/lib/root:$LD_LIBRARY_PATH
ENV PYTHONPATH /usr/local/physics/lib:/usr/local/lib/root:$PYTHONPATH

USER root

RUN apt-get update && apt-get install -y python3-pip python3-numpy \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --upgrade cython tqdm \
    && rm -rf /home/builder/.cache

RUN rm -f /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python

WORKDIR /usr/local/src

COPY . /usr/local/src/physics_code
RUN mkdir compile \
    && cd compile \
    && cmake /usr/local/src/physics_code/current -DCMAKE_INSTALL_PREFIX=/usr/local/physics \
    && make -j$(nproc) \
    && make install

WORKDIR /data
ENTRYPOINT /bin/bash

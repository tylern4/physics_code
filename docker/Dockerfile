FROM rootproject/root-ubuntu16

ENV PATH /usr/local/physics/bin:$PATH
ENV LD_LIBRARY_PATH /usr/local/physics/lib:/usr/local/lib/root:$LD_LIBRARY_PATH
ENV PYTHONPATH /usr/local/physics/lib:/usr/local/lib/root:$PYTHONPATH

USER root
COPY current /usr/local/src

RUN apt-get update && apt-get install -y python-pip \ 
    && rm -rf /var/lib/apt/lists/* \
    && pip install --upgrade pip cython tqdm \
    && rm -rf /home/builder/.cache

WORKDIR /usr/local/src
RUN mkdir compile \
    && cd compile \
    && cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/physics \
    && make -j$(nproc) \
    && sudo make install

WORKDIR /data
ENTRYPOINT /bin/bash

FROM rootproject/root-ubuntu16

USER root
COPY current /usr/local/src

RUN apt-get update && apt-get install -y python-pip \ 
    && rm -rf /var/lib/apt/lists/* \
    && pip install --upgrade pip cython \
    && rm -rf /home/builder/.cache

WORKDIR /usr/local/src
RUN mkdir compile \
    && cd compile \
    && cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/physics \
    && make -j$(nproc) \
    && sudo make install

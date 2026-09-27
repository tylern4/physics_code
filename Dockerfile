# syntax=docker/dockerfile:1
#
# Container for the full CLAS E1D p(pi+) chain:
#
#   ROOT files (exp + sim)  ->  csv_maker / csv_maker_mc  ->  CSVs
#                          ->  crossSectionsCSV.py       ->  cross sections
#                          ->  analysis/coarse_xsec.py   ->  coarse xsec + plots
#
# ROOT is built from source rather than taken from a prebuilt image, so this
# file is architecture independent: it works on linux/amd64 and linux/arm64.
# The only arch-dependent thing is how long the ROOT build takes.
#
#   podman build -t e1d:latest .
#   podman buildx build --platform linux/amd64,linux/arm64 -t e1d:latest .
#
# Then:
#   scripts/run_full_chain.sh --exp ... --sim ... --momcorr ... --out ...
#
# Why ROOT from source: the cheaper distribution channels do not cover both
# architectures. Ubuntu 24.04 packages ROOT for amd64 only, and the prebuilt
# rootproject/root images are amd64-only, so an arm64 build has to compile ROOT
# regardless. Building it here rather than pulling it means one Dockerfile
# covers both.
#
# Only the ROOT components the analysis actually needs are built. The code uses
# TChain/TTree/TFile (Core, RIO), TH1D/TH2D (Hist), TLorentzVector (Physics),
# TF1/TMath (Maths), TCanvas/TGraph (Grafpad) and TChain::SetURL (Net). The
# expensive parts of a default ROOT build -- Cling (the interactive C++ and
# Fortran interpreters, which drag in LLVM) and PyROOT -- are disabled, which is
# where most of the build time goes.

# ===========================================================================
# Stage 1: build ROOT from source
# ===========================================================================
FROM docker.io/library/ubuntu:24.04 AS rootbuild

ARG ROOT_VERSION=6.34.00
ARG BUILD_JOBS=4
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl xz-utils bzip2 \
        cmake make g++ gfortran pkg-config \
        libcurl4-openssl-dev libssl-dev zlib1g-dev libexpat1-dev \
        libpng-dev libjpeg-dev libtiff-dev libfreetype6-dev \
        liblz4-dev libzstd-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# The source tarball is used instead of a git clone: no history, and the fetch
# is a single compressed download. ROOT tags the releases v6-34-00, i.e. with
# dashes, so the version above is translated when forming the tag and the
# directory the tarball unpacks into.
#
# --http1.1 and the retry flags are not optional in practice: fetching the
# ~90 MB ROOT tarball over HTTP/2 fails intermittently with
# "curl: (92) HTTP/2 stream ... was not closed cleanly".
RUN root_tag="v$(echo "${ROOT_VERSION}" | tr '.' '-')" \
    && echo "fetching ROOT ${ROOT_VERSION} (${root_tag})" \
    && curl --http1.1 -fsSL --retry 5 --retry-delay 3 --retry-connrefused \
        -o root.tar.gz \
        "https://github.com/root-project/root/archive/refs/tags/${root_tag}.tar.gz" \
    && tar xzf root.tar.gz \
    && mv "root-$(echo "${ROOT_VERSION}" | tr '.' '-')" src \
    && rm root.tar.gz \
    && test -f src/CMakeLists.txt

# cxx/for=OFF drops Cling. python=OFF drops PyROOT. The rest are optional
# features this analysis never touches, disabled so their dependencies are not
# required and their code is not compiled.
RUN cmake -S src -B build \
      -DCMAKE_INSTALL_PREFIX=/opt/root \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_SHARED_LIBS=ON \
      -Droot6=ON \
      -Dbuild-dirsingle=ON \
      -DCMAKE_CXX_STANDARD=17 \
      -Dcore=ON -Dmath=ON -Dhist=ON -Drio=ON -Dnet=ON \
      -Dgrafpad=ON -Dphysics=ON \
      -Dcxx=OFF -Dfor=OFF -Dpython=OFF \
      -Djava=OFF -Djava-cling=OFF \
      -Dx11=OFF -Dglx=OFF -Dchirp=OFF -Dasimage=OFF \
      -Dqt=OFF -Dweb=OFF -Deve=OFF -Dtmtd=OFF \
      -Dmysql=OFF -Dpgsql=OFF -Dsqlite=OFF -Dodbc=OFF -Dxroot=OFF \
      -Dkrb5=OFF -Dldap=OFF -Dgserv=OFF \
      -Dgfal=OFF -Ddcache=OFF -Dgl2d=OFF -Dgl3d=OFF -Demboss=OFF \
      -Dclhep=OFF -Dbuiltin_gsl=OFF -Dgsl=OFF \
      -Dfftw3=ON -Dfftw3_fft=OFF \
      -Dminizip=ON -Dzlib=ON -Dlz4=ON -Dzstd=ON \
      -Drpath=OFF \
      -Dbenchmark=OFF
# NB: do not add -DTEST=OFF here. ROOT's CMakeLists uses the bare word TEST as a
# condition, and defining a cache variable of that name shadows CMake's
# if(TEST) command form; configure then dies at CMakeLists.txt:289 with
# 'if given arguments: "(" "var" "MATCHES" ... Unknown arguments specified'.
    && cmake --build build -j "${BUILD_JOBS}" \
    && cmake --install build \
    && /opt/root/bin/root-config --version

ENV ROOTSYS=/opt/root
ENV PATH=/opt/root/bin:$PATH

# ===========================================================================
# Stage 2: build the analysis code, the MAID Cython module, and the venv
# ===========================================================================
FROM docker.io/library/ubuntu:24.04 AS build

ARG DEBIAN_FRONTEND=noninteractive

# gfortran is needed both for src/*.f (fid_cuts.f, vertex_correction.f) and for
# the MAID2007 Fortran library. libgfortran5-dev supplies the runtime that
# maid_interface.so links against; without it the module builds but fails to
# dlopen on undefined __gfortran_* symbols.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl xz-utils bzip2 \
        cmake make g++ gfortran git \
        python3-dev python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Reuse the ROOT just built instead of installing it again.
COPY --from=rootbuild /opt/root /opt/root
ENV ROOTSYS=/opt/root
ENV PATH=/opt/root/bin:$PATH

RUN gfortran --version | head -1 && cmake --version | head -1 && root-config --version

# --- MAID2007 parameter tables --------------------------------------------
# Not redistributable from this repo (see .gitignore), so fetch them here. About
# 118 MB, needed at run time by the phi bin-centering correction. Pinned to a
# commit so the image is reproducible.
ARG SPP_TBL_URL=https://github.com/tylern4/onepigen/raw/7ad5c99aa5b526dd73311500d680458efb793494/spp_tbl.tar.gz
RUN mkdir -p /src/physics_code/exe/maid_xsection \
    && curl --http1.1 -fsSL --retry 5 --retry-delay 3 --retry-connrefused \
        "${SPP_TBL_URL}" -o /tmp/spp_tbl.tar.gz \
    && tar xzf /tmp/spp_tbl.tar.gz -C /src/physics_code/exe/maid_xsection \
    && rm /tmp/spp_tbl.tar.gz \
    && ls /src/physics_code/exe/maid_xsection/spp_tbl/ | wc -l

WORKDIR /src/physics_code

# .dockerignore keeps this to source: no .venv, no build trees, no ROOT files.
COPY . /src/physics_code

# The MAID Cython extension needs the aao_rad_lund submodule, which is absent
# from a plain COPY of a fresh clone. lz4xx and zstr are not needed: their
# add_subdirectory() calls are commented out in the top-level CMakeLists.txt.
RUN if [ -d .git ] && [ ! -f external/aao_rad_lund/CMakeLists.txt ]; then \
        git submodule update --init external/aao_rad_lund; \
    fi; \
    test -f external/aao_rad_lund/CMakeLists.txt \
    || { echo "ERROR: external/aao_rad_lund missing."; \
         echo "Run: git submodule update --init external/aao_rad_lund"; exit 1; }

# --- python environment ----------------------------------------------------
# A venv on the system interpreter, so the interpreter and the Cython module
# agree on the libpython they were built against.
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/venv/bin/pip install --no-cache-dir \
        numpy pandas pyarrow boost-histogram scipy lmfit matplotlib tqdm \
        cython setuptools

# --- C++ analysis code -----------------------------------------------------
# NUM_THREADS is only a default; csv_maker reads it from the environment.
ENV NUM_THREADS=8
RUN cmake -S /src/physics_code -B /src/physics_code/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DPython3_EXECUTABLE=/opt/venv/bin/python \
        -DROOT_DIR=/opt/root \
    && cmake --build /src/physics_code/build -j "$(nproc)" \
    && ls -la /src/physics_code/build/csv_maker /src/physics_code/build/csv_maker_mc

# --- MAID Cython extension -------------------------------------------------
# CYTHON_EXECUTABLE must be explicit: with Cython only inside the venv,
# find_package(Cython) reports not-found and the .pyx files are skipped with no
# error at all, so the build looks successful but produces no module.
RUN cmake -S external/aao_rad_lund -B external/aao_rad_lund/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DPython3_EXECUTABLE=/opt/venv/bin/python \
        -DCYTHON_EXECUTABLE=/opt/venv/bin/cython \
    && cmake --build external/aao_rad_lund/build -j "$(nproc)" \
    && ls -la external/aao_rad_lund/build/maid_interface.so

# --- smoke test the extension before shipping the stage --------------------
# A bad dlopen only shows at the first MAID call, which would otherwise surface
# halfway through a multi-hour run.
RUN cd /src/physics_code/exe/maid_xsection \
    && PYTHONPATH=/src/physics_code/external/aao_rad_lund/build \
       CLAS_PARMS=/src/physics_code/exe/maid_xsection \
       /opt/venv/bin/python -c "\
from maid_interface import maid_2007_Npi as m; \
s = m(4.81726, 1.23, 0.5, 0.0, 90.0); \
print('maid_interface OK, sigma =', s); \
assert s == 18.7191162109375, f'unexpected value {s}'"

# ===========================================================================
# Stage 3: runtime
# ===========================================================================
FROM docker.io/library/ubuntu:24.04 AS runtime

ARG DEBIAN_FRONTEND=noninteractive

# Same reasoning as the builder: libgfortran is needed at run time by the
# executables and maid_interface.so.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        gfortran libgfortran5-dev \
        python3 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=rootbuild /opt/root /opt/root

# Default uid/gid. 1001 is the owner of the NFS export, so files written back
# to the share keep the right ownership instead of appearing as root. Override
# with --user.
ARG RUN_UID=1001
ARG RUN_GID=1001
RUN groupadd -g "${RUN_GID}" e1d 2>/dev/null || true \
    && useradd -u "${RUN_UID}" -g "${RUN_GID}" -m -s /bin/bash e1d 2>/dev/null || true

ENV ROOTSYS=/opt/root \
    HOME=/home/e1d \
    TMPDIR=/tmp \
    NUM_THREADS=8 \
    PATH=/opt/e1d/bin:/opt/venv/bin:/opt/root/bin:$PATH \
    PYTHONPATH=/opt/e1d/lib:/opt/e1d/python \
    CLAS_PARMS=/opt/e1d

WORKDIR /work

# Executables on PATH so the chain script can just call `csv_maker`.
COPY --from=build /src/physics_code/build/ /opt/e1d/bin/
# The MAID Cython module, where calc_xsections expects to find it.
COPY --from=build /src/physics_code/external/aao_rad_lund/build/maid_interface.so \
     /opt/e1d/lib/maid_interface.so
# Parameter tables, and the analysis package that resolves them via CLAS_PARMS.
COPY --from=build /src/physics_code/exe/maid_xsection/spp_tbl/ /opt/e1d/spp_tbl/
COPY --from=build /src/physics_code/python/ /opt/e1d/python/
COPY --from=build /src/physics_code/analysis/ /opt/e1d/analysis/
COPY --from=build /opt/venv /opt/venv

COPY scripts/run_full_chain.sh /opt/e1d/bin/run_full_chain.sh
RUN chmod +x /opt/e1d/bin/run_full_chain.sh

# Root files are large; put the working area somewhere with room.
VOLUME ["/work"]

USER ${RUN_UID}:${RUN_GID}

# Prove the image is wired up: the executables and the Python module respond and
# MAID actually loads, so a failure surfaces at build time not mid-analysis.
RUN set -eux; \
    root-config --version; \
    csv_maker --help 2>&1 | head -1; \
    /opt/venv/bin/python -c "import calc_xsections as c; assert c.maid is not None, 'MAID not loaded'; print('calc_xsections OK, MAID loaded')"

CMD ["/opt/e1d/bin/run_full_chain.sh", "--help"]

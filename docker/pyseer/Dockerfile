# syntax=docker/dockerfile:1
FROM ubuntu:22.04
WORKDIR /opt/pyseer

ARG PYTHON_VERSION=3.11
RUN apt update -q -y \
&& apt install -y \
  python${PYTHON_VERSION}-dev \
  python3-pip \
  git \
  gfortran \
  bwa \
  bedtools \
  bedops \
&& rm -rf /var/lib/apt/lists/* \
&& ln -sf "python${PYTHON_VERSION}" /usr/bin/python3 \
&& ln -sf python3 /usr/bin/python

RUN python -m pip install pybedtools

ARG PYSEER_VERSION
RUN git clone --depth 1 --branch "$PYSEER_VERSION" https://github.com/mgalardini/pyseer . \
&& rm -rf .git \
&& python -m pip install . \
&& chmod 755 scripts/* 

RUN find . -name '*.py' -exec \
  sed -i "1i #!/usr/bin/env python" {} +

ENV PATH="/opt/pyseer/scripts:${PATH}"
FROM ubuntu:16.04

MAINTAINER Bastian Rieck

RUN apt-get update && apt-get install -y \
  build-essential                        \ # Just to be on the safe side
  cmake                                  \
  git                                    \
  libboost-dev                           \
  libboost-regex-dev                     \
  libeigen3-dev                          \
  libflann-dev                           \
  libgomp1                               \
  libtinyxml2-dev                        \
  python3

WORKDIR /tmp

# Somewhat superfluous because the client most likely downloaded this
# file already by cloning the repo.
RUN git clone https://github.com/Submanifold/Aleph.git \
  && cd Aleph                                          \
  && mkdir build                                       \
  && cd build                                          \
  && cmake ../ && make                                 \
  && make install

# Cleanup to reduce docker image size
RUN rm -rf /var/lib/apt/lists/*

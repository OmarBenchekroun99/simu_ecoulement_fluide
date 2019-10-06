# Dockerfile for gitlab continuous integration
FROM ubuntu:bionic
MAINTAINER Jean-Baptiste.Keck@imag.fr

# upgrade initial image
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get full-upgrade -y

# get build tools and required libraries
RUN apt-get install -y make
RUN apt-get install -y scilab
RUN apt-get install -y ffmpeg
RUN apt-get install -y texlive-full
RUN apt-get install -y openjdk-8-jdk
RUN apt-get install -y evince

# clean cached packages
RUN rm -rf /var/lib/apt/lists/*
RUN rm -rf /tmp/*

CMD ["/bin/bash"]

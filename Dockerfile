FROM ubuntu:20.04
LABEL author="Rodrigo Martin <rodrigo.martin@bsc.es>"

# Install dependencies
RUN apt-get update && apt-get upgrade -y && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    python3 \
    python3-pip \
    build-essential \
    libz-dev \
    libglib2.0-dev \
    libbz2-dev \
    liblzma-dev \
    autoconf \
    samtools

# Python dependencies
RUN pip install pysam variant-extractor


# Copy scripts from src
COPY src/ /genome-mosaic-maker/
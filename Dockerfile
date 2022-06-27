FROM ubuntu:bionic
MAINTAINER T.N. Wylie <twylie@wustl.edu>

###############################################################################
#                                  ViroSearch                                 #
###############################################################################

LABEL \
description = "ViroSearch: EBV variant review pipeline."

RUN apt-get update -y && apt-get upgrade -y && apt-get install -y \
    build-essential \
    bzip2 \
    cmake \
    git \
    libnss-sss \
    libtbb2 \
    libtbb-dev \
    liblzma-dev \
    ncurses-dev \
    nodejs \
    python3.7-dev \
    python3-pip \
    unzip \
    wget \
    zlib1g \
    zlib1g-dev \
    curl \
    zsh \
    libbz2-dev \
    autoconf \
    automake \
    zile \
    default-jre \
    varscan \
    bcftools \
    && apt-get clean

# Install BWA for sequence alignments.

RUN \
git clone https://github.com/lh3/bwa.git virosearch_install/bwa && \
cd virosearch_install/bwa && \
make && \
cp bwa /usr/bin

# Install SAMTOOLS for SAM/BAM handling.

RUN \
curl -L https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 > virosearch_install/samtools-1.9.tar.bz2 && \
bzip2 -d virosearch_install/samtools-1.9.tar.bz2 && \
tar xvf virosearch_install/samtools-1.9.tar -C virosearch_install/ && \
cd virosearch_install/samtools-1.9 && \
./configure --without-curses --disable-lzma && \
make && \
cp samtools /usr/bin

# Install VSEARCH for low-complexity filter methods.

RUN \
git clone https://github.com/torognes/vsearch.git virosearch_install/vsearch && \
cd virosearch_install/vsearch && \
./autogen.sh && \
./configure && \
make && \
make install

# BBMap suite installation.

RUN cd virosearch_install/ && \
wget 'https://sourceforge.net/projects/bbmap/files/BBMap_38.96.tar.gz/download' && \
tar xvfz download && \
cp -pr bbmap/* /usr/bin

# File clean-up.

RUN \
cd / && \
\rm -rf virosearch_install

# Install Snakemake for pipeline development.

RUN python3.7 -m pip install --upgrade setuptools && python3.7 -m pip install --upgrade pip

RUN \
python3.7 -m pip install snakemake==7.8.3

# Python installations.

RUN \
python3.7 -m pip install pandas && \
python3.7 -m pip install pyyaml && \
python3.7 -m pip install pysam && \
python3.7 -m pip install vcfpy && \
python3.7 -m pip install cutadapt

# Python installs.

RUN python3.7 -m pip install biopython

# Installing Python code from local instance.

COPY ./bin/filter_pidv_pe.py /usr/bin/
COPY ./bin/virosearch.py /usr/bin/virosearch
COPY ./virosearch/viromatch /usr/lib/python3.7/viromatch
COPY ./virosearch/lib /usr/lib/python3.7/virosearch
COPY ./virosearch/recipes /usr/lib/python3.7/virosearch/recipes
COPY ./virosearch/washu /usr/lib/python3.7/virosearch/washu

# __END__

# T.N. Wylie  <twylie@wustl.edu>
# Last Update: Mon Jun 27 09:39:44 CDT 2022


#https://hub.docker.com/r/jsschrepping/velocyto-docker/dockerfile

FROM debian:buster
LABEL \
    author="Shixiang Wang" \
    maintainer="Shixiang Wang and Yankun Zhao" \
    email="w_shixiang@163.com" \
    description="Modified Docker Image for MetaFusion" \
    org.label-schema.license="Academic Free License v.3.0"
ENV DEBIAN_FRONTEND noninteractive
ENV PATH /opt/conda/bin:$PATH

# System deps
RUN sed -i 's/deb.debian.org/mirrors.aliyun.com/g' /etc/apt/sources.list && \
    apt-get update --yes && \
    apt-get install -y git sudo apache2 python2.7 python-pip wget && \
    apt-get install --yes --no-install-recommends \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    gcc \
    mono-mcs \
    build-essential \
    cmake \
    less \
    libbamtools-dev \
    libboost-dev \
    libboost-iostreams-dev \
    libboost-log-dev \
    libboost-system-dev \
    libboost-test-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libz-dev \
    nano \
    pandoc && \
    rm -rf /var/lib/apt/lists/* && \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/*

# Install Miniconda2
RUN wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda2-4.7.12-Linux-x86_64.sh && \
    bash Miniconda2-4.7.12-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda2-4.7.12-Linux-x86_64.sh && \
    echo "channels:\n  - defaults\nshow_channel_urls: true\ndefault_channels:\n  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main\ncustom_channels:\n  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud\n  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud" > ~/.condarc && \
    cat ~/.condarc

# Install softwares with Conda, Pip and R
RUN conda update -n base conda && \
    conda install -c conda-forge -y r-base=4.1.0 && \
    conda install -c bioconda -c conda-forge -y bedtools==2.27 && \
    R -e 'options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")); install.packages("BiocManager"); BiocManager::install("RBGL"); BiocManager::install("data.table")' && \
    pip config set global.index-url https://mirrors.aliyun.com/pypi/simple/ && \
    pip install pybedtools==0.8.1 && \
    pip install pandas==0.24.2

# Perl
COPY cpanmin.us.txt /opt/
RUN perl --version && PERL_MM_USE_DEFAULT=1 perl -MCPAN -e 'CPAN::HandleConfig->edit("urllist", "unshift", "https://mirrors.tuna.tsinghua.edu.cn/CPAN/"); mkmyconfig' && \
    cat /opt/cpanmin.us.txt | perl - App::cpanminus && \
    cpanm Set::IntervalTree JSON::XS && \
    apt autoremove -y && apt clean -y && apt purge -y && rm -rf /tmp/* /var/tmp/* && \
    conda clean -ay

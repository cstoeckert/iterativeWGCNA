FROM rocker/r-ver:3.3.3

RUN apt-get update && apt-get install -y build-essential wget unzip curl python python-dev python-matplotlib libicu-dev libssl-dev libffi-dev libxml2-dev libxslt1-dev zlib1g-dev libreadline-dev && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

ENV LD_LIBRARY_PATH="/lib/x86_64-linux-gnu:/lib/x86_64-linux-gnu/:/usr/local/lib/R/lib"
ENV LDFLAGS="-L/lib/x86_64-linux-gnu/:/usr/lib/x86_64-linux-gnu/:/usr/local/lib/R/lib"

RUN R -e "source('http://bioconductor.org/biocLite.R');biocLite(c('GO.db', 'preprocessCore', 'impute', 'AnnotationDbi'));install.packages(c('data.table','matrixStats', 'checkmate', 'htmlTable', 'Hmisc', 'WGCNA'))"

COPY . /usr/local/iterativeWGCNA

# For some reason the sym links for these are missing so adding ....
WORKDIR /lib/x86_64-linux-gnu
RUN ln -s libpcre.so.3 libpcre.so; ln -s liblzma.so.5 liblzma.so; ln -s libbz2.so.1 libbz2.so

WORKDIR /usr/local/iterativeWGCNA
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
RUN python get-pip.py
RUN pip install rpy2==2.7.9 --force-reinstall
RUN pip install iterativeWGCNA

WORKDIR /home/docker

ENTRYPOINT ["iterativeWGCNA"]


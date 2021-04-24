FROM r-base:latest

RUN apt-get update && apt-get install -y build-essential wget unzip python3 python python-setuptools python3-pip python3-rpy2 && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('WGCNA')"

COPY . /usr/local/iterativeWGCNA

WORKDIR /usr/local/iterativeWGCNA

RUN pip install iterativeWGCNA

WORKDIR /home/docker

ENTRYPOINT ["iterativeWGCNA"]



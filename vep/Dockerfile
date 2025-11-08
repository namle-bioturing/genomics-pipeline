FROM ensemblorg/ensembl-vep:release_115.2

USER root

RUN apt-get update -y && apt install -y less git bedtools bcftools wget python3 python3-pip python3-setuptools python3-dev

WORKDIR /workspace
FROM continuumio/miniconda3
  
RUN conda update -n base -c defaults conda && \
    conda install -n base -c omnia -c conda-forge openmm -y && \
    conda clean --yes --all

RUN git clone --branch v1.7-docker https://github.com/ozzzzz/pdbfixer.git && \
    cd /pdbfixer && \
    python setup.py install && \
    cd / && \
    rm -rf /pdbfixer

ENTRYPOINT pdbfixer
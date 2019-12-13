FROM continuumio/miniconda3

RUN conda update -n base -c defaults conda

RUN conda install -c omnia -c conda-forge openmm -y

RUN git clone https://github.com/pandegroup/pdbfixer.git

WORKDIR /pdbfixer

RUN python setup.py install

ENTRYPOINT pdbfixer

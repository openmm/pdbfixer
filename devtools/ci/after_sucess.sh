#!/bin/bash

echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" == "true" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi


if [[ "2.7 3.3" =~ "$python" ]]; then
    conda install --yes binstar jinja2
    binstar -t $BINSTAR_TOKEN upload --force -u omnia -p pdbfixer-dev $HOME/miniconda/conda-bld/*/pdbfixer-dev-*.tar.bz2
fi

if [[ "$python" != "2.7" ]]; then
    echo "No deploy on PYTHON_VERSION=${python}"; exit 0
fi

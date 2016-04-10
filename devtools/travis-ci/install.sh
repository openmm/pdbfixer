# Temporarily change directory to $HOME to install software
pushd .
cd $HOME

# Install Miniconda
MINICONDA=Miniconda2-latest-Linux-x86_64.sh
MINICONDA_HOME=$HOME/miniconda
MINICONDA_MD5=$(curl -s http://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget -q http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b -p $MINICONDA_HOME

# Configure miniconda
export PIP_ARGS="-U"
export PATH=$MINICONDA_HOME/bin:$PATH
conda update --yes conda
conda install --yes conda-build jinja2 anaconda-client pip

# Restore original directory
popd

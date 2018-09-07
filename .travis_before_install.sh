#!/bin/bash
if [ $TRAVIS_OS_NAME == 'linux' ]; then
  echo "Installing deps for linux"
  sudo add-apt-repository ppa:nschloe/swig-backports -y
  sudo apt-get -qq update
  sudo apt-get install -y swig3.0

elif [ $TRAVIS_OS_NAME == 'osx' ]; then
  echo "Installing deps for OSX"
  if [ $PYTHON_VERSION == "2.7" ]; then
    CONDA_VER='2'
    CONDA_URL="https://repo.continuum.io/miniconda/Miniconda${CONDA_VER}-latest-MacOSX-x86_64.sh"
  elif [ $PYTHON_VERSION == "3.7" ]; then
    CONDA_VER='3'
    CONDA_URL="https://repo.continuum.io/miniconda/Miniconda${CONDA_VER}-latest-MacOSX-x86_64.sh"
  # NOTE: Miniconda stopped supporting 3.6, but pysam is late...
  elif [ $PYTHON_VERSION == "3.6" ]; then
    CONDA_VER='3'
    CONDA_URL='https://repo.continuum.io/miniconda/Miniconda3-4.5.4-MacOSX-x86_64.sh'
  else
    echo "Conda only supports 2.7 and 3.6"
  fi
  curl "${CONDA_URL}" -o "miniconda.sh"
  bash "miniconda.sh" -b -p $HOME/miniconda
  echo "$PATH"
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate
  # Use pip from conda
  conda install -y pip
  pip --version

else
  echo "OS not recognized: $TRAVIS_OS_NAME"
  exit 1
fi

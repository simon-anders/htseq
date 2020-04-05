#!/bin/bash
if [ $TRAVIS_OS_NAME == 'linux' ]; then
  echo "Installing deps for linux"
  sudo add-apt-repository ppa:nschloe/swig-backports -y
  sudo apt-get -qq update
  sudo apt-get install -y swig3.0

elif [ $TRAVIS_OS_NAME == 'osx' ]; then
  echo "Find out OSX version"
  osx_version=$(sw_vers -productVersion)
  echo "OSX version: $osx_version"

  echo "Installing deps for OSX"
  CONDA_URL="https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
  curl "${CONDA_URL}" -o "miniconda.sh"
  bash "miniconda.sh" -b -p $HOME/miniconda
  echo "$PATH"
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate

  # Make conda environment and activate
  conda create -n travis python=$CONDA_PY
  conda activate travis

  # Use pip from conda
  conda install -y pip
  pip --version

else
  echo "OS not recognized: $TRAVIS_OS_NAME"
  exit 1
fi

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
  osx_version1=$(echo $osx_version | cut -d. -f1)
  osx_version2=$(echo $osx_version | cut -d. -f2)
  if [ $osx_version1 == "10" ] && [ $osx_version2 -ge 14 ]; then
    echo "Installing C headers for OSX Mojave and later"
    #open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg
    sudo installer -pkg /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg -target /
  fi

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

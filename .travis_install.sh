#!/bin/bash
# try and make wheels
if [ $DOCKER_IMAGE ]; then
  docker run --rm -v `pwd`:/io $DOCKER_IMAGE /io/buildwheels.sh 
  if [ $? != 0 ]; then
      exit 1
  fi
  ls wheelhouse/
  if [ $? != 0 ]; then
      exit 1
  fi

# compile normally
else
  if [ $TRAVIS_OS_NAME == 'osx' ]; then
    export PATH="$HOME/miniconda/bin:$PATH"
    source $HOME/miniconda/bin/activate
  fi

  # setuptools < 18.0 has issues with Cython as a dependency
  pip install "$CYTHON_INSTALL"
  
  if [ $TRAVIS_OS_NAME == 'linux' ]; then
    sed -i "s|pysam>=0.9.0|$PYSAM_VERSION|" requirements.txt
  elif [ $TRAVIS_OS_NAME == 'osx' ]; then
    sed -i "" "s|pysam>=0.9.0|$PYSAM_VERSION|" requirements.txt
  else
    echo "OS not recognized: $TRAVIS_OS_NAME"
    exit 1
  fi
  if [ $? != 0 ]; then
      exit 1
  fi

  pip install -r requirements.txt
  
  pip install $PYPI HTSeq
  if [ $? != 0 ]; then
      exit 1
  fi
fi

## OSX makes wheels as well
#if [ $TRAVIS_OS_NAME == 'osx' ]; then
#  mkdir wheelhouse
#  pip wheel . -w wheelhouse/
#  if [ $? != 0 ]; then
#      exit 1
#  fi
#fi

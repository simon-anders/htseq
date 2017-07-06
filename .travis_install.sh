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
  # setuptools < 18.0 has issues with Cython as a dependency
  pip install Cython
  if [ $? != 0 ]; then
      exit 1
  fi
  
  if [ $TRAVIS_OS_NAME == 'linux' ]; then
    sed -i "s|pysam>=0.9.0|$PYSAM_VERSION|" setup.py
  elif [ $TRAVIS_OS_NAME == 'osx' ]; then
    sed -i "" "s|pysam>=0.9.0|$PYSAM_VERSION|" setup.py
  else
    echo "OS not recognized: $TRAVIS_OS_NAME"
    exit 1
  fi
  if [ $? != 0 ]; then
      exit 1
  fi
  
  # old setuptools also has a bug for extras, but it still compiles
  pip install -v '.[htseq-qa]'
  if [ $? != 0 ]; then
      exit 1
  fi
fi

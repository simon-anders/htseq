#!/bin/bash
if [ $TRAVIS_OS_NAME == 'linux' ]; then
  echo "Installing deps for linux"
  sudo add-apt-repository ppa:nschloe/swig-backports -y
  sudo apt-get -qq update
  sudo apt-get install -y swig3.0
elif [ $TRAVIS_OS_NAME == 'osx' ]; then
  echo "Installing deps for OSX"
  #TODO
else
 echo "OS not recognized: $TRAVIS_OS_NAME"
 exit 1
fi

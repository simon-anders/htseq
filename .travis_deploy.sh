#!/bin/bash
# only deploy builds for a release_<sematic-version>_RC?? tag to testpypi
if [ -z $TRAVIS_TAG ]; then
  exit 0
fi
TAG1=$(echo $TRAVIS_TAG | cut -f1 -d_)
TAG2=$(echo $TRAVIS_TAG | cut -f2 -d_)
TAG3=$(echo $TRAVIS_TAG | cut -f3 -d_)
if [ -z $TAG2 ] || [ -z $TAG3 ]; then
  exit 0;
fi
if [ $TAG1 != 'release' ] || [ $TAG2 != $(cat VERSION) ]; then
  exit 0;
fi

# do not deploy outside of manylinux1
if [ -z $DOCKER_IMAGE ]; then
  exit 0
fi


# deploy onto pypitest unless you have no RC
if [ ${TAG3:0:2} == 'RC' ]; then
  TWINE_REPOSITORY='https://testpypi.python.org/pypi'
  echo 'Deploying to testpypi'
else
  #FIXME
  #TWINE_REPOSITORY='https://pypi.python.org/pypi'
  echo 'Deploying to production pypi'
fi
   
# Wheels are already tested in docker image
docker run -e TWINE_REPOSITORY -e TWINE_USERNAME -e TWINE_PASSWORD --rm -v $(pwd):/io $DOCKER_IMAGE /io/deloywheels.sh 

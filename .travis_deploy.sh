#!/bin/bash
# only deploy builds for a release_<sematic-version>_RC?? tag to testpypi
if [ -z $TRAVIS_TAG ]; then
  echo 'No TRAVIS_TAG, exit'
  exit 0
fi
TAG1=$(echo $TRAVIS_TAG | cut -f1 -d_)
TAG2=$(echo $TRAVIS_TAG | cut -f2 -d_)
TAG3=$(echo $TRAVIS_TAG | cut -f3 -d_)
if [ -z $TAG2 ]; then
  echo 'No TAG2, exit'
  exit 0;
fi
if [ $TAG1 != 'release' ] || [ $TAG2 != $(cat VERSION) ]; then
  echo 'No release tag or wrong version, exit'
  exit 0;
fi

# do not deploy on linux outside of manylinux
if [ -z $DOCKER_IMAGE ] && [ $TRAVIS_OS_NAME != 'osx' ]; then
  echo 'Not inside manylinux docker image and not OSX, exit'
  exit 0
fi

# OSX only deploys on latest pysam
if [ $TRAVIS_OS_NAME == 'osx' ] && [ $PYSAM_VERSION != 'pysam>=0.13.0' ]; then
  echo "OSX on older pysam, exit"
  exit 0
fi

# deploy onto pypitest unless you have no RC
if [ -z $TAG3 ]; then
  TWINE_PASSWORD=${TWINE_PASSWORD_PYPI}
  TWINE_REPOSITORY='https://upload.pypi.org/legacy/'
  echo 'Deploying to production pypi'
elif [ ${TAG3:0:2} == 'RC' ]; then
  TWINE_PASSWORD=${TWINE_PASSWORD_TESTPYPI}
  TWINE_REPOSITORY='https://test.pypi.org/legacy/'
  echo 'Deploying to testpypi'
else
  echo "Tag not recognized: $TRAVIS_TAG"
  exit 1
fi
   
if [ $DOCKER_IMAGE ]; then
  # Wheels are already tested in docker image
  docker run -e TWINE_REPOSITORY="$TWINE_REPOSITORY" -e TWINE_USERNAME="$TWINE_USERNAME" -e TWINE_PASSWORD="$TWINE_PASSWORD" --rm -v $(pwd):/io $DOCKER_IMAGE /io/deploywheels.sh 

elif [ $TRAVIS_OS_NAME == 'osx' ]; then
  # OSX deployment
  echo "Deploying for OSX"

  # Only deploy on 10.14 to ensure 10.9+ compatibility and Mojave header/linker changes
  osx_version=$(sw_vers -productVersion)
  echo "OSX version: $osx_version"
  osx_ver1=$(echo $osx_version | cut -d. -f1)
  osx_ver2=$(echo $osx_version | cut -d. -f2)
  if [ $osx_ver1 -lt 10 ] || [ $osx_ver2 -lt 14 ]; then
    echo "OSX version not for deployment: $osx_version"
    exit 1
  fi

  HTSEQ_VERSION=$(cat VERSION)
  echo "TWINE_REPOSITORY=$TWINE_REPOSITORY"
  echo "TWINE_USERNAME=$TWINE_USERNAME"
  echo "TWINE_PASSWORD=$TWINE_PASSWORD"
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate
  conda activate travis

  pip --version
  pip install twine

  # Figure out architecture string
  PYVER=$(echo $CONDA_PY | sed 's/\.//')
  PYARCH=cp${PYVER}-cp${PYVER}m

  echo "Contents of wheelhouse:"
  ls wheelhouse
  TWINE_WHEEL=$(ls wheelhouse/HTSeq-${HTSEQ_VERSION}-${PYARCH}*.whl)
  echo "TWINE_WHEEL=$TWINE_WHEEL"

  echo "Uploading..."
  twine upload  --repository-url "${TWINE_REPOSITORY}" -u "${TWINE_USERNAME}" -p "${TWINE_PASSWORD}" "${TWINE_WHEEL}"
  if [ $? != 0 ]; then
    echo "Upload failed" 
    exit 1
  fi
  echo "Upload complete"

else
  echo "No DOCKER_IMAGE and not OSX, we should not be here!"
  exit 1
fi

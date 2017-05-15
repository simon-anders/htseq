#!/bin/bash
# Wheels are already tested in docker image
if [ $DOCKER_IMAGE ]; then
  docker run --rm -v $(pwd):/io $DOCKER_IMAGE /io/testwheels.sh
  if [ $? != 0 ]; then
      exit 1
  fi
fi

PYTHON=${PYTHON:-python}
echo "python: ${PYTHON}"

py_fdn="python$(python -c 'from __future__ import print_function; import sys; print(sys.version_info[0])')/"
echo "py_fdn: ${py_fdn}"

echo 'Running tests...'

echo 'Doctests...'
${PYTHON} "${py_fdn}test/test.py"
if [ $? != 0 ]; then
    exit 1
fi
echo 'done!'

echo 'htseq-count...'
${PYTHON} "${py_fdn}test/test_htseq-count.py"
if [ $? != 0 ]; then
    exit 1
fi
echo 'done!'

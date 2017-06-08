#!/bin/bash
# Wheels are already tested in docker image
if [ $DOCKER_IMAGE ]; then
  docker run --rm -v $(pwd):/io $DOCKER_IMAGE /io/testwheels.sh
  exit $?
fi

PYTHON=${PYTHON:-python}
echo "python: ${PYTHON}"

py_fdn="python$(${PYTHON} --version 2>&1 | cut -f2 -d ' ' | cut -f1 -d'.')/"
echo "py_fdn: ${py_fdn}"

echo 'Running tests...'

echo 'General tests...'
${PYTHON} "${py_fdn}test/test_other.py"
if [ $? != 0 ]; then
    exit 1
fi
echo 'done!'

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

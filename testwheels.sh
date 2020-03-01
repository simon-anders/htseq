#!/bin/bash
#
# Build manylinux wheels for HTSeq. Based on the example at
# <https://github.com/pypa/python-manylinux-demo>
#
# It is best to run this in a fresh clone of the repository!
#
# Run this within the repository root:
#   docker run --rm -v $(pwd):/io quay.io/pypa/manylinux2010_x86_64 /io/buildwheels.sh
#
# The wheels will be put into the wheelhouse/ subdirectory.
#
# For interactive tests:
#   docker run -it -v $(pwd):/io quay.io/pypa/manylinux2010_x86_64 /bin/bash

set -xeuo pipefail

# For convenience, if this script is called from outside of a docker container,
# it starts a container and runs itself inside of it.
if ! grep -q docker /proc/1/cgroup; then
  # We are not inside a container
  exec docker run --rm -v $(pwd):/io quay.io/pypa/manylinux2010_x86_64 /io/$0
fi

# Python 2.6 is not supported
rm -rf /opt/python/cp26*
rm -rf /opt/python/cpython-2.6*

# Python 2.7 is deprecated
rm -rf /opt/python/cp27*
rm -rf /opt/python/cpython-2.7*

# Python 3.3-4 is not supported:
rm -rf /opt/python/cp33*
rm -rf /opt/python/cp34*

# Install packages and test them
PYBINS="/opt/python/*/bin"
for PYBIN in ${PYBINS}; do
    ${PYBIN}/pip install HTSeq --no-index -f /io/wheelhouse
    (cd /io; ls; DOCKER_IMAGE='' PYTHON=${PYBIN}/python PATH=${PYBIN}:${PATH} ./.travis_test.sh)
    if [ $? != 0 ]; then
        exit 1
    fi
done


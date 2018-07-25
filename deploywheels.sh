#!/bin/bash
#
# Build manylinux1 wheels for HTSeq. Based on the example at
# <https://github.com/pypa/python-manylinux-demo>
#
# It is best to run this in a fresh clone of the repository!
#
# Run this within the repository root:
#   docker run --rm -v $(pwd):/io quay.io/pypa/manylinux1_x86_64 /io/buildwheels.sh
#
# The wheels will be put into the wheelhouse/ subdirectory.
#
# For interactive tests:
#   docker run -it -v $(pwd):/io quay.io/pypa/manylinux1_x86_64 /bin/bash

set -xeuo pipefail

# For convenience, if this script is called from outside of a docker container,
# it starts a container and runs itself inside of it.
if ! grep -q docker /proc/1/cgroup; then
  # We are not inside a container
  exec docker run --rm -v $(pwd):/io quay.io/pypa/manylinux1_x86_64 /io/$0
fi

# Python 2.6 is not supported
rm -rf /opt/python/cp26*

# Python 3.3 is not supported:
rm -rf /opt/python/cp33*

# Python 3.7 is not supported YET:
rm -rf /opt/python/cp37*

# Deploy binary packages
HTSEQ_VERSION=$(cat /io/VERSION)
PYBINS="/opt/python/*/bin"
ERRS=0
for PYBIN in ${PYBINS}; do
  PYVER=$(basename $(dirname ${PYBIN}))
  echo "PYVER=$PYVER"
  echo "TWINE_REPOSITORY=$TWINE_REPOSITORY"
  echo "TWINE_USERNAME=$TWINE_USERNAME"
  echo "TWINE_PASSWORD=$TWINE_PASSWORD"
  ${PYBIN}/pip install twine
  ${PYBIN}/twine upload --repository-url "${TWINE_REPOSITORY}" -u "${TWINE_USERNAME}" -p "${TWINE_PASSWORD}" /io/wheelhouse/HTSeq-${HTSEQ_VERSION}-${PYVER}-manylinux1_x86_64.whl
  if [ $? != 0 ]; then
    ERRS=1
  fi
done

# Deploy source code
${PYBIN}/twine upload --repository-url "${TWINE_REPOSITORY}" -u "${TWINE_USERNAME}" -p "${TWINE_PASSWORD}" /io/wheelhouse/HTSeq-${HTSEQ_VERSION}.tar.gz
if [ $? != 0 ]; then
  ERRS=1
fi
exit $ERRS

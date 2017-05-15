#!/bin/bash
py_fdn="python$(python -c 'from __future__ import print_function; import sys; print(sys.version_info[0])')"
echo "py_fdn: ${py_fdn}"

echo 'Running tests...'

echo 'Doctests...'
python "${py_fdn}test/test.py"
if [ $? != 0 ]; then
    exit 1
fi
echo 'done!'

echo 'htseq-count...'
python "${py_fdn}test/test_htseq-count.py"
if [ $? != 0 ]; then
    exit 1
fi
echo 'done!'

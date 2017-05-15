#!/bin/bash
py_fdn=python$(python --version | cut -f2 -d' ' | cut -f1 -d'.')
echo "py_fdn: ${py_fdn}/"

echo 'Running tests...'

echo 'Doctests...'
python "${py_fdn}test/test.py"
echo 'done!'

echo 'htseq-count...'
python "${py_fdn}test/test_htseq-count.py"
echo 'done!'

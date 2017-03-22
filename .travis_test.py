#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import subprocess

py_fdn = 'python'+str(sys.version_info[0])+'/'
print('py_fdn:', py_fdn)
print('Symlink tests and docs...')
for fdn in ['doc', 'test']:
    if os.path.islink(fdn):
        os.unlink(fdn)
    os.symlink(py_fdn+fdn, fdn)
print('Running tests...')
subprocess.call('python test/test.py', shell=True)

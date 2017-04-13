#!/usr/bin/env python
from __future__ import print_function
import sys
import subprocess

py_fdn = 'python'+str(sys.version_info[0])+'/'
print('py_fdn:', py_fdn)
print('Running tests...')
subprocess.check_call('python '+py_fdn+'test/test.py', shell=True)
print('done!')

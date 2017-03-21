#!/usr/bin/env python
import sys
import subprocess

py_fdn = 'python'+str(sys.version_info[0])+'/'
subprocess.call('python '+py_fdn+'test/test.py')

#!/usr/bin/env python
import sys
import subprocess

#subprocess.call('ls', shell=True)
py_fdn = 'python'+str(sys.version_info[0])+'/'
print('py_fdn:', py_fdn)
subprocess.call('python '+py_fdn+'test/test.py', shell=True)

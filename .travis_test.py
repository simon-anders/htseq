#!/usr/bin/env python
import sys
import subprocess

subprocess.call('ls', shell=True)
py_fdn = 'python'+str(sys.version_info[0])+'/'
subprocess.call('python '+py_fdn+'test/test.py', shell=True)

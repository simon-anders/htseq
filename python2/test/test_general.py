from __future__ import print_function
import sys
import os
import glob
import distutils.util
import doctest

build_dir = "build/lib.%s-%s" % (distutils.util.get_platform(), sys.version[0:3])

sys.path.insert(0, os.path.join(os.getcwd(), build_dir))
import HTSeq

py_fdn = 'python'+str(sys.version_info[0])


def test_fasta_parser():
    print("Test Fasta parser")
    for seq in HTSeq.FastaReader('example_data/fastaExLong.fa'):
        pass
    print("Test passed")


def test_pickle():
    import pickle

    print('Test pickling and inpickling')
    fn = 'example_data/pickle_test.pickle'
    pickles = [
            {'name': 'HTSeq.Sequence',
             'object': HTSeq.Sequence(b'ACTG', 'sequence')},
            ]

    for pic in pickles:
        with open(fn, 'w') as f:
            print('Pickling '+pic['name'])
            pickle.dump(pic['object'], f)
            print('Done')

        with open(fn, 'r') as f:
            print('Unpickling '+pic['name'])
            pickle.load(f)
            print('Done')

    os.remove(fn)


if len(sys.argv) == 1:
    test_fasta_parser()
    test_pickle()
else:
    print("Wrong usage")
    print("Call without arguments to run all tests")

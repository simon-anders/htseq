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
    pickles = [
            {'name': 'HTSeq.Sequence',
             'object': HTSeq.Sequence(b'ACTG', 'sequence'),
             'assert_properties': ('seq', 'name', 'descr')},
            ]

    for pic in pickles:
        print('Pickling '+pic['name'])
        pickled = pickle.dumps(pic['object'])
        print('Done')

        print('Unpickling '+pic['name'])
        unpick = pickle.loads(pickled)
        print('Done')

        if 'assert_properties' in pic:
            print('Checking serialized/deserialized')
            for prop in pic['assert_properties']:
                assert getattr(pic['object'], prop) == getattr(unpick, prop)
            print('Done')


if len(sys.argv) == 1:
    test_fasta_parser()
    test_pickle()
else:
    print("Wrong usage")
    print("Call without arguments to run all tests")

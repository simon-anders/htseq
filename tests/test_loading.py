#!/usr/bin/env python

import unittest


class TestFusion(unittest.TestCase):
    def test_01(self):
        import HTSeq
        self.assertTrue(len(str(HTSeq.__version__)) > 0)


def main():
    unittest.main()


if __name__ == '__main__':
    main()

#!/usr/bin/env python

import unittest
import os
import shutil
import sys
import subprocess
import ihm.reader
import RMF

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))


class Tests(unittest.TestCase):

    def test_simple(self):
        """Test model building"""
        os.chdir(os.path.join(TOPDIR, 'modeling_scripts'))
        if os.path.exists('./only_xl'):
            shutil.rmtree('./only_xl')
        p = subprocess.check_call(["python", "modeling_test.py", "--test"])

        try:
            if os.path.exists('./only_xl/run0/output/initial.0.rmf3'):
                print('rmf3 is present')
        except FileNotFoundError:
                print('rmf3 is not present')

        rmf_fn_1 = os.path.join(TOPDIR, 'modeling_scripts' , 'only_xl', 'run0', 'output', 'rmfs', '0.rmf3')
        assert os.path.isfile(rmf_fn_1)
        rh1 = RMF.open_rmf_file_read_only(rmf_fn_1)
        print(rh1.get_number_of_frames())
        assert rh1.get_number_of_frames() == 20

        if os.path.exists(os.path.join(TOPDIR, 'modeling_scripts' , 'only_xl')):
            shutil.rmtree(os.path.join(TOPDIR, 'modeling_scripts' , 'only_xl'))


if __name__ == '__main__':
    unittest.main()

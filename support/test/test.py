#!/usr/bin/env python

import unittest
import os
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


class Tests(unittest.TestCase):
    def test_modeling_script(self):
        """Test the main modeling script"""
        subprocess.check_call(["python", "modeling.py", "--test"],
                              cwd=os.path.join(TOPDIR, 'modeling'))
        # todo: assert that it generated correct outputs


if __name__ == '__main__':
    unittest.main()

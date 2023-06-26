#!/usr/bin/env python

import unittest
import os
import subprocess
import IMP
import pickle


TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


class Tests(unittest.TestCase):
    def test_modeling_script(self):
        """Test the main modeling script"""
        subprocess.check_call(["python", "modeling.py", "--test"],
                              cwd=os.path.join(TOPDIR, 'modeling'))
        # todo: assert that it generated correct outputs

    def test_pickle(self):
        """Test that pickled ReplicaExchange object works"""
        # Set up modeling but don't run sampling
        os.chdir(os.path.join(TOPDIR, 'modeling'))
        with open('modeling.py') as fh:
            contents = fh.read().replace('rex.execute_macro()', '')
        g = {}
        exec(contents, g)
        rex = g['rex']
        del g
        rex.vars['number_of_frames'] = 2

        dump = pickle.dumps((rex.model, rex))

        # Run the original ReplicaExchange and get the final score
        IMP.random_number_generator.seed(99)
        rex.execute_macro()
        rs = IMP.pmi.tools.get_restraint_set(rex.model)
        original_score = rs.evaluate(False)
        del rex, rs

        # With the same random seed, we should get the exact same trajectory
        # with the pickled object
        newm, newrex = pickle.loads(dump)
        IMP.random_number_generator.seed(99)
        newrex.execute_macro()
        rs = IMP.pmi.tools.get_restraint_set(newrex.model)
        new_score = rs.evaluate(False)
        self.assertAlmostEqual(original_score, new_score, delta=1e-4)


if __name__ == '__main__':
    unittest.main()

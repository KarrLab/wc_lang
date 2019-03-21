""" Test SBML export transform

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2019-03-18
:Copyright: 2019, Karr Lab
:License: MIT
"""
from wc_lang.core import Model
from wc_lang.transform.prep_for_sbml import PrepForSbmlTransform
import os
import unittest
import wc_lang.io


class PrepForSbmlTransformTestCase(unittest.TestCase):
    def test_run(self):
        filename = os.path.join(os.path.dirname(__file__), '..', 'fixtures', 'sbml-io.xlsx')
        model = wc_lang.io.Reader().run(filename)[Model][0]
        transform = PrepForSbmlTransform()
        transform.run(model)

        filename = os.path.join(os.path.dirname(__file__), '..', 'fixtures', 'sbml-io-transformed.xlsx')
        sbml_model = wc_lang.io.Reader().run(filename)[Model][0]

        self.assertTrue(model.is_equal(sbml_model))

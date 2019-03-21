""" Tests API

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-03-12
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_lang
import wc_lang.io
import wc_lang.config
import wc_lang.sbml
import wc_lang.sbml.io
import types
import unittest


class ApiTestCase(unittest.TestCase):
    def test(self):
        self.assertIsInstance(wc_lang, types.ModuleType)
        self.assertIsInstance(wc_lang.Model, type)
        self.assertIsInstance(wc_lang.io, types.ModuleType)
        self.assertIsInstance(wc_lang.config, types.ModuleType)
        self.assertIsInstance(wc_lang.config.get_config, types.FunctionType)
        self.assertIsInstance(wc_lang.sbml, types.ModuleType)
        self.assertIsInstance(wc_lang.sbml.io, types.ModuleType)

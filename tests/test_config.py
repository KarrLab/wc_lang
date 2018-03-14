""" Tests of config module

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-03-14
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_lang import config
import unittest


class ConfigTestCase(unittest.TestCase):
    def test_get_config(self):
        config.get_config()

    def test_get_debug_logs_config(self):
        config.get_debug_logs_config()

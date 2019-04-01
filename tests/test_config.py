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

    def test_validate_config(self):
        config.get_config(extra={'wc_lang': {
            'dfba': {
                'flux_bounds': {
                    'min_reversible': -1e3,
                    'min_irreversible': -1e3,
                    'max': 1e3,
                },
            },
        },
        })

        with self.assertRaisesRegex(ValueError, 'dFBA reversible flux bound must be greater'):
            config.get_config(extra={'wc_lang': {
                'dfba': {
                    'flux_bounds': {
                        'min_reversible': 1e4,
                        'min_irreversible': -1e3,
                        'max': 1e3,
                    },
                },
            },
            })
        with self.assertRaisesRegex(ValueError, 'dFBA irreversible flux bound must be greater'):
            config.get_config(extra={'wc_lang': {
                'dfba': {
                    'flux_bounds': {
                        'min_reversible': -1e3,
                        'min_irreversible': 1e4,
                        'max': 1e3,
                    },
                },
            },
            })

    def test_get_debug_logs_config(self):
        config.get_debug_logs_config()

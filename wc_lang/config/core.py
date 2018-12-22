""" Configuration

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-10-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import configobj
import os
import pkg_resources
import wc_utils.config
import wc_utils.debug_logs.config


def get_config(extra=None):
    """ Get configuration

    Args:
        extra (:obj:`dict`, optional): additional configuration to override

    Returns:
        :obj:`configobj.ConfigObj`: nested dictionary with the configuration settings loaded from the configuration source(s).
    """
    paths = wc_utils.config.ConfigPaths(
        default=pkg_resources.resource_filename('wc_lang', 'config/core.default.cfg'),
        schema=pkg_resources.resource_filename('wc_lang', 'config/core.schema.cfg'),
        user=(
            'wc_lang.cfg',
            os.path.expanduser('~/.wc/wc_lang.cfg'),
        ),
    )

    config = wc_utils.config.ConfigManager(paths).get_config(extra=extra)
    validate_config(config)
    return config


def validate_config(config):
    """ Validate configuration

    * Check that flux_min_bound >= 0 and flux_max_bound >= flux_min_bound

    Args:
        config (:obj:`configobj.ConfigObj`): nested dictionary with the configuration settings

    Raises:
        :obj:`ValueError`: if minimum dFBA flux bound is negative or the
            maximum dFBA flux bound is less than the minimum dFBA flux bound
    """
    flux_min_bound_reversible = config['wc_lang']['dfba']['flux_min_bound_reversible']
    flux_min_bound_irreversible = config['wc_lang']['dfba']['flux_min_bound_irreversible']
    flux_max_bound = config['wc_lang']['dfba']['flux_max_bound']

    if flux_max_bound < flux_min_bound_reversible:
        raise ValueError(("minimum dFBA reversible flux bound must be greater than or equal to "
                          "the maximum bound:\n"
                          "  flux_min_bound_reversible={}\n"
                          "  flux_max_bound={}").format(
            flux_min_bound_reversible,
            flux_max_bound))

    if flux_max_bound < flux_min_bound_irreversible:
        raise ValueError(("minimum dFBA irreversible flux bound must be greater than or equal to "
                          "the maximum bound:\n"
                          "  flux_min_bound_irreversible={}\n"
                          "  flux_max_bound={}").format(
            flux_min_bound_irreversible,
            flux_max_bound))


def get_debug_logs_config(extra=None):
    """ Get debug logs configuration

    Args:
        extra (:obj:`dict`, optional): additional configuration to override

    Returns:
        :obj:`configobj.ConfigObj`: nested dictionary with the configuration settings loaded from the configuration source(s).
    """
    paths = wc_utils.debug_logs.config.paths.deepcopy()
    paths.default = pkg_resources.resource_filename('wc_lang', 'config/debug.default.cfg')
    paths.user = (
        'wc_lang.debug.cfg',
        os.path.expanduser('~/.wc/wc_lang.debug.cfg'),
    )
    return wc_utils.config.ConfigManager(paths).get_config(extra=extra)

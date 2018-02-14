""" Configuration

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-10-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

import configobj
import os
import pkg_resources
import wc_utils.config.core
import wc_utils.debug_logs.config


def get_config(extra=None):
    """ Get configuration

    Args:
        extra (:obj:`dict`, optional): additional configuration to override

    Returns:
        :obj:`configobj.ConfigObj`: nested dictionary with the configuration settings loaded from the configuration source(s).
    """
    paths = wc_utils.config.core.ConfigPaths(
        default=pkg_resources.resource_filename('wc_lang', 'config/core.default.cfg'),
        schema=pkg_resources.resource_filename('wc_lang', 'config/core.schema.cfg'),
        user=(
            'wc_lang.core.cfg',
            os.path.expanduser('~/.wc/wc_lang.core.cfg'),
        ),
    )

    return wc_utils.config.core.ConfigManager(paths).get_config(extra=extra)


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
    return wc_utils.config.core.ConfigManager(paths).get_config(extra=extra)

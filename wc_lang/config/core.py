""" Configuration

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2017-10-22
:Copyright: 2017, Karr Lab
:License: MIT
"""

from pkg_resources import resource_filename
from wc_utils.config.core import ConfigPaths
from wc_utils.debug_logs.config import paths as debug_logs_default_paths
import os


core = ConfigPaths(
    default=resource_filename('wc_lang', 'config/core.default.cfg'),
    schema=resource_filename('wc_lang', 'config/core.schema.cfg'),
    user=(
        'wc_lang.core.cfg',
        os.path.expanduser('~/.wc/wc_lang.core.cfg'),
    ),
)

debug_logs = debug_logs_default_paths.deepcopy()
debug_logs.default = resource_filename('wc_lang', 'config/debug.default.cfg')
debug_logs.user = (
    'wc_lang.debug.cfg',
    os.path.expanduser('~/.wc/wc_lang.debug.cfg'),
)

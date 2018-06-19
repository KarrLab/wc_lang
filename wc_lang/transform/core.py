""" Transform models.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from abc import ABCMeta, abstractmethod
from six import with_metaclass
from wc_lang import Model
import sys


def get_transforms():
    """ Get dictionary of available transform classes

    Returns:
        :obj:`dict` of `str`: `class`: dictionary of available transform classes
    """
    module = sys.modules[__name__.rpartition('.')[0]]
    transforms = {}
    for attr_name in dir(module):
        attr = getattr(module, attr_name)
        if isinstance(attr, type) and issubclass(attr, Transform) and attr is not Transform:
            transforms[attr.Meta.id] = attr

    return transforms


class Transform(with_metaclass(ABCMeta, object)):

    @abstractmethod
    def run(self, model):
        """ Transform a model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: transformed model
        """
        pass  # pragma: no cover

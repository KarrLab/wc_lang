""" Transform models.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from abc import ABC, abstractmethod
from wc_lang.core import Model
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


class Transform(ABC):

    @abstractmethod
    def run(self, model):
        """ Transform a model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: transformed model
        """
        pass  # pragma: no cover


class CompositeTransform(Transform):
    """ Transform which is composed of multiple transforms 

    Attributes:
        transforms (:obj:`list` of :obj:`Transform`): transforms
        DEFAULT_TRANSFORMS (:obj:`tuple` of :obj:`Transform`): default transforms
    """
    DEFAULT_TRANSFORMS = ()

    def __init__(self, transforms=None):
        """
        Args:
            transforms (:obj:`list` of :obj:`Transform`, optional): list of transforms
        """
        self.transforms = transforms or list(self.DEFAULT_TRANSFORMS)

    def run(self, model):
        """ Transform model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: same model, but transformed
        """
        for transform in self.transforms:
            transform.run(model)

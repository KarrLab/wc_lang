""" Generate a reduced model by nullifying specific attributes

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-03-20
:Copyright: 2019, Karr Lab
:License: MIT
"""

from .core import Transform


class NullifyAttrsTransform(Transform):
    """ Generate a reduced model by nullifying specific attributes

    Attributes:
        attrs_to_keep (:obj:`str`): name of meta attribute which encodes lists of the
            attributes to keep (not nullify)
    """

    class Meta(object):
        id = 'NullifyAttrs'
        label = 'Generate a reduced model by nullifying specific attributes'

    def __init__(self, attrs_to_keep):
        self.attrs_to_keep = attrs_to_keep

    def run(self, model):
        """ Transform model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: same model, but transformed
        """
        for obj in model.get_related():
            Meta = obj.__class__.Meta
            attrs_to_keep = Meta.child_attrs[self.attrs_to_keep]

            assert not set(attrs_to_keep).difference(set(Meta.attributes.keys()))

            discard_attr_names = set(Meta.attributes.keys()).difference(set(attrs_to_keep))
            for attr_name in discard_attr_names:
                attr = Meta.attributes[attr_name]
                none_value = attr.get_none_value()
                setattr(obj, attr_name, none_value)

        return model

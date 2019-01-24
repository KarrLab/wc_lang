""" Change a value of an attribute of a model

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .core import Transform
from wc_lang.core import Model
import json


class ChangeValueTransform(Transform):
    """ Change a value of an attribute of a model

    Attributes:
        attr_path (:obj:`list` of :obj:`list` of :obj:`str`): list that
            represents the path to an attribute or nested attribute of a model

            Examples:

            * model.name --> 'name'
            * model.reactions.get_one(id='rxn_1').reversible --> (('reactions', {'id': 'rxn_1'}), 'reversible')
            * model.parameters.get_one(id='param_1').value --> (('parameters', {'id': 'param_1'}), 'value')
            * model \
                .reactions.get_one(id='rxn_1') \
                .rate_laws.get_one(direction=RateLawDirection.forward) \
                .expression \
                .parameters.get_one(id='param_1') \
                .value
              --> (
                    ('reactions', {'id': 'rxn_1'},
                    ('rate_laws', {'direction': RateLawDirection.forward}),
                    'expression',
                    ('parameters', {'id': 'param_1'}),
                    'value',
                  )
        value (:obj:`object`): new value
    """

    class Meta(object):
        id = 'ChangeValue'
        label = 'Change a value of an attribute of a model'

    def __init__(self, attr_path=None, value=None):
        """
        Args:
            attr_path (:obj:`list` of :obj:`list` of :obj:`str`, optional): list that
                represents the path to an attribute or nested attribute of a model
            value (:obj:`object`, optional): new value
        """
        self.attr_path = attr_path
        self.value = value

    def run(self, model):
        """ Change a value of an attribute of a model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: same model, but with a different value of an attribute
        """
        model.set_nested_attr_val(self.attr_path, self.value)
        return model

    def attr_path_to_str(self):
        """ Generate a string representation of `attr_path`

        Returns:
            :obj:`str`: string representation of `attr_path`R
        """
        return json.dumps(self.attr_path)

    @staticmethod
    def attr_path_from_str(str):
        """ Generate `attr_path` from its string representation
        Args:
            str (:obj:`str`): string representation of `attr_path`

        Returns:
            :obj:`Object`: `attr_path`
        """
        return json.loads(str)

    def __eq__(self, other, tol=0.):
        """ Compare two :obj:`ChangeValueTransform` objects

        Args:
            other (:obj:`Object`): other object
            tol (:obj:`float`, optional): equality tolerance

        Returns:
            :obj:`bool`: true if :obj:`ChangeValueTransform` objects are semantically equal
        """
        if other.__class__ is not self.__class__:
            return False

        if self.attr_path != other.attr_path:
            return False

        attr = Model.get_nested_attr(self.attr_path)
        if not attr.value_equal(self.value, other.value, tol=tol):
            return False

        return True

""" Change a value of an attribute of a model

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .core import Transform
from wc_lang.core import Model


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

    def __init__(self, attr_path, value):
        """
        Args:
            attr_path (:obj:`list` of :obj:`list` of :obj:`str`): list that 
                represents the path to an attribute or nested attribute of a model
            value (:obj:`object`): new value
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
        model.set_nested_attr(self.attr_path, self.value)
        return model

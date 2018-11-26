""" Change a value of an attribute of a model

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from .core import Transform
from wc_lang import Compartment, Function, Parameter, RateLawDirection, Reaction, Species
import copy
import itertools


class ChangeValueTransform(Transform):
    """ Change a value of an attribute of a model

    Attributes:
        target_type (:obj:`type`): type of the target to change, E.g.
            :obj:`wc_lang.Compartment`, :obj:`wc_lang.Function`
        target_id (:obj:`str`): id of the target to change
        target_attr (:obj:`list`): list of names of the nested attribute to change
        value (:obj:`object`): new value

    Supports the following combinations of :obj:`target_type`, :obj:`target_id` and :obj:`target_attr`

        +-------------+-----------------------------------+-----------------------------------------------------------+
        | target_type | target_id                         | target_attr                                               |
        +=============+===================================+===========================================================+
        | Compartment | Compartment.id                    | ['initial_volume']                                        |
        +-------------+-----------------------------------+-----------------------------------------------------------+
        | Function    | Function.id                       | ['expression', 'expression']                              |
        +-------------+-----------------------------------+-----------------------------------------------------------+
        | Parameter   | Parameter.id                      | ['value']                                                 |
        | Parameter   | Parameter.id                      | ['units']                                                 |
        +-------------+-----------------------------------+-----------------------------------------------------------+
        | Reaction    | Reaction.id                       | ['reversible']                                            |
        | Reaction    | Reaction.id                       | ['min_flux']                                              |
        | Reaction    | Reaction.id                       | ['max_flux']                                              |
        | Reaction    | Reaction.id                       | ['rate_laws', RateLawDirection, 'expression', 'expression'] |
        +-------------+-----------------------------------+-----------------------------------------------------------+
        | Species     | SpeciesType.id [ Compartment.id ] | ['concentration', 'value']                                |
        | Species     | SpeciesType.id [ Compartment.id ] | ['concentration', 'units']                                |
        +-------------+-----------------------------------+-----------------------------------------------------------+
    """

    class Meta(object):
        id = 'ChangeValue'
        label = 'Change a value of an attribute of a model'

    def __init__(self, target_type, target_id, target_attr, value):
        """
        Args:
            target_type (:obj:`type`): type of the target to change, E.g.
                :obj:`wc_lang.Compartment`, :obj:`wc_lang.Function`
            target_id (:obj:`str`): id of the target to change
            target_attr (:obj:`list`): list of names of the nested attribute to change
            value (:obj:`object`): new value
        """
        self.target_type = target_type
        self.target_id = target_id
        self.target_attr = target_attr
        self.value = value

    def run(self, model):
        """ Change a value of an attribute of a model

        Args:
            model (:obj:`Model`): model definition

        Returns:
            :obj:`Model`: same model, but with a different value of an attribute
        """
        # get object with id :obj:`target_id` of type :obj:`target_type`

        if self.target_type in [Reaction]:
            for submodel in model.submodels:
                target_obj = submodel.reactions.get_one(id=self.target_id)
                if target_obj:
                    break
        elif self.target_type in [Species]:
            species_type_id, _, compartment_id = self.target_id[0:-1].partition('[')
            species_type = model.species_types.get_one(id=species_type_id)
            compartment = model.compartments.get_one(id=compartment_id)
            target_obj = species_type.species.get_one(compartment=compartment)
        else:
            target_objs = getattr(model, self.target_type.Meta.attributes['model'].related_name)
            target_obj = target_objs.get_one(id=self.target_id)

        # get attribute with path :obj:`target_attr`
        if self.target_type == Reaction and self.target_attr[0] == 'rate_laws':
            target_attr = target_obj.rate_laws.get_one(direction=RateLawDirection[self.target_attr[1]])
            attr_names = self.target_attr[2:-1]
        else:
            target_attr = target_obj
            attr_names = self.target_attr[0:-1]

        for attr_name in attr_names:
            target_attr = getattr(target_attr, attr_name)

        # change value
        setattr(target_attr, self.target_attr[-1], self.value)

        # return model
        return model

""" Data model to represent biochemical models.

This module defines classes that represent the schema of a biochemical model:

* :obj:`Taxon`
* :obj:`Model`
* :obj:`Submodel`
* :obj:`DfbaObjective`
* :obj:`Compartment`
* :obj:`SpeciesType`
* :obj:`Species`
* :obj:`Concentration`
* :obj:`Reaction`
* :obj:`SpeciesCoefficient`
* :obj:`RateLaw`
* :obj:`RateLawExpression`
* :obj:`BiomassComponent`
* :obj:`BiomassReaction`
* :obj:`Parameter`
* :obj:`Reference`
* :obj:`DatabaseReference`

These are all instances of `obj_model.Model`, an alias for `obj_model.Model`.
A biochemical model may contain a list of instances of each of these classes, interlinked
by object references. For example, a :obj:`Reaction` will reference its constituent
:obj:`SpeciesCoefficient` instances, and the :obj:`RateLaw` that describes the reaction's rate.

This module also defines numerous classes that serve as attributes of these classes.

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016-2017, Karr Lab
:License: MIT
"""
# TODO: for determinism, replace remaining list(set(list1)) expressions with det_dedupe(list1) in other packages

from enum import Enum, EnumMeta
from itertools import chain
from math import ceil, floor, exp, log, log10, isnan
from natsort import natsorted, ns
from six import with_metaclass, string_types
import collections
import datetime
import networkx
import pkg_resources
import re
import six
import stringcase
import sys
import token

from obj_model import (BooleanAttribute, EnumAttribute,
                       FloatAttribute, PositiveFloatAttribute,
                       IntegerAttribute, PositiveIntegerAttribute,
                       RegexAttribute, SlugAttribute, StringAttribute, LongStringAttribute, UrlAttribute,
                       DateTimeAttribute,
                       OneToOneAttribute, ManyToOneAttribute, ManyToManyAttribute, OneToManyAttribute,
                       InvalidModel, InvalidObject, InvalidAttribute, TabularOrientation)
import obj_model
from wc_utils.util.enumerate import CaseInsensitiveEnum, CaseInsensitiveEnumMeta
from wc_utils.util.list import det_dedupe
from wc_lang.sbml.util import (wrap_libsbml, str_to_xmlstr, LibSBMLError,
                               init_sbml_model, create_sbml_parameter, add_sbml_unit)
from wc_lang.expression_utils import (ParsedExpression, ParsedExpressionError,
                                      LinearExpressionVerifier, WcLangToken, TokCodes)

with open(pkg_resources.resource_filename('wc_lang', 'VERSION'), 'r') as file:
    wc_lang_version = file.read().strip()

# wc_lang generates obj_model SchemaWarning warnings because some Models lack primary attributes.
# These models include :obj:`RateLaw`, :obj:`SpeciesCoefficient`, :obj:`RateLawExpression`, and :obj:`Species`.
# However, these are not needed by the workbook and delimiter-separated representations of
# models on disk. Therefore, suppress the warnings.
import warnings
warnings.filterwarnings('ignore', '', obj_model.SchemaWarning, 'obj_model')

# configuration
import wc_lang.config.core
config_wc_lang = wc_lang.config.core.get_config()['wc_lang']

EXTRACELLULAR_COMPARTMENT_ID = config_wc_lang['EXTRACELLULAR_COMPARTMENT_ID']


class TaxonRankMeta(CaseInsensitiveEnumMeta):

    def __getattr__(cls, name):
        """ Get value by name

        Args:
            name (:obj:`str`): attribute name

        Returns:
            :obj:`TaxonRank`: taxonomic rank
        """
        if name.lower() == 'class':
            name = 'classis'
        return super(TaxonRankMeta, cls).__getattr__(name)

    def __getitem__(cls, name):
        """ Get value by name

        Args:
            name (:obj:`str`): attribute name

        Returns:
            :obj:`TaxonRank`: taxonomic rank
        """
        lower_name = name.lower()

        if lower_name in ['varietas', 'strain']:
            name = 'variety'
        elif lower_name == 'tribus':
            name = 'tribe'
        elif lower_name == 'familia':
            name = 'family'
        elif lower_name == 'ordo':
            name = 'order'
        elif lower_name == 'class':
            name = 'classis'
        elif lower_name in ['division', 'divisio']:
            name = 'phylum'
        elif lower_name == 'regnum':
            name = 'kingdom'
        return super(TaxonRankMeta, cls).__getitem__(name)


class TaxonRank(with_metaclass(TaxonRankMeta, int, Enum)):
    """ Taxonomic ranks """
    domain = 1
    kingdom = 2
    phylum = 3
    classis = 4
    order = 5
    family = 6
    tribe = 7
    genus = 8
    species = 9
    variety = 10


class CompartmentType(int, CaseInsensitiveEnum):
    """ Compartment types """
    abstract = 1
    # physical_1d = 2
    # physical_2d = 2
    physical_3d = 3


class SubmodelAlgorithm(int, CaseInsensitiveEnum):
    """ Submodel algorithms """
    dfba = 1
    ode = 2
    ssa = 3


class SpeciesTypeType(int, CaseInsensitiveEnum):
    """ Types of species types """
    metabolite = 1
    protein = 2
    dna = 3
    rna = 4
    pseudo_species = 5


ConcentrationUnit = Enum('ConcentrationUnit', type=int, names=[
    ('molecules', 1),
    ('M', 2),
    ('mM', 3),
    ('uM', 4),
    ('nM', 5),
    ('pM', 6),
    ('fM', 7),
    ('aM', 8),
    ('mol dm^-2', 9),
])
ConcentrationUnit.Meta = {
    ConcentrationUnit['molecules']: {
        'xml_id': 'molecules',
        'substance_units': {'kind': 'item', 'exponent': 1, 'scale': 0},
        'volume_units': None,
    },
    ConcentrationUnit['M']: {
        'xml_id': 'M',
        'substance_units': {'kind': 'mole', 'exponent': 1, 'scale': 0},
        'volume_units': {'kind': 'litre', 'exponent': -1, 'scale': 0},
    },
    ConcentrationUnit['mM']: {
        'xml_id': 'mM',
        'substance_units': {'kind': 'mole', 'exponent': 1, 'scale': -3},
        'volume_units': {'kind': 'litre', 'exponent': -1, 'scale': 0},
    },
    ConcentrationUnit['uM']: {
        'xml_id': 'uM',
        'substance_units': {'kind': 'mole', 'exponent': 1, 'scale': -6},
        'volume_units': {'kind': 'litre', 'exponent': -1, 'scale': 0},
    },
    ConcentrationUnit['nM']: {
        'xml_id': 'nM',
        'substance_units': {'kind': 'mole', 'exponent': 1, 'scale': -9},
        'volume_units': {'kind': 'litre', 'exponent': -1, 'scale': 0},
    },
    ConcentrationUnit['pM']: {
        'xml_id': 'pM',
        'substance_units': {'kind': 'mole', 'exponent': 1, 'scale': -12},
        'volume_units': {'kind': 'litre', 'exponent': -1, 'scale': 0},
    },
    ConcentrationUnit['fM']: {
        'xml_id': 'fM',
        'substance_units': {'kind': 'mole', 'exponent': 1, 'scale': -15},
        'volume_units': {'kind': 'litre', 'exponent': -1, 'scale': 0},
    },
    ConcentrationUnit['aM']: {
        'xml_id': 'aM',
        'substance_units': {'kind': 'mole', 'exponent': 1, 'scale': -18},
        'volume_units': {'kind': 'litre', 'exponent': -1, 'scale': 0},
    },
    ConcentrationUnit['mol dm^-2']: {
        'xml_id': 'mol_per_dm_2',
        'substance_units': {'kind': 'mole', 'exponent': 1, 'scale': 0},
        'volume_units': {'kind': 'metre', 'exponent': -2, 'scale': -1},
    },
}


class RateLawDirection(int, CaseInsensitiveEnum):
    """ Rate law directions """
    backward = -1
    forward = 1


class RateLawType(int, CaseInsensitiveEnum):
    """ SBO rate law types """
    hill = 192
    mass_action = 12
    michaelis_menten = 29
    modular = 527
    other = 1


RateLawType = CaseInsensitiveEnum('RateLawType', type=int, names=[
    ('hill', 192),
    ('mass-action', 12),
    ('michaelis-menten', 29),
    ('modular', 527),
    ('other', 1),
])

RateLawUnits = Enum('RateLawUnits', type=int, names=[
    ('s^-1', 1),
    ('M s^-1', 2),
])


class ParameterType(int, Enum):
    """ SBO parameter types """
    k_cat = 25
    v_max = 186
    K_m = 27
    K_i = 261
    other = 2


class ReferenceType(int, CaseInsensitiveEnum):
    """ Reference types """
    article = 1
    book = 2
    online = 3
    proceedings = 4
    thesis = 5

    inbook = 6
    incollection = 7
    inproceedings = 8

    misc = 9


class ReactionParticipantAttribute(ManyToManyAttribute):
    """ Reaction participants """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(ReactionParticipantAttribute, self).__init__('SpeciesCoefficient', related_name=related_name,
                                                           min_related=1,
                                                           verbose_name=verbose_name,
                                                           verbose_related_name=verbose_related_name,
                                                           help=help)

    def serialize(self, participants, encoded=None):
        """ Serialize related object

        Args:
            participants (:obj:`list` of :obj:`SpeciesCoefficient`): Python representation of reaction participants
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: simple Python representation
        """
        if not participants:
            return ''

        comps = set([part.species.compartment for part in participants])
        if len(comps) == 1:
            global_comp = comps.pop()
        else:
            global_comp = None

        if global_comp:
            participants = natsorted(participants, lambda part: part.species.species_type.id, alg=ns.IGNORECASE)
        else:
            participants = natsorted(participants, lambda part: (
                part.species.species_type.id, part.species.compartment.id), alg=ns.IGNORECASE)

        lhs = []
        rhs = []
        for part in participants:
            if part.coefficient < 0:
                lhs.append(part.serialize(show_compartment=global_comp is None, show_coefficient_sign=False))
            elif part.coefficient > 0:
                rhs.append(part.serialize(show_compartment=global_comp is None, show_coefficient_sign=False))

        if global_comp:
            return '[{}]: {} ==> {}'.format(global_comp.get_primary_attribute(), ' + '.join(lhs), ' + '.join(rhs))
        else:
            return '{} ==> {}'.format(' + '.join(lhs), ' + '.join(rhs))

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of `list` of `SpeciesCoefficient`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        errors = []

        id = r'[a-z][a-z0-9_]*'
        stoch = r'\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\)'
        gbl_part = r'({} *)*({})'.format(stoch, id)
        lcl_part = r'({} *)*({}\[{}\])'.format(stoch, id, id)
        gbl_side = r'{}( *\+ *{})*'.format(gbl_part, gbl_part)
        lcl_side = r'{}( *\+ *{})*'.format(lcl_part, lcl_part)
        gbl_pattern = r'^\[({})\]: *({}|) *==> *({}|)$'.format(id, gbl_side, gbl_side)
        lcl_pattern = r'^({}|) *==> *({}|)$'.format(lcl_side, lcl_side)

        value = value.strip(' ')
        global_match = re.match(gbl_pattern, value, flags=re.I)
        local_match = re.match(lcl_pattern, value, flags=re.I)

        if global_match:
            if global_match.group(1) in objects[Compartment]:
                global_comp = objects[Compartment][global_match.group(1)]
            else:
                global_comp = None
                errors.append('Undefined compartment "{}"'.format(global_match.group(1)))
            lhs = global_match.group(2)
            rhs = global_match.group(14)

        elif local_match:
            global_comp = None
            lhs = local_match.group(1)
            rhs = local_match.group(13)

        else:
            return (None, InvalidAttribute(self, ['Incorrectly formatted participants: {}'.format(value)]))

        lhs_parts, lhs_errors = self.deserialize_side(-1., lhs, objects, global_comp)
        rhs_parts, rhs_errors = self.deserialize_side(1., rhs, objects, global_comp)

        parts = lhs_parts + rhs_parts
        errors.extend(lhs_errors)
        errors.extend(rhs_errors)

        if errors:
            return (None, InvalidAttribute(self, errors))
        return (parts, None)

    def deserialize_side(self, direction, value, objects, global_comp):
        """ Deserialize the LHS or RHS of a reaction expression

        Args:
            direction (:obj:`float`): -1. indicates LHS, +1. indicates RHS
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            global_comp (:obj:`Compartment`): global compartment of the reaction

        Returns:
            :obj:`tuple`:

                * :obj:`list` of :obj:`SpeciesCoefficient`: list of species coefficients
                * :obj:`list` of :obj:`Exception`: list of errors
        """
        parts_str = re.findall(r'(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*)(\[([a-z][a-z0-9_]*)\])*', value, flags=re.I)

        if global_comp:
            temp = [part[4] for part in parts_str]
        else:
            temp = [part[4] + '[' + part[6] + ']' for part in parts_str]
        repeated_parts = [item for item, count in collections.Counter(temp).items() if count > 1]
        if repeated_parts:
            return ([], ['Participants are repeated\n  {}'.format('\n  '.join(repeated_parts))])

        parts = []
        errors = []
        for part in parts_str:
            part_errors = []

            if part[4] in objects[SpeciesType]:
                species_type = objects[SpeciesType][part[4]]
            else:
                part_errors.append('Undefined species type "{}"'.format(part[4]))

            if global_comp:
                compartment = global_comp
            elif part[6] in objects[Compartment]:
                compartment = objects[Compartment][part[6]]
            else:
                part_errors.append('Undefined compartment "{}"'.format(part[6]))

            coefficient = direction * float(part[1] or 1.)

            if part_errors:
                errors += part_errors
            else:
                species_id = Species.gen_id(species_type.id, compartment.id)
                species, error = Species.deserialize(species_id, objects)
                if error:
                    raise ValueError('Invalid species "{}"'.format(species_id)
                                     )  # pragma: no cover; unreachable due to above error checking of species types and compartments

                if coefficient != 0:
                    if SpeciesCoefficient not in objects:
                        objects[SpeciesCoefficient] = {}
                    serialized_value = SpeciesCoefficient._serialize(species, coefficient)
                    if serialized_value in objects[SpeciesCoefficient]:
                        rxn_part = objects[SpeciesCoefficient][serialized_value]
                    else:
                        rxn_part = SpeciesCoefficient(species=species, coefficient=coefficient)
                        objects[SpeciesCoefficient][serialized_value] = rxn_part
                    parts.append(rxn_part)

        return (parts, errors)

    def validate(self, obj, value):
        """ Determine if `value` is a valid value of the attribute

        Args:
            obj (:obj:`Reaction`): object being validated
            value (:obj:`list` of :obj:`SpeciesCoefficient`): value of attribute to validate

        Returns:
            :obj:`InvalidAttribute` or None: None if attribute is valid, other return list of errors as an instance of `InvalidAttribute`
        """
        error = super(ReactionParticipantAttribute, self).validate(obj, value)
        if error:
            return error

        # check that LHS and RHS are different
        net_coeffs = {}
        for spec_coeff in value:
            net_coeffs[spec_coeff.species] = \
                net_coeffs.get(spec_coeff.species, 0) + \
                spec_coeff.coefficient
            if net_coeffs[spec_coeff.species] == 0:
                net_coeffs.pop(spec_coeff.species)
        if not net_coeffs:
            return InvalidAttribute(self, ['LHS and RHS must be different'])
        return None


class RateLawExpressionAttribute(ManyToOneAttribute):
    """ Rate law expression """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(RateLawExpressionAttribute, self).__init__('RateLawExpression',
                                                         related_name=related_name, min_related=1, min_related_rev=1,
                                                         verbose_name=verbose_name, verbose_related_name=verbose_related_name, help=help)

    def serialize(self, rate_law_expression, encoded=None):
        """ Serialize related object

        Args:
            rate_law_expression (:obj:`RateLawExpression`): the related `RateLawExpression`
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: simple Python representation of the rate law expression
        """
        return rate_law_expression.serialize()

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of `object`, `InvalidAttribute` or `None`: tuple of cleaned value and cleaning error
        """
        return RateLawExpression.deserialize(value, objects)


class ExpressionAttribute(OneToOneAttribute):
    """ Expression attribute """

    def serialize(self, expression, encoded=None):
        """ Serialize related object

        Args:
            expression (:obj:`obj_model.Model`): the referenced Expression
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: simple Python representation
        """
        return expression.serialize()

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of `object`, `InvalidAttribute` or `None`: tuple of cleaned value and cleaning error
        """
        return self.related_class.deserialize(value, objects)


class Model(obj_model.Model):
    """ Model

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        version (:obj:`str`): version of the model
        url (:obj:`str`): url of the model Git repository
        branch (:obj:`str`): branch of the model Git repository
        revision (:obj:`str`): revision of the model Git repository
        wc_lang_version (:obj:`str`): version of ``wc_lang``
        author (:obj:`str`): author(s)
        author_organization (:obj:`str`): author organization(s)
        author_email (:obj:`str`): author emails(s)
        comments (:obj:`str`): comments
        created (:obj:`datetime`): date created
        updated (:obj:`datetime`): date updated

    Related attributes:
        taxon (:obj:`Taxon`): taxon
        submodels (:obj:`list` of :obj:`Submodel`): submodels
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species types
        species (:obj:`list` of :obj:`Species`): species
        concentrations (:obj:`list` of :obj:`Concentration`): concentrations
        observables (:obj:`list` of :obj:`Observable`): observables
        functions (:obj:`list` of :obj:`Function`): functions
        dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        biomass_reactions (:obj:`list` of :obj:`BiomassReaction`): biomass reactions
        parameters (:obj:`list` of :obj:`Parameter`): parameters
        stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
        references (:obj:`list` of :obj:`Reference`): references
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
    """
    id = SlugAttribute()
    name = StringAttribute()
    version = RegexAttribute(min_length=1, pattern=r'^[0-9]+\.[0-9+]\.[0-9]+', flags=re.I)
    url = obj_model.core.StringAttribute(verbose_name='URL')
    branch = obj_model.core.StringAttribute()
    revision = obj_model.core.StringAttribute()
    wc_lang_version = RegexAttribute(min_length=1, pattern=r'^[0-9]+\.[0-9+]\.[0-9]+', flags=re.I,
                                     default=wc_lang_version, verbose_name='wc_lang version')
    author = LongStringAttribute()
    author_organization = LongStringAttribute()
    author_email = LongStringAttribute()
    comments = LongStringAttribute()
    created = DateTimeAttribute()
    updated = DateTimeAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'version',
                           'url', 'branch', 'revision',
                           'wc_lang_version',
                           'author', 'author_organization', 'author_email',
                           'comments',
                           'created', 'updated')
        tabular_orientation = TabularOrientation.column

    def __init__(self, **kwargs):
        """
        Args:
            **kwargs (:obj:`dict`, optional): dictionary of keyword arguments with keys equal to the names of the model attributes

        Raises:
            :obj:`TypeError`: if keyword argument is not a defined attribute
        """
        super(Model, self).__init__(**kwargs)

        if 'created' not in kwargs:
            self.created = datetime.datetime.now().replace(microsecond=0)
        if 'updated' not in kwargs:
            self.updated = datetime.datetime.now().replace(microsecond=0)

    def validate(self):
        """ Determine if the model is valid

        * Networks of observables and functions are acyclic

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(Model, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        cyclic_deps = self._get_cyclic_deps()
        for model_type, metadata in cyclic_deps.items():
            errors.append(InvalidAttribute(metadata['attribute'], [
                'The following instances of {} cannot have cyclic depencencies:\n  {}'.format(
                    model_type.__name__, '\n  '.join(metadata['models']))]))

        if errors:
            return InvalidObject(self, errors)
        return None

    def _get_cyclic_deps(self):
        """ Verify that the networks of depencencies for observables and functions are acyclic

        Returns:
            :obj:`dict`: dictionary of dictionary of lists of objects with cyclic dependencies,
                keyed by type
        """
        cyclic_deps = {}
        for model_type in (Observable, Function):
            # get name of attribute that contains instances of model_type
            for attr_name, attr in self.__class__.Meta.related_attributes.items():
                if attr.primary_class == model_type:
                    break

            # get all instances of type in model
            models = getattr(self, attr_name)

            # get name of self-referential attribute, if any
            expression_type = model_type.Meta.expression_model
            for self_ref_attr_name, self_ref_attr in expression_type.Meta.attributes.items():
                if isinstance(self_ref_attr, obj_model.RelatedAttribute) and self_ref_attr.related_class == model_type:
                    break

            # find cyclic dependencies
            digraph = networkx.DiGraph()
            for model in models:
                for other_model in getattr(model.expression, self_ref_attr_name):
                    digraph.add_edge(model.id, other_model.id)
            cycles = list(networkx.simple_cycles(digraph))
            if cycles:
                cyclic_deps[model_type] = {
                    'attribute': attr,
                    'attribute_name': attr_name,
                    'models': set(),
                }
            for cycle in cycles:
                cyclic_deps[model_type]['models'].update(set(cycle))

        return cyclic_deps

    def get_submodels(self, __type=None, **kwargs):
        """ Get all submodels

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`Submodel`: submodels
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.submodels.get(__type=__type, **kwargs)

    def get_compartments(self, __type=None, **kwargs):
        """ Get all compartments

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`Compartment`: compartments
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.compartments.get(__type=__type, **kwargs)

    def get_species_types(self, __type=None, **kwargs):
        """ Get all species types

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`SpeciesType`: species types
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.species_types.get(__type=__type, **kwargs)

    def get_species(self, __type=None, **kwargs):
        """ Get all species from submodels

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`Species`: species
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.species.get(__type=__type, **kwargs)

    def get_concentrations(self, __type=None, **kwargs):
        """ Get all concentrations from species types

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`Concentration`: concentations
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.concentations.get(__type=__type, **kwargs)

    def get_observables(self, __type=None, **kwargs):
        """ Get all observables

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`Observables`: observables
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.observables.get(__type=__type, **kwargs)

    def get_functions(self, __type=None, **kwargs):
        """ Get all functions

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`Function`: functions
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.functions.get(__type=__type, **kwargs)

    def get_dfba_objs(self, __type=None, **kwargs):
        """ Get all dFBA objectives

        Returns:
            :obj:`list` of :obj:`DfbaObjective`: dFBA objectives
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.dfba_objs.get(__type=__type, **kwargs)

    def get_reactions(self, __type=None, **kwargs):
        """ Get all reactions from submodels

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`Reaction`: reactions
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.reactions.get(__type=__type, **kwargs)

    def get_rate_laws(self, __type=None, **kwargs):
        """ Get all rate laws from reactions

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`RateLaw`: rate laws
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.rate_laws.get(__type=__type, **kwargs)

    def get_biomass_reactions(self, __type=None, **kwargs):
        """ Get all biomass reactions used by submodels

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`BiomassReaction`: biomass reactions
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.biomass_reactions.get(__type=__type, **kwargs)

    def get_parameters(self, __type=None, **kwargs):
        """ Get all parameters from model and submodels

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`Parameter`: parameters
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.parameters.get(__type=__type, **kwargs)

    def get_stop_conditions(self, __type=None, **kwargs):
        """ Get all stop conditions

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`StopCondition`: stop conditions
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.stop_conditions.get(__type=__type, **kwargs)

    def get_references(self, __type=None, **kwargs):
        """ Get all references from model and children

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`Reference`: references
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.references.get(__type=__type, **kwargs)

    def get_components(self, __type=None, **kwargs):
        """ Find model component of `type` with `id`

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`obj_model.Model`: component with `id`, or `None` if there is no component with `id`=`id`
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        if __type:
            type_names = [stringcase.snakecase(__type.__name__)]
        else:
            type_names = [
                'submodel', 'compartment', 'species_type', 'species',
                'concentration', 'observable', 'function',
                'dfba_obj', 'reaction', 'rate_law', 'biomass_reaction',
                'parameter', 'stop_condition', 'reference',
            ]

        components = []
        for type_name in type_names:
            get_func = getattr(self, 'get_' + type_name + 's')
            components.extend(get_func(__type=__type, **kwargs))

        return components


class Taxon(obj_model.Model):
    """ Biological taxon (e.g. family, genus, species, strain, etc.)

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        rank (:obj:`TaxonRank`): rank
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = OneToOneAttribute(Model, related_name='taxon')
    rank = EnumAttribute(TaxonRank, default=TaxonRank.species)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='taxa')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'rank',
                           'comments', 'references')
        tabular_orientation = TabularOrientation.column


class Submodel(obj_model.Model):
    """ Submodel

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        algorithm (:obj:`SubmodelAlgorithm`): algorithm
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        dfba_obj (:obj:`DfbaObjective`): objective function for a dFBA submodel;
            if not initialized, then `biomass_reaction` is used as the objective function
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        biomass_reactions (:obj:`list` of :obj:`BiomassReaction`): the growth reaction for a dFBA submodel
        parameters (:obj:`list` of :obj:`Parameter`): parameters
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='submodels')
    algorithm = EnumAttribute(SubmodelAlgorithm, default=SubmodelAlgorithm.ssa)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='submodels')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'algorithm', 'comments', 'references')
        indexed_attrs_tuples = (('id',), )

    def validate(self):
        """ Determine if the submodel is valid

        * dFBA submodel has an objective

        .. todo :: Check that the submodel uses water consistently -- either in all compartments or in none

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(Submodel, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        if self.algorithm == SubmodelAlgorithm.dfba:
            if not self.dfba_obj:
                errors.append(InvalidAttribute(self.Meta.related_attributes['dfba_obj'],
                                               ['dFBA submodel must have an objective']))

        if errors:
            return InvalidObject(self, errors)
        return None

    def get_compartments(self):
        """ Get compartments in submodel

        Returns:
            :obj:`list` of :obj:`Compartment`: compartments in submodel
        """
        compartments = []
        for species in self.get_species():
            compartments.append(species.compartment)
        return det_dedupe(compartments)

    def get_specices_types(self):
        """ Get species types in submodel

        Returns:
            :obj:`list` of :obj:`SpeciesType`: species types in submodel
        """
        species_types = []
        for species in self.get_species():
            species_types.append(species.species_type)
        return det_dedupe(species_types)

    def get_species(self):
        """ Get species in submodel

        Returns:
            :obj:`list` of :obj:`Species`: species in submodel
        """
        species = []
        for reaction in self.get_reactions():
            species.extend(reaction.get_species())
        for biomass_reaction in self.get_biomass_reactions():
            for biomass_component in biomass_reaction.biomass_components:
                species.append(biomass_component.species)
        for observable in self.get_observables():
            species.extend(observable.expression.species)
        for function in self.get_functions():
            species.extend(function.expression.species)
        for rate_law in self.get_rate_laws():
            species.extend(rate_law.expression.modifiers)
        return det_dedupe(species)

    def get_observables(self):
        """ Get observables in submodel

        Returns:
            :obj:`list` of :obj:`Observable`: observables in submodel
        """
        obs = []
        for function in self.get_functions():
            obs.extend(function.expression.observables)
        for rate_law in self.get_rate_laws():
            obs.extend(rate_law.expression.observables)
        obs = det_dedupe(obs)
        obs_to_flats = list(obs)
        while obs_to_flats:
            obs_to_flat = obs_to_flats.pop()
            obs.extend(obs_to_flat.expression.observables)
            obs_to_flats.extend(obs_to_flat.expression.observables)
        return det_dedupe(obs)

    def get_functions(self):
        """ Get functions in submodel

        Returns:
            :obj:`list` of :obj:`Function`: functions in submodel
        """
        funcs = []
        for rate_law in self.get_rate_laws():
            funcs.extend(rate_law.expression.functions)
        funcs = det_dedupe(funcs)
        funcs_to_flats = list(funcs)
        while funcs_to_flats:
            funcs_to_flat = funcs_to_flats.pop()
            funcs.extend(funcs_to_flat.expression.functions)
            funcs_to_flats.extend(funcs_to_flat.expression.functions)
        return det_dedupe(funcs)

    def get_dfba_objs(self):
        """ Get dFBA objectives in submodel

        Returns:
            :obj:`list` of :obj:`DfbaObjective`: dFBA objectives in submodel
        """
        if self.dfba_obj:
            return [self.dfba_obj]
        return []

    def get_reactions(self):
        """ Get reactions in submodel

        Returns:
            :obj:`list` of :obj:`Reaction`: reactions in submodel
        """
        reactions = list(self.reactions)
        for dfba_obj in self.get_dfba_objs():
            reactions.extend(dfba_obj.expression.reactions)
        return det_dedupe(reactions)

    def get_rate_laws(self):
        """ Get rate laws in submodel

        Returns:
            :obj:`list` of :obj:`RateLaw`: rate laws in submodel
        """
        rate_laws = []
        for reaction in self.get_reactions():
            rate_laws.extend(reaction.rate_laws)
        return det_dedupe(rate_laws)

    def get_biomass_reactions(self):
        """ Get biomass reactions in submodel

        Returns:
            :obj:`list` of :obj:`BiomassReaction`: biomass
                reactions in submodel
        """
        bm_rxns = list(self.biomass_reactions)
        for dfba_obj in self.get_dfba_objs():
            bm_rxns.extend(dfba_obj.expression.biomass_reactions)
        return det_dedupe(bm_rxns)

    def get_parameters(self):
        """ Get parameters in submodel

        Returns:
            :obj:`list` of :obj:`Parameter`: parameters in submodel
        """
        parameters = []
        for rate_law in self.get_rate_laws():
            parameters.extend(rate_law.expression.parameters)
        for function in self.get_functions():
            parameters.extend(function.expression.parameters)
        return det_dedupe(parameters)

    def get_references(self):
        """ Get references of submodel

        Returns:
            :obj:`list` of :obj:`Reference`: references in submodel
        """
        types = [
            'compartment',
            'species_type',
            'specie',
            'observable',
            'function',
            'dfba_obj',
            'reaction',
            'rate_law',
            'biomass_reaction',
            'parameter',
        ]
        references = []
        for type in types:
            get_func = getattr(self, 'get_' + type + 's')
            for obj in get_func():
                references.extend(obj.references)
        return references

    def get_components(self):
        """ Get components of submodel

        Returns:
            :obj:`list` of :obj:`obj_model.Model`: components in submodel
        """
        return self.get_compartments() + \
            self.get_species_types() + \
            self.get_species() + \
            self.get_observables() + \
            self.get_functions() + \
            self.get_dfba_objs() + \
            self.get_reactions() + \
            self.get_rate_laws() + \
            self.get_biomass_reactions() + \
            self.get_parameters() + \
            self.get_references()

    def add_to_sbml_doc(self, sbml_document):
        """ Add this Submodel to a libsbml SBML document as a `libsbml.model`.

        Args:
             sbml_document (:obj:`obj`): a `libsbml` SBMLDocument

        Returns:
            :obj:`libsbml.model`: the libsbml model

        Raises:
            :obj:`LibSBMLError`: if calling `libsbml` raises an error
        """
        sbml_model = wrap_libsbml(sbml_document.getModel)
        wrap_libsbml(sbml_model.setIdAttribute, self.id)
        if self.name:
            wrap_libsbml(sbml_model.setName, self.name)
        # compartment, dfba_obj, and parameters are created separately
        if self.comments:
            wrap_libsbml(sbml_model.appendNotes, str_to_xmlstr(self.comments))
        return sbml_model


class DfbaObjectiveExpression(obj_model.Model):
    """ A mathematical expression of Reactions and BiomassReactions

    The expression used by a :obj:`DfbaObjective`.

    Attributes:
        expression (:obj:`str`): mathematical expression
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        reactions (:obj:`list` of :obj:`Reaction`): Reactions used by this expression
        biomass_reactions (:obj:`list` of :obj:`Species`): Biomass reactions used by this expression

    Related attributes:
        dfba_obj (:obj:`DfbaObjective`): dFBA objective
    """

    expression = LongStringAttribute(primary=True, unique=True, default='')
    reactions = OneToManyAttribute('Reaction', related_name='dfba_obj_expression')
    biomass_reactions = OneToManyAttribute('BiomassReaction', related_name='dfba_obj_expression')

    class Meta(obj_model.Model.Meta):
        """
        Attributes:
            valid_functions (:obj:`tuple` of `builtin_function_or_method`): tuple of functions that
                can be used in this `Function`s `expression`
            valid_models (:obj:`tuple` of `str`): names of `obj_model.Model`s in this module that a
                `DfbaObjective` is allowed to reference in its `expression`
        """
        tabular_orientation = TabularOrientation.inline
        valid_functions = ()
        valid_models = ('Reaction', 'BiomassReaction')

    def validate(self):
        """ Determine if the dFBA objective expression is valid

        * Check that the expression is a linear function
        * Check if expression is a function of at least one reaction or biomass reaction
        * Check that the reactions and biomass reactions belong to the same submodel

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors in an `InvalidObject` instance
        """
        invalid_obj = super(DfbaObjectiveExpression, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        expr_errors = []
        if not self.reactions and not self.biomass_reactions:
            expr_errors.append('Expression must be a function of at least one reaction or biomass reaction')

        if self.dfba_obj and self.dfba_obj.submodel:
            missing_rxns = set(self.reactions).difference(set(self.dfba_obj.submodel.reactions))
            missing_bm_rxns = set(self.biomass_reactions).difference(set(self.dfba_obj.submodel.biomass_reactions))

            if missing_rxns:
                expr_errors.append(('dFBA submodel {} must contain the following reactions '
                                    'that are in its objective function:\n  {}').format(
                    self.dfba_obj.submodel.id,
                    '\n  '.join(rxn.id for rxn in missing_rxns)))
            if missing_bm_rxns:
                expr_errors.append(('dFBA submodel {} must contain the following biomass reactions '
                                    'that are in its objective function:\n  {}').format(
                    self.dfba_obj.submodel.id,
                    '\n  '.join(rxn.id for rxn in missing_bm_rxns)))

        if expr_errors:
            errors.append(InvalidAttribute(self.Meta.attributes['expression'], expr_errors))

        if errors:
            return InvalidObject(self, errors)
        return ExpressionMethods.validate(self)

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return ExpressionMethods.serialize(self)

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of :obj:`DfbaObjectiveExpression`, `InvalidAttribute` or `None`: tuple
                of cleaned value and cleaning error
        """
        return ExpressionMethods.deserialize(cls, value, objects)


class DfbaObjective(obj_model.Model):
    """ dFBA objective function

    Attributes:
        id (:obj:`str`): identifier equal to `dfba-obj-{submodel.id}`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        submodel (:obj:`Submodel`): the `Submodel` which uses this `DfbaObjective`
        expression (:obj:`DfbaObjectiveExpression`): mathematical expression of the objective function
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='dfba_objs')
    submodel = OneToOneAttribute(Submodel, related_name='dfba_obj', min_related=1)
    expression = ExpressionAttribute('DfbaObjectiveExpression', related_name='dfba_obj',
                                     min_related=1, min_related_rev=1)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='dfba_objs')

    class Meta(obj_model.Model.Meta):
        verbose_name = 'dFBA objective'
        attribute_order = ('id', 'name', 'submodel', 'expression', 'comments', 'references')
        expression_model = DfbaObjectiveExpression

    @staticmethod
    def gen_id(submodel_id):
        """ Generate identifier

        Args:
            submodel_id (:obj:`str`): submodel id

        Returns:
            :obj:`str`: identifier
        """
        return 'dfba-obj-{}'.format(submodel_id)

    def validate(self):
        """ Validate that the dFBA objective is valid

        * Check if the identifier is equal to `dfba-obj-{submodel.id}]`

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(DfbaObjective, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        if self.submodel and self.id != self.gen_id(self.submodel.id):
            errors.append(InvalidAttribute(self.Meta.attributes['id'],
                                           ['Id must be {}'.format(self.gen_id(self.submodel.id))]))

        if errors:
            return InvalidObject(self, errors)
        return None

    ACTIVE_OBJECTIVE = 'active_objective'

    def add_to_sbml_doc(self, sbml_document):
        """ Add this DfbaObjective to a libsbml SBML document in a `libsbml.model.ListOfObjectives`.

        This uses version 2 of the 'Flux Balance Constraints' extension. SBML assumes that an
        DfbaObjective is a linear combination of reaction fluxes.

        Args:
             sbml_document (:obj:`obj`): a `libsbml` SBMLDocument

        Returns:
            :obj:`libsbml.Objective`: the libsbml Objective that's created

        Raises:
            :obj:`LibSBMLError`: if calling `libsbml` raises an error
        """
        # issue warning if objective function not linear
        if not self.expression._parsed_expression.is_linear:
            warnings.warn("submodel '{}' can't add non-linear objective function to SBML FBC model".format(
                self.submodel.id), UserWarning)
            return
        sbml_model = wrap_libsbml(sbml_document.getModel)
        fbc_model_plugin = wrap_libsbml(sbml_model.getPlugin, 'fbc')
        sbml_objective = wrap_libsbml(fbc_model_plugin.createObjective)
        wrap_libsbml(sbml_objective.setType, 'maximize')
        # In SBML 3 FBC 2, the 'activeObjective' attribute must be set on ListOfObjectives.
        # Since a submodel has only one Objective, it must be the active one.
        wrap_libsbml(sbml_objective.setIdAttribute, DfbaObjective.ACTIVE_OBJECTIVE)
        list_of_objectives = wrap_libsbml(fbc_model_plugin.getListOfObjectives)
        wrap_libsbml(list_of_objectives.setActiveObjective, DfbaObjective.ACTIVE_OBJECTIVE)
        for idx, reaction in enumerate(self.expression.reactions):
            sbml_flux_objective = wrap_libsbml(sbml_objective.createFluxObjective)
            wrap_libsbml(sbml_flux_objective.setReaction, reaction.id)
            wrap_libsbml(sbml_flux_objective.setCoefficient,
                         self.expression._parsed_expression.lin_coeffs[Reaction][reaction])
        for idx, biomass_reaction in enumerate(self.expression.biomass_reactions):
            sbml_flux_objective = wrap_libsbml(sbml_objective.createFluxObjective)
            wrap_libsbml(sbml_flux_objective.setReaction, biomass_reaction.id)
            wrap_libsbml(sbml_flux_objective.setCoefficient,
                         self.expression._parsed_expression.lin_coeffs[BiomassReaction][biomass_reaction])

        return sbml_objective

    def get_products(self, __type=None, **kwargs):
        """ Get the species produced by this objective function

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`Species`: species produced by this objective function
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        products = []
        for reaction in self.expression.reactions:
            if reaction.reversible:
                for part in reaction.participants:
                    if part.species.has_attr_vals(__type=__type, **kwargs):
                        products.append(part.species)
            else:
                for part in reaction.participants:
                    if 0 < part.coefficient:
                        if part.species.has_attr_vals(__type=__type, **kwargs):
                            products.append(part.species)

        tmp_species_ids = []
        for biomass_reaction in self.expression.biomass_reactions:
            for biomass_component in biomass_reaction.biomass_components:
                if 0 < biomass_component.coefficient:
                    tmp_species_ids.append(biomass_component.species.id)
        with self.submodel as submodel:
            tmp_species = Species.get(tmp_species_ids, self.submodel.get_species())
            for tmp_specie_id, tmp_specie in zip(tmp_species_ids, tmp_species):
                if tmp_specie.has_attr_vals(__type=__type, **kwargs):
                    products.append(tmp_specie)
        return det_dedupe(products)


class Compartment(obj_model.Model):
    """ Compartment

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        type (:obj:`CompartmentType`): type
        initial_volume (:obj:`float`): initial volume (L)
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        species (:obj:`list` of :obj:`Species`): species in this compartment
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='compartments')
    type = EnumAttribute(CompartmentType, default=CompartmentType.physical_3d)
    initial_volume = FloatAttribute(min=0)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='compartments')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'type', 'initial_volume',
                           'comments', 'references')

    def add_to_sbml_doc(self, sbml_document):
        """ Add this Compartment to a libsbml SBML document.

        Args:
             sbml_document (:obj:`obj`): a `libsbml` SBMLDocument

        Returns:
            :obj:`libsbml.compartment`: the libsbml compartment that's created

        Raises:
            :obj:`LibSBMLError`: if calling `libsbml` raises an error
        """
        sbml_model = wrap_libsbml(sbml_document.getModel)
        sbml_compartment = wrap_libsbml(sbml_model.createCompartment)
        wrap_libsbml(sbml_compartment.setIdAttribute, self.id)
        wrap_libsbml(sbml_compartment.setName, self.name)
        wrap_libsbml(sbml_compartment.setSpatialDimensions, 3)
        wrap_libsbml(sbml_compartment.setSize, self.initial_volume)
        wrap_libsbml(sbml_compartment.setConstant, False)
        if self.comments:
            wrap_libsbml(sbml_compartment.setNotes, self.comments, True)
        return sbml_compartment


class SpeciesType(obj_model.Model):
    """ Species type

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        structure (:obj:`str`): structure (InChI for metabolites; sequence for DNA, RNA, proteins)
        empirical_formula (:obj:`str`): empirical formula
        molecular_weight (:obj:`float`): molecular weight
        charge (:obj:`int`): charge
        type (:obj:`SpeciesTypeType`): type
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        species (:obj:`list` of :obj:`Species`): species
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
        concentrations (:obj:`list` of :obj:`Concentration`): concentrations
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='species_types')
    structure = LongStringAttribute()
    empirical_formula = RegexAttribute(pattern=r'^([A-Z][a-z]?\d*)*$')
    molecular_weight = PositiveFloatAttribute()
    charge = IntegerAttribute()
    type = EnumAttribute(SpeciesTypeType, default=SpeciesTypeType.metabolite)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='species_types')

    class Meta(obj_model.Model.Meta):
        verbose_name = 'Species type'
        attribute_order = ('id', 'name', 'structure', 'empirical_formula',
                           'molecular_weight', 'charge', 'type', 'comments', 'references')

        indexed_attrs_tuples = (('id',), )

    def is_carbon_containing(self):
        """ Returns `True` is species contains at least one carbon atom.

        Returns:
            :obj:`bool`: `True` is species contains at least one carbon atom.
        """
        # todo: move to compiled model
        return re.match('C[1-9]', self.empirical_formula) is not None


class Species(obj_model.Model):
    """ Species (tuple of species type, compartment)

    Attributes:
        id (:obj:`str`): identifier equal to `{species_type.id}[{compartment.id}]`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        species_type (:obj:`SpeciesType`): species type
        compartment (:obj:`Compartment`): compartment
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        concentration (:obj:`Concentration`): concentration
        species_coefficients (:obj:`list` of :obj:`SpeciesCoefficient`): participations in reactions and observables
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        observable_expressions (:obj:`list` of :obj:`ObservableExpression`): observable expressions
        function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        biomass_components (:obj:`list` of :obj:`BiomassComponent`): biomass components
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='species')
    species_type = ManyToOneAttribute(SpeciesType, related_name='species', min_related=1)
    compartment = ManyToOneAttribute(Compartment, related_name='species', min_related=1)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='species')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'species_type', 'compartment', 'comments', 'references')
        frozen_columns = 1
        # unique_together = (('species_type', 'compartment', ), )
        ordering = ('id',)
        indexed_attrs_tuples = (('species_type', 'compartment'), )
        token_pattern = (token.NAME, token.LSQB, token.NAME, token.RSQB)

    @staticmethod
    def gen_id(species_type_id, compartment_id):
        """ Generate identifier

        Args:
            species_type_id (:obj:`str`): species type id
            compartment_id (:obj:`str`): species type id

        Returns:
            :obj:`str`: identifier
        """
        return '{}[{}]'.format(species_type_id, compartment_id)

    def validate(self):
        """ Check that the species is valid

        * Check if the identifier is equal to `{species_type.id}[{compartment.id}]`

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(Species, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        if self.id != self.gen_id(self.species_type.id, self.compartment.id):
            errors.append(InvalidAttribute(self.Meta.attributes['id'],
                                           ['Id must be {}'.format(self.gen_id(self.species_type.id, self.compartment.id))]))

        if errors:
            return InvalidObject(self, errors)
        return None

    @staticmethod
    def get(ids, species_iterator):
        """ Find some Species instances

        Args:
            ids (:obj:`Iterator` of `str`): an iterator over some species identifiers
            species_iterator (:obj:`Iterator`): an iterator over some species

        Returns:
            :obj:`list` of :obj:`Species` or `None`: each element of the `list` corresponds to an element
                of `ids` and contains either a `Species` with `id()` equal to the element in `ids`,
                or `None` indicating that `species_iterator` does not contain a matching `Species`
        """
        # TODO: this costs O(|ids||species_iterator|); replace with O(|ids|) operation using obj_model.Manager.get()
        rv = []
        for id in ids:
            s = None
            for specie in species_iterator:
                if specie.id == id:
                    s = specie
                    break
            rv.append(s)
        return rv

    def gen_sbml_id(self):
        """ Make a Species id that satisfies the SBML string id syntax.

        Replaces the '[' and ']' in Species.id with double-underscores '__'.
        See Finney and Hucka, "Systems Biology Markup Language (SBML) Level 2: Structures and
        Facilities for Model Definitions", 2003, section 3.4.

        Returns:
            :obj:`str`: SBML id
        """
        return '{}__{}__'.format(self.species_type.id, self.compartment.id)

    @staticmethod
    def sbml_id_to_id(sbml_id):
        """ Convert an `sbml_id` to its species id.

        Returns:
            :obj:`str`: species id
        """
        return sbml_id.replace('__', '[', 1).replace('__', ']', 1)

    def add_to_sbml_doc(self, sbml_document):
        """ Add this Species to a libsbml SBML document.

        Args:
             sbml_document (:obj:`obj`): a `libsbml` SBMLDocument

        Returns:
            :obj:`libsbml.species`: the libsbml species that's created

        Raises:
            :obj:`LibSBMLError`: if calling `libsbml` raises an error
        """
        sbml_model = wrap_libsbml(sbml_document.getModel)
        sbml_species = wrap_libsbml(sbml_model.createSpecies)
        # initDefaults() isn't wrapped in wrap_libsbml because it returns None
        sbml_species.initDefaults()
        wrap_libsbml(sbml_species.setIdAttribute, self.gen_sbml_id())

        # add some SpeciesType data
        wrap_libsbml(sbml_species.setName, self.species_type.name)
        if self.species_type.comments:
            wrap_libsbml(sbml_species.setNotes, self.species_type.comments, True)

        # set Compartment, which must already be in the SBML document
        wrap_libsbml(sbml_species.setCompartment, self.compartment.id)

        # set the Initial Concentration
        wrap_libsbml(sbml_species.setInitialConcentration, self.concentration.value)

        # set units
        unit_xml_id = ConcentrationUnit.Meta[self.concentration.units]['xml_id']
        wrap_libsbml(sbml_species.setSubstanceUnits, unit_xml_id)

        return sbml_species


class Concentration(obj_model.Model):
    """ Species concentration

    Attributes:
        id (:obj:`str`): identifier equal to `conc-{species.id}`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        species (:obj:`Species`): species
        value (:obj:`float`): value
        units (:obj:`ConcentrationUnit`): units; default units is `M`
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='concentrations')
    species = OneToOneAttribute(Species, related_name='concentration')
    value = FloatAttribute(min=0)
    units = EnumAttribute(ConcentrationUnit, default=ConcentrationUnit.M)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='concentrations')

    class Meta(obj_model.Model.Meta):
        # unique_together = (('species', ), )
        attribute_order = ('id', 'name', 'species', 'value', 'units', 'comments', 'references')

        frozen_columns = 1
        ordering = ('id',)

    @staticmethod
    def gen_id(species_id):
        """ Generate string representation

        Args:
            species_id (:obj:`str`): species id

        Returns:
            :obj:`str`: string representation
        """
        return 'conc-{}'.format(species_id)

    def validate(self):
        """ Check that the concentration is valid

        * Validate that identifier is equal to `conc-{species.id}]`

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(Concentration, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        if self.id != self.gen_id(self.species.id):
            errors.append(InvalidAttribute(self.Meta.attributes['id'],
                                           ['Id must be {}'.format(self.gen_id(self.species.id))]))

        if errors:
            return InvalidObject(self, errors)
        return None


class ExpressionMethods(object):
    """ Generic methods for mathematical expressions
    """

    @staticmethod
    def serialize(model_obj):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return model_obj.expression

    @classmethod
    def deserialize(cls, model_cls, value, objects):
        """ Deserialize expression

        Args:
            model_cls (:obj:`type`): expression class
            value (:obj:`str`): string representation of the mathematical expression, in a
                Python expression
            objects (:obj:`dict`): dictionary of objects which can be used in `expression`, grouped by model

        Returns:
            :obj:`tuple`: on error return (`None`, `InvalidAttribute`),
                otherwise return (object in this class with instantiated `_parsed_expression`, `None`)
        """
        # objects must contain all objects types in valid_models
        value = value or ''

        used_model_types = []
        for used_model in model_cls.Meta.valid_models:
            used_model_type = globals()[used_model]
            used_model_types.append(used_model_type)
        expr_field = 'expression'
        try:
            _parsed_expression = ParsedExpression(model_cls, expr_field, value, objects)
        except ParsedExpressionError as e:
            attr = model_cls.Meta.attributes['expression']
            return (None, InvalidAttribute(attr, [str(e)]))
        rv = _parsed_expression.tokenize()
        if rv[0] is None:
            attr = model_cls.Meta.attributes['expression']
            errors = rv[1]
            return (None, InvalidAttribute(attr, errors))
        _, used_objects = rv
        if model_cls not in objects:
            objects[model_cls] = {}
        if value in objects[model_cls]:
            obj = objects[model_cls][value]
        else:
            obj = model_cls(expression=value)
            objects[model_cls][value] = obj

            for attr_name, attr in model_cls.Meta.attributes.items():
                if isinstance(attr, obj_model.RelatedAttribute) and \
                        attr.related_class.__name__ in model_cls.Meta.valid_models:
                    attr_value = list(used_objects.get(attr.related_class, {}).values())
                    setattr(obj, attr_name, attr_value)
        obj._parsed_expression = _parsed_expression

        # check expression is linear
        obj._parsed_expression.is_linear, _ = LinearExpressionVerifier().validate(
            obj._parsed_expression.wc_tokens)
        cls.set_lin_coeffs(obj)

        return (obj, None)

    @classmethod
    def set_lin_coeffs(cls, obj):
        """ Set the linear coefficients for the related objects

        Args:
            obj (:obj:`obj_model.Model`): expression object
        """
        model_cls = obj.__class__
        parsed_expr = obj._parsed_expression
        tokens = parsed_expr.wc_tokens
        is_linear = parsed_expr.is_linear

        if is_linear:
            default_val = 0.
        else:
            default_val = float('nan')

        parsed_expr.lin_coeffs = lin_coeffs = {}
        for attr_name, attr in model_cls.Meta.attributes.items():
            if isinstance(attr, obj_model.RelatedAttribute) and \
                    attr.related_class.__name__ in model_cls.Meta.valid_models:
                lin_coeffs[attr.related_class] = {}

        for related_class, related_objs in parsed_expr.related_objects.items():
            for related_obj in related_objs.values():
                lin_coeffs[related_class][related_obj] = default_val

        if not is_linear:
            return

        sense = 1.
        cur_coeff = 1.
        for token in tokens:
            if token.tok_code == TokCodes.op and token.token_string == '+':
                sense = 1.
                cur_coeff = 1.
            elif token.tok_code == TokCodes.op and token.token_string == '-':
                sense = -1.
                cur_coeff = 1.
            elif token.tok_code == TokCodes.number:
                cur_coeff = float(token.token_string)
            elif token.tok_code == TokCodes.wc_lang_obj_id:
                lin_coeffs[token.model_type][token.model] += sense * cur_coeff

    @classmethod
    def validate(cls, model_obj, return_type=None, check_linear=False):
        """ Determine whether an expression model is valid by eval'ing its deserialized expression

        Args:
            model_obj (`Expression`): expression object
            return_type (:obj:`type`, optional): if provided, an expression's required return type
            check_linear (:obj:`bool`, optional): if :obj:`True`, validate that the expression is a
                linear function

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors in an `InvalidObject` instance
        """
        model_cls = model_obj.__class__

        # generate _parsed_expression
        objects = {}
        for related_attr_name, related_attr in model_cls.Meta.attributes.items():
            if isinstance(related_attr, obj_model.RelatedAttribute):
                objects[related_attr.related_class] = {
                    m.get_primary_attribute(): m for m in getattr(model_obj, related_attr_name)
                }
        try:
            model_obj._parsed_expression = ParsedExpression(model_obj.__class__, 'expression',
                                                            model_obj.expression, objects)
        except ParsedExpressionError as e:
            attr = model_cls.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, [str(e)])
            return InvalidObject(model_obj, [attr_err])

        is_valid, errors = model_obj._parsed_expression.tokenize()
        if is_valid is None:
            attr = model_cls.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, errors)
            return InvalidObject(model_obj, [attr_err])
        model_obj._parsed_expression.is_linear, _ = LinearExpressionVerifier().validate(
            model_obj._parsed_expression.wc_tokens)
        cls.set_lin_coeffs(model_obj)

        # check related objects matches the tokens of the _parsed_expression
        related_objs = {}
        for related_attr_name, related_attr in model_cls.Meta.attributes.items():
            if isinstance(related_attr, obj_model.RelatedAttribute):
                related_model_objs = getattr(model_obj, related_attr_name)
                if related_model_objs:
                    related_objs[related_attr.related_class] = set(related_model_objs)

        token_objs = {}
        token_obj_ids = {}
        for token in model_obj._parsed_expression.wc_tokens:
            if token.model_type is not None:
                if token.model_type not in token_objs:
                    token_objs[token.model_type] = set()
                    token_obj_ids[token.model_type] = set()
                token_objs[token.model_type].add(token.model)
                token_obj_ids[token.model_type].add(token.token_string)

        if related_objs != token_objs:
            attr = model_cls.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, ['Related objects must match the tokens of the analyzed expression'])
            return InvalidObject(model_obj, [attr_err])

        # check expression is valid
        try:
            rv = model_obj._parsed_expression.test_eval()
            if return_type is not None:
                if not isinstance(rv, return_type):
                    attr = model_cls.Meta.attributes['expression']
                    attr_err = InvalidAttribute(attr,
                                                ["Evaluating '{}', a {} expression, should return a {} but it returns a {}".format(
                                                    model_obj.expression, model_obj.__class__.__name__,
                                                    return_type.__name__, type(rv).__name__)])
                    return InvalidObject(model_obj, [attr_err])
        except ParsedExpressionError as e:
            attr = model_cls.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, [str(e)])
            return InvalidObject(model_obj, [attr_err])

        # check expression is linear
        if check_linear and not model_obj._parsed_expression.is_linear:
            attr = model_cls.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, ['Expression must be linear'])
            return InvalidObject(model_obj, [attr_err])

        # return `None` to indicate valid object
        return None

    @staticmethod
    def make_expression_obj(model_type, expression, objects):
        """ Make an expression object

        Args:
            model_type (:obj:`type`): an `obj_model.Model` that uses a mathemetical expression, like
                `Function` and `Observable`
            expression (:obj:`str`): the expression used by the `model_type` being created
            objects (:obj:`dict` of `dict`): all objects that are referenced in `expression`

        Returns:
            :obj:`tuple`: if successful, (`obj_model.Model`, `None`) containing a new instance of
                `model_type`'s expression helper class; otherwise, (`None`, `InvalidAttribute`)
                reporting the error
        """
        expr_model_type = model_type.Meta.expression_model
        return expr_model_type.deserialize(expression, objects)

    @classmethod
    def make_obj(cls, model, model_type, id, expression, objects, allow_invalid_objects=False):
        """ Make a model that contains an expression by using its expression helper class

        For example, this uses `FunctionExpression` to make a `Function`.

        Args:
            model (:obj:`obj_model.Model`): a `wc_lang.core.Model` which is the root model
            model_type (:obj:`type`): an `obj_model.Model` that uses a mathemetical expression, like
                `Function` and `Observable`
            id (:obj:`str`): the id of the `model_type` being created
            expression (:obj:`str`): the expression used by the `model_type` being created
            objects (:obj:`dict` of `dict`): all objects that are referenced in `expression`
            allow_invalid_objects (:obj:`bool`, optional): if set, return object - not error - if
                the expression object does not validate

        Returns:
            :obj:`obj_model.Model` or `InvalidAttribute`: a new instance of `model_type`, or,
                if an error occurs, an `InvalidAttribute` reporting the error
        """
        expr_model_obj, error = cls.make_expression_obj(model_type, expression, objects)
        if error:
            return error
        error_or_none = expr_model_obj.validate()
        if error_or_none is not None and not allow_invalid_objects:
            return error_or_none
        related_name = model_type.Meta.attributes['model'].related_name
        related_in_model = getattr(model, related_name)
        new_obj = related_in_model.create(id=id, expression=expr_model_obj)
        return new_obj


class ObservableExpression(obj_model.Model):
    """ A mathematical expression of Observables and Species

    The expression used by a `Observable`.

    Attributes:
        expression (:obj:`str`): mathematical expression for an Observable
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        species (:obj:`list` of :obj:`Species`): Species used by this Observable expression
        observables (:obj:`list` of :obj:`Observable`): other Observables used by this Observable expression
    """

    expression = LongStringAttribute(primary=True, unique=True, default='')
    species = ManyToManyAttribute(Species, related_name='observable_expressions')
    observables = ManyToManyAttribute('Observable', related_name='observable_expressions')

    class Meta(obj_model.Model.Meta):
        """
        Attributes:
            valid_Observables (:obj:`tuple` of `builtin_Observable_or_method`): tuple of Observables that
                can be used in this `Observable`s `expression`
            valid_models (:obj:`tuple` of `str`): names of `obj_model.Model`s in this module that a
                `Observable` is allowed to reference in its `expression`
        """
        tabular_orientation = TabularOrientation.inline
        valid_models = ('Species', 'Observable')

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return ExpressionMethods.serialize(self)

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of :obj:`ObservableExpression`, `InvalidAttribute` or `None`:
                tuple of cleaned value and cleaning error
        """
        return ExpressionMethods.deserialize(cls, value, objects)

    def validate(self):
        """ Check that the observable is valid

        * Check that the expression is a linear function

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        return ExpressionMethods.validate(self, check_linear=True)


class Observable(obj_model.Model):
    """ Observable: a linear function of other Observbles and Species

    Attributes:
        id (:obj:`str`): unique id
        name (:obj:`str`): name
        model (:obj:`Model`): model
        expression (:obj:`ObservableExpression`): mathematical expression for an Observable
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        observable_expressions (:obj:`list` of :obj:`ObservableExpression`): observable expressions
        function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        stop_condition_expressions (:obj:`list` of :obj:`StopConditionExpression`): stop condition expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='observables')
    expression = ExpressionAttribute('ObservableExpression', related_name='observable')
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='observables')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'expression', 'comments', 'references')
        expression_model = ObservableExpression


class FunctionExpression(obj_model.Model):
    """ A mathematical expression of Functions, Observbles, Parameters and Python functions

    The expression used by a :obj:`Function`.

    Attributes:
        expression (:obj:`str`): mathematical expression for a Function
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        species (:obj:`list` of :obj:`Species`): Species used by this function expression
        observables (:obj:`list` of :obj:`Observable`): Observables used by this function expression
        parameters (:obj:`list` of :obj:`Parameter`): Parameters used by this function expression
        functions (:obj:`list` of :obj:`Function`): other Functions used by this function expression

    Related attributes:
        function (:obj:`Function`): function
    """
    expression = LongStringAttribute(primary=True, unique=True, default='')
    species = ManyToManyAttribute(Species, related_name='function_expressions')
    observables = ManyToManyAttribute(Observable, related_name='function_expressions')
    parameters = ManyToManyAttribute('Parameter', related_name='function_expressions')
    functions = ManyToManyAttribute('Function', related_name='function_expressions')

    class Meta(obj_model.Model.Meta):
        """
        Attributes:
            valid_functions (:obj:`tuple` of `builtin_function_or_method`): tuple of functions that
                can be used in this `Function`s `expression`
            valid_models (:obj:`tuple` of `str`): names of `obj_model.Model`s in this module that a
                `Function` is allowed to reference in its `expression`
        """
        tabular_orientation = TabularOrientation.inline
        valid_functions = (ceil, floor, exp, pow, log, log10, min, max)
        valid_models = ('Parameter', 'Species', 'Observable', 'Function')

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return ExpressionMethods.serialize(self)

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of :obj:`FunctionExpression`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        return ExpressionMethods.deserialize(cls, value, objects)

    def validate(self):
        """ Check that the function is valid

        * Check that the expression is a valid Python function

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        return ExpressionMethods.validate(self)


class Function(obj_model.Model):
    """ Function: a mathematical expression of Functions, Observbles, Parameters and Python functions

    Attributes:
        id (:obj:`str`): unique id
        name (:obj:`str`): name
        model (:obj:`Model`): model
        expression (:obj:`FunctionExpression`): mathematical expression for a Function
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        stop_condition_expressions (:obj:`list` of :obj:`StopConditionExpression`): stop condition expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='functions')
    expression = ExpressionAttribute('FunctionExpression', related_name='function')
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='functions')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'expression', 'comments', 'references')
        expression_model = FunctionExpression


class StopConditionExpression(obj_model.Model):
    """ A mathematical expression of Functions, Observables, Parameters and Python functions

    The expression used by a :obj:`StopCondition`.

    Attributes:
        expression (:obj:`str`): mathematical expression for a StopCondition
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        observables (:obj:`list` of :obj:`Observable`): Observables used by this stop condition expression
        parameters (:obj:`list` of :obj:`Parameter`): Parameters used by this stop condition expression
        functions (:obj:`list` of :obj:`Function`): Functions used by this stop condition expression
    """

    expression = LongStringAttribute(primary=True, unique=True, default='')
    observables = ManyToManyAttribute(Observable, related_name='stop_condition_expressions')
    parameters = ManyToManyAttribute('Parameter', related_name='stop_condition_expressions')
    functions = ManyToManyAttribute('Function', related_name='stop_condition_expressions')

    class Meta(obj_model.Model.Meta):
        """
        Attributes:
            valid_functions (:obj:`tuple` of `builtin_function_or_method`): tuple of functions that
                can be used in this `StopCondition`s `expression`
            valid_models (:obj:`tuple` of `str`): names of `obj_model.Model`s in this module that a
                `StopCondition` is allowed to reference in its `expression`
        """
        tabular_orientation = TabularOrientation.inline
        valid_functions = (ceil, floor, exp, pow, log, log10, min, max)
        valid_models = ('Parameter', 'Observable', 'Function')

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return ExpressionMethods.serialize(self)

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of :obj:`StopConditionExpression`, `InvalidAttribute` or `None`: tuple of
                cleaned value and cleaning error
        """
        return ExpressionMethods.deserialize(cls, value, objects)

    def validate(self):
        """ Check that the stop condition is valid

        * Check that the expression is a Boolean function

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        return ExpressionMethods.validate(self, return_type=bool)


class StopCondition(obj_model.Model):
    """ StopCondition: Simulation of a model terminates when its StopCondition is true.

    A mathematical expression of Functions, Observbles, Parameters and Python functions `StopCondition`s
    are optional. It must return a Boolean.

    Attributes:
        id (:obj:`str`): unique id
        name (:obj:`str`): name
        model (:obj:`Model`): model
        expression (:obj:`StopConditionExpression`): mathematical expression for a StopCondition
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        expressions (:obj:`Expressions`): expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='stop_conditions')
    expression = ExpressionAttribute('StopConditionExpression', related_name='stop_condition')
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='stop_conditions')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'expression', 'comments', 'references')
        expression_model = StopConditionExpression


class Reaction(obj_model.Model):
    """ Reaction

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        submodel (:obj:`Submodel`): submodel that reaction belongs to
        participants (:obj:`list` of :obj:`SpeciesCoefficient`): participants
        reversible (:obj:`bool`): indicates if reaction is thermodynamically reversible
        min_flux (:obj:`float`): minimum flux bound for solving an FBA model; negative for reversible reactions
        max_flux (:obj:`float`): maximum flux bound for solving an FBA model
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws; if present, rate_laws[0] is the forward
            rate law, and rate_laws[0] is the backward rate law
        dfba_obj_expression (:obj:`DfbaObjectiveExpression`): dFBA objectie expression
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='reactions')
    submodel = ManyToOneAttribute(Submodel, related_name='reactions')
    participants = ReactionParticipantAttribute(related_name='reactions')
    reversible = BooleanAttribute()
    min_flux = FloatAttribute(nan=True)
    max_flux = FloatAttribute(min=0, nan=True)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='reactions')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'submodel', 'participants', 'reversible', 'min_flux', 'max_flux', 'comments', 'references')
        indexed_attrs_tuples = (('id',), )

    def validate(self):
        """ Check if the reaction is valid

        * If the submodel is ODE or SSA, check that the reaction has a forward rate law
        * If the submodel is ODE or SSA and the reaction is reversible, check that the reaction has a
          backward rate law
        * Check that `min_flux` <= `max_flux`

        .. todo :: Check reaction is mass-balanced
        .. todo :: Check reaction is charge-balanced

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors in an `InvalidObject` instance
        """
        invalid_obj = super(Reaction, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        # check that rate laws are defined as needed for ODE, SSA submodels
        rl_errors = []

        for_rl = self.rate_laws.get_one(direction=RateLawDirection.forward)
        rev_rl = self.rate_laws.get_one(direction=RateLawDirection.backward)
        if self.submodel and self.submodel.algorithm in [SubmodelAlgorithm.ode, SubmodelAlgorithm.ssa]:
            if not for_rl:
                rl_errors.append('Reaction in {} submodel must have a forward rate law'.format(
                    self.submodel.algorithm.name))
            if self.reversible and not rev_rl:
                rl_errors.append('Reversible reaction in {} submodel must have a backward rate law'.format(
                    self.submodel.algorithm.name))
        if not self.reversible and rev_rl:
            rl_errors.append('Irreversible reaction in {} submodel cannot have a backward rate law'.format(
                self.submodel.algorithm.name))

        if rl_errors:
            errors.append(InvalidAttribute(self.Meta.related_attributes['rate_laws'], rl_errors))

        # check min, max fluxes
        if self.submodel and self.submodel.algorithm is not SubmodelAlgorithm.dfba:
            if not isnan(self.min_flux):
                errors.append(InvalidAttribute(self.Meta.attributes['min_flux'],
                                               ['Minimum flux should be NaN for reactions in non-dFBA submodels']))
            if not isnan(self.max_flux):
                errors.append(InvalidAttribute(self.Meta.attributes['max_flux'],
                                               ['Maximum flux should be NaN for reactions in non-dFBA submodels']))

        if not isnan(self.min_flux) and not isnan(self.min_flux) and self.min_flux > self.max_flux:
            errors.append(InvalidAttribute(self.Meta.attributes['max_flux'],
                                           ['Maximum flux must be least the minimum flux']))

        if self.reversible and not isnan(self.min_flux) and self.min_flux >= 0:
            errors.append(InvalidAttribute(self.Meta.attributes['min_flux'],
                                           ['Minimum flux for reversible reaction should be negative or NaN']))
        if not self.reversible and not isnan(self.min_flux) and self.min_flux < 0:
            errors.append(InvalidAttribute(self.Meta.attributes['min_flux'],
                                           ['Minimum flux for irreversible reaction should be non-negative']))

        # return errors
        if errors:
            return InvalidObject(self, errors)
        return None

    def get_species(self, __type=None, **kwargs):
        """ Get species

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list`: list of `Species`
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        species = []

        for part in self.participants:
            if part.species.has_attr_vals(__type=__type, **kwargs):
                species.append(part.species)

        for rate_law in self.rate_laws:
            if rate_law.expression:
                species.extend(rate_law.expression.modifiers.get(__type=__type, **kwargs))

        return det_dedupe(species)

    def add_to_sbml_doc(self, sbml_document):
        """ Add this Reaction to a libsbml SBML document.

        Args:
             sbml_document (:obj:`obj`): a `libsbml` SBMLDocument

        Returns:
            :obj:`libsbml.reaction`: the libsbml reaction that's created

        Raises:
            :obj:`LibSBMLError`: if calling `libsbml` raises an error
        """
        sbml_model = wrap_libsbml(sbml_document.getModel)

        # create SBML reaction in SBML document
        sbml_reaction = wrap_libsbml(sbml_model.createReaction)
        wrap_libsbml(sbml_reaction.setIdAttribute, self.id)
        wrap_libsbml(sbml_reaction.setName, self.name)
        wrap_libsbml(sbml_reaction.setReversible, self.reversible)
        wrap_libsbml(sbml_reaction.setFast, False)
        if self.comments:
            wrap_libsbml(sbml_reaction.setNotes, self.comments, True)

        # write reaction participants to SBML document
        for participant in self.participants:
            if participant.coefficient < 0:
                species_reference = wrap_libsbml(sbml_reaction.createReactant)
                wrap_libsbml(species_reference.setStoichiometry, -participant.coefficient)
            elif 0 < participant.coefficient:
                species_reference = wrap_libsbml(sbml_reaction.createProduct)
                wrap_libsbml(species_reference.setStoichiometry, participant.coefficient)
            wrap_libsbml(species_reference.setSpecies, participant.species.gen_sbml_id())
            wrap_libsbml(species_reference.setConstant, True)

        # for dFBA submodels, write flux bounds to SBML document
        # uses version 2 of the 'Flux Balance Constraints' extension
        if self.submodel.algorithm == SubmodelAlgorithm.dfba:
            fbc_reaction_plugin = wrap_libsbml(sbml_reaction.getPlugin, 'fbc')
            for bound in ['lower', 'upper']:
                # make a unique ID for each flux bound parameter
                # ids for wc_lang Parameters all start with 'parameter'
                param_id = "_reaction_{}_{}_bound".format(self.id, bound)
                param = create_sbml_parameter(sbml_model, id=param_id, value=self.min_flux,
                                              units='mmol_per_gDW_per_hr')
                if bound == 'lower':
                    wrap_libsbml(param.setValue, self.min_flux)
                    wrap_libsbml(fbc_reaction_plugin.setLowerFluxBound, param_id)
                if bound == 'upper':
                    wrap_libsbml(param.setValue, self.max_flux)
                    wrap_libsbml(fbc_reaction_plugin.setUpperFluxBound, param_id)
        return sbml_reaction


class SpeciesCoefficient(obj_model.Model):
    """ A tuple of a species and a coefficient

    Attributes:
        species (:obj:`Species`): species
        coefficient (:obj:`float`): coefficient

    Related attributes:
        reaction (:obj:`Reaction`): reaction
        observables (:obj:`Observable`): observables
    """
    species = ManyToOneAttribute(Species, related_name='species_coefficients')
    coefficient = FloatAttribute(nan=False)

    class Meta(obj_model.Model.Meta):
        unique_together = (('species', 'coefficient'),)
        attribute_order = ('species', 'coefficient')
        frozen_columns = 1
        tabular_orientation = TabularOrientation.inline
        ordering = ('species',)

    def serialize(self, show_compartment=True, show_coefficient_sign=True):
        """ Serialize related object

        Args:
            show_compartment (:obj:`bool`, optional): if true, show compartment
            show_coefficient_sign (:obj:`bool`, optional): if true, show coefficient sign

        Returns:
            :obj:`str`: simple Python representation
        """
        return self._serialize(self.species, self.coefficient,
                               show_compartment=show_compartment, show_coefficient_sign=show_coefficient_sign)

    @staticmethod
    def _serialize(species, coefficient, show_compartment=True, show_coefficient_sign=True):
        """ Serialize values

        Args:
            species (:obj:`Species`): species
            coefficient (:obj:`float`): coefficient
            show_compartment (:obj:`bool`, optional): if true, show compartment
            show_coefficient_sign (:obj:`bool`, optional): if true, show coefficient sign

        Returns:
            :obj:`str`: simple Python representation
        """
        coefficient = float(coefficient)

        if not show_coefficient_sign:
            coefficient = abs(coefficient)

        if coefficient == 1:
            coefficient_str = ''
        elif coefficient % 1 == 0 and abs(coefficient) < 1000:
            coefficient_str = '({:.0f}) '.format(coefficient)
        else:
            coefficient_str = '({:e}) '.format(coefficient)

        if show_compartment:
            return '{}{}'.format(coefficient_str, species.serialize())
        else:
            return '{}{}'.format(coefficient_str, species.species_type.get_primary_attribute())

    @classmethod
    def deserialize(cls, value, objects, compartment=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            compartment (:obj:`Compartment`, optional): compartment

        Returns:
            :obj:`tuple` of :obj:`SpeciesCoefficient`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        errors = []

        if compartment:
            pattern = r'^(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*)$'
        else:
            pattern = r'^(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*\[[a-z][a-z0-9_]*\])$'

        match = re.match(pattern, value, flags=re.I)
        if match:
            errors = []

            coefficient = float(match.group(2) or 1.)

            if compartment:
                species_id = Species.gen_id(match.group(5), compartment.get_primary_attribute())
            else:
                species_id = match.group(5)

            species, error = Species.deserialize(species_id, objects)
            if error:
                return (None, error)

            serialized_val = cls._serialize(species, coefficient)
            if cls not in objects:
                objects[cls] = {}
            if serialized_val in objects[cls]:
                obj = objects[cls][serialized_val]
            else:
                obj = cls(species=species, coefficient=coefficient)
                objects[cls][serialized_val] = obj
            return (obj, None)

        else:
            attr = cls.Meta.attributes['species']
            return (None, InvalidAttribute(attr, ['Invalid species coefficient']))


class RateLaw(obj_model.Model):
    """ Rate law

    Attributes:
        id (:obj:`str`): identifier equal to `{reaction.id}-{direction.name}`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        reaction (:obj:`Reaction`): reaction
        direction (:obj:`RateLawDirection`): direction
        type (:obj:`RateLawType`): type
        expression (:obj:`RateLawExpression`): expression
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='rate_laws')
    reaction = ManyToOneAttribute(Reaction, related_name='rate_laws')
    direction = EnumAttribute(RateLawDirection, default=RateLawDirection.forward)
    type = EnumAttribute(RateLawType, default=RateLawType.other)
    expression = RateLawExpressionAttribute(related_name='rate_laws')
    units = EnumAttribute(RateLawUnits, default=RateLawUnits['s^-1'])
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='rate_laws')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'reaction', 'direction', 'type',
                           'expression', 'units',
                           'comments', 'references')
        # unique_together = (('reaction', 'direction'), )
        ordering = ('id',)

    @staticmethod
    def gen_id(reaction_id, direction_name):
        """ Generate identifier

        Args:
            reaction_id (:obj:`str`): reaction id
            direction_name (:obj:`str`): direction name

        Returns:
            :obj:`str`: identifier
        """
        return '{}-{}'.format(reaction_id, direction_name)

    def validate(self):
        """ Determine whether this `RateLaw` is valid

        * Check if identifier equal to `{reaction.id}-{direction.name}`

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors in an `InvalidObject` instance
        """
        invalid_obj = super(RateLaw, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        if self.reaction and self.direction is not None and \
                self.id != self.gen_id(self.reaction.id, self.direction.name):
            errors.append(InvalidAttribute(
                self.Meta.attributes['id'],
                ['Id must be {}'.format(self.gen_id(self.reaction.id, self.direction.name))]))

        """ return errors or `None` to indicate valid object """
        if errors:
            return InvalidObject(self, errors)
        return None


class RateLawExpression(obj_model.Model):
    """ Rate law expression

    Attributes:
        expression (:obj:`str`): mathematical expression of the rate law
        modifiers (:obj:`list` of :obj:`Species`): species whose concentrations are used in the rate law
        parameters (:obj:`list` of :obj:`Parameter`): parameters whose values are used in the rate law

    Related attributes:
        rate_law (:obj:`RateLaw`): the `RateLaw` which uses this `RateLawExpression`
    """
    expression = LongStringAttribute(primary=True, unique=True, default='')
    modifiers = ManyToManyAttribute(Species, related_name='rate_law_expressions')
    parameters = ManyToManyAttribute('Parameter', related_name='rate_law_expressions')
    observables = ManyToManyAttribute('Observable', related_name='rate_law_expressions')
    functions = ManyToManyAttribute('Function', related_name='rate_law_expressions')

    class Meta(obj_model.Model.Meta):
        """
        Attributes:
            valid_functions (:obj:`tuple` of `builtin_function_or_method`): tuple of functions that
                can be used in a `RateLawExpression`s `expression`
            valid_models (:obj:`tuple` of `str`): names of `obj_model.Model`s in this module that a
                `RateLawExpression` is allowed to reference in its `expression`
        """
        attribute_order = ('expression', 'modifiers', 'parameters')
        tabular_orientation = TabularOrientation.inline
        ordering = ('rate_law',)
        valid_functions = (ceil, floor, exp, pow, log, log10, min, max)
        valid_models = ('Parameter', 'Species', 'Observable', 'Function')

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return ExpressionMethods.serialize(self)

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of :obj:`RateLawExpression`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        return ExpressionMethods.deserialize(cls, value, objects)

    def validate(self):
        """ Determine whether a `RateLawExpression` is valid

        * Check that all of the modifiers and parameters contribute to the expression
        * Check that the modifiers and parameters encompass of the named entities in the expression

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors in an `InvalidObject` instance
        """
        return ExpressionMethods.validate(self)


class BiomassComponent(obj_model.Model):
    """ BiomassComponent

    A biomass reaction contains a list of BiomassComponent instances. Distinct BiomassComponents
    enable separate comments and references for each one.

    Attributes:
        id (:obj:`str`): unique identifier per BiomassComponent
        name (:obj:`str`): name
        biomass_reaction (:obj:`BiomassReaction`): the biomass reaction that uses the biomass component
        coefficient (:obj:`float`): the specie's reaction coefficient
        species (:obj:`Species`): species
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = SlugAttribute()
    name = StringAttribute()
    biomass_reaction = ManyToOneAttribute('BiomassReaction', related_name='biomass_components')
    coefficient = FloatAttribute()
    species = ManyToOneAttribute(Species, related_name='biomass_components')
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='biomass_components')

    class Meta(obj_model.Model.Meta):
        unique_together = (('biomass_reaction', 'species'), )
        attribute_order = ('id', 'name', 'biomass_reaction',
                           'coefficient', 'species',
                           'comments', 'references')


class BiomassReaction(obj_model.Model):
    """ A pseudo-reaction used to represent the interface between metabolism and other
    cell processes.

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        submodel (:obj:`Submodel`): submodel that uses this reaction
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        dfba_obj_expression (:obj:`DfbaObjectiveExpression`): dFBA objectie expression
        biomass_components (:obj:`list` of :obj:`BiomassComponent`): the components of this biomass reaction
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='biomass_reactions')
    submodel = ManyToOneAttribute('Submodel', related_name='biomass_reactions')
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='biomass_reactions')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'submodel', 'comments', 'references')
        indexed_attrs_tuples = (('id',), )

    def add_to_sbml_doc(self, sbml_document):
        """ Add a BiomassReaction to a libsbml SBML document.

        BiomassReactions are added to the SBML document because they can be used in a dFBA submodel's
        objective function. In fact the default objective function is the submodel's biomass reaction.
        Since SBML does not define BiomassReaction as a separate class, BiomassReactions are added
        to the SBML model as SBML reactions.
        CheckModel ensures that wc_lang BiomassReactions and Reactions have distinct ids.

        Args:
             sbml_document (:obj:`obj`): a `libsbml` SBMLDocument

        Returns:
            :obj:`libsbml.reaction`: the libsbml reaction that's created

        Raises:
            :obj:`LibSBMLError`: if calling `libsbml` raises an error
        """
        sbml_model = wrap_libsbml(sbml_document.getModel)

        # create SBML reaction in SBML document
        sbml_reaction = wrap_libsbml(sbml_model.createReaction)
        wrap_libsbml(sbml_reaction.setIdAttribute, self.id)
        wrap_libsbml(sbml_reaction.setName, self.name)
        wrap_libsbml(sbml_reaction.setReversible, False)
        wrap_libsbml(sbml_reaction.setFast, False)
        if self.comments:
            wrap_libsbml(sbml_reaction.setNotes, self.comments, True)

        # write biomass reaction participants to SBML document
        for biomass_component in self.biomass_components:
            if biomass_component.coefficient < 0:
                species_reference = wrap_libsbml(sbml_reaction.createReactant)
                wrap_libsbml(species_reference.setStoichiometry, -biomass_component.coefficient)
            elif 0 < biomass_component.coefficient:
                species_reference = wrap_libsbml(sbml_reaction.createProduct)
                wrap_libsbml(species_reference.setStoichiometry, biomass_component.coefficient)
            id = biomass_component.species.gen_sbml_id()
            wrap_libsbml(species_reference.setSpecies, id)
            wrap_libsbml(species_reference.setConstant, True)

        # the biomass reaction does not constrain the optimization, so set its bounds to 0 and INF
        fbc_reaction_plugin = wrap_libsbml(sbml_reaction.getPlugin, 'fbc')
        for bound in ['lower', 'upper']:
            # make a unique ID for each flux bound parameter
            # ids for wc_lang Parameters all start with 'parameter'
            param_id = "_biomass_reaction_{}_{}_bound".format(self.id, bound)
            param = create_sbml_parameter(sbml_model, id=param_id, value=0,
                                          units='mmol_per_gDW_per_hr')
            if bound == 'lower':
                wrap_libsbml(param.setValue, 0)
                wrap_libsbml(fbc_reaction_plugin.setLowerFluxBound, param_id)
            if bound == 'upper':
                wrap_libsbml(param.setValue, float('inf'))
                wrap_libsbml(fbc_reaction_plugin.setUpperFluxBound, param_id)
        return sbml_reaction


class Parameter(obj_model.Model):
    """ Parameter

    Attributes:
        id (:obj:`str`): unique identifier per model/submodel
        name (:obj:`str`): name
        model (:obj:`Model`): model
        type (:obj:`ParameterType`): parameter type
        value (:obj:`float`): value
        units (:obj:`str`): units of value
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        observable_expressions (:obj:`list` of :obj:`ObservableExpression`): observable expressions
        function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        stop_condition_expressions (:obj:`list` of :obj:`StopConditionExpression`): stop condition expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='parameters')
    type = EnumAttribute(ParameterType, default=ParameterType.other)
    value = FloatAttribute(min=0)
    units = StringAttribute()
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='parameters')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'type',
                           'value', 'units',
                           'comments', 'references')

    def add_to_sbml_doc(self, sbml_document):
        """ Add this Parameter to a libsbml SBML document.

        Args:
             sbml_document (:obj:`obj`): a `libsbml` SBMLDocument

        Returns:
            :obj:`libsbml.Parameter`: the libsbml Parameter that's created

        Raises:
            :obj:`LibSBMLError`: if calling `libsbml` raises an error
        """
        sbml_model = wrap_libsbml(sbml_document.getModel)
        # prefix id with 'parameter' so ids for wc_lang Parameters don't collide with ids for other libsbml parameters
        sbml_id = "parameter_{}".format(self.id)
        # TODO: use a standard unit ontology to map self.units to SBML model units
        if self.units == 'dimensionless':
            sbml_parameter = create_sbml_parameter(sbml_model, sbml_id, self.value, 'dimensionless_ud',
                                                   name=self.name)
        elif self.units == 's':
            sbml_parameter = create_sbml_parameter(sbml_model, sbml_id, self.value, 'second',
                                                   name=self.name)
        elif self.units == 'mmol/gDCW/h':
            sbml_parameter = create_sbml_parameter(sbml_model, sbml_id, self.value, 'mmol_per_gDW_per_hr',
                                                   name=self.name)
        else:
            sbml_parameter = create_sbml_parameter(sbml_model, sbml_id, self.value, 'dimensionless_ud',
                                                   name=self.name)

        return sbml_parameter


class Reference(obj_model.Model):
    """ Reference

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        title (:obj:`str`): title
        author (:obj:`str`): author(s)
        editor (:obj:`str`): editor(s)
        year (:obj:`int`): year
        type (:obj:`ReferenceType`): type
        publication (:obj:`str`): publication title
        publisher (:obj:`str`): publisher
        series (:obj:`str`): series
        volume (:obj:`str`): volume
        number (:obj:`str`): number
        issue (:obj:`str`): issue
        edition (:obj:`str`): edition
        chapter (:obj:`str`): chapter
        pages (:obj:`str`): page range
        comments (:obj:`str`): comments

    Related attributes:
        database_references (:obj:`list` of :obj:`DatabaseReference`): database references
        taxa (:obj:`list` of :obj:`Taxon`): taxa
        submodels (:obj:`list` of :obj:`Submodel`): submodels
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species types
        species (:obj:`list` of :obj:`Species`): species
        concentrations (:obj:`list` of :obj:`Concentration`): concentrations
        observables (:obj:`list` of :obj:`Observable`): observables
        functions (:obj:`list` of :obj:`Function`): functions
        dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        biomass_components (:obj:`list` of :obj:`BiomassComponent`): biomass components
        parameters (:obj:`list` of :obj:`Parameter`): parameters
        stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='references')
    title = StringAttribute()
    author = StringAttribute()
    editor = StringAttribute()
    year = PositiveIntegerAttribute()
    type = EnumAttribute(ReferenceType)
    publication = StringAttribute()
    publisher = StringAttribute()
    series = StringAttribute()
    volume = StringAttribute()
    number = StringAttribute()
    issue = StringAttribute()
    edition = StringAttribute()
    chapter = StringAttribute()
    pages = StringAttribute()
    comments = LongStringAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'title', 'author', 'editor', 'year', 'type', 'publication', 'publisher',
                           'series', 'volume', 'number', 'issue', 'edition', 'chapter', 'pages',
                           'comments')


class DatabaseReference(obj_model.Model):
    """ Reference to a source database entry

    Attributes:
        database (:obj:`str`): database name
        id (:obj:`str`): id of database entry
        url (:obj:`str`): URL of database entry
        model (:obj:`Model`): model
        taxon (:obj:`Taxon`): taxon
        submodel (:obj:`Submodel`): submodel
        species_type (:obj:`SpeciesType`): species type
        reaction (:obj:`Reaction`): reaction
        reference (:obj:`Reference`): reference
    """

    database = StringAttribute(min_length=1)
    id = StringAttribute(verbose_name='ID', min_length=1)
    url = UrlAttribute(verbose_name='URL')
    model = ManyToOneAttribute(Model, related_name='database_references')
    taxon = ManyToOneAttribute(Taxon, related_name='database_references')
    submodel = ManyToOneAttribute(Submodel, related_name='database_references')
    compartment = ManyToOneAttribute(Compartment, related_name='database_references')
    species_type = ManyToOneAttribute(SpeciesType, related_name='database_references')
    reaction = ManyToOneAttribute(Reaction, related_name='database_references')
    reference = ManyToOneAttribute(Reference, related_name='database_references')

    class Meta(obj_model.Model.Meta):
        unique_together = (('database', 'id', ), )
        attribute_order = ('database', 'id', 'url',
                           'model', 'taxon', 'submodel', 'compartment', 'species_type', 'reaction', 'reference')
        frozen_columns = 2
        ordering = ('database', 'id', )

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return '{}: {}'.format(self.database, self.id)

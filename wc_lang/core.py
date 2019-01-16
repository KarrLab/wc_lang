""" Data model to represent composite, multi-algorithmic biochemical models.

This module defines classes that represent the schema of a biochemical model:

* :obj:`Taxon`
* :obj:`Model`
* :obj:`Submodel`
* :obj:`DfbaObjective`
* :obj:`Compartment`
* :obj:`SpeciesType`
* :obj:`Species`
* :obj:`DistributionInitConcentration`
* :obj:`Reaction`
* :obj:`SpeciesCoefficient`
* :obj:`RateLaw`
* :obj:`RateLawExpression`
* :obj:`DfbaObjSpecies`
* :obj:`DfbaObjReaction`
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

from enum import Enum, EnumMeta
from math import ceil, floor, exp, log, log10, isnan
from natsort import natsorted, ns
from obj_model import (BooleanAttribute, EnumAttribute,
                       FloatAttribute,
                       IntegerAttribute, PositiveIntegerAttribute,
                       RegexAttribute, SlugAttribute, StringAttribute, LongStringAttribute, UrlAttribute,
                       DateTimeAttribute,
                       OneToOneAttribute, ManyToOneAttribute, ManyToManyAttribute, OneToManyAttribute,
                       InvalidObject, InvalidAttribute, TabularOrientation)
from obj_model.expression import (ExpressionOneToOneAttribute, ExpressionManyToOneAttribute,
                                  ExpressionStaticTermMeta, ExpressionDynamicTermMeta,
                                  ExpressionExpressionTermMeta, Expression,
                                  ParsedExpression, ParsedExpressionError)
from obj_model.ontology import OntologyAttribute
from six import with_metaclass
from wc_lang.sbml.util import (wrap_libsbml, str_to_xmlstr, LibSBMLError,
                               create_sbml_parameter)
from wc_utils.util.chem import EmpiricalFormula
from wc_utils.util.enumerate import CaseInsensitiveEnum, CaseInsensitiveEnumMeta
from wc_utils.util.list import det_dedupe
from wc_utils.util.ontology import wcm_ontology
from wc_utils.util.units import unit_registry
import collections
import datetime
import networkx
import obj_model
import obj_model.chem
import pkg_resources
import pronto.term
import re
import six
import stringcase
import token

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


class Unit(int, Enum):
    """ Base class for units """
    pass


class TimeUnit(Unit):
    """ Time units """
    s = 1


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


class TemperatureUnit(Unit):
    """ Temperature units """
    C = 1


class PhUnit(Unit):
    """ pH units """
    dimensionless = 1


class MassUnit(Unit):
    """ Mass units """
    g = 1


VolumeUnit = Unit('VolumeUnit', names=[
    ('l', 1),
    # ('dm^2', 2),
])


DensityUnit = Unit('DensityUnit', names=[
    ('g l^-1', 1),
    # ('g dm^-2', 2),
])


class ObservableCoefficientUnit(Unit):
    """ Observable coefficient units """
    dimensionless = 1


class MoleculeCountUnit(Unit):
    """ Units of molecule counts """
    molecule = 1


ConcentrationUnit = Unit('ConcentrationUnit', names=[
    ('molecule', 1),
    ('M', 2),
    ('mM', 3),
    ('uM', 4),
    ('nM', 5),
    ('pM', 6),
    ('fM', 7),
    ('aM', 8),
    # ('mol dm^-2', 9),
])
ConcentrationUnit.Meta = {
    ConcentrationUnit['molecule']: {
        'xml_id': 'molecule',
        'substance_units': {'kind': 'item', 'exponent': 1, 'scale': 0},
        'volume_units': {'kind': 'item', 'exponent': 1, 'scale': 0},
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
    # ConcentrationUnit['mol dm^-2']: {
    #    'xml_id': 'mol_per_dm_2',
    #    'substance_units': {'kind': 'mole', 'exponent': 1, 'scale': 0},
    #    'volume_units': {'kind': 'metre', 'exponent': -2, 'scale': -1},
    #},
}


class ReactionParticipantUnit(Unit):
    """ Units of reaction participants """
    dimensionless = 1


class RateLawDirection(int, CaseInsensitiveEnum):
    """ Rate law directions """
    backward = -1
    forward = 1


ReactionRateUnit = Unit('ReactionRateUnit', names=[
    ('s^-1', 1),
])

ReactionFluxBoundUnit = Unit('ReactionFluxBoundUnit', names=[
    ('M s^-1', 1),
])

DfbaObjectiveUnit = Unit('DfbaObjectiveUnit', names=[
    ('dimensionless', 1),
])


class DfbaCellSizeUnit(Unit):
    """ dFBA cell size units """
    l = 1
    gDCW = 2


DfbaObjectiveCoefficientUnit = Unit('DfbaObjectiveCoefficientUnit', names=[
    ('s', 1),
])

DfbaObjSpeciesUnit = Unit('DfbaObjSpeciesUnit', names=[
    ('M s^-1', 1),
    ('mol gDCW^-1 s^-1', 2),
])


class StopConditionUnit(Unit):
    """ Stop condition units """
    dimensionless = 1


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
                species_id = Species.gen_id_static(species_type.id, compartment.id)
                species, error = Species.deserialize(species_id, objects)
                if error:
                    errors.extend(error.messages)

                elif coefficient != 0:
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

        # check element and charge balance
        validate_element_charge_balance = wc_lang.config.core.get_config()[
            'wc_lang']['validation']['validate_element_charge_balance']
        if validate_element_charge_balance:
            delta_formula = EmpiricalFormula()
            delta_charge = 0.
            errors = []

            for part in value:
                if part.species.species_type.empirical_formula:
                    delta_formula += part.species.species_type.empirical_formula * part.coefficient

                if part.species.species_type.charge is None:
                    errors.append('Charge must be defined for {}'.format(part.species.species_type.id))
                else:
                    delta_charge += part.species.species_type.charge * part.coefficient

            if not errors:
                if delta_formula:
                    errors.append('Reaction is element imbalanced: {}'.format(delta_formula))
                if delta_charge != 0.:
                    errors.append('Reaction is charge imbalanced: {}'.format(delta_charge))

            if errors:
                return InvalidAttribute(self, errors)

        # return None
        return None


class DatabaseReferenceOneToManyAttribute(OneToManyAttribute):
    def __init__(self, related_name='', verbose_name='Database references', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(DatabaseReferenceOneToManyAttribute, self).__init__('DatabaseReference', related_name=related_name,
                                                                  verbose_name=verbose_name,
                                                                  verbose_related_name=verbose_related_name,
                                                                  help=help)

    def serialize(self, db_refs, encoded=None):
        """ Serialize related object

        Args:
            db_refs (:obj:`list` of :obj:`DatabaseReference`): Python representation of database references
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: string representation
        """
        sorted_db_refs = sorted(db_refs, key=lambda db_ref: (db_ref.database, db_ref.id))
        return ', '.join(db_ref.serialize() for db_ref in sorted_db_refs)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of `list` of `DatabaseReference`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        value = value or ''
        value = value.strip()
        if not value:
            return ([], None)

        db_refs = set()
        errors = []
        for val in value.split(','):
            db_ref, invalid = DatabaseReference.deserialize(val, objects)
            if invalid:
                errors.extend(invalid.messages)
            else:
                db_refs.add(db_ref)

        if errors:
            return (None, InvalidAttribute(self, errors))
        else:
            return (det_dedupe(db_refs), None)


class DatabaseReferenceManyToManyAttribute(ManyToManyAttribute):
    def __init__(self, related_name='', verbose_name='Database references', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(DatabaseReferenceManyToManyAttribute, self).__init__('DatabaseReference', related_name=related_name,
                                                                   verbose_name=verbose_name,
                                                                   verbose_related_name=verbose_related_name,
                                                                   help=help)

    def serialize(self, db_refs, encoded=None):
        """ Serialize related object

        Args:
            db_refs (:obj:`list` of :obj:`DatabaseReference`): Python representation of database references
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: string representation
        """
        sorted_db_refs = sorted(db_refs, key=lambda db_ref: (db_ref.database, db_ref.id))
        return ', '.join(db_ref.serialize() for db_ref in sorted_db_refs)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`tuple` of `list` of `DatabaseReference`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        value = value or ''
        value = value.strip()
        if not value:
            return ([], None)

        db_refs = set()
        errors = []
        for val in value.split(','):
            db_ref, invalid = DatabaseReference.deserialize(val, objects)
            if invalid:
                errors.extend(invalid.messages)
            else:
                db_refs.add(db_ref)

        if errors:
            return (None, InvalidAttribute(self, errors))
        else:
            return (det_dedupe(db_refs), None)


class CommentAttribute(obj_model.LongStringAttribute):
    """ Comment attribute """

    SEPARATOR = '\n\n'

    def merge(self, left, right, right_objs_in_left, left_objs_in_right):
        """ Merge an attribute of elements of two models

        Args:
            left (:obj:`Model`): an element in a model to merge
            right (:obj:`Model`): an element in a second model to merge
            right_objs_in_left (:obj:`dict`): mapping from objects in right model to objects in left model
            left_objs_in_right (:obj:`dict`): mapping from objects in left model to objects in right model

        Raises:
            :obj:`ValueError`: if the attributes of the elements of the models are different
        """
        left_val = getattr(left, self.name)
        right_val = getattr(right, self.name)
        if left_val != right_val:
            setattr(left, self.name, left_val + self.SEPARATOR + right_val)


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
        time_units (:obj:`TimeUnit`): time units
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        comments (:obj:`str`): comments
        created (:obj:`datetime`): date created
        updated (:obj:`datetime`): date updated

    Related attributes:
        taxon (:obj:`Taxon`): taxon
        env (:obj:`Environment`): environment
        submodels (:obj:`list` of :obj:`Submodel`): submodels
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species types
        species (:obj:`list` of :obj:`Species`): species
        distribution_init_concentrations (:obj:`list` of :obj:`DistributionInitConcentration`):
            distributions of initial concentrations of species at the beginning of
            each cell cycle
        observables (:obj:`list` of :obj:`Observable`): observables
        functions (:obj:`list` of :obj:`Function`): functions
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        dfba_obj_reactions (:obj:`list` of :obj:`DfbaObjReaction`): dFBA objective reactions
        dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
        stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
        parameters (:obj:`list` of :obj:`Parameter`): parameters
        evidences (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = SlugAttribute()
    name = StringAttribute()
    version = RegexAttribute(min_length=1, pattern=r'^[0-9]+\.[0-9+]\.[0-9]+', flags=re.I)
    url = UrlAttribute(verbose_name='URL')
    branch = StringAttribute()
    revision = StringAttribute()
    wc_lang_version = RegexAttribute(min_length=1, pattern=r'^[0-9]+\.[0-9+]\.[0-9]+', flags=re.I,
                                     default=wc_lang_version, verbose_name='wc_lang version')
    author = LongStringAttribute()
    author_organization = LongStringAttribute()
    author_email = LongStringAttribute()
    time_units = EnumAttribute(TimeUnit, default=TimeUnit.s)
    db_refs = DatabaseReferenceOneToManyAttribute(related_name='model')
    comments = CommentAttribute()
    created = DateTimeAttribute()
    updated = DateTimeAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'version',
                           'url', 'branch', 'revision',
                           'wc_lang_version',
                           'author', 'author_organization', 'author_email',
                           'time_units',
                           'db_refs', 'comments',
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

        * Network of compartments is rooted and acyclic
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

        # Network of compartments is rooted and acyclic
        digraph = networkx.DiGraph()
        for comp in self.compartments:
            for sub_comp in comp.sub_compartments:
                digraph.add_edge(comp.id, sub_comp.id)
        if list(networkx.simple_cycles(digraph)):
            errors.append(InvalidAttribute(self.Meta.related_attributes['compartments'],
                                           ['Compartment parent/child relations cannot be cyclic']))

        # Networks of observables and functions are acyclic
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
            expression_type = model_type.Meta.expression_term_model
            for self_ref_attr_name, self_ref_attr in expression_type.Meta.attributes.items():
                if isinstance(self_ref_attr, obj_model.RelatedAttribute) and self_ref_attr.related_class == model_type:
                    break

            # find cyclic dependencies
            digraph = networkx.DiGraph()
            for model in models:
                if model.expression:
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

    def get_root_compartments(self, __type=None, **kwargs):
        """ Get the root compartments (compartments that either have no parent compartment or whose
        parent compartment is not cellular)

        Returns:
            :obj:`list` of :obj:`Compartment`: root compartments
        """
        roots = []
        for comp in self.get_compartments(__type=__type, **kwargs):
            if comp.parent_compartment is None \
                    or comp.parent_compartment.biological_type != wcm_ontology['WCM:cellular_compartment']:
                roots.append(comp)
        return roots

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

    def get_distribution_init_concentrations(self, __type=None, **kwargs):
        """ Get all initial distributions of concentrations of species at the
        beginning of each cell cycle

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`DistributionInitConcentration`: initial distributions
                of concentrations of species at the beginning of each cell cycle
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.distribution_init_concentrations.get(__type=__type, **kwargs)

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

    def get_dfba_objs(self, __type=None, **kwargs):
        """ Get all dFBA objectives

        Returns:
            :obj:`list` of :obj:`DfbaObjective`: dFBA objectives
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.dfba_objs.get(__type=__type, **kwargs)

    def get_dfba_obj_reactions(self, __type=None, **kwargs):
        """ Get all dFBA objective reactions used by submodels

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            **kwargs (:obj:`dict` of `str`:`object`): dictionary of attribute name/value pairs to find matching
                objects

        Returns:
            :obj:`list` of :obj:`DfbaObjReaction`: dFBA objective reactions
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.dfba_obj_reactions.get(__type=__type, **kwargs)

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

    def get_evidence(self, __type=None, **kwargs):
        """ Get all evidence for model

        Returns:
            :obj:`list` of :obj:`Evidence`: evidence for model
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.evidences.get(__type=__type, **kwargs)

    def get_interpretations(self, __type=None, **kwargs):
        """ Get all interpretations for model

        Returns:
            :obj:`list` of :obj:`Interpretation`: interpretations for model
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.interpretations.get(__type=__type, **kwargs)

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
            type_names = [stringcase.snakecase(__type.__name__) + 's']
        else:
            type_names = [
                'submodels', 'compartments', 'species_types', 'species',
                'distribution_init_concentrations', 'observables', 'functions',
                'dfba_objs', 'reactions', 'rate_laws', 'dfba_obj_reactions',
                'stop_conditions', 'parameters', 'evidence', 'interpretations',
                'references',
            ]

        components = []
        for type_name in type_names:
            get_func = getattr(self, 'get_' + type_name)
            components.extend(get_func(__type=__type, **kwargs))

        return components

    def merge_attrs(self, other, other_objs_in_self, self_objs_in_other):
        """ Merge attributes of two objects

        Args:
            other (:obj:`Model`): other model
            other_objs_in_self (:obj:`dict`): dictionary that maps instances of objects in another model to objects
                in a model
            self_objs_in_other (:obj:`dict`): dictionary that maps instances of objects in a model to objects
                in another model
        """
        self.Meta.attributes['comments'].merge(self, other, other_objs_in_self, self_objs_in_other)
        for attr in self.Meta.attributes.values():
            if isinstance(attr, obj_model.RelatedAttribute):
                attr.merge(self, other, other_objs_in_self, self_objs_in_other)


class Taxon(obj_model.Model):
    """ Biological taxon (e.g. family, genus, species, strain, etc.)

    Attributes:
        id (:obj:`str`): unique identifier equal to 'taxon'
        name (:obj:`str`): name
        model (:obj:`Model`): model
        rank (:obj:`TaxonRank`): rank
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = RegexAttribute(pattern=r'^taxon$', primary=True, unique=True)
    name = StringAttribute()
    model = OneToOneAttribute(Model, related_name='taxon')
    rank = EnumAttribute(TaxonRank, default=TaxonRank.species)
    db_refs = DatabaseReferenceOneToManyAttribute(related_name='taxon')
    comments = CommentAttribute()
    references = OneToManyAttribute('Reference', related_name='taxon')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'rank',
                           'db_refs', 'comments', 'references')
        tabular_orientation = TabularOrientation.column


class Environment(obj_model.Model):
    """ Environment

    Attributes:
        id (:obj:`str`): unique identifier equal to 'env'
        name (:obj:`str`): name
        model (:obj:`Model`): model
        temp (:obj:`float`): temperature
        temp_units (:obj:`TemperatureUnit`): temperature units
        ph (:obj:`float`): pH
        ph_units (:obj:`PhUnit`): pH units
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = RegexAttribute(pattern=r'^env$', primary=True, unique=True)
    name = StringAttribute()
    model = OneToOneAttribute(Model, related_name='env')
    temp = FloatAttribute(verbose_name='Temperature')
    temp_units = EnumAttribute(TemperatureUnit, default=TemperatureUnit.C, verbose_name='Temperature units')
    ph = FloatAttribute(verbose_name='pH')
    ph_units = EnumAttribute(PhUnit, default=PhUnit.dimensionless, verbose_name='pH units')
    db_refs = DatabaseReferenceOneToManyAttribute(related_name='env')
    comments = CommentAttribute()
    references = OneToManyAttribute('Reference', related_name='env')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'temp', 'temp_units', 'ph', 'ph_units',
                           'db_refs', 'comments', 'references')
        tabular_orientation = TabularOrientation.column


class Submodel(obj_model.Model):
    """ Submodel

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        framework (:obj:`pronto.term.Term`): modeling framework (e.g. dynamic flux balance analysis)
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        dfba_obj (:obj:`DfbaObjective`): objective function for a dFBA submodel;
            if not initialized, then `dfba_obj_reaction` is used as the objective function
        dfba_obj_reactions (:obj:`list` of :obj:`DfbaObjReaction`): the growth reaction for a dFBA submodel
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='submodels')
    framework = OntologyAttribute(wcm_ontology,
                                  namespace='WCM',
                                  terms=wcm_ontology['WCM:modeling_framework'].rchildren(),
                                  default=wcm_ontology['WCM:stochastic_simulation_algorithm'],
                                  none=False)
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='submodels')
    evidence = ManyToManyAttribute('Evidence', related_name='submodels')
    interpretations = ManyToManyAttribute('Interpretation', related_name='submodels')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='submodels')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'framework',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        indexed_attrs_tuples = (('id',), )
        merge = obj_model.ModelMerge.append

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

        if self.framework == wcm_ontology['WCM:dynamic_flux_balance_analysis']:
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

    def get_species_types(self):
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
        for dfba_obj_reaction in self.get_dfba_obj_reactions():
            for dfba_obj_species in dfba_obj_reaction.dfba_obj_species:
                species.append(dfba_obj_species.species)
        for observable in self.get_observables():
            species.extend(observable.expression.species)
        for function in self.get_functions():
            species.extend(function.expression.species)
        for rate_law in self.get_rate_laws():
            species.extend(rate_law.expression.species)
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

    def get_reactions(self):
        """ Get reactions in submodel

        Returns:
            :obj:`list` of :obj:`Reaction`: reactions in submodel
        """
        reactions = list(self.reactions)
        for dfba_obj in self.get_dfba_objs():
            if dfba_obj.expression:
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

    def get_dfba_objs(self):
        """ Get dFBA objectives in submodel

        Returns:
            :obj:`list` of :obj:`DfbaObjective`: dFBA objectives in submodel
        """
        if self.dfba_obj:
            return [self.dfba_obj]
        return []

    def get_dfba_obj_reactions(self):
        """ Get dFBA objective reactions in submodel

        Returns:
            :obj:`list` of :obj:`DfbaObjReaction`: dFBA objective reactions in submodel
        """
        rxns = list(self.dfba_obj_reactions)
        for dfba_obj in self.get_dfba_objs():
            if dfba_obj.expression:
                rxns.extend(dfba_obj.expression.dfba_obj_reactions)
        return det_dedupe(rxns)

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

    def get_evidence(self):
        """ Get evidence of submodel

        Returns:
            :obj:`list` of :obj:`Evidence`: evidence for submodel
        """
        types = [
            'compartments',
            'species_types',
            'species',
            'observables',
            'functions',
            'dfba_objs',
            'reactions',
            'rate_laws',
            'dfba_obj_reactions',
            'parameters',
            'interpretations',
        ]
        evidence = list(self.evidence)
        for type in types:
            get_func = getattr(self, 'get_' + type)
            for obj in get_func():
                evidence.extend(obj.evidence)
        return det_dedupe(evidence)

    def get_interpretations(self):
        """ Get interpretations of submodel

        Returns:
            :obj:`list` of :obj:`Interpretation`: interpretations for submodel
        """
        types = [
            'compartments',
            'species_types',
            'species',
            'observables',
            'functions',
            'dfba_objs',
            'reactions',
            'rate_laws',
            'dfba_obj_reactions',
            'parameters',
        ]
        interpretations = list(self.interpretations)
        for type in types:
            get_func = getattr(self, 'get_' + type)
            for obj in get_func():
                interpretations.extend(obj.interpretations)
        return det_dedupe(interpretations)

    def get_references(self):
        """ Get references of submodel

        Returns:
            :obj:`list` of :obj:`Reference`: references in submodel
        """
        types = [
            'compartments',
            'species_types',
            'species',
            'observables',
            'functions',
            'dfba_objs',
            'reactions',
            'rate_laws',
            'dfba_obj_reactions',
            'parameters',
            'evidence',
            'interpretations',
        ]
        references = list(self.references)
        for type in types:
            get_func = getattr(self, 'get_' + type)
            for obj in get_func():
                references.extend(obj.references)
        return det_dedupe(references)

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
            self.get_reactions() + \
            self.get_rate_laws() + \
            self.get_dfba_objs() + \
            self.get_dfba_obj_reactions() + \
            self.get_parameters() + \
            self.get_evidence() + \
            self.get_interpretations() + \
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


class DfbaObjectiveExpression(obj_model.Model, Expression):
    """ A mathematical expression of Reactions and DfbaObjReactions

    The expression used by a :obj:`DfbaObjective`.

    Attributes:
        expression (:obj:`str`): mathematical expression
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        reactions (:obj:`list` of :obj:`Reaction`): reactions used by this expression
        dfba_obj_reactions (:obj:`list` of :obj:`Species`): dFBA objective reactions used by this expression

    Related attributes:
        dfba_obj (:obj:`DfbaObjective`): dFBA objective
    """

    expression = LongStringAttribute(primary=True, unique=True, default='')
    reactions = OneToManyAttribute('Reaction', related_name='dfba_obj_expression',
                                   verbose_related_name='dFBA objective expression')
    dfba_obj_reactions = OneToManyAttribute('DfbaObjReaction', related_name='dfba_obj_expression',
                                            verbose_name='dFBA objective reactions', verbose_related_name='dFBA objective expression')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        tabular_orientation = TabularOrientation.inline
        expression_valid_functions = ()
        expression_term_models = ('Reaction', 'DfbaObjReaction')
        verbose_name = 'dFBA objective expression'
        merge = obj_model.ModelMerge.append

    def validate(self):
        """ Determine if the dFBA objective expression is valid

        * Check that the expression is a linear function
        * Check if expression is a function of at least one reaction or dFBA objective reaction
        * Check that the reactions and dFBA objective reactions belong to the same submodel

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
        if not self.reactions and not self.dfba_obj_reactions:
            expr_errors.append('Expression must be a function of at least one reaction or dFBA objective reaction')

        if self.dfba_obj and self.dfba_obj.submodel:
            missing_rxns = set(self.reactions).difference(set(self.dfba_obj.submodel.reactions))
            missing_dfba_obj_rxns = set(self.dfba_obj_reactions).difference(set(self.dfba_obj.submodel.dfba_obj_reactions))

            if missing_rxns:
                expr_errors.append(('dFBA submodel {} must contain the following reactions '
                                    'that are in its objective function:\n  {}').format(
                    self.dfba_obj.submodel.id,
                    '\n  '.join(rxn.id for rxn in missing_rxns)))
            if missing_dfba_obj_rxns:
                expr_errors.append(('dFBA submodel {} must contain the following dFBA objective reactions '
                                    'that are in its objective function:\n  {}').format(
                    self.dfba_obj.submodel.id,
                    '\n  '.join(rxn.id for rxn in missing_dfba_obj_rxns)))

        if expr_errors:
            errors.append(InvalidAttribute(self.Meta.attributes['expression'], expr_errors))

        if errors:
            return InvalidObject(self, errors)
        return Expression.validate(self, self.dfba_obj)

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return Expression.serialize(self)

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
        return Expression.deserialize(cls, value, objects)


class DfbaObjective(obj_model.Model):
    """ dFBA objective function

    Attributes:
        id (:obj:`str`): identifier equal to `dfba-obj-{submodel.id}`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        submodel (:obj:`Submodel`): the `Submodel` which uses this `DfbaObjective`
        expression (:obj:`DfbaObjectiveExpression`): mathematical expression of the objective function
        units (:obj:`DfbaObjectiveUnit`): units
        reaction_rate_units (:obj:`ReactionRateUnit`): reaction rate units
        coefficient_units (:obj:`DfbaObjectiveCoefficientUnit`): coefficient units
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='dfba_objs', verbose_related_name='dFBA objectives')
    submodel = OneToOneAttribute(Submodel, related_name='dfba_obj', min_related=1, verbose_related_name='dFBA objective')
    expression = ExpressionOneToOneAttribute(DfbaObjectiveExpression, related_name='dfba_obj',
                                             min_related=1, min_related_rev=1, verbose_related_name='dFBA objective')
    units = EnumAttribute(DfbaObjectiveUnit, default=DfbaObjectiveUnit['dimensionless'])
    reaction_rate_units = EnumAttribute(ReactionRateUnit, default=ReactionRateUnit['s^-1'])
    coefficient_units = EnumAttribute(DfbaObjectiveCoefficientUnit, default=DfbaObjectiveCoefficientUnit['s'])
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='dfba_objs', verbose_related_name='dFBA objectives')
    evidence = ManyToManyAttribute('Evidence', related_name='dfba_objs', verbose_related_name='dFBA objectives')
    interpretations = ManyToManyAttribute('Interpretation', related_name='dfba_objs', verbose_related_name='dFBA objectives')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='dfba_objs', verbose_related_name='dFBA objectives')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        verbose_name = 'dFBA objective'
        attribute_order = ('id', 'name', 'submodel', 'expression', 'units', 'reaction_rate_units', 'coefficient_units',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        expression_term_model = DfbaObjectiveExpression
        expression_term_units = 'units'
        merge = obj_model.ModelMerge.append

    def gen_id(self):
        """ Generate identifier

        Returns:
            :obj:`str`: identifier
        """
        return 'dfba-obj-{}'.format(self.submodel.id)

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

        if self.submodel and self.id != self.gen_id():
            errors.append(InvalidAttribute(self.Meta.attributes['id'],
                                           ['Id must be {}'.format(self.gen_id())]))
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
        for idx, dfba_obj_reaction in enumerate(self.expression.dfba_obj_reactions):
            sbml_flux_objective = wrap_libsbml(sbml_objective.createFluxObjective)
            wrap_libsbml(sbml_flux_objective.setReaction, dfba_obj_reaction.id)
            wrap_libsbml(sbml_flux_objective.setCoefficient,
                         self.expression._parsed_expression.lin_coeffs[DfbaObjReaction][dfba_obj_reaction])

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

        # products of reactions
        for reaction in self.expression.reactions:
            if reaction.reversible:
                for part in reaction.participants:
                    if part.species.has_attr_vals(__type=__type, **kwargs):
                        products.append(part.species)
            else:
                for part in reaction.participants:
                    if part.coefficient > 0:
                        if part.species.has_attr_vals(__type=__type, **kwargs):
                            products.append(part.species)

        # products of dFBA objective reactions
        for dfba_obj_reaction in self.expression.dfba_obj_reactions:
            for dfba_obj_species in dfba_obj_reaction.dfba_obj_species:
                if dfba_obj_species.value > 0 and dfba_obj_species.species.has_attr_vals(__type=__type, **kwargs):
                    products.append(dfba_obj_species.species)

        # return unique list
        return det_dedupe(products)


class Compartment(obj_model.Model):
    """ Compartment

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        biological_type (:obj:`pronto.term.Term`): biological type
        physical_type (:obj:`pronto.term.Term`): physical type
        geometry (:obj:`pronto.term.Term`): geometry
        parent_compartment (:obj:`Compartment`): parent compartment
        mass_units (:obj:`MassUnit`): mass units
        mean_init_volume (:obj:`float`): mean initial volume
        std_init_volume (:obj:`float`): standard  deviation of the mean initial volume
        init_volume_units (:obj:`VolumeUnit`): units of volume
        init_density (:obj:`Parameter`): function that calculates the density during the initialization of
            each simulation
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        sub_compartments (:obj:`list` of :obj:`Compartment`): compartments contained in this compartment
        species (:obj:`list` of :obj:`Species`): species in this compartment
        function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        stop_condition_expressions (:obj:`list` of :obj:`StopConditionExpression`): stop condition expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='compartments')
    biological_type = OntologyAttribute(wcm_ontology,
                                        namespace='WCM',
                                        terms=wcm_ontology['WCM:biological_compartment'].rchildren(),
                                        default=wcm_ontology['WCM:cellular_compartment'],
                                        none=True)
    physical_type = OntologyAttribute(wcm_ontology,
                                      namespace='WCM',
                                      terms=wcm_ontology['WCM:physical_compartment'].rchildren(),
                                      default=wcm_ontology['WCM:fluid_compartment'],
                                      none=True)
    geometry = OntologyAttribute(wcm_ontology,
                                 namespace='WCM',
                                 terms=wcm_ontology['WCM:geometric_compartment'].rchildren(),
                                 default=wcm_ontology['WCM:3D_compartment'],
                                 none=True)
    parent_compartment = ManyToOneAttribute('Compartment', related_name='sub_compartments')
    distribution_init_volume = OntologyAttribute(wcm_ontology,
                                                 namespace='WCM',
                                                 terms=wcm_ontology['WCM:random_distribution'].rchildren(),
                                                 default=wcm_ontology['WCM:normal_distribution'],
                                                 verbose_name='Initial volume distribution')
    mass_units = EnumAttribute(MassUnit, default=MassUnit.g)
    mean_init_volume = FloatAttribute(min=0, verbose_name='Initial volume mean')
    std_init_volume = FloatAttribute(min=0, verbose_name='Initial volume standard deviation')
    init_volume_units = EnumAttribute(VolumeUnit, default=VolumeUnit.l, verbose_name='Initial volume units')
    init_density = OneToOneAttribute('Parameter', related_name='density_compartment', verbose_name='Initial density')
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='compartments')
    evidence = ManyToManyAttribute('Evidence', related_name='compartments')
    interpretations = ManyToManyAttribute('Interpretation', related_name='compartments')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='compartments')

    class Meta(obj_model.Model.Meta, ExpressionDynamicTermMeta):
        attribute_order = ('id', 'name',
                           'biological_type', 'physical_type', 'geometry', 'parent_compartment',
                           'mass_units',
                           'distribution_init_volume', 'mean_init_volume', 'std_init_volume', 'init_volume_units',
                           'init_density',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        expression_term_units = 'mass_units'

    def validate(self):
        """ Check that the compartment is valid

        * Check that the units of density are `g l^-1`

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(Compartment, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        if self.geometry == wcm_ontology['WCM:3D_compartment']:
            if not self.init_density:
                errors.append(InvalidAttribute(self.Meta.attributes['init_density'],
                                               ['Initial density must be defined for 3D compartments']))
            elif unit_registry.parse_expression(self.init_density.units).to_base_units().units != \
                    unit_registry.parse_expression(DensityUnit['g l^-1'].name).to_base_units().units:
                errors.append(InvalidAttribute(self.Meta.attributes['init_density'],
                                               ['Initial density of 3D compartment must have units `{}`'.format(
                                                DensityUnit['g l^-1'].name)]))

        if errors:
            return InvalidObject(self, errors)
        return None

    def get_sub_compartments(self, nested=False):
        """ Get sub-compartments, optionally including all nested sub-compartments

        Returns:
            :obj:`list` of :obj:`Compartment`: list of sub-compartments
        """
        if nested and self.sub_compartments:
            nested_compartments = list(self.sub_compartments)
            to_expand = list(self.sub_compartments)
            while to_expand:
                comp = to_expand.pop()
                for sub_comp in comp.sub_compartments:
                    if sub_comp not in nested_compartments:
                        nested_compartments.append(sub_comp)
                        to_expand.append(sub_comp)
            return nested_compartments
        else:
            return self.sub_compartments

    def get_tot_mean_init_volume(self):
        """ Get total mean initial volume of compartment and nested
        sub-compartments

        Returns:
            :obj:`float`: total mean initial volume of compartment and nested
                sub-compartments
        """
        tot = self.mean_init_volume
        for comp in self.get_sub_compartments(nested=True):
            tot += comp.mean_init_volume
        return tot

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
        wrap_libsbml(sbml_compartment.setSize, 1.0)  # todo: set based on calculated concentrations and density
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
        empirical_formula (:obj:`EmpiricalFormula`): empirical formula
        molecular_weight (:obj:`float`): molecular weight
        charge (:obj:`int`): charge
        type (:obj:`pronto.term.Term`): type
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        species (:obj:`list` of :obj:`Species`): species
        distribution_init_concentrations (:obj:`list` of :obj:`DistributionInitConcentration`):
            distribution of initial concentrations of species at the beginning of
            each cell cycle
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='species_types')
    structure = LongStringAttribute()
    empirical_formula = obj_model.chem.EmpiricalFormulaAttribute()
    molecular_weight = FloatAttribute(min=0)
    charge = IntegerAttribute()
    type = OntologyAttribute(wcm_ontology,
                             namespace='WCM',
                             terms=wcm_ontology['WCM:species_type'].rchildren(),
                             default=wcm_ontology['WCM:metabolite'],
                             none=True)
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='species_types', verbose_related_name='species types')
    evidence = ManyToManyAttribute('Evidence', related_name='species_types')
    interpretations = ManyToManyAttribute('Interpretation', related_name='species_types')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='species_types')

    class Meta(obj_model.Model.Meta):
        verbose_name = 'Species type'
        attribute_order = ('id', 'name', 'structure', 'empirical_formula',
                           'molecular_weight', 'charge', 'type',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        indexed_attrs_tuples = (('id',), )

    def has_carbon(self):
        """ Returns `True` is species contains at least one carbon atom.

        Returns:
            :obj:`bool`: `True` is species contains at least one carbon atom.
        """
        return self.empirical_formula and self.empirical_formula['C'] > 0


class Species(obj_model.Model):
    """ Species (tuple of species type, compartment)

    Attributes:
        id (:obj:`str`): identifier equal to `{species_type.id}[{compartment.id}]`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        species_type (:obj:`SpeciesType`): species type
        compartment (:obj:`Compartment`): compartment
        units (:obj:`MoleculeCountUnit`): units of counts
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        distribution_init_concentration (:obj:`DistributionInitConcentration`):
            distribution of initial concentration
        species_coefficients (:obj:`list` of :obj:`SpeciesCoefficient`): participations in reactions and observables
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        observable_expressions (:obj:`list` of :obj:`ObservableExpression`): observable expressions
        stop_condition_expressions (:obj:`list` of :obj:`StopConditionExpression`): stop condition expressions
        function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='species')
    species_type = ManyToOneAttribute(SpeciesType, related_name='species', min_related=1)
    compartment = ManyToOneAttribute(Compartment, related_name='species', min_related=1)
    units = EnumAttribute(MoleculeCountUnit, default=MoleculeCountUnit['molecule'])
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='species')
    evidence = ManyToManyAttribute('Evidence', related_name='species')
    interpretations = ManyToManyAttribute('Interpretation', related_name='species')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='species')

    class Meta(obj_model.Model.Meta, ExpressionDynamicTermMeta):
        attribute_order = ('id', 'name', 'species_type', 'compartment', 'units',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        frozen_columns = 1
        # unique_together = (('species_type', 'compartment', ), )
        indexed_attrs_tuples = (('species_type', 'compartment'), )
        expression_term_token_pattern = (token.NAME, token.LSQB, token.NAME, token.RSQB)
        expression_term_units = 'units'

    def gen_id(self):
        """ Generate identifier

        Returns:
            :obj:`str`: identifier
        """
        return self.gen_id_static(self.species_type.id, self.compartment.id)

    @staticmethod
    def gen_id_static(species_type_id, compartment_id):
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

        if self.id != self.gen_id():
            errors.append(InvalidAttribute(self.Meta.attributes['id'],
                                           ['Id must be {}'.format(self.gen_id())]))

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
        # TODO: this costs O(|ids||species_iterator|); replace with O(|ids|) operation
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

        # set the initial concentration
        wrap_libsbml(sbml_species.setInitialConcentration, self.distribution_init_concentration.mean)

        # set units
        unit_xml_id = ConcentrationUnit.Meta[self.distribution_init_concentration.units]['xml_id']
        wrap_libsbml(sbml_species.setSubstanceUnits, unit_xml_id)

        return sbml_species


class DistributionInitConcentration(obj_model.Model):
    """ Distribution of the initial concentration of a species
    at the beginning of each cell cycle

    Attributes:
        id (:obj:`str`): identifier equal to `dist-init-conc-{species.id}`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        species (:obj:`Species`): species
        distribution (:obj:`pronto.term.Term`): distribution
        mean (:obj:`float`): mean concentration in a population of single cells at the
            beginning of each cell cycle
        std (:obj:`float`): standard deviation of the concentration in a population of
            single cells at the beginning of each cell cycle
        units (:obj:`ConcentrationUnit`): units; default units is `M`
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='distribution_init_concentrations')
    species = OneToOneAttribute(Species, min_related=1, related_name='distribution_init_concentration')
    distribution = OntologyAttribute(wcm_ontology,
                                     namespace='WCM',
                                     terms=wcm_ontology['WCM:random_distribution'].rchildren(),
                                     default=wcm_ontology['WCM:normal_distribution'])
    mean = FloatAttribute(min=0)
    std = FloatAttribute(min=0, verbose_name='Standard deviation')
    units = EnumAttribute(ConcentrationUnit, default=ConcentrationUnit.M)
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='distribution_init_concentrations')
    evidence = ManyToManyAttribute('Evidence', related_name='distribution_init_concentrations')
    interpretations = ManyToManyAttribute('Interpretation', related_name='distribution_init_concentrations')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='distribution_init_concentrations')

    class Meta(obj_model.Model.Meta):
        # unique_together = (('species', ), )
        attribute_order = ('id', 'name', 'species',
                           'distribution', 'mean', 'std', 'units',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        verbose_name = 'Initial species concentration'
        frozen_columns = 1

    def gen_id(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return 'dist-init-conc-{}'.format(self.species.id)

    def validate(self):
        """ Check that the distribution of initial concentrations
        at the beginning of each cell cycle is valid

        * Validate that identifier is equal to `dist-init-conc-{species.id}]`

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(DistributionInitConcentration, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        if self.species and self.id != self.gen_id():
            errors.append(InvalidAttribute(self.Meta.attributes['id'],
                                           ['Id must be {}'.format(self.gen_id())]))

        if errors:
            return InvalidObject(self, errors)
        return None


class ObservableExpression(obj_model.Model, Expression):
    """ A mathematical expression of Observables and Species

    The expression used by a `Observable`.

    Attributes:
        expression (:obj:`str`): mathematical expression for an Observable
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        species (:obj:`list` of :obj:`Species`): Species used by this Observable expression
        observables (:obj:`list` of :obj:`Observable`): other Observables used by this Observable expression

    Related attributes:
        observable (:obj:`Observable`): observable
    """

    expression = LongStringAttribute(primary=True, unique=True, default='')
    species = ManyToManyAttribute(Species, related_name='observable_expressions')
    observables = ManyToManyAttribute('Observable', related_name='observable_expressions')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        tabular_orientation = TabularOrientation.inline
        expression_term_models = ('Species', 'Observable')
        expression_is_linear = True

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return Expression.serialize(self)

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
        return Expression.deserialize(cls, value, objects)

    def validate(self):
        """ Check that the observable is valid

        * Check that the expression is a linear function

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        return Expression.validate(self, self.observable)


class Observable(obj_model.Model):
    """ Observable: a linear function of other Observbles and Species

    Attributes:
        id (:obj:`str`): unique id
        name (:obj:`str`): name
        model (:obj:`Model`): model
        expression (:obj:`ObservableExpression`): mathematical expression for an Observable
        units (:obj:`MoleculeCountUnit`): units of expression
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
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
    expression = ExpressionOneToOneAttribute(ObservableExpression, related_name='observable',
                                             min_related=1, min_related_rev=1)
    units = EnumAttribute(MoleculeCountUnit, default=MoleculeCountUnit['molecule'])
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='observables')
    evidence = ManyToManyAttribute('Evidence', related_name='observables')
    interpretations = ManyToManyAttribute('Interpretation', related_name='observables')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='observables')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'name', 'expression', 'units',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        expression_term_model = ObservableExpression
        expression_term_units = 'units'


class FunctionExpression(obj_model.Model, Expression):
    """ A mathematical expression of Functions, Observbles, Parameters and Python functions

    The expression used by a :obj:`Function`.

    Attributes:
        expression (:obj:`str`): mathematical expression for a Function
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        species (:obj:`list` of :obj:`Species`): Species used by this function expression
        observables (:obj:`list` of :obj:`Observable`): Observables used by this function expression
        parameters (:obj:`list` of :obj:`Parameter`): Parameters used by this function expression
        functions (:obj:`list` of :obj:`Function`): other Functions used by this function expression
        compartments (:obj:`list` of :obj:`Compartment`): Compartments used by this stop condition expression

    Related attributes:
        function (:obj:`Function`): function
    """
    expression = LongStringAttribute(primary=True, unique=True, default='')
    parameters = ManyToManyAttribute('Parameter', related_name='function_expressions')
    species = ManyToManyAttribute(Species, related_name='function_expressions')
    observables = ManyToManyAttribute(Observable, related_name='function_expressions')
    functions = ManyToManyAttribute('Function', related_name='function_expressions')
    compartments = ManyToManyAttribute(Compartment, related_name='function_expressions')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        tabular_orientation = TabularOrientation.inline
        expression_term_models = ('Parameter', 'Species', 'Observable', 'Function', 'Compartment')

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return Expression.serialize(self)

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
        return Expression.deserialize(cls, value, objects)

    def validate(self):
        """ Check that the function is valid

        * Check that the expression is a valid Python function

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        return Expression.validate(self, self.function)


class Function(obj_model.Model):
    """ Function: a mathematical expression of Functions, Observbles, Parameters and Python functions

    Attributes:
        id (:obj:`str`): unique id
        name (:obj:`str`): name
        model (:obj:`Model`): model
        expression (:obj:`FunctionExpression`): mathematical expression for a Function
        units (:obj:`str`): units
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
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
    expression = ExpressionOneToOneAttribute(FunctionExpression, related_name='function',
                                             min_related=1, min_related_rev=1)
    units = LongStringAttribute()
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='functions')
    evidence = ManyToManyAttribute('Evidence', related_name='functions')
    interpretations = ManyToManyAttribute('Interpretation', related_name='functions')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='functions')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'name', 'expression', 'units',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        expression_term_model = FunctionExpression
        expression_term_units = 'units'

    def validate(self):
        """ Check that the Function is valid

        * Check that `expression` has units `units`

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(Function, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        # check that units are valid
        if self.expression and hasattr(self.expression, '_parsed_expression') and self.expression._parsed_expression:
            exp_units = unit_registry.parse_expression(self.units).to_base_units().units
            try:
                calc = self.expression._parsed_expression.test_eval(with_units=True)
            except ParsedExpressionError as error:
                errors.append(InvalidAttribute(self.Meta.attributes['units'], [str(error)]))
            else:
                if hasattr(calc, 'units'):
                    calc_units = calc.to_base_units().units
                else:
                    calc_units = unit_registry.parse_expression('dimensionless')

                if calc_units != exp_units:
                    errors.append(InvalidAttribute(self.Meta.attributes['units'],
                                                   ['Units of "{}" should be "{}" not "{}"'.format(
                                                    self.expression.expression, exp_units, calc_units)]))

        # return errors
        if errors:
            return InvalidObject(self, errors)
        return None


class StopConditionExpression(obj_model.Model, Expression):
    """ A mathematical expression of Functions, Observables, Parameters and Python functions

    The expression used by a :obj:`StopCondition`.

    Attributes:
        expression (:obj:`str`): mathematical expression for a StopCondition
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        observables (:obj:`list` of :obj:`Observable`): Observables used by this stop condition expression
        parameters (:obj:`list` of :obj:`Parameter`): Parameters used by this stop condition expression
        functions (:obj:`list` of :obj:`Function`): Functions used by this stop condition expression
        compartments (:obj:`list` of :obj:`Compartment`): Compartments used by this stop condition expression

    Related attributes:
        stop_condition (:obj:`StopCondition`): stop condition
    """

    expression = LongStringAttribute(primary=True, unique=True, default='')
    parameters = ManyToManyAttribute('Parameter', related_name='stop_condition_expressions')
    species = ManyToManyAttribute(Species, related_name='stop_condition_expressions')
    observables = ManyToManyAttribute(Observable, related_name='stop_condition_expressions')
    functions = ManyToManyAttribute(Function, related_name='stop_condition_expressions')
    compartments = ManyToManyAttribute(Compartment, related_name='stop_condition_expressions')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        tabular_orientation = TabularOrientation.inline
        expression_term_models = ('Parameter', 'Species', 'Observable', 'Function', 'Compartment')
        expression_type = bool

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return Expression.serialize(self)

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
        return Expression.deserialize(cls, value, objects)

    def validate(self):
        """ Check that the stop condition is valid

        * Check that the expression is a Boolean function

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        return Expression.validate(self, self.stop_condition)


class StopCondition(obj_model.Model):
    """ StopCondition: Simulation of a model terminates when its StopCondition is true.

    A mathematical expression of Functions, Observbles, Parameters and Python functions `StopCondition`s
    are optional. It must return a Boolean.

    Attributes:
        id (:obj:`str`): unique id
        name (:obj:`str`): name
        model (:obj:`Model`): model
        expression (:obj:`StopConditionExpression`): mathematical expression for a StopCondition
        units (:obj:`StopConditionUnit`): units
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        expressions (:obj:`Expressions`): expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='stop_conditions')
    expression = ExpressionOneToOneAttribute(StopConditionExpression, related_name='stop_condition',
                                             min_related=1, min_related_rev=1)
    units = EnumAttribute(StopConditionUnit, default=StopConditionUnit.dimensionless)
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='stop_conditions')
    evidence = ManyToManyAttribute('Evidence', related_name='stop_conditions')
    interpretations = ManyToManyAttribute('Interpretation', related_name='stop_conditions')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='stop_conditions')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'name', 'expression', 'units',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        expression_term_model = StopConditionExpression
        expression_term_units = 'units'

    def validate(self):
        """ Check that the stop condition is valid

        * Check that `expression` has units `units`

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(StopCondition, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        # check that units are valid
        if self.expression \
                and hasattr(self.expression, '_parsed_expression') \
                and self.expression._parsed_expression:
            try:
                test_val = self.expression._parsed_expression.test_eval(with_units=True)
            except ParsedExpressionError as error:
                errors.append(InvalidAttribute(self.Meta.attributes['units'], [str(error)]))
            else:
                if hasattr(test_val, 'units'):
                    errors.append(InvalidAttribute(self.Meta.attributes['units'], ['Units must be dimensionless']))
        else:
            if self.expression:
                expression = self.expression.expression
            else:
                expression = None
            errors.append(InvalidAttribute(self.Meta.attributes['expression'],
                                           ['Expression for {} could not be parsed: {}'.format(
                                            self.id, expression)]))

        if self.units != StopConditionUnit.dimensionless:
            errors.append(InvalidAttribute(self.Meta.attributes['units'], ['Units must be dimensionless']))

        # return errors
        if errors:
            return InvalidObject(self, errors)
        return None


class Reaction(obj_model.Model):
    """ Reaction

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        submodel (:obj:`Submodel`): submodel that reaction belongs to
        participants (:obj:`list` of :obj:`SpeciesCoefficient`): participants
        reversible (:obj:`bool`): indicates if reaction is thermodynamically reversible
        flux_min (:obj:`float`): minimum flux bound for solving an FBA model; negative for reversible reactions
        flux_max (:obj:`float`): maximum flux bound for solving an FBA model
        flux_bound_units (:obj:`ReactionFluxBoundUnit`): units for the minimum and maximum fluxes
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws; if present, rate_laws[0] is the forward
            rate law, and rate_laws[0] is the backward rate law
        dfba_obj_expression (:obj:`DfbaObjectiveExpression`): dFBA objective expression
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='reactions')
    submodel = ManyToOneAttribute(Submodel, related_name='reactions')
    participants = ReactionParticipantAttribute(related_name='reactions')
    reversible = BooleanAttribute()
    rate_units = EnumAttribute(ReactionRateUnit, default=ReactionRateUnit['s^-1'])
    flux_min = FloatAttribute(nan=True)
    flux_max = FloatAttribute(min=0, nan=True)
    flux_bound_units = EnumAttribute(ReactionFluxBoundUnit, default=None, none=True)
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='reactions')
    evidence = ManyToManyAttribute('Evidence', related_name='reactions')
    interpretations = ManyToManyAttribute('Interpretation', related_name='reactions')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='reactions')

    class Meta(obj_model.Model.Meta, ExpressionDynamicTermMeta):
        attribute_order = ('id', 'name', 'submodel',
                           'participants', 'reversible',
                           'rate_units', 'flux_min', 'flux_max', 'flux_bound_units',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        indexed_attrs_tuples = (('id',), )
        expression_term_units = 'rate_units'
        merge = obj_model.ModelMerge.append

    def validate(self):
        """ Check if the reaction is valid

        * If the submodel is ODE or SSA, check that the reaction has a forward rate law
        * If the submodel is ODE or SSA and the reaction is reversible, check that the reaction has a
          backward rate law
        * Check flux units are not None if flux_min or flux_max is defined
        * Check that `flux_min` <= `flux_max`
        * Check that reaction is element and charge balanced

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
        if self.submodel and isinstance(self.submodel.framework, pronto.term.Term) and self.submodel.framework.id in [
            'WCM:ordinary_differential_equations',
            'WCM:stochastic_simulation_algorithm',
        ]:
            if not for_rl:
                rl_errors.append('Reaction in {} submodel must have a forward rate law'.format(
                    self.submodel.framework.name))
            if self.reversible and not rev_rl:
                rl_errors.append('Reversible reaction in {} submodel must have a backward rate law'.format(
                    self.submodel.framework.name))
        if not self.reversible and rev_rl:
            rl_errors.append('Irreversible reaction cannot have a backward rate law')

        if rl_errors:
            errors.append(InvalidAttribute(self.Meta.related_attributes['rate_laws'], rl_errors))

        # check min, max fluxes
        if (not isnan(self.flux_min) or not isnan(self.flux_max)) and self.flux_bound_units is None:
            errors.append(InvalidAttribute(self.Meta.attributes['flux_bound_units'],
                                           ['Units must be defined for the flux bounds']))

        if self.submodel and self.submodel.framework != wcm_ontology['WCM:dynamic_flux_balance_analysis']:
            if not isnan(self.flux_min):
                errors.append(InvalidAttribute(self.Meta.attributes['flux_min'],
                                               ['Minimum flux should be NaN for reactions in non-dFBA submodels']))
            if not isnan(self.flux_max):
                errors.append(InvalidAttribute(self.Meta.attributes['flux_max'],
                                               ['Maximum flux should be NaN for reactions in non-dFBA submodels']))

        if not isnan(self.flux_min) and not isnan(self.flux_min) and self.flux_min > self.flux_max:
            errors.append(InvalidAttribute(self.Meta.attributes['flux_max'],
                                           ['Maximum flux must be least the minimum flux']))

        if self.reversible and not isnan(self.flux_min) and self.flux_min >= 0:
            errors.append(InvalidAttribute(self.Meta.attributes['flux_min'],
                                           ['Minimum flux for reversible reaction should be negative or NaN']))
        if not self.reversible and not isnan(self.flux_min) and self.flux_min < 0:
            errors.append(InvalidAttribute(self.Meta.attributes['flux_min'],
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
                species.extend(rate_law.expression.species.get(__type=__type, **kwargs))

        return det_dedupe(species)

    def get_reactants(self):
        """ Get the species consumed by the reaction

        Returns:
            :obj:`list` of :obj:`Species`: reactant species
        """
        species = []
        for part in self.participants:
            if part.coefficient < 0:
                species.append(part.species)
        return det_dedupe(species)

    def get_products(self):
        """ Get the species produced by the reaction

        Returns:
            :obj:`list` of :obj:`Species`: product species
        """
        species = []
        for part in self.participants:
            if part.coefficient > 0:
                species.append(part.species)
        return det_dedupe(species)

    def get_modifiers(self, direction=RateLawDirection.forward):
        """ Get species in the expression that are not reactants in the reaction

        Returns:
            :obj:`list` of :obj:`Species`: species in the expression that are not reactants in the reaction
        """
        return list(set(self.rate_laws.get_one(direction=direction).expression.species) - set(self.get_reactants()))

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
        if self.submodel.framework == wcm_ontology['WCM:dynamic_flux_balance_analysis']:
            fbc_reaction_plugin = wrap_libsbml(sbml_reaction.getPlugin, 'fbc')
            for bound in ['lower', 'upper']:
                # make a unique ID for each flux bound parameter
                # ids for wc_lang Parameters all start with 'parameter'
                param_id = "_reaction_{}_{}_bound".format(self.id, bound)
                param = create_sbml_parameter(sbml_model, id=param_id, value=self.flux_min,
                                              units='mmol_per_gDW_per_hr')
                if bound == 'lower':
                    wrap_libsbml(param.setValue, self.flux_min)
                    wrap_libsbml(fbc_reaction_plugin.setLowerFluxBound, param_id)
                if bound == 'upper':
                    wrap_libsbml(param.setValue, self.flux_max)
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
        ordering = ('species', 'coefficient')

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
                species_id = Species.gen_id_static(match.group(5), compartment.get_primary_attribute())
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


class RateLawExpression(obj_model.Model, Expression):
    """ Rate law expression

    Attributes:
        expression (:obj:`str`): mathematical expression of the rate law
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        species (:obj:`list` of :obj:`Species`): species whose dynamic concentrations are used in the rate law
        parameters (:obj:`list` of :obj:`Parameter`): parameters whose values are used in the rate law
        compartments (:obj:`list` of :obj:`Compartment`): Compartments used by this stop condition expression

    Related attributes:
        rate_law (:obj:`RateLaw`): the `RateLaw` which uses this `RateLawExpression`
    """
    expression = LongStringAttribute(primary=True, unique=True, default='')
    parameters = ManyToManyAttribute('Parameter', related_name='rate_law_expressions')
    species = ManyToManyAttribute(Species, related_name='rate_law_expressions')
    observables = ManyToManyAttribute(Observable, related_name='rate_law_expressions')
    functions = ManyToManyAttribute(Function, related_name='rate_law_expressions')
    compartments = ManyToManyAttribute(Compartment, related_name='rate_law_expressions')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        attribute_order = ('expression', 'species', 'parameters')
        tabular_orientation = TabularOrientation.inline
        ordering = ('expression',)
        expression_term_models = ('Parameter', 'Species', 'Observable', 'Function', 'Compartment')

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return Expression.serialize(self)

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
        return Expression.deserialize(cls, value, objects)

    def validate(self):
        """ Determine whether a `RateLawExpression` is valid

        * Check that all of the species and parameters contribute to the expression
        * Check that the species and parameters encompass of the named entities in the expression

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors in an `InvalidObject` instance
        """
        return Expression.validate(self, self.rate_laws[0])


class RateLaw(obj_model.Model):
    """ Rate law

    Attributes:
        id (:obj:`str`): identifier equal to `{reaction.id}-{direction.name}`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        reaction (:obj:`Reaction`): reaction
        direction (:obj:`RateLawDirection`): direction
        type (:obj:`pronto.term.Term`): type
        expression (:obj:`RateLawExpression`): expression
        units (:obj:`ReactionRateUnit`): units
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='rate_laws')
    reaction = ManyToOneAttribute(Reaction, related_name='rate_laws')
    direction = EnumAttribute(RateLawDirection, default=RateLawDirection.forward)
    type = OntologyAttribute(wcm_ontology,
                             namespace='WCM',
                             terms=wcm_ontology['WCM:rate_law'].rchildren(),
                             default=None, none=True)
    expression = ExpressionManyToOneAttribute(RateLawExpression, min_related=1, min_related_rev=1, related_name='rate_laws')
    units = EnumAttribute(ReactionRateUnit, default=ReactionRateUnit['s^-1'])
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='rate_laws')
    evidence = ManyToManyAttribute('Evidence', related_name='rate_laws')
    interpretations = ManyToManyAttribute('Interpretation', related_name='rate_laws')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='rate_laws')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'name', 'reaction', 'direction', 'type',
                           'expression', 'units',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        # unique_together = (('reaction', 'direction'), )
        expression_term_model = RateLawExpression
        expression_term_units = 'units'

    def gen_id(self):
        """ Generate identifier

        Returns:
            :obj:`str`: identifier
        """
        return '{}-{}'.format(self.reaction.id, self.direction.name)

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

        # check ID
        if self.reaction and self.direction is not None and \
                self.id != self.gen_id():
            errors.append(InvalidAttribute(
                self.Meta.attributes['id'],
                ['Id must be {}'.format(self.gen_id())]))

        # check that units are valid
        if self.expression and self.reaction and self.units != self.reaction.rate_units:
            errors.append(InvalidAttribute(self.Meta.attributes['units'],
                                           ['Units must the same as reaction rate units']))

        if self.expression \
                and hasattr(self.expression, '_parsed_expression') \
                and self.expression._parsed_expression:
            exp_units = unit_registry.parse_expression(self.units.name).to_base_units().units
            try:
                calc = self.expression._parsed_expression.test_eval(with_units=True)
            except ParsedExpressionError as error:
                errors.append(InvalidAttribute(self.Meta.attributes['units'], [str(error)]))
            else:
                if hasattr(calc, 'units'):
                    calc_units = calc.to_base_units().units
                else:
                    calc_units = unit_registry.parse_expression('dimensionless')

                if calc_units != exp_units:
                    errors.append(InvalidAttribute(self.Meta.attributes['units'],
                                                   ['Units of "{}" should be "{}" not "{}"'.format(
                                                    self.expression.expression, exp_units, calc_units)]))
        else:
            if self.expression:
                expression = self.expression.expression
            else:
                expression = None
            errors.append(InvalidAttribute(self.Meta.attributes['expression'],
                                           ['Expression for {} could not be parsed: {}'.format(
                                            self.id, expression)]))

        """ return errors or `None` to indicate valid object """
        if errors:
            return InvalidObject(self, errors)
        return None


class DfbaObjSpecies(obj_model.Model):
    """ DfbaObjSpecies

    A dFBA objective reaction contains a list of DfbaObjSpecies instances. Distinct DfbaObjSpecies
    enable separate comments and references for each one.

    Attributes:
        id (:obj:`str`): unique identifier per DfbaObjSpecies equal to
            `dfba-net-species-{dfba_obj_reaction.id}-{species.id}`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        dfba_obj_reaction (:obj:`DfbaObjReaction`): the dFBA objective reaction that uses the dFBA objective species
        species (:obj:`Species`): species
        value (:obj:`float`): the specie's reaction coefficient
        units (:obj:`DfbaObjSpeciesUnit`): units of the value
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='dfba_obj_species')
    dfba_obj_reaction = ManyToOneAttribute('DfbaObjReaction', min_related=1, related_name='dfba_obj_species',
                                           verbose_name='dFBA objective reaction',
                                           verbose_related_name='dFBA objective species')
    species = ManyToOneAttribute(Species, min_related=1, related_name='dfba_obj_species',
                                 verbose_related_name='dFBA objective species')
    value = FloatAttribute()
    units = EnumAttribute(DfbaObjSpeciesUnit, default=DfbaObjSpeciesUnit['M s^-1'])
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='dfba_obj_species',
                                                   verbose_related_name='dFBA objective species')
    evidence = ManyToManyAttribute('Evidence', related_name='dfba_obj_species',
                                   verbose_related_name='dFBA objective species')
    interpretations = ManyToManyAttribute('Interpretation', related_name='dfba_obj_species',
                                          verbose_related_name='dFBA objective species')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='dfba_obj_species',
                                     verbose_related_name='dFBA objective species')

    class Meta(obj_model.Model.Meta):
        # unique_together = (('dfba_obj_reaction', 'species'), )
        attribute_order = ('id', 'name', 'dfba_obj_reaction',
                           'species', 'value', 'units',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        verbose_name = 'dFBA objective species'
        verbose_name_plural = 'dFBA objective species'
        merge = obj_model.ModelMerge.append

    def gen_id(self):
        """ Generate identifier equal to
        `dfba-net-species-{dfba_obj_reaction.id}-{species.id}`

        Returns:
            :obj:`str`: identifier
        """
        return 'dfba-net-species-{}-{}'.format(self.dfba_obj_reaction.id, self.species.id)

    def validate(self):
        """ Validate that the dFBA objective is valid

        * Check if the identifier is equal to
          `dfba-net-species-{dfba_obj_reaction.id}-{species.id}`
        * Units consistent with units of cell size

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(DfbaObjSpecies, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        # id
        if self.dfba_obj_reaction and self.species and self.id != self.gen_id():
            errors.append(InvalidAttribute(self.Meta.attributes['id'], ['Id must be {}'.format(
                self.gen_id())]))

        # units consistent with units of cell size
        if self.dfba_obj_reaction and \
                ((self.units == DfbaObjSpeciesUnit['M s^-1'] and
                    self.dfba_obj_reaction.cell_size_units != DfbaCellSizeUnit.l) or
                 (self.units == DfbaObjSpeciesUnit['mol gDCW^-1 s^-1'] and
                    self.dfba_obj_reaction.cell_size_units != DfbaCellSizeUnit.gDCW)):
            errors.append(InvalidAttribute(
                self.Meta.attributes['units'],
                ['Units {} are not consistent with cell size units {}'.format(
                    self.units.name, self.dfba_obj_reaction.cell_size_units.name)]))

        if errors:
            return InvalidObject(self, errors)
        return None


class DfbaObjReaction(obj_model.Model):
    """ A pseudo-reaction used to represent the interface between metabolism and other
    cell processes.

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        submodel (:obj:`Submodel`): submodel that uses this reaction
        units (:obj:`ReactionRateUnit`): rate units
        cell_size_units (:obj:`DfbaCellSizeUnit`): cell size units
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        dfba_obj_expression (:obj:`DfbaObjectiveExpression`): dFBA objectie expression
        dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): the components of this dFBA objective reaction
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='dfba_obj_reactions', verbose_related_name='dFBA objective reactions')
    submodel = ManyToOneAttribute(Submodel, related_name='dfba_obj_reactions', verbose_related_name='dFBA objective reactions')
    units = EnumAttribute(ReactionRateUnit, default=ReactionRateUnit['s^-1'])
    cell_size_units = EnumAttribute(DfbaCellSizeUnit, default=DfbaCellSizeUnit.l)
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='dfba_obj_reactions',
                                                   verbose_related_name='dFBA objective reactions')
    evidence = ManyToManyAttribute('Evidence', related_name='dfba_obj_reactions',
                                   verbose_related_name='dFBA objective reactions')
    interpretations = ManyToManyAttribute('Interpretation', related_name='dfba_obj_reactions',
                                          verbose_related_name='dFBA objective reactions')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='dfba_obj_reactions', verbose_related_name='dFBA objective reactions')

    class Meta(obj_model.Model.Meta, ExpressionDynamicTermMeta):
        attribute_order = ('id', 'name', 'submodel', 'units', 'cell_size_units',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        indexed_attrs_tuples = (('id',), )
        verbose_name = 'dFBA objective reaction'
        expression_term_units = 'units'
        merge = obj_model.ModelMerge.append

    def add_to_sbml_doc(self, sbml_document):
        """ Add a DfbaObjReaction to a libsbml SBML document.

        DfbaObjReactions are added to the SBML document because they can be used in a dFBA submodel's
        objective function. In fact the default objective function is the submodel's dFBA objective reaction.
        Since SBML does not define DfbaObjReaction as a separate class, DfbaObjReactions are added
        to the SBML model as SBML reactions.
        CheckModel ensures that wc_lang DfbaObjReactions and Reactions have distinct ids.

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

        # write dFBA objective reaction participants to SBML document
        for dfba_obj_species in self.dfba_obj_species:
            if dfba_obj_species.value < 0:
                species_reference = wrap_libsbml(sbml_reaction.createReactant)
                wrap_libsbml(species_reference.setStoichiometry, -dfba_obj_species.value)
            elif 0 < dfba_obj_species.value:
                species_reference = wrap_libsbml(sbml_reaction.createProduct)
                wrap_libsbml(species_reference.setStoichiometry, dfba_obj_species.value)
            id = dfba_obj_species.species.gen_sbml_id()
            wrap_libsbml(species_reference.setSpecies, id)
            wrap_libsbml(species_reference.setConstant, True)

        # the dFBA objective reaction does not constrain the optimization, so set its bounds to 0 and INF
        fbc_reaction_plugin = wrap_libsbml(sbml_reaction.getPlugin, 'fbc')
        for bound in ['lower', 'upper']:
            # make a unique ID for each flux bound parameter
            # ids for wc_lang Parameters all start with 'parameter'
            param_id = "_dfba_obj_reaction_{}_{}_bound".format(self.id, bound)
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
        type (:obj:`pronto.term.Term`): parameter type
        value (:obj:`float`): value
        std (:obj:`float`): standard error of the value
        units (:obj:`str`): units of the value and standard error
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        density_compartment (:obj:`Compartment`): compartments whose density is represented by the parameter
        observable_expressions (:obj:`list` of :obj:`ObservableExpression`): observable expressions
        function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        stop_condition_expressions (:obj:`list` of :obj:`StopConditionExpression`): stop condition expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='parameters')
    type = OntologyAttribute(wcm_ontology,
                             namespace='WCM',
                             terms=wcm_ontology['WCM:parameter'].rchildren(),
                             default=None, none=True)
    value = FloatAttribute()
    std = FloatAttribute(min=0, verbose_name='Standard error')
    units = StringAttribute(min_length=1)
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='parameters')
    evidence = ManyToManyAttribute('Evidence', related_name='parameters')
    interpretations = ManyToManyAttribute('Interpretation', related_name='parameters')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='parameters')

    class Meta(obj_model.Model.Meta, ExpressionStaticTermMeta):
        attribute_order = ('id', 'name', 'type',
                           'value', 'std', 'units',
                           'db_refs', 'evidence', 'interpretations', 'comments', 'references')
        expression_term_value = 'value'
        expression_term_units = 'units'

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


class Evidence(obj_model.Model):
    """ Evidence

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        value (:obj:`str`): value
        std (:obj:`str`): standard error of the value
        units (:obj:`str`): units
        type (:obj:`pronto.term.Term`): type
        taxon (:obj:`str`): taxon in which the evidence was observed
        genetic_variant (:obj:`str`): genetic variant in which the evidence was observed
        temp (:obj:`float`): temperature at which the evidence was observed
        temp_units (:obj:`TemperatureUnit`): temperature units
        ph (:obj:`float`): pH at which the evidence was observed
        ph_units (:obj:`PhUnit): pH units
        growth_media (:obj:`str`): growth media at which the evidence was observed
        condition (:obj:`str`): experimental conditions (e.g. control)
        experiment_type (:obj:`str`): type of experiment (e.g. RNA-seq)
        experiment_design (:obj:`str`): experimental design
        measurement_method (:obj:`str`): method used to measure data (e.g. deep sequencing)
        analysis_method (:obj:`str`): method used to analyze data (e.g. Cufflinks)
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence underlying reduced evidence
            (e.g. individual observations underlying an average)
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        submodels (:obj:`list` of :obj:`Submodel`): submodels
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species types
        species (:obj:`list` of :obj:`Species`): species
        distribution_init_concentrations (:obj:`list` of :obj:`DistributionInitConcentration`):
            distributions of initial concentrations of species at the beginning of each
            cell cycle
        observables (:obj:`list` of :obj:`Observable`): observables
        functions (:obj:`list` of :obj:`Function`): functions
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        dfba_obj_reactions (:obj:`list` of :obj:`DfbaObjReaction`): dFBA objective reactions
        dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
        stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
        parameters (:obj:`list` of :obj:`Parameter`): parameters
        interpretations (:obj:`list` of :obj:`Interpretation`): interpretations
        reduced_evidences (:obj:`list` of :obj:`Evidence`): reduced evidence that the evidence
            supports (e.g. averages supported by this and other evidence)
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='evidences')
    value = StringAttribute()
    std = StringAttribute(verbose_name='Standard error')
    units = StringAttribute()
    type = OntologyAttribute(wcm_ontology,
                             namespace='WCM',
                             terms=wcm_ontology['WCM:evidence'].rchildren(),
                             default=None, none=True)
    taxon = StringAttribute()
    genetic_variant = StringAttribute()
    temp = FloatAttribute(nan=True, verbose_name='Temperature')
    temp_units = EnumAttribute(TemperatureUnit, none=True, verbose_name='Temperature units')
    ph = FloatAttribute(nan=True, verbose_name='pH')
    ph_units = EnumAttribute(PhUnit, none=True, verbose_name='pH units')
    growth_media = LongStringAttribute()
    condition = LongStringAttribute()
    experiment_type = LongStringAttribute()
    experiment_design = LongStringAttribute()
    measurement_method = LongStringAttribute()
    analysis_method = LongStringAttribute()
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='evidences')
    evidence = ManyToManyAttribute('Evidence', related_name='reduced_evidences')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='evidences')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'value', 'std', 'units',
                           'type',
                           'taxon', 'genetic_variant',
                           'temp', 'temp_units', 'ph', 'ph_units', 'growth_media', 'condition',
                           'experiment_type', 'experiment_design', 'measurement_method', 'analysis_method',
                           'db_refs', 'evidence', 'comments', 'references')
        verbose_name_plural = 'Evidence'

    def validate(self):
        """ Determine if the evidence is valid

        * temperature units are defined if the temperature is not None
        * pH units are defined if the pH is not None

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(Evidence, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        if self.temp is not None and not isnan(self.temp) \
                and not isinstance(self.temp_units, TemperatureUnit):
            errors.append(InvalidAttribute(self.Meta.attributes['temp_units'],
                                           ['Temperature units must be defined']))
        if self.ph is not None and not isnan(self.ph) \
                and not isinstance(self.ph_units, PhUnit):
            errors.append(InvalidAttribute(self.Meta.attributes['ph_units'],
                                           ['pH units must be defined']))

        if errors:
            return InvalidObject(self, errors)
        return None


class Interpretation(obj_model.Model):
    """ Interpretation of evidence

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        value (:obj:`str`): value
        std (:obj:`str`): standard error of the value
        units (:obj:`str`): units
        type (:obj:`pronto.term.Term`): type
        method (:obj:`str`): procedure which produced the interpretation
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        evidence (:obj:`list` of :obj:`Evidence`): evidence underlying reduced evidence
            (e.g. individual observations underlying an average)
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:
        submodels (:obj:`list` of :obj:`Submodel`): submodels
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species types
        species (:obj:`list` of :obj:`Species`): species
        distribution_init_concentrations (:obj:`list` of :obj:`DistributionInitConcentration`):
            distributions of initial concentrations of species at the beginning of each
            cell cycle
        observables (:obj:`list` of :obj:`Observable`): observables
        functions (:obj:`list` of :obj:`Function`): functions
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        dfba_obj_reactions (:obj:`list` of :obj:`DfbaObjReaction`): dFBA objective reactions
        dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
        stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
        parameters (:obj:`list` of :obj:`Parameter`): parameters
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='interpretations')
    value = StringAttribute()
    std = StringAttribute(verbose_name='Standard error')
    units = StringAttribute()
    type = OntologyAttribute(wcm_ontology,
                             namespace='WCM',
                             terms=wcm_ontology['WCM:interpretation'].rchildren(),
                             default=None, none=True)
    method = LongStringAttribute()
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='interpretations')
    evidence = ManyToManyAttribute('Evidence', related_name='interpretations')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='interpretations')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'value', 'std', 'units',
                           'type', 'method',
                           'db_refs', 'evidence', 'comments', 'references')


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
        type (:obj:`proto.term.Term`): type
        publication (:obj:`str`): publication title
        publisher (:obj:`str`): publisher
        series (:obj:`str`): series
        volume (:obj:`str`): volume
        number (:obj:`str`): number
        issue (:obj:`str`): issue
        edition (:obj:`str`): edition
        chapter (:obj:`str`): chapter
        pages (:obj:`str`): page range
        db_refs (:obj:`list` of :obj:`DatabaseReference`): database references
        comments (:obj:`str`): comments

    Related attributes:
        taxon (:obj:`Taxon`): taxon
        env (:obj:`Environment`): environment
        submodels (:obj:`list` of :obj:`Submodel`): submodels
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species types
        species (:obj:`list` of :obj:`Species`): species
        distribution_init_concentrations (:obj:`list` of :obj:`DistributionInitConcentration`):
            distributions of initial concentrations of species at the beginning of
            each cell cycle
        observables (:obj:`list` of :obj:`Observable`): observables
        functions (:obj:`list` of :obj:`Function`): functions
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
        stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
        parameters (:obj:`list` of :obj:`Parameter`): parameters
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='references')
    title = StringAttribute()
    author = StringAttribute()
    editor = StringAttribute()
    year = PositiveIntegerAttribute()
    type = OntologyAttribute(wcm_ontology,
                             namespace='WCM',
                             terms=wcm_ontology['WCM:reference'].rchildren(),
                             default=None, none=True)
    publication = StringAttribute()
    publisher = StringAttribute()
    series = StringAttribute()
    volume = StringAttribute()
    number = StringAttribute()
    issue = StringAttribute()
    edition = StringAttribute()
    chapter = StringAttribute()
    pages = StringAttribute()
    db_refs = DatabaseReferenceManyToManyAttribute(related_name='references')
    comments = CommentAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'title', 'author', 'editor', 'year', 'type', 'publication', 'publisher',
                           'series', 'volume', 'number', 'issue', 'edition', 'chapter', 'pages',
                           'db_refs', 'comments')


class DatabaseReference(obj_model.Model):
    """ Reference to a source database entry

    Attributes:
        database (:obj:`str`): database name
        id (:obj:`str`): id of database entry

    Related attributes:
        model (:obj:`Model`): model
        taxon (:obj:`Taxon`): taxon
        env (:obj:`Environment`): environment
        submodels (:obj:`list` of :obj:`Submodel`): submodels
        compartments (:obj:`list` of :obj:`Compartment`): compartments
        species_types (:obj:`list` of :obj:`SpeciesType`): species types
        species (:obj:`list` of :obj:`Species`): species
        distribution_init_concentrations (:obj:`list` of :obj:`DistributionInitConcentration`):
            distributions of initial concentrations of species at the beginning of
            each cell cycle
        observables (:obj:`list` of :obj:`Observable`): observables
        functions (:obj:`list` of :obj:`Function`): functions
        reactions (:obj:`list` of :obj:`Reaction`): reactions
        rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        dfba_obj_reactions (:obj:`list` of :obj:`DfbaObjReaction`): dFBA objective reactions
        dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
        stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
        parameters (:obj:`list` of :obj:`Parameter`): parameters
        references (:obj:`list` of :obj:`Reference`): references
    """

    database = StringAttribute(min_length=1)
    id = StringAttribute(min_length=1)

    class Meta(obj_model.Model.Meta):
        unique_together = (('database', 'id', ), )
        tabular_orientation = TabularOrientation.inline
        attribute_order = ('database', 'id')
        frozen_columns = 2
        ordering = ('database', 'id', )

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return '{}: {}'.format(self.database, self.id)

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of :obj:`DatabaseReference`, `InvalidAttribute` or `None`: tuple
                of cleaned value and cleaning error
        """
        if ': ' not in value:
            return (None, InvalidAttribute(cls.Meta.attributes['id'], ['Invalid format']))

        database, _, id = value.strip().partition(': ')
        db_ref = cls(database=database.strip(), id=id.strip())

        if DatabaseReference not in objects:
            objects[DatabaseReference] = {}

        serialized_val = db_ref.serialize()
        if serialized_val in objects[DatabaseReference]:
            db_ref = objects[DatabaseReference][serialized_val]
        else:
            objects[DatabaseReference][serialized_val] = db_ref

        return (db_ref, None)


class Validator(obj_model.Validator):
    def run(self, model, get_related=True):
        """ Validate a model and return its errors

        Args:
            model (:obj:`Model`): model
            get_related (:obj:`bool`, optional): if true, get all related objects

        Returns:
            :obj:`InvalidObjectSet` or `None`: list of invalid objects/models and their errors
        """
        return super(Validator, self).run(model, get_related=get_related)

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
* :obj:`Identifier`

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
from math import ceil, floor, exp, log, log10, isinf, isnan
from natsort import natsorted, ns
from obj_model import (BooleanAttribute, EnumAttribute,
                       FloatAttribute,
                       IntegerAttribute, PositiveIntegerAttribute,
                       RegexAttribute, SlugAttribute, StringAttribute, LongStringAttribute,
                       UrlAttribute, EmailAttribute, DateTimeAttribute,
                       OneToOneAttribute, ManyToOneAttribute, ManyToManyAttribute, OneToManyAttribute,
                       ManyToOneRelatedManager,
                       InvalidObject, InvalidAttribute, TabularOrientation)
from obj_model.expression import (ExpressionOneToOneAttribute, ExpressionManyToOneAttribute,
                                  ExpressionStaticTermMeta, ExpressionDynamicTermMeta,
                                  ExpressionExpressionTermMeta, Expression,
                                  ParsedExpression, ParsedExpressionError)
from obj_model.ontology import OntologyAttribute
from obj_model.units import UnitAttribute
from wc_lang.sbml.util import SbmlModelMixin, SbmlAssignmentRuleMixin, LibSbmlInterface
from wc_utils.util.chem import EmpiricalFormula, OpenBabelUtils
from wc_utils.util.enumerate import CaseInsensitiveEnum, CaseInsensitiveEnumMeta
from wc_utils.util.list import det_dedupe
from wc_onto import onto
from wc_utils.util.ontology import are_terms_equivalent
from wc_utils.util.units import unit_registry, are_units_equivalent
from wc_utils.workbook.core import get_column_letter
import bpforms
import collections
import datetime
import libsbml
import networkx
import obj_model
import obj_model.chem
import openbabel
import pkg_resources
import pronto.term
import re
import scipy.constants
import stringcase
import token
import wc_lang.config.core
import warnings

with open(pkg_resources.resource_filename('wc_lang', 'VERSION'), 'r') as file:
    wc_lang_version = file.read().strip()

# wc_lang generates obj_model SchemaWarning warnings because some Models lack primary attributes.
# These models include :obj:`RateLaw`, :obj:`SpeciesCoefficient`, :obj:`RateLawExpression`, and :obj:`Species`.
# However, these are not needed by the workbook and delimiter-separated representations of
# models on disk. Therefore, suppress the warnings.
warnings.filterwarnings('ignore', '', obj_model.SchemaWarning, 'obj_model')

# configuration

call_libsbml = LibSbmlInterface.call_libsbml


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


class TaxonRank(int, Enum, metaclass=TaxonRankMeta):
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


class RateLawDirection(int, CaseInsensitiveEnum):
    """ Rate law directions """
    backward = -1
    forward = 1


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
            :obj:`list` of :obj:`SpeciesCoefficient`: cleaned value
            :obj:`InvalidAttribute`: cleaning error
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
            :obj:`list` of :obj:`SpeciesCoefficient`: list of species coefficients
            :obj:`list` of :obj:`Exception`: list of errors
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
                species_id = Species._gen_id(species_type.id, compartment.id)
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
                if part.species.species_type.structure and part.species.species_type.structure.empirical_formula:
                    delta_formula += part.species.species_type.structure.empirical_formula * part.coefficient

                if part.species.species_type.structure is None or part.species.species_type.structure.charge is None:
                    errors.append('Charge must be defined for {}'.format(part.species.species_type.id))
                else:
                    delta_charge += part.species.species_type.structure.charge * part.coefficient

            if not errors:
                if delta_formula:
                    errors.append('Reaction is element imbalanced: {}'.format(delta_formula))
                if delta_charge != 0.:
                    errors.append('Reaction is charge imbalanced: {}'.format(delta_charge))

            if errors:
                return InvalidAttribute(self, errors)

        # return None
        return None

    def get_excel_validation(self):
        """ Get Excel validation

        Returns:
            :obj:`wc_utils.workbook.io.FieldValidation`: validation
        """
        validation = super(ManyToManyAttribute, self).get_excel_validation()

        related_class = Species
        related_ws = related_class.Meta.verbose_name_plural
        related_col = get_column_letter(related_class.get_attr_index(related_class.Meta.primary_attribute) + 1)
        source = '{}:{}'.format(related_ws, related_col)

        validation.ignore_blank = False
        input_message = ['Enter a reaction string using species from "{}".'.format(source)]
        error_message = ['Value must be a reaction string using species from "{}".'.format(source)]

        input_message.append(('Examples:\n'
                              '* [c]: atp + h2o ==> adp + pi + h\n'
                              '* glc[e] + atp[c] + h2o[c] ==> glc[e] + adp[c] + pi[c] + h[c]\n'
                              '* (3) Na[c] + (2) K[e] ==> (3) Na[e] + (2) K[c]'))

        error_message.append(('Examples:\n'
                              '* [c]: atp + h2o ==> adp + pi + h\n'
                              '* glc[e] + atp[c] + h2o[c] ==> glc[e] + adp[c] + pi[c] + h[c]\n'
                              '* (3) Na[c] + (2) K[e] ==> (3) Na[e] + (2) K[c]'))

        validation.input_message = '\n\n'.join(input_message)
        validation.error_message = '\n\n'.join(error_message)

        return validation


class EvidenceManyToManyAttribute(ManyToManyAttribute):
    """ Many to many attribute for evidence """

    def serialize(self, evidence, encoded=None):
        """ Serialize related object

        Args:
            evidence (:obj:`list` of :obj:`Evidence`): Python representation of evidence
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: simple Python representation
        """
        return '; '.join(ev.serialize() for ev in evidence)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`list` of :obj:`Evidence`: cleaned value
            :obj:`InvalidAttribute`: cleaning error
        """
        value = (value or '').strip()
        if not value:
            return ([], None)

        objs = []
        errors = []
        for v in value.split(';'):
            obj, error = Evidence.deserialize(v.strip(), objects)
            if error:
                errors.extend(error.messages)
            else:
                objs.append(obj)

        if errors:
            return (None, InvalidAttribute(self, errors))
        return (objs, None)

    def get_excel_validation(self):
        """ Get Excel validation

        Returns:
            :obj:`wc_utils.workbook.io.FieldValidation`: validation
        """
        validation = super(ManyToManyAttribute, self).get_excel_validation()

        validation.ignore_blank = False

        related_class = Observation
        related_ws = related_class.Meta.verbose_name_plural
        related_col = get_column_letter(related_class.get_attr_index(related_class.Meta.primary_attribute) + 1)
        source = '{}:{}'.format(related_ws, related_col)

        input_message = ['Enter a list of observations from "{}".'.format(source)]
        error_message = ['Value must be a list of observations from "{}".'.format(source)]

        input_message.append(('Examples:\n'
                              '* Obs1(+, s=50, q=100); Obs2(-)\n'
                              '* Obs3(~, s=90)\n'
                              '* Obs4(+, q=30)'))

        error_message.append(('Examples:\n'
                              '* Obs1(+, s=50, q=100); Obs2(-)\n'
                              '* Obs3(~, s=90)\n'
                              '* Obs4(+, q=30)'))

        validation.input_message = '\n\n'.join(input_message)
        validation.error_message = '\n\n'.join(error_message)

        return validation


class IdentifierOneToManyAttribute(OneToManyAttribute):
    def __init__(self, related_name='', verbose_name='Identifiers', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(IdentifierOneToManyAttribute, self).__init__('Identifier', related_name=related_name,
                                                           verbose_name=verbose_name,
                                                           verbose_related_name=verbose_related_name,
                                                           help=help)

    def serialize(self, identifiers, encoded=None):
        """ Serialize related object

        Args:
            identifiers (:obj:`list` of :obj:`Identifier`): Python representation of identifiers
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: string representation
        """
        sorted_ids = sorted(identifiers, key=lambda id: (id.namespace, id.id))
        return ', '.join(id.serialize() for id in sorted_ids)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`list` of :obj:`Identifier`: cleaned value
            :obj:`InvalidAttribute`: cleaning error
        """
        value = value or ''
        value = value.strip()
        if not value:
            return ([], None)

        identifiers = set()
        errors = []
        for val in value.split(','):
            id, invalid = Identifier.deserialize(val, objects)
            if invalid:
                errors.extend(invalid.messages)
            else:
                identifiers.add(id)

        if errors:
            return (None, InvalidAttribute(self, errors))
        else:
            return (det_dedupe(identifiers), None)

    def get_excel_validation(self):
        """ Get Excel validation

        Returns:
            :obj:`wc_utils.workbook.io.FieldValidation`: validation
        """
        validation = super(OneToManyAttribute, self).get_excel_validation()

        validation.ignore_blank = True
        input_message = ['Enter a comma-separated list of identifiers in external namespaces.']
        error_message = ['Value must be a comma-separated list of identifiers in external namespaces.']

        input_message.append(('Examples:\n'
                              '* doi: 10.1016/j.tcb.2015.09.004\n'
                              '* chebi: CHEBI:15377, kegg.compound: C00001'))

        error_message.append(('Examples:\n'
                              '* doi: 10.1016/j.tcb.2015.09.004\n'
                              '* chebi: CHEBI:15377, kegg.compound: C00001'))

        validation.input_message = '\n\n'.join(input_message)
        validation.error_message = '\n\n'.join(error_message)

        return validation


class IdentifierManyToManyAttribute(ManyToManyAttribute):
    def __init__(self, related_name='', verbose_name='Identifiers', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(IdentifierManyToManyAttribute, self).__init__('Identifier', related_name=related_name,
                                                            verbose_name=verbose_name,
                                                            verbose_related_name=verbose_related_name,
                                                            help=help)

    def serialize(self, identifiers, encoded=None):
        """ Serialize related object

        Args:
            identifiers (:obj:`list` of :obj:`Identifier`): Python representation of identifiers
            encoded (:obj:`dict`, optional): dictionary of objects that have already been encoded

        Returns:
            :obj:`str`: string representation
        """
        sorted_ids = sorted(identifiers, key=lambda id: (id.namespace, id.id))
        return ', '.join(id.serialize() for id in sorted_ids)

    def deserialize(self, value, objects, decoded=None):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            decoded (:obj:`dict`, optional): dictionary of objects that have already been decoded

        Returns:
            :obj:`list` of :obj:`Identifier`: cleaned value
            :obj:`InvalidAttribute`: cleaning error
        """
        value = value or ''
        value = value.strip()
        if not value:
            return ([], None)

        identifiers = set()
        errors = []
        for val in value.split(','):
            id, invalid = Identifier.deserialize(val, objects)
            if invalid:
                errors.extend(invalid.messages)
            else:
                identifiers.add(id)

        if errors:
            return (None, InvalidAttribute(self, errors))
        else:
            return (det_dedupe(identifiers), None)

    def get_excel_validation(self):
        """ Get Excel validation

        Returns:
            :obj:`wc_utils.workbook.io.FieldValidation`: validation
        """
        validation = super(ManyToManyAttribute, self).get_excel_validation()

        validation.ignore_blank = True
        input_message = ['Enter a comma-separated list of identifiers in external namespaces.']
        error_message = ['Value must be a comma-separated list of identifiers in external namespaces.']

        input_message.append(('Examples:\n'
                              '* doi: 10.1016/j.tcb.2015.09.004\n'
                              '* chebi: CHEBI:15377, kegg.compound: C00001'))

        error_message.append(('Examples:\n'
                              '* doi: 10.1016/j.tcb.2015.09.004\n'
                              '* chebi: CHEBI:15377, kegg.compound: C00001'))

        validation.input_message = '\n\n'.join(input_message)
        validation.error_message = '\n\n'.join(error_message)

        return validation


class SubmodelsToModelRelatedManager(ManyToOneRelatedManager):
    """ Submodels to model related manager """

    def gen_models(self):
        """ Generate separate models for each submodel, as well as a "core" model which contains the
        model integration framework (compartments, species types, species, observables, functions,
        stop conditions, etc.) and all model objects that are not associated with at least 1 submodel.

        Returns:
            :obj:`Model`: model with objects that (a) form the model integration framework or (b) are
                not associated with any submodel

            :obj:`list` of :obj:`Model`: one model for each submodel
        """
        model = self.object

        # if 1 or fewer submodels, return model
        if len(model.submodels) <= 1:
            return (model.copy(), [])

        # cut submodels from model
        submodels = self.cut(kind='submodel')

        # get objects that won't be exported with at least one submodel to export as "core" model
        core_model = model.copy()
        objs_in_core = set(core_model.get_related())
        for submodel in core_model.submodels:
            objs_in_core.discard(submodel)
            objs_in_core.difference_update(set(submodel.get_children(kind='submodel')))
        objs_in_core.add(core_model)

        # add additional objects to "core" model that define model integration framework
        # * compartments
        # * species types
        # * species
        # * observables
        # * functions
        # * stop conditions
        for obj in list(objs_in_core):
            objs_in_core.update(obj.get_children(kind='core_model'))

        # generate "core" model which has
        # * the model integration framework and
        # * other objects not exported with at least one 1 submodel
        for obj in objs_in_core:
            obj.cut_relations(objs_in_core)

        # return core model and separate models for each submodel
        return (core_model, [submodel.model for submodel in submodels])


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
        """
        left_val = getattr(left, self.name)
        right_val = getattr(right, self.name)
        if left_val != right_val:
            setattr(left, self.name, left_val + self.SEPARATOR + right_val)


class Model(obj_model.Model, SbmlModelMixin):
    """ Model

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        version (:obj:`str`): version of the model
        url (:obj:`str`): url of the model Git repository
        branch (:obj:`str`): branch of the model Git repository
        revision (:obj:`str`): hash for git commit of the wc_lang version used by a model
        wc_lang_version (:obj:`str`): version of ``wc_lang``
        time_units (:obj:`unit_registry.Unit`): time units
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        comments (:obj:`str`): comments
        created (:obj:`datetime.datetime`): date created
        updated (:obj:`datetime.datetime`): date updated

    Related attributes:

        * taxon (:obj:`Taxon`): taxon
        * env (:obj:`Environment`): environment
        * submodels (:obj:`list` of :obj:`Submodel`): submodels
        * compartments (:obj:`list` of :obj:`Compartment`): compartments
        * species_types (:obj:`list` of :obj:`SpeciesType`): species types
        * species (:obj:`list` of :obj:`Species`): species
        * distribution_init_concentrations (:obj:`list` of :obj:`DistributionInitConcentration`):
          distributions of initial concentrations of species at the beginning of
          each cell cycle
        * observables (:obj:`list` of :obj:`Observable`): observables
        * functions (:obj:`list` of :obj:`Function`): functions
        * reactions (:obj:`list` of :obj:`Reaction`): reactions
        * rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        * dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        * dfba_obj_reactions (:obj:`list` of :obj:`DfbaObjReaction`): dFBA objective reactions
        * dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
        * stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
        * parameters (:obj:`list` of :obj:`Parameter`): parameters
        * observations (:obj:`list` of :obj:`Observation`): observations
        * conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        * references (:obj:`list` of :obj:`Reference`): references
        * authors (:obj:`list` of :obj:`Author`): authors
        * changes (:obj:`list` of :obj:`Change`): changes
    """
    id = SlugAttribute()
    name = StringAttribute()
    version = RegexAttribute(pattern=r'^[0-9]+(\.[0-9]+(\.[0-9a-z]+)?)?$', flags=re.I)
    url = UrlAttribute(verbose_name='URL')
    branch = StringAttribute()
    revision = StringAttribute()
    wc_lang_version = RegexAttribute(pattern=r'^[0-9]+(\.[0-9]+(\.[0-9a-z]+)?)?$', flags=re.I,
                                     default=wc_lang_version, verbose_name='wc_lang version')
    time_units = UnitAttribute(unit_registry,
                               choices=[unit_registry.parse_units('s')],
                               default=unit_registry.parse_units('s'))
    identifiers = IdentifierOneToManyAttribute(related_name='model')
    comments = CommentAttribute()
    created = DateTimeAttribute()
    updated = DateTimeAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'version',
                           'url', 'branch', 'revision',
                           'wc_lang_version',
                           'time_units',
                           'identifiers', 'comments',
                           'created', 'updated')
        tabular_orientation = TabularOrientation.column
        children = {
            'submodel': ('taxon', 'env'),
            'core_model': ('taxon', 'env',
                           'compartments', 'species_types', 'species', 'distribution_init_concentrations',
                           'observables', 'functions', 'stop_conditions'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'version', 'url', 'branch', 'revision',
                     'wc_lang_version', 'time_units', 'identifiers', 'comments', 'updated', 'created'),
            'wc_sim': ('id', 'version', 'wc_lang_version', 'time_units'),
        }

    def __init__(self, **kwargs):
        """
        Args:
            kwargs (:obj:`dict`, optional): dictionary of keyword arguments with keys equal to the names of the model attributes
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
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching
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
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

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
                    or not are_terms_equivalent(comp.parent_compartment.biological_type,
                                                onto['WC:cellular_compartment']):
                roots.append(comp)
        return roots

    def get_species_types(self, __type=None, **kwargs):
        """ Get all species types

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

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
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

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
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

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
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

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
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

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
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

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
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

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
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

        Returns:
            :obj:`list` of :obj:`DfbaObjReaction`: dFBA objective reactions
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.dfba_obj_reactions.get(__type=__type, **kwargs)

    def get_dfba_obj_species(self, __type=None, **kwargs):
        """ Get all dFBA objective species used by submodels

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

        Returns:
            :obj:`list` of :obj:`DfbaObjSpecies`: dFBA objective species
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.dfba_obj_species.get(__type=__type, **kwargs)

    def get_parameters(self, __type=None, **kwargs):
        """ Get all parameters from model and submodels

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

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
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

        Returns:
            :obj:`list` of :obj:`StopCondition`: stop conditions
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.stop_conditions.get(__type=__type, **kwargs)

    def get_observations(self, __type=None, **kwargs):
        """ Get all observations for model

        Returns:
            :obj:`list` of :obj:`Observation`: observations for model
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.observations.get(__type=__type, **kwargs)

    def get_observation_sets(self, __type=None, **kwargs):
        """ Get all observation sets for model

        Returns:
            :obj:`list` of :obj:`ObservationSet`: observation sets for model
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.observation_sets.get(__type=__type, **kwargs)

    def get_evidence(self, __type=None, **kwargs):
        """ Get all evidence for model

        Returns:
            :obj:`list` of :obj:`Evidence`: evidence for model
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        evidence = []
        for obs in self.observations:
            evidence.extend(obs.evidence.get(__type=__type, **kwargs))
        return evidence

    def get_conclusions(self, __type=None, **kwargs):
        """ Get all conclusions for model

        Returns:
            :obj:`list` of :obj:`Conclusion`: conclusions for model
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.conclusions.get(__type=__type, **kwargs)

    def get_references(self, __type=None, **kwargs):
        """ Get all references from model and children

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

        Returns:
            :obj:`list` of :obj:`Reference`: references
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.references.get(__type=__type, **kwargs)

    def get_authors(self, __type=None, **kwargs):
        """ Get all authors from model and children

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

        Returns:
            :obj:`list` of :obj:`Author`: authors
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.authors.get(__type=__type, **kwargs)

    def get_changes(self, __type=None, **kwargs):
        """ Get all changes from model and children

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

        Returns:
            :obj:`list` of :obj:`Change`: changes
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')

        return self.changes.get(__type=__type, **kwargs)

    def get_components(self, __type=None, **kwargs):
        """ Find model component of `type` with `id`

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching objects

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
                'dfba_objs', 'reactions', 'rate_laws', 'dfba_obj_reactions', 'dfba_obj_species',
                'stop_conditions', 'parameters', 'observations', 'observation_sets', 'evidence', 'conclusions',
                'references', 'authors', 'changes',
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

    def export_to_sbml(self, sbml_model):
        """ Add this metadata about this model to a SBML model.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`libsbml.Model`: SBML model
        """
        if self.submodels:
            return

        call_libsbml(sbml_model.setIdAttribute, self.gen_sbml_id())
        call_libsbml(sbml_model.setName, self.name)
        LibSbmlInterface.set_commments(self, sbml_model)

        return sbml_model

    def export_relations_to_sbml(self, sbml_model, sbml):
        """ Add relationships to/from object to SBML model.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            sbml (:obj:`libsbml.Model`): SBML model
        """
        if self.submodels:
            return

        annots = []

        annots.extend(['version', 'url', 'branch', 'revision', 'wc_lang_version',
                       'identifiers',
                       'updated', 'created'])

        if self.taxon:
            annots.extend(['taxon.id', 'taxon.name', 'taxon.rank', 'taxon.identifiers', 'taxon.comments'])

        if self.env:
            annots.extend(['env.id', 'env.name', 'env.temp', 'env.temp_units',
                           'env.identifiers', 'env.comments'])

        xml_annotation = '<annotation><wcLang:annotation>' \
                         + LibSbmlInterface.gen_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml) \
                         + LibSbmlInterface.gen_authors_annotation(self) \
                         + '</wcLang:annotation></annotation>'
        call_libsbml(sbml.setAnnotation, xml_annotation)

    def import_from_sbml(self, sbml):
        """ Load from SBML model

        Args:
            sbml (:obj:`libsbml.Model`): SBML model
        """
        self.id = self.parse_sbml_id(call_libsbml(sbml.getIdAttribute))
        self.name = call_libsbml(sbml.getName)
        LibSbmlInterface.get_commments(self, sbml)

    def import_relations_from_sbml(self, sbml, objs):
        """ Load relationships from SBML model

        Args:
            sbml (:obj:`libsbml.Model`): SBML model
            objs (:obj:`dict`): dictionary that maps WC-Lang types to dictionaries that
                map the ids of WC-Lang objects to WC-Lang objects
        """
        parsed_annots = LibSbmlInterface.parse_annotations(sbml)
        annots = []

        # identifiers
        annots.extend(['version', 'url', 'branch', 'revision', 'wc_lang_version', 'identifiers', 'updated', 'created'])

        if 'taxon.id' in parsed_annots:
            self.taxon = Taxon()
            annots.extend(['taxon.id', 'taxon.name', 'taxon.rank', 'taxon.identifiers', 'taxon.comments'])

        if 'env.id' in parsed_annots:
            self.env = Environment()
            annots.extend(['env.id', 'env.name', 'env.temp', 'env.temp_units', 'env.identifiers', 'env.comments'])

        LibSbmlInterface.get_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml, objs)

        # authors
        LibSbmlInterface.get_authors_annotation(self, sbml, objs)


class Taxon(obj_model.Model, SbmlModelMixin):
    """ Biological taxon (e.g. family, genus, species, strain, etc.)

    Attributes:
        id (:obj:`str`): unique identifier equal to 'taxon'
        name (:obj:`str`): name
        model (:obj:`Model`): model
        rank (:obj:`TaxonRank`): rank
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = RegexAttribute(pattern=r'^taxon$', primary=True, unique=True)
    name = StringAttribute()
    model = OneToOneAttribute(Model, related_name='taxon')
    rank = EnumAttribute(TaxonRank, default=TaxonRank.species)
    identifiers = IdentifierOneToManyAttribute(related_name='taxon')
    comments = CommentAttribute()
    references = OneToManyAttribute('Reference', related_name='taxon')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'rank',
                           'identifiers', 'comments', 'references')
        tabular_orientation = TabularOrientation.column
        children = {
            'submodel': ('identifiers', 'references'),
            'core_model': ('identifiers', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'rank', 'identifiers', 'comments'),
            'wc_sim': (),
        }


class Environment(obj_model.Model, SbmlModelMixin):
    """ Environment

    Attributes:
        id (:obj:`str`): unique identifier equal to 'env'
        name (:obj:`str`): name
        model (:obj:`Model`): model
        temp (:obj:`float`): temperature
        temp_units (:obj:`unit_registry.Unit`): temperature units
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = RegexAttribute(pattern=r'^env$', primary=True, unique=True)
    name = StringAttribute()
    model = OneToOneAttribute(Model, related_name='env')
    temp = FloatAttribute(verbose_name='Temperature')
    temp_units = UnitAttribute(unit_registry,
                               choices=(unit_registry.parse_units('celsius'),),
                               default=unit_registry.parse_units('celsius'),
                               verbose_name='Temperature units')
    identifiers = IdentifierOneToManyAttribute(related_name='env')
    comments = CommentAttribute()
    references = OneToManyAttribute('Reference', related_name='env')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'temp', 'temp_units',
                           'identifiers', 'comments', 'references')
        tabular_orientation = TabularOrientation.column
        children = {
            'submodel': ('identifiers', 'references'),
            'core_model': ('identifiers', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'temp', 'temp_units', 'identifiers', 'comments'),
            'wc_sim': (),
        }


class Submodel(obj_model.Model, SbmlModelMixin):
    """ Submodel

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        framework (:obj:`pronto.term.Term`): modeling framework (e.g. dynamic flux balance analysis)
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:

        * reactions (:obj:`list` of :obj:`Reaction`): reactions
        * dfba_obj (:obj:`DfbaObjective`): objective function for a dFBA submodel;
          if not initialized, then `dfba_obj_reaction` is used as the objective function
        * dfba_obj_reactions (:obj:`list` of :obj:`DfbaObjReaction`): the growth reaction for a dFBA submodel
        * changes (:obj:`list` of :obj:`Change`): changes
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='submodels', related_manager=SubmodelsToModelRelatedManager)
    framework = OntologyAttribute(onto,
                                  namespace='WC',
                                  terms=onto['WC:modeling_framework'].rchildren(),
                                  default=onto['WC:stochastic_simulation_algorithm'],
                                  none=False)
    identifiers = IdentifierManyToManyAttribute(related_name='submodels')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='submodels')
    conclusions = ManyToManyAttribute('Conclusion', related_name='submodels')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='submodels')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name', 'framework',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        indexed_attrs_tuples = (('id',), )
        merge = obj_model.ModelMerge.append
        children = {
            'submodel': ('model', 'reactions', 'dfba_obj', 'dfba_obj_reactions',
                         'identifiers', 'evidence', 'conclusions', 'references', 'changes'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'framework', 'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'framework'),
        }

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

        if are_terms_equivalent(self.framework, onto['WC:dynamic_flux_balance_analysis']):
            if not self.dfba_obj:
                errors.append(InvalidAttribute(self.Meta.related_attributes['dfba_obj'],
                                               ['dFBA submodel must have an objective']))

        if errors:
            return InvalidObject(self, errors)
        return None

    def get_children(self, kind=None, __type=None, recursive=True, __include_stop_conditions=True,
                     **kwargs):
        """ Get a kind of children.

        If :obj:`kind` is :obj:`None`, children are defined to be the values of the related attributes defined
        in each class.

        Args:
            kind (:obj:`str`, optional): kind of children to get
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            recursive (:obj:`bool`, optional): if :obj:`True`, get children recursively
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs

        Returns:
            :obj:`list` of :obj:`Model`: children
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')
        if '__include_stop_conditions' in kwargs:
            __include_stop_conditions = kwargs.pop('__include_stop_conditions')

        children = self.get_immediate_children(kind=kind, __include_stop_conditions=__include_stop_conditions)

        # get recursive children
        if recursive:
            objs_to_explore = children
            children = set(children)
            while objs_to_explore:
                obj_to_explore = objs_to_explore.pop()
                for child in obj_to_explore.get_immediate_children(kind=kind):
                    if child not in children:
                        children.add(child)
                        objs_to_explore.append(child)
            children = list(children)

        # filter by type/attributes
        matches = []
        for child in children:
            if child.has_attr_vals(__type=__type, **kwargs):
                matches.append(child)
        children = matches

        # return children
        return children

    def get_immediate_children(self, kind=None, __type=None, __include_stop_conditions=True,
                               **kwargs):
        """ Get a kind of children.

        If :obj:`kind` is :obj:`None`, children are defined to be the values of the related attributes defined
        in each class.

        Args:
            kind (:obj:`str`, optional): kind of children to get
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs

        Returns:
            :obj:`list` of :obj:`Model`: children
        """
        if '__type' in kwargs:
            __type = kwargs.pop('__type')
        if '__include_stop_conditions' in kwargs:
            __include_stop_conditions = kwargs.pop('__include_stop_conditions')

        immediate_children = super(Submodel, self).get_immediate_children(kind=kind)

        if kind == 'submodel' and __include_stop_conditions:
            all_children = set(self.get_children(kind=kind, __include_stop_conditions=False))

            for stop_condition in self.model.stop_conditions:
                stop_condition_children = stop_condition.get_children(kind='submodel')
                stop_condition_species = set()
                for child in stop_condition_children:
                    if isinstance(child, Species):
                        stop_condition_species.add(child)
                if not stop_condition_species.difference(all_children):
                    immediate_children.append(stop_condition)

            immediate_children = det_dedupe(immediate_children)

        # filter by type/attributes
        matches = []
        for child in immediate_children:
            if child.has_attr_vals(__type=__type, **kwargs):
                matches.append(child)
        immediate_children = matches

        return immediate_children

    def export_to_sbml(self, sbml_model):
        """ Add metadata about the submodel to a SBML model.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`libsbml.Model`: SBML model
        """
        call_libsbml(sbml_model.setIdAttribute, self.gen_sbml_id())
        call_libsbml(sbml_model.setName, self.name)
        LibSbmlInterface.set_commments(self, sbml_model)

        return sbml_model

    def export_relations_to_sbml(self, sbml_model, sbml):
        """ Add relationships to/from object to SBML model.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            sbml (:obj:`libsbml.Model`): SBML model
        """
        annots = []

        annots.extend(['framework', 'identifiers'])

        annots.extend(['model.id', 'model.name',  'model.version', 'model.url', 'model.branch', 'model.revision', 'model.wc_lang_version',
                       'model.identifiers',
                       'model.comments', 'model.updated', 'model.created'])

        if self.model.taxon:
            annots.extend(['model.taxon.id', 'model.taxon.name', 'model.taxon.rank', 'model.taxon.identifiers', 'model.taxon.comments'])

        if self.model.env:
            annots.extend(['model.env.id', 'model.env.name',
                           'model.env.temp', 'model.env.temp_units',
                           'model.env.identifiers', 'model.env.comments'])

        xml_annotation = '<annotation><wcLang:annotation>' \
            + LibSbmlInterface.gen_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml) \
            + LibSbmlInterface.gen_authors_annotation(self.model) \
            + '</wcLang:annotation></annotation>'
        call_libsbml(sbml.setAnnotation, xml_annotation)

    def import_from_sbml(self, sbml):
        """ Load from SBML model

        Args:
            sbml (:obj:`libsbml.Model`): SBML model
        """
        self.id = self.parse_sbml_id(call_libsbml(sbml.getIdAttribute))
        self.name = call_libsbml(sbml.getName)
        LibSbmlInterface.get_commments(self, sbml)

    def import_relations_from_sbml(self, sbml, objs):
        """ Load relationships from SBML model

        Args:
            sbml (:obj:`libsbml.Model`): SBML model
            objs (:obj:`dict`): dictionary that maps WC-Lang types to dictionaries that
                map the ids of WC-Lang objects to WC-Lang objects
        """
        parsed_annots = LibSbmlInterface.parse_annotations(sbml)
        annots = []

        # identifiers
        annots.extend(['framework', 'identifiers'])

        annots.extend(['model.id', 'model.name',  'model.version', 'model.url', 'model.branch', 'model.revision', 'model.wc_lang_version',
                       'model.identifiers', 'model.comments', 'model.updated', 'model.created'])

        if 'model.taxon.id' in parsed_annots:
            self.model.taxon = Taxon()
            annots.extend(['model.taxon.id', 'model.taxon.name', 'model.taxon.rank', 'model.taxon.identifiers', 'model.taxon.comments'])
        if 'model.env.id' in parsed_annots:
            self.model.env = Environment()
            annots.extend(['model.env.id', 'model.env.name',
                           'model.env.temp', 'model.env.temp_units',
                           'model.env.identifiers', 'model.env.comments'])

        LibSbmlInterface.get_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml, objs)

        # authors
        LibSbmlInterface.get_authors_annotation(self.model, sbml, objs)


class DfbaObjectiveExpression(obj_model.Model, Expression, SbmlModelMixin):
    """ A mathematical expression of Reactions and DfbaObjReactions

    The expression used by a :obj:`DfbaObjective`.

    Attributes:
        expression (:obj:`str`): mathematical expression
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        reactions (:obj:`list` of :obj:`Reaction`): reactions used by this expression
        dfba_obj_reactions (:obj:`list` of :obj:`Species`): dFBA objective reactions used by this expression

    Related attributes:

        * dfba_obj (:obj:`DfbaObjective`): dFBA objective
    """

    expression = LongStringAttribute(primary=True, unique=True, default='')
    reactions = OneToManyAttribute('Reaction', related_name='dfba_obj_expression',
                                   verbose_related_name='dFBA objective expression')
    dfba_obj_reactions = OneToManyAttribute('DfbaObjReaction', related_name='dfba_obj_expression',
                                            verbose_name='dFBA objective reactions', verbose_related_name='dFBA objective expression')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        tabular_orientation = TabularOrientation.cell
        expression_valid_functions = ()
        expression_term_models = ('Reaction', 'DfbaObjReaction')
        verbose_name = 'dFBA objective expression'
        merge = obj_model.ModelMerge.append
        expression_unit_registry = unit_registry
        children = {
            'submodel': ('reactions', 'dfba_obj_reactions'),
            'core_model': ('reactions', 'dfba_obj_reactions'),
        }
        child_attrs = {
            'sbml': ('expression', 'reactions', 'dfba_obj_reactions'),
            'wc_sim': ('expression', 'reactions', 'dfba_obj_reactions'),
        }

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

    def merge_attrs(self, other, other_objs_in_self, self_objs_in_other):
        """ Merge attributes of two objects

        Args:
            other (:obj:`obj_model.Model`): other model
            other_objs_in_self (:obj:`dict`): dictionary that maps instances of objects in another model to objects
                in a model
            self_objs_in_other (:obj:`dict`): dictionary that maps instances of objects in a model to objects
                in another model
        """
        super(DfbaObjectiveExpression, self).merge_attrs(other, other_objs_in_self, self_objs_in_other)
        Expression.merge_attrs(self, other, other_objs_in_self, self_objs_in_other)


class DfbaObjective(obj_model.Model, SbmlModelMixin):
    """ dFBA objective function

    Attributes:
        id (:obj:`str`): identifier equal to `dfba-obj-{submodel.id}`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        submodel (:obj:`Submodel`): the `Submodel` which uses this `DfbaObjective`
        expression (:obj:`DfbaObjectiveExpression`): mathematical expression of the objective function
        units (:obj:`unit_registry.Unit`): units
        reaction_rate_units (:obj:`unit_registry.Unit`): reaction rate units
        coefficient_units (:obj:`unit_registry.Unit`): coefficient units
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='dfba_objs', verbose_related_name='dFBA objectives')
    submodel = OneToOneAttribute(Submodel, related_name='dfba_obj', min_related=1, verbose_related_name='dFBA objective')
    expression = ExpressionOneToOneAttribute(DfbaObjectiveExpression, related_name='dfba_obj',
                                             min_related=1, min_related_rev=1, verbose_related_name='dFBA objective')
    units = UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('dimensionless'),),
                          default=unit_registry.parse_units('dimensionless'))
    reaction_rate_units = UnitAttribute(unit_registry,
                                        choices=(unit_registry.parse_units('s^-1'),),
                                        default=unit_registry.parse_units('s^-1'))
    coefficient_units = UnitAttribute(unit_registry,
                                      choices=(unit_registry.parse_units('s'),),
                                      default=unit_registry.parse_units('s'))
    identifiers = IdentifierManyToManyAttribute(related_name='dfba_objs', verbose_related_name='dFBA objectives')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='dfba_objs', verbose_related_name='dFBA objectives')
    conclusions = ManyToManyAttribute('Conclusion', related_name='dfba_objs', verbose_related_name='dFBA objectives')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='dfba_objs', verbose_related_name='dFBA objectives')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        verbose_name = 'dFBA objective'
        attribute_order = ('id', 'name', 'submodel', 'expression', 'units',
                           'reaction_rate_units', 'coefficient_units',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        expression_term_model = DfbaObjectiveExpression
        expression_term_units = 'units'
        merge = obj_model.ModelMerge.append
        children = {
            'submodel': ('expression', 'identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('expression', 'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'submodel', 'expression', 'units',
                     'reaction_rate_units', 'coefficient_units', 'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'submodel', 'expression', 'units',
                       'reaction_rate_units', 'coefficient_units'),
        }

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

    def export_to_sbml(self, sbml_model):
        """ Add this dFBA objective to a SBML model.

        This uses version 2 of the 'Flux Balance Constraints' extension. SBML assumes that an
        DfbaObjective is a linear combination of reaction fluxes.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`libsbml.Objective`: SBML objective
        """
        # warn if objective is not linear
        if not self.expression._parsed_expression.is_linear:
            warnings.warn("SBML export doesn't support the non-linear objective for submodel '{}'".format(
                self.submodel.id), WcLangWarning)
            return

        sbml_plugin = call_libsbml(sbml_model.getPlugin, 'fbc')
        sbml = call_libsbml(sbml_plugin.createObjective)

        # sense
        call_libsbml(sbml.setType, 'maximize')

        # set objective as active
        # In SBML 3 FBC 2, the 'activeObjective' attribute must be set on ListOfObjectives.
        # Since a submodel has only one Objective, it must be the active one.
        list_of_objectives = call_libsbml(sbml_plugin.getListOfObjectives)
        call_libsbml(list_of_objectives.setActiveObjective, self.gen_sbml_id())

        # id, name
        call_libsbml(sbml.setIdAttribute, self.gen_sbml_id())
        call_libsbml(sbml.setName, self.name)

        # comments
        LibSbmlInterface.set_commments(self, sbml)

        # return SBML objective
        return sbml

    def export_relations_to_sbml(self, sbml_model, sbml):
        """ Add relationships to/from object to SBML objective.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            sbml (:obj:`libsbml.Objective`): SBML objective
        """
        # expression
        for rxn in self.expression.reactions:
            sbml_rxn_id = rxn.gen_sbml_id()
            coeff = self.expression._parsed_expression.lin_coeffs[Reaction][rxn]
            sbml_flux_obj = call_libsbml(sbml.createFluxObjective)
            call_libsbml(sbml_flux_obj.setReaction, sbml_rxn_id)
            call_libsbml(sbml_flux_obj.setCoefficient, coeff)

        for dfba_obj_rxn in self.expression.dfba_obj_reactions:
            sbml_rxn_id = dfba_obj_rxn.gen_sbml_id()
            coeff = self.expression._parsed_expression.lin_coeffs[DfbaObjReaction][dfba_obj_rxn]
            sbml_flux_obj = call_libsbml(sbml.createFluxObjective)
            call_libsbml(sbml_flux_obj.setReaction, sbml_rxn_id)
            call_libsbml(sbml_flux_obj.setCoefficient, coeff)

        # units, identifiers
        annots = ['units', 'reaction_rate_units', 'coefficient_units', 'identifiers']
        LibSbmlInterface.set_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml)

    def import_from_sbml(self, sbml):
        """ Load from SBML objective

        Args:
            sbml (:obj:`libsbml.Objective`): SBML objective
        """
        # id, name
        self.id = self.parse_sbml_id(call_libsbml(sbml.getIdAttribute))
        self.name = call_libsbml(sbml.getName)

        # comments
        LibSbmlInterface.get_commments(self, sbml)

    def import_relations_from_sbml(self, sbml, objs):
        """ Load relationships from SBML objective

        Args:
            sbml (:obj:`libsbml.Objective`): SBML objective
            objs (:obj:`dict`): dictionary that maps WC-Lang types to dictionaries that
                map the ids of WC-Lang objects to WC-Lang objects
        """
        # submodel
        self.submodel = self.model.submodels[0]

        # expression
        expression = []
        for i_flux_obj in range(call_libsbml(sbml.getNumFluxObjectives, returns_int=True)):
            sbml_flux_obj = call_libsbml(sbml.getFluxObjective, i_flux_obj)
            sbml_rxn_id = call_libsbml(sbml_flux_obj.getReaction)
            coeff = call_libsbml(sbml_flux_obj.getCoefficient)

            if sbml_rxn_id.startswith('Reaction__'):
                wc_lang_rxn_id = Reaction.parse_sbml_id(sbml_rxn_id)
            else:
                wc_lang_rxn_id = DfbaObjReaction.parse_sbml_id(sbml_rxn_id)

            expression.append('{} * {}'.format(coeff, wc_lang_rxn_id))

        expression.sort()
        self.expression, error = DfbaObjectiveExpression.deserialize(' + '.join(expression), objs)
        assert error is None, str(error)

        # units, identifiers
        annots = ['units', 'reaction_rate_units', 'coefficient_units', 'identifiers']
        LibSbmlInterface.get_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml, objs)

    def get_products(self, __type=None, **kwargs):
        """ Get the species produced by this objective function

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching
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


class InitVolume(obj_model.Model, SbmlModelMixin):
    """ Initial volume of a compartment

    Attributes:
        distribution (:obj:`proto.Term`): distribution
        mean (:obj:`float`): mean initial volume
        std (:obj:`float`): standard  deviation of the mean initial volume
        units (:obj:`unit_registry.Unit`): units of volume

    Related attributes:

        * compartments (:obj:`list` of :obj:`Compartment`): compartment
    """
    distribution = OntologyAttribute(onto,
                                     namespace='WC',
                                     terms=onto['WC:random_distribution'].rchildren(),
                                     default=onto['WC:normal_distribution'])
    mean = FloatAttribute(min=0)
    std = FloatAttribute(min=0, verbose_name='Standard deviation')
    units = UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('l'),),
                          default=unit_registry.parse_units('l'))

    class Meta(obj_model.Model.Meta):
        tabular_orientation = TabularOrientation.multiple_cells
        unique_together = (('distribution', 'mean', 'std', 'units'), )
        attribute_order = ('distribution', 'mean', 'std', 'units')
        children = {
            'submodel': (),
            'core_model': (),
        }
        child_attrs = {
            'sbml': ('distribution', 'mean', 'std', 'units'),
            'wc_sim': ('distribution', 'mean', 'std', 'units'),
        }

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return '__'.join([str(self.distribution), str(self.mean), str(self.std), str(self.units)])


class Ph(obj_model.Model, SbmlModelMixin):
    """ pH of a compartment

    Attributes:
        distribution (:obj:`proto.Term`): distribution
        mean (:obj:`float`): mean initial volume
        std (:obj:`float`): standard  deviation of the mean initial volume
        units (:obj:`unit_registry.Unit`): units of volume

    Related attributes:

        * compartments (:obj:`list` of :obj:`Compartment`): compartment
    """
    distribution = OntologyAttribute(onto,
                                     namespace='WC',
                                     terms=onto['WC:random_distribution'].rchildren(),
                                     default=onto['WC:normal_distribution'])
    mean = FloatAttribute()
    std = FloatAttribute(min=0, verbose_name='Standard deviation')
    units = UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('dimensionless'),),
                          default=unit_registry.parse_units('dimensionless'))

    class Meta(obj_model.Model.Meta):
        tabular_orientation = TabularOrientation.multiple_cells
        unique_together = (('distribution', 'mean', 'std', 'units'), )
        attribute_order = ('distribution', 'mean', 'std', 'units')
        children = {
            'submodel': (),
            'core_model': (),
        }
        child_attrs = {
            'sbml': ('distribution', 'mean', 'std', 'units'),
            'wc_sim': (),
        }

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return '__'.join([str(self.distribution), str(self.mean), str(self.std), str(self.units)])


class Compartment(obj_model.Model, SbmlModelMixin):
    """ Compartment

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        biological_type (:obj:`pronto.term.Term`): biological type
        physical_type (:obj:`pronto.term.Term`): physical type
        geometry (:obj:`pronto.term.Term`): geometry
        parent_compartment (:obj:`Compartment`): parent compartment
        mass_units (:obj:`unit_registry.Unit`): mass units
        init_volume (:obj:`InitVolume`): initial volume
        init_density (:obj:`Parameter`): function that calculates the density during the initialization of
            each simulation
        ph (:obj:`Ph`): pH
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:

        * sub_compartments (:obj:`list` of :obj:`Compartment`): compartments contained in this compartment
        * species (:obj:`list` of :obj:`Species`): species in this compartment
        * function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        * rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        * stop_condition_expressions (:obj:`list` of :obj:`StopConditionExpression`): stop condition expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='compartments')
    biological_type = OntologyAttribute(onto,
                                        namespace='WC',
                                        terms=onto['WC:biological_compartment'].rchildren(),
                                        default=onto['WC:cellular_compartment'],
                                        none=True)
    physical_type = OntologyAttribute(onto,
                                      namespace='WC',
                                      terms=onto['WC:physical_compartment'].rchildren(),
                                      default=onto['WC:fluid_compartment'],
                                      none=True)
    geometry = OntologyAttribute(onto,
                                 namespace='WC',
                                 terms=onto['WC:geometric_compartment'].rchildren(),
                                 default=onto['WC:3D_compartment'],
                                 none=True)
    parent_compartment = ManyToOneAttribute('Compartment', related_name='sub_compartments')
    mass_units = UnitAttribute(unit_registry,
                               choices=(unit_registry.parse_units('g'),),
                               default=unit_registry.parse_units('g'))
    init_volume = ManyToOneAttribute(InitVolume, related_name='compartments', verbose_name='Initial volume')
    init_density = OneToOneAttribute('Parameter', related_name='density_compartment', verbose_name='Initial density')
    ph = ManyToOneAttribute(Ph, related_name='compartments', verbose_name='pH')
    identifiers = IdentifierManyToManyAttribute(related_name='compartments')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='compartments')
    conclusions = ManyToManyAttribute('Conclusion', related_name='compartments')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='compartments')

    class Meta(obj_model.Model.Meta, ExpressionDynamicTermMeta):
        attribute_order = ('id', 'name',
                           'biological_type', 'physical_type', 'geometry', 'parent_compartment',
                           'mass_units', 'init_volume', 'init_density', 'ph',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        expression_term_units = 'mass_units'
        children = {
            'submodel': (  # 'parent_compartment', 'sub_compartments',
                'init_volume', 'init_density', 'ph', 'identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': (
                'parent_compartment', 'sub_compartments', 'init_volume', 'init_density', 'ph',
                'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'biological_type', 'physical_type', 'geometry',
                     'parent_compartment', 'mass_units', 'init_volume', 'init_density', 'ph',
                     'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'mass_units', 'init_volume', 'init_density'),
        }

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

        if are_terms_equivalent(self.geometry, onto['WC:3D_compartment']):
            if not self.init_density:
                errors.append(InvalidAttribute(self.Meta.attributes['init_density'],
                                               ['Initial density must be defined for 3D compartments']))
            elif not are_units_equivalent(self.init_density.units, unit_registry.parse_units('g l^-1'),
                                          check_same_magnitude=True):
                errors.append(InvalidAttribute(self.Meta.attributes['init_density'],
                                               ['Initial density of 3D compartment must have units `{}`'.format(
                                                'g l^-1')]))

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
        tot = self.init_volume.mean
        for comp in self.get_sub_compartments(nested=True):
            tot += comp.init_volume.mean
        return tot

    def export_to_sbml(self, sbml_model):
        """ Add this compartment to a SBML model.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`libsbml.Compartment`: SBML compartment

        Raises:
            :obj:`ValueError`: if the geometry cannot be exported to SBML
        """
        sbml = call_libsbml(sbml_model.createCompartment)

        call_libsbml(sbml.setIdAttribute, self.gen_sbml_id())
        call_libsbml(sbml.setName, self.name)

        if self.geometry == onto['WC:3D_compartment']:
            call_libsbml(sbml.setSpatialDimensions, 3)
        else:
            raise ValueError('Unsupported geometry {}'.format(str(self.geometry)))

        if self.init_volume:
            call_libsbml(sbml.setSize, self.init_volume.mean)
            LibSbmlInterface.set_unit(sbml.setUnits, self.init_volume.units)
        call_libsbml(sbml.setConstant, False)

        if self.init_density:
            param_id = '__mass__' + self.gen_sbml_id()
            LibSbmlInterface.create_parameter(sbml_model, param_id, self.init_volume.mean *
                                              self.init_density.value, self.mass_units, constant=False)

            rule = LibSbmlInterface.call_libsbml(sbml_model.createAssignmentRule)
            rule_id = '__mass_rule__' + self.gen_sbml_id()
            LibSbmlInterface.call_libsbml(rule.setIdAttribute, rule_id)
            LibSbmlInterface.call_libsbml(rule.setVariable, param_id)
            formula_parts = []
            for species in self.species:
                formula_part = '{} * {} gram'.format(species.gen_sbml_id(), species.species_type.structure.molecular_weight)
                formula_parts.append(formula_part)
            str_formula = '({}) / {} {}'.format(' + '.join(formula_parts), scipy.constants.Avogadro,
                                                LibSbmlInterface.gen_unit_id(Species.Meta.attributes['units'].choices[0]))
            sbml_formula = call_libsbml(libsbml.parseL3Formula, str_formula)
            call_libsbml(rule.setMath, sbml_formula)

        LibSbmlInterface.set_commments(self, sbml)

        return sbml

    def export_relations_to_sbml(self, sbml_model, sbml):
        """ Add relationships to/from object to SBML compartment.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            sbml (:obj:`libsbml.Compartment`): SBML compartment
        """
        annots = ['biological_type', 'physical_type',
                  'parent_compartment',
                  'init_volume.distribution', 'init_volume.std', 'init_density',
                  'identifiers']
        if self.ph:
            annots.extend(['ph.distribution', 'ph.mean', 'ph.std', 'ph.units'])
        LibSbmlInterface.set_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml)

    def import_from_sbml(self, sbml):
        """ Load from SBML compartment

        Args:
            sbml (:obj:`libsbml.Compartment`): SBML compartment

        Raises:
            :obj:`ValueError`: if the geometry cannot be imported from SBML
        """
        self.id = self.parse_sbml_id(call_libsbml(sbml.getIdAttribute))
        self.name = call_libsbml(sbml.getName)

        dims = call_libsbml(sbml.getSpatialDimensions, returns_int=True)
        if dims == 3:
            self.geometry == onto['WC:3D_compartment']
        else:
            raise ValueError('Unsupported spatial dimensions {}'.format(dims))

        if call_libsbml(sbml.isSetSize):
            mean = call_libsbml(sbml.getSize)
            units = LibSbmlInterface.get_unit(sbml.getUnits)
            init_volume = InitVolume(mean=mean, units=units)
            for comp in self.model.compartments:
                if comp.init_volume and comp.init_volume.mean == init_volume.mean:
                    init_volume = comp.init_volume
                    break
            self.init_volume = init_volume

        param_id = '__mass__' + self.gen_sbml_id()
        sbml_model = call_libsbml(sbml.getModel)
        sbml_param = sbml_model.getParameter(param_id)  # not wrapped in `call_libsbml` because `None` return is valid
        if sbml_param:
            _, _, _, self.mass_units = LibSbmlInterface.parse_parameter(sbml_param)

        LibSbmlInterface.get_commments(self, sbml)

    def import_relations_from_sbml(self, sbml, objs):
        """ Load relationships from SBML compartment

        Args:
            sbml (:obj:`libsbml.Compartment`): SBML compartment
            objs (:obj:`dict`): dictionary that maps WC-Lang types to dictionaries that
                map the ids of WC-Lang objects to WC-Lang objects
        """
        parsed_annots = LibSbmlInterface.parse_annotations(sbml)
        annots = ['biological_type', 'physical_type',
                  'init_volume.distribution', 'init_volume.std',
                  'parent_compartment', 'init_density',
                  'identifiers']
        if 'ph.distribution' in parsed_annots:
            self.ph = Ph()
            annots.extend(['ph.distribution', 'ph.mean', 'ph.std', 'ph.units'])
        LibSbmlInterface.get_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml, objs)

        if 'ph.distribution' in parsed_annots:
            for comp in self.model.compartments:
                if comp.ph and comp.ph.serialize() == self.ph.serialize():
                    self.ph = comp.ph
                    break


class ChemicalStructureFormat(int, CaseInsensitiveEnum):
    """ Format of a chemical structure """
    SMILES = 0
    BpForms = 1


ChemicalStructureAlphabet = CaseInsensitiveEnum('ChemicalStructureAlphabet', list(bpforms.util.get_alphabets().keys()), type=int)
# :obj:`CaseInsensitiveEnum`: Alphabet of a BpForms-encoded chemical structure


class ChemicalStructure(obj_model.Model, SbmlModelMixin):
    """ Structure of a chemical compound

    Attributes:
        value (:obj:`str`)
        format (:obj:`ChemicalStructureFormat`): format of the structure
        alphabet (:obj:`ChemicalStructureAlphabet`): alphabet of BpForms-encoded structure

        empirical_formula (:obj:`EmpiricalFormula`): empirical formula
        molecular_weight (:obj:`float`): molecular weight
        charge (:obj:`int`): charge
    """
    value = LongStringAttribute()
    format = EnumAttribute(ChemicalStructureFormat, none=True)
    alphabet = EnumAttribute(ChemicalStructureAlphabet, none=True)

    empirical_formula = obj_model.chem.EmpiricalFormulaAttribute()
    molecular_weight = FloatAttribute(min=0)
    charge = IntegerAttribute()

    class Meta(obj_model.Model.Meta):
        tabular_orientation = TabularOrientation.multiple_cells
        unique_together = (('value', 'format', 'alphabet',
                            'empirical_formula', 'molecular_weight', 'charge',), )
        attribute_order = ('value', 'format', 'alphabet',
                           'empirical_formula', 'molecular_weight', 'charge',)
        children = {
            'submodel': (),
            'core_model': (),
        }
        child_attrs = {
            'sbml': ('value', 'format', 'alphabet', 'empirical_formula', 'molecular_weight', 'charge',),
            'wc_sim': ('molecular_weight', 'charge'),
        }

    def get_structure(self):
        """ Get structure

        Returns:
            :obj:`openbabel.OBMol` of :obj:`bpforms.BpForm`: structure

        Raises:
            :obj:`ValueError`: if the structure cannot be parsed
        """
        if self.format is None:
            return None

        if self.format == ChemicalStructureFormat.SMILES:
            mol = openbabel.OBMol()
            conv = openbabel.OBConversion()
            assert conv.SetInFormat('smi')
            conv.ReadString(mol, self.value)
            return mol

        if self.format == ChemicalStructureFormat.BpForms and self.alphabet is not None:
            return bpforms.util.get_form(self.alphabet.name)().from_str(self.value)

        raise ValueError('Unsupported format {}'.format(str(self.format)))

    def validate(self):
        """ Check that the structure is valid

        * Format provided when structure is not None
        * Value provided when format is not None
        * Alphabet provided for BpForms-encoded structures
        * Empirical formula, molecular weight, charge match structure (when)

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(ChemicalStructure, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        if self.value and self.format is None:
            errors.append(InvalidAttribute(self.Meta.attributes['format'],
                                           ['A format must be defined for the structure']))

        if not self.value and self.format is not None:
            errors.append(InvalidAttribute(self.Meta.attributes['value'],
                                           ['The format should be None for None structures']))

        if self.format == ChemicalStructureFormat.BpForms and self.alphabet is None:
            errors.append(InvalidAttribute(self.Meta.attributes['alphabet'], [
                          'An alphabet must be defined for BpForms-encoded structures']))

        if self.format != ChemicalStructureFormat.BpForms and self.alphabet is not None:
            errors.append(InvalidAttribute(self.Meta.attributes['alphabet'], [
                          'An alphabet can only be defined for BpForms-encoded structures']))

        if self.value:
            exp_formula = None
            exp_charge = None
            if self.format == ChemicalStructureFormat.SMILES:
                mol = self.get_structure()
                exp_formula = OpenBabelUtils.get_formula(mol)
                exp_charge = mol.GetTotalCharge()
            elif self.format == ChemicalStructureFormat.BpForms and self.alphabet is not None:
                form = self.get_structure()
                exp_formula = form.get_formula()
                exp_charge = form.get_charge()

            if self.empirical_formula is not None and self.empirical_formula != exp_formula:
                errors.append(InvalidAttribute(self.Meta.attributes['empirical_formula'],
                                               ['Empirical formula does not match structure {} != {}'.format(
                                                str(self.empirical_formula), str(exp_formula))]))
            if self.charge is not None and self.charge != exp_charge:
                errors.append(InvalidAttribute(self.Meta.attributes['charge'],
                                               ['Charge does not match structure {} != {}'.format(
                                                self.charge, exp_charge)]))

            if not errors:
                self.empirical_formula = exp_formula
                self.charge = exp_charge

        if self.empirical_formula:
            exp_mol_wt = self.empirical_formula.get_molecular_weight()
            if self.molecular_weight is not None \
                    and not isnan(self.molecular_weight) \
                    and abs(self.molecular_weight - exp_mol_wt) / exp_mol_wt > 1e-3:
                errors.append(InvalidAttribute(self.Meta.attributes['molecular_weight'],
                                               ['Molecular weight does not match structure {} != {}'.format(
                                                self.molecular_weight, exp_mol_wt)]))
            else:
                self.molecular_weight = exp_mol_wt

        if errors:
            return InvalidObject(self, errors)
        return None

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return '__'.join([str(self.value), str(self.format), str(self.alphabet),
                          str(self.empirical_formula), str(self.molecular_weight), str(self.charge)])

    def has_carbon(self):
        """ Returns `True` is species contains at least one carbon atom.

        Returns:
            :obj:`bool`: `True` is species contains at least one carbon atom.
        """
        return self.empirical_formula and self.empirical_formula['C'] > 0


class SpeciesType(obj_model.Model, SbmlModelMixin):
    """ Species type

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        structure (:obj:`ChemicalStructure`): structure (InChI for metabolites; sequence for DNA, RNA, proteins)        
        type (:obj:`pronto.term.Term`): type
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:

        * species (:obj:`list` of :obj:`Species`): species
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='species_types')
    structure = ManyToOneAttribute(ChemicalStructure, related_name='species_types')
    type = OntologyAttribute(onto,
                             namespace='WC',
                             terms=onto['WC:species_type'].rchildren(),
                             default=onto['WC:metabolite'],
                             none=True)
    identifiers = IdentifierManyToManyAttribute(related_name='species_types', verbose_related_name='species types')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='species_types')
    conclusions = ManyToManyAttribute('Conclusion', related_name='species_types')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='species_types')

    class Meta(obj_model.Model.Meta):
        verbose_name = 'Species type'
        attribute_order = ('id', 'name', 'structure', 'type',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        indexed_attrs_tuples = (('id',), )
        children = {
            'submodel': ('structure', 'identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('structure', 'species', 'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'structure', 'type',
                     'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'structure'),
        }


class Species(obj_model.Model, SbmlModelMixin):
    """ Species (tuple of species type, compartment)

    Attributes:
        id (:obj:`str`): identifier equal to `{species_type.id}[{compartment.id}]`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        species_type (:obj:`SpeciesType`): species type
        compartment (:obj:`Compartment`): compartment
        units (:obj:`unit_registry.Unit`): units of counts
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:

        * distribution_init_concentration (:obj:`DistributionInitConcentration`):
          distribution of initial concentration
        * species_coefficients (:obj:`list` of :obj:`SpeciesCoefficient`): participations in reactions and observables
        * rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        * observable_expressions (:obj:`list` of :obj:`ObservableExpression`): observable expressions
        * stop_condition_expressions (:obj:`list` of :obj:`StopConditionExpression`): stop condition expressions
        * function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        * dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='species')
    species_type = ManyToOneAttribute(SpeciesType, related_name='species', min_related=1)
    compartment = ManyToOneAttribute(Compartment, related_name='species', min_related=1)
    units = UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('molecule'),),
                          default=unit_registry.parse_units('molecule'))
    identifiers = IdentifierManyToManyAttribute(related_name='species')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='species')
    conclusions = ManyToManyAttribute('Conclusion', related_name='species')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='species')

    class Meta(obj_model.Model.Meta, ExpressionDynamicTermMeta):
        attribute_order = ('id', 'name', 'species_type', 'compartment', 'units',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        frozen_columns = 1
        # unique_together = (('species_type', 'compartment', ), )
        indexed_attrs_tuples = (('species_type', 'compartment'), )
        expression_term_token_pattern = (token.NAME, token.LSQB, token.NAME, token.RSQB)
        expression_term_units = 'units'
        children = {
            'submodel': ('species_type', 'compartment', 'distribution_init_concentration',
                         'identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('species_type', 'compartment', 'distribution_init_concentration',
                           'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'species_type', 'compartment', 'units',
                     'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'species_type', 'compartment', 'units'),
        }

    def gen_id(self):
        """ Generate identifier

        Returns:
            :obj:`str`: identifier
        """
        return self._gen_id(self.species_type.id, self.compartment.id)

    @staticmethod
    def _gen_id(species_type_id, compartment_id):
        """ Generate identifier

        Args:
            species_type_id (:obj:`str`): species type id
            compartment_id (:obj:`str`): compartment id

        Returns:
            :obj:`str`: identifier
        """
        return '{}[{}]'.format(species_type_id, compartment_id)

    @classmethod
    def parse_id(cls, id):
        """
        Args:
            id (:obj:`str`): identifier

        Returns:
            :obj:`str`: species type id
            :obj:`str`: compartment id

        Raises:
            :obj:`ValueError`: if the id does not have the format `{species.id}[{compartment.id}]`
        """
        st = cls.species_type.related_class.id.pattern[1:-1]
        comp = cls.compartment.related_class.id.pattern[1:-1]
        match = re.match(r'^(' + st + r')\[(' + comp + r')\]$', id)
        if not match:
            raise ValueError('{} is not a valid id')

        return (match.group(1), match.group(2))

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

    def export_to_sbml(self, sbml_model):
        """ Add this species to a SBML model.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`libsbml.Species`: SBML species
        """
        sbml = call_libsbml(sbml_model.createSpecies)
        sbml.initDefaults()  # isn't wrapped in call_libsbml because it returns None

        # id, name
        call_libsbml(sbml.setIdAttribute, self.gen_sbml_id())
        call_libsbml(sbml.setName, self.name)

        # initial concentration
        if self.distribution_init_concentration:
            if are_units_equivalent(self.distribution_init_concentration.units, unit_registry.parse_units('molecule'),
                                    check_same_magnitude=True):
                init_amount = self.distribution_init_concentration.mean
            else:
                prefix = ((1. * self.distribution_init_concentration.units) /
                          (1. * unit_registry.parse_units('M'))).to_base_units().magnitude
                init_amount = self.distribution_init_concentration.mean \
                    * self.compartment.init_volume.mean \
                    * scipy.constants.Avogadro \
                    * prefix
            call_libsbml(sbml.setInitialAmount, init_amount)

        # units
        LibSbmlInterface.set_unit(sbml.setSubstanceUnits, self.units)
        LibSbmlInterface.call_libsbml(sbml.setHasOnlySubstanceUnits, True)

        # comments
        LibSbmlInterface.set_commments(self, sbml)

        return sbml

    def export_relations_to_sbml(self, sbml_model, sbml):
        """ Add relationships to/from object to SBML species.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            sbml (:obj:`libsbml.Species`): SBML species
        """
        # compartment
        call_libsbml(sbml.setCompartment, self.compartment.gen_sbml_id())

        # species type, initial concentration, identifiers
        annots = ['species_type.id', 'species_type.name',
                  'species_type.type', 'species_type.identifiers',
                  'species_type.comments',
                  'identifiers']

        if self.species_type.structure:
            annots.extend(['species_type.structure.value', 'species_type.structure.format',
                           'species_type.structure.alphabet', 'species_type.structure.empirical_formula',
                           'species_type.structure.molecular_weight', 'species_type.structure.charge'])

        if self.distribution_init_concentration:
            annots.extend(['distribution_init_concentration.id',
                           'distribution_init_concentration.name',
                           'distribution_init_concentration.distribution',
                           'distribution_init_concentration.mean',
                           'distribution_init_concentration.std',
                           'distribution_init_concentration.units',
                           'distribution_init_concentration.identifiers',
                           'distribution_init_concentration.comments'])

        LibSbmlInterface.set_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml)

    def import_from_sbml(self, sbml):
        """ Load from SBML species

        Args:
            sbml (:obj:`libsbml.Species`): SBML species
        """
        # id, name
        self.id = self.parse_sbml_id(call_libsbml(sbml.getIdAttribute))
        self.name = call_libsbml(sbml.getName)

        # units
        self.units = LibSbmlInterface.get_unit(sbml.getSubstanceUnits)

        # comments
        LibSbmlInterface.get_commments(self, sbml)

    def import_relations_from_sbml(self, sbml, objs):
        """ Load relationships from SBML species

        Args:
            sbml (:obj:`libsbml.Species`): SBML species
            objs (:obj:`dict`): dictionary that maps WC-Lang types to dictionaries that
                map the ids of WC-Lang objects to WC-Lang objects
        """
        self.compartment = objs[Compartment][Compartment.parse_sbml_id(call_libsbml(sbml.getCompartment, ))]

        # identifiers
        parsed_annots = LibSbmlInterface.parse_annotations(sbml)
        annots = []

        annots.extend(['identifiers', 'species_type.identifiers'])

        # species type
        self.species_type = self.model.species_types.get_or_create(
            id=parsed_annots['species_type.id'])
        annots.extend(['species_type.name',
                       'species_type.type', 'species_type.comments'])

        structure_annots = ['species_type.structure.value', 'species_type.structure.format',
                            'species_type.structure.alphabet', 'species_type.structure.empirical_formula',
                            'species_type.structure.molecular_weight', 'species_type.structure.charge']
        if set(parsed_annots).intersection(set(structure_annots)):
            structure = self.species_type.structure = ChemicalStructure()
            annots.extend(structure_annots)
        else:
            structure = None

        # initial concentration
        if call_libsbml(sbml.isSetInitialAmount):
            self.distribution_init_concentration = self.model.distribution_init_concentrations.create()
            annots.extend(['distribution_init_concentration.id',
                           'distribution_init_concentration.name',
                           'distribution_init_concentration.distribution',
                           'distribution_init_concentration.mean',
                           'distribution_init_concentration.std',
                           'distribution_init_concentration.units',
                           'distribution_init_concentration.identifiers',
                           'distribution_init_concentration.comments'])

        LibSbmlInterface.get_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml, objs)

        if structure:
            for st in self.model.species_types:
                if st != self.species_type and st.structure and st.structure.serialize() == structure.serialize():
                    self.species_type.structure = st.structure
                    break


class DistributionInitConcentration(obj_model.Model, SbmlModelMixin):
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
        units (:obj:`unit_registry.Unit`): units; default units is `M`
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='distribution_init_concentrations',
                               verbose_related_name='Initial species concentrations')
    species = OneToOneAttribute(Species, min_related=1, related_name='distribution_init_concentration',
                                verbose_related_name='Initial species concentration')
    distribution = OntologyAttribute(onto,
                                     namespace='WC',
                                     terms=onto['WC:random_distribution'].rchildren(),
                                     default=onto['WC:normal_distribution'])
    mean = FloatAttribute(min=0)
    std = FloatAttribute(min=0, verbose_name='Standard deviation')
    units = UnitAttribute(unit_registry,
                          choices=(
                              unit_registry.parse_units('molecule'),
                              unit_registry.parse_units('M'),
                              unit_registry.parse_units('mM'),
                              unit_registry.parse_units('uM'),
                              unit_registry.parse_units('nM'),
                              unit_registry.parse_units('pM'),
                              unit_registry.parse_units('fM'),
                              unit_registry.parse_units('aM'),
                          ),
                          default=unit_registry.parse_units('M'))
    identifiers = IdentifierManyToManyAttribute(related_name='distribution_init_concentrations',
                                                verbose_related_name='Initial species concentrations')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='distribution_init_concentrations',
                                           verbose_related_name='Initial species concentrations')
    conclusions = ManyToManyAttribute('Conclusion', related_name='distribution_init_concentrations',
                                      verbose_related_name='Initial species concentrations')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='distribution_init_concentrations',
                                     verbose_related_name='Initial species concentrations')

    class Meta(obj_model.Model.Meta):
        # unique_together = (('species', ), )
        attribute_order = ('id', 'name', 'species',
                           'distribution', 'mean', 'std', 'units',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        verbose_name = 'Initial species concentration'
        frozen_columns = 1
        children = {
            'submodel': ('identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('species',
                           'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'species', 'distribution', 'mean', 'std', 'units', 'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'species', 'distribution', 'mean', 'std', 'units'),
        }

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


class ObservableExpression(obj_model.Model, Expression, SbmlModelMixin):
    """ A mathematical expression of Observables and Species

    The expression used by a `Observable`.

    Attributes:
        expression (:obj:`str`): mathematical expression for an Observable
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        species (:obj:`list` of :obj:`Species`): Species used by this Observable expression
        observables (:obj:`list` of :obj:`Observable`): other Observables used by this Observable expression

    Related attributes:

        * observable (:obj:`Observable`): observable
    """

    expression = LongStringAttribute(primary=True, unique=True, default='')
    species = ManyToManyAttribute(Species, related_name='observable_expressions')
    observables = ManyToManyAttribute('Observable', related_name='observable_expressions')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        tabular_orientation = TabularOrientation.cell
        expression_term_models = ('Species', 'Observable')
        expression_is_linear = True
        expression_unit_registry = unit_registry
        children = {
            'submodel': ('species', 'observables'),
            'core_model': ('species', 'observables'),
        }
        child_attrs = {
            'sbml': ('expression', 'species', 'observables'),
            'wc_sim': ('expression', 'species', 'observables'),
        }

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

    def merge_attrs(self, other, other_objs_in_self, self_objs_in_other):
        """ Merge attributes of two objects

        Args:
            other (:obj:`obj_model.Model`): other model
            other_objs_in_self (:obj:`dict`): dictionary that maps instances of objects in another model to objects
                in a model
            self_objs_in_other (:obj:`dict`): dictionary that maps instances of objects in a model to objects
                in another model
        """
        super(ObservableExpression, self).merge_attrs(other, other_objs_in_self, self_objs_in_other)
        Expression.merge_attrs(self, other, other_objs_in_self, self_objs_in_other)


class Observable(obj_model.Model, SbmlAssignmentRuleMixin):
    """ Observable: a linear function of other Observbles and Species

    Attributes:
        id (:obj:`str`): unique id
        name (:obj:`str`): name
        model (:obj:`Model`): model
        expression (:obj:`ObservableExpression`): mathematical expression for an Observable
        units (:obj:`unit_registry.Unit`): units of expression
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:

        * observable_expressions (:obj:`list` of :obj:`ObservableExpression`): observable expressions
        * function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        * rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        * stop_condition_expressions (:obj:`list` of :obj:`StopConditionExpression`): stop condition expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='observables')
    expression = ExpressionOneToOneAttribute(ObservableExpression, related_name='observable',
                                             min_related=1, min_related_rev=1)
    units = UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('molecule'),),
                          default=unit_registry.parse_units('molecule'))
    identifiers = IdentifierManyToManyAttribute(related_name='observables')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='observables')
    conclusions = ManyToManyAttribute('Conclusion', related_name='observables')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='observables')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'name', 'expression', 'units',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        expression_term_model = ObservableExpression
        expression_term_units = 'units'
        children = {
            'submodel': ('expression', 'identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('expression', 'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'expression', 'units', 'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'expression', 'units'),
        }


class FunctionExpression(obj_model.Model, Expression, SbmlModelMixin):
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

        * function (:obj:`Function`): function
    """
    expression = LongStringAttribute(primary=True, unique=True, default='')
    parameters = ManyToManyAttribute('Parameter', related_name='function_expressions')
    species = ManyToManyAttribute(Species, related_name='function_expressions')
    observables = ManyToManyAttribute(Observable, related_name='function_expressions')
    functions = ManyToManyAttribute('Function', related_name='function_expressions')
    compartments = ManyToManyAttribute(Compartment, related_name='function_expressions')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        tabular_orientation = TabularOrientation.cell
        expression_term_models = ('Parameter', 'Species', 'Observable', 'Function', 'Compartment')
        expression_unit_registry = unit_registry
        children = {
            'submodel': ('parameters', 'species', 'observables', 'functions', 'compartments'),
            'core_model': ('parameters', 'species', 'observables', 'functions', 'compartments'),
        }
        child_attrs = {
            'sbml': ('expression', 'parameters', 'species', 'observables', 'functions', 'compartments'),
            'wc_sim': ('expression', 'parameters', 'species', 'observables', 'functions', 'compartments'),
        }

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

    def merge_attrs(self, other, other_objs_in_self, self_objs_in_other):
        """ Merge attributes of two objects

        Args:
            other (:obj:`obj_model.Model`): other model
            other_objs_in_self (:obj:`dict`): dictionary that maps instances of objects in another model to objects
                in a model
            self_objs_in_other (:obj:`dict`): dictionary that maps instances of objects in a model to objects
                in another model
        """
        super(FunctionExpression, self).merge_attrs(other, other_objs_in_self, self_objs_in_other)
        Expression.merge_attrs(self, other, other_objs_in_self, self_objs_in_other)


class Function(obj_model.Model, SbmlAssignmentRuleMixin):
    """ Function: a mathematical expression of Functions, Observbles, Parameters and Python functions

    Attributes:
        id (:obj:`str`): unique id
        name (:obj:`str`): name
        model (:obj:`Model`): model
        expression (:obj:`FunctionExpression`): mathematical expression for a Function
        units (:obj:`unit_registry.Unit`): units
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:

        * function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        * rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        * stop_condition_expressions (:obj:`list` of :obj:`StopConditionExpression`): stop condition expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='functions')
    expression = ExpressionOneToOneAttribute(FunctionExpression, related_name='function',
                                             min_related=1, min_related_rev=1)
    units = UnitAttribute(unit_registry)
    identifiers = IdentifierManyToManyAttribute(related_name='functions')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='functions')
    conclusions = ManyToManyAttribute('Conclusion', related_name='functions')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='functions')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'name', 'expression', 'units',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        expression_term_model = FunctionExpression
        expression_term_units = 'units'
        children = {
            'submodel': ('expression', 'identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('expression', 'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'expression', 'units', 'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'expression', 'units'),
        }

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
            try:
                calc = self.expression._parsed_expression.test_eval(with_units=True)
            except ParsedExpressionError as error:
                errors.append(InvalidAttribute(self.Meta.attributes['units'], [str(error)]))
            else:
                if hasattr(calc, 'units'):
                    calc_units = calc.units or unit_registry.parse_units('dimensionless')
                else:
                    calc_units = unit_registry.parse_units('dimensionless')

                exp_units = self.units
                if not are_units_equivalent(exp_units, calc_units, check_same_magnitude=True):
                    errors.append(InvalidAttribute(self.Meta.attributes['units'],
                                                   ['Units of "{}" should be "{}" not "{}"'.format(
                                                    self.expression.expression, str(exp_units), str(calc_units))]))

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
        species (:obj:`list` of :obj:`Species`): Species used by this stop condition expression
        observables (:obj:`list` of :obj:`Observable`): Observables used by this stop condition expression
        parameters (:obj:`list` of :obj:`Parameter`): Parameters used by this stop condition expression
        functions (:obj:`list` of :obj:`Function`): Functions used by this stop condition expression
        compartments (:obj:`list` of :obj:`Compartment`): Compartments used by this stop condition expression

    Related attributes:

        * stop_condition (:obj:`StopCondition`): stop condition
    """

    expression = LongStringAttribute(primary=True, unique=True, default='')
    parameters = ManyToManyAttribute('Parameter', related_name='stop_condition_expressions')
    species = ManyToManyAttribute(Species, related_name='stop_condition_expressions')
    observables = ManyToManyAttribute(Observable, related_name='stop_condition_expressions')
    functions = ManyToManyAttribute(Function, related_name='stop_condition_expressions')
    compartments = ManyToManyAttribute(Compartment, related_name='stop_condition_expressions')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        tabular_orientation = TabularOrientation.cell
        expression_term_models = ('Parameter', 'Species', 'Observable', 'Function', 'Compartment')
        expression_type = bool
        expression_unit_registry = unit_registry
        children = {
            'submodel': ('parameters', 'species', 'observables', 'functions', 'compartments'),
            'core_model': ('parameters', 'species', 'observables', 'functions', 'compartments'),
        }
        child_attrs = {
            'sbml': (),
            'wc_sim': ('expression', 'parameters', 'species', 'observables', 'functions', 'compartments')
        }

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

    def merge_attrs(self, other, other_objs_in_self, self_objs_in_other):
        """ Merge attributes of two objects

        Args:
            other (:obj:`obj_model.Model`): other model
            other_objs_in_self (:obj:`dict`): dictionary that maps instances of objects in another model to objects
                in a model
            self_objs_in_other (:obj:`dict`): dictionary that maps instances of objects in a model to objects
                in another model
        """
        super(StopConditionExpression, self).merge_attrs(other, other_objs_in_self, self_objs_in_other)
        Expression.merge_attrs(self, other, other_objs_in_self, self_objs_in_other)


class StopCondition(obj_model.Model):
    """ StopCondition: Simulation of a model terminates when one of its stop conditions is true.

    A Boolean expression of Functions, Observbles, Parameters and Python functions. Stop conditions
    are optional.

    Attributes:
        id (:obj:`str`): unique id
        name (:obj:`str`): name
        model (:obj:`Model`): model
        expression (:obj:`StopConditionExpression`): mathematical expression for a StopCondition
        units (:obj:`unit_registry.Unit`): units
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:

        * expressions (:obj:`Expressions`): expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='stop_conditions')
    expression = ExpressionOneToOneAttribute(StopConditionExpression, related_name='stop_condition',
                                             min_related=1, min_related_rev=1)
    units = UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('dimensionless'),),
                          default=unit_registry.parse_units('dimensionless'))
    identifiers = IdentifierManyToManyAttribute(related_name='stop_conditions')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='stop_conditions')
    conclusions = ManyToManyAttribute('Conclusion', related_name='stop_conditions')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='stop_conditions')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'name', 'expression', 'units',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        expression_term_model = StopConditionExpression
        expression_term_units = 'units'
        children = {
            'submodel': ('expression', 'identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('expression', 'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': (),
            'wc_sim': ('id', 'model', 'expression', 'units')
        }

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

        # return errors
        if errors:
            return InvalidObject(self, errors)
        return None


class FluxBounds(obj_model.Model, SbmlModelMixin):
    """ Flux bounds 

    Attributes:
        min (:obj:`float`): minimum flux bound for solving an FBA model; negative for reversible reactions
        max (:obj:`float`): maximum flux bound for solving an FBA model
        units (:obj:`unit_registry.Unit`): units for the minimum and maximum fluxes    

    Related attributes:

        * reactions (:obj:`list` of :obj:`Reaction`): reactions
    """
    min = FloatAttribute(nan=True, verbose_name='Minimum')
    max = FloatAttribute(min=0, nan=True, verbose_name='Maximum')
    units = UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('M s^-1'),),
                          default=None, none=True, verbose_name='Units')

    class Meta(obj_model.Model.Meta):
        tabular_orientation = TabularOrientation.multiple_cells
        unique_together = (('min', 'max', 'units'), )
        attribute_order = ('min', 'max', 'units')
        children = {
            'submodel': (),
            'core_model': (),
        }
        child_attrs = {
            'sbml': ('min', 'max', 'units'),
            'wc_sim': ('min', 'max', 'units'),
        }

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return '__'.join([str(self.min), str(self.max), str(self.units)])


class Reaction(obj_model.Model, SbmlModelMixin):
    """ Reaction

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        submodel (:obj:`Submodel`): submodel that reaction belongs to
        participants (:obj:`list` of :obj:`SpeciesCoefficient`): participants
        reversible (:obj:`bool`): indicates if reaction is thermodynamically reversible
        flux_bounds (:obj:`FluxBounds`): flux bounds
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:

        * rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws; if present, rate_laws[0] is the forward
          rate law, and rate_laws[0] is the backward rate law
        * dfba_obj_expression (:obj:`DfbaObjectiveExpression`): dFBA objective expression
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='reactions')
    submodel = ManyToOneAttribute(Submodel, related_name='reactions')
    participants = ReactionParticipantAttribute(related_name='reactions')
    reversible = BooleanAttribute()
    rate_units = UnitAttribute(unit_registry,
                               choices=(unit_registry.parse_units('s^-1'),),
                               default=unit_registry.parse_units('s^-1'))
    flux_bounds = ManyToOneAttribute(FluxBounds, related_name='reactions')
    identifiers = IdentifierManyToManyAttribute(related_name='reactions')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='reactions')
    conclusions = ManyToManyAttribute('Conclusion', related_name='reactions')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='reactions')

    class Meta(obj_model.Model.Meta, ExpressionDynamicTermMeta):
        attribute_order = ('id', 'name', 'submodel',
                           'participants', 'reversible',
                           'rate_units', 'flux_bounds',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        indexed_attrs_tuples = (('id',), )
        expression_term_units = 'rate_units'
        merge = obj_model.ModelMerge.append
        children = {
            'submodel': ('participants', 'rate_laws', 'flux_bounds',
                         'identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('participants', 'rate_laws', 'flux_bounds',
                           'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'submodel', 'participants', 'reversible',
                     'rate_units', 'flux_bounds',
                     'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'submodel', 'participants', 'reversible',
                       'rate_units', 'flux_bounds'),
        }

    def validate(self):
        """ Check if the reaction is valid

        * If the submodel is ODE or SSA, check that the reaction has a forward rate law
        * If the submodel is ODE or SSA and the reaction is reversible, check that the reaction has a
          backward rate law
        * Check flux units are not None if flux_bounds.min or flux_bounds.max is defined
        * Check that `flux_bounds.min` <= `flux_bounds.max`
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
            'WC:ordinary_differential_equations',
            'WC:stochastic_simulation_algorithm',
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
        if self.flux_bounds \
                and (not isnan(self.flux_bounds.min) or not isnan(self.flux_bounds.max)) \
                and self.flux_bounds.units is None:
            errors.append(InvalidAttribute(self.Meta.attributes['flux_bounds'],
                                           ['Units must be defined for the flux bounds']))

        if self.submodel and not are_terms_equivalent(self.submodel.framework, onto['WC:dynamic_flux_balance_analysis']):
            if self.flux_bounds is not None:
                errors.append(InvalidAttribute(self.Meta.attributes['flux_bounds'],
                                               ['Flux bounds should be None reactions in non-dFBA submodels']))

        if self.flux_bounds \
                and not isnan(self.flux_bounds.min) \
                and not isnan(self.flux_bounds.min) \
                and self.flux_bounds.min > self.flux_bounds.max:
            errors.append(InvalidAttribute(self.Meta.attributes['flux_bounds'],
                                           ['Maximum flux must be least the minimum flux']))

        if self.reversible and self.flux_bounds and not isnan(self.flux_bounds.min) and self.flux_bounds.min >= 0:
            errors.append(InvalidAttribute(self.Meta.attributes['flux_bounds'],
                                           ['Minimum flux for reversible reaction should be negative or NaN']))
        if not self.reversible and self.flux_bounds and not isnan(self.flux_bounds.min) and self.flux_bounds.min < 0:
            errors.append(InvalidAttribute(self.Meta.attributes['flux_bounds'],
                                           ['Minimum flux for irreversible reaction should be non-negative']))

        # return errors
        if errors:
            return InvalidObject(self, errors)
        return None

    def get_species(self, __type=None, **kwargs):
        """ Get species

        Args:
            __type (:obj:`types.TypeType` or :obj:`tuple` of :obj:`types.TypeType`): subclass(es) of :obj:`Model`
            kwargs (:obj:`dict` of :obj:`str` --> :obj:`object`): dictionary of attribute name/value pairs to find matching
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

    def export_to_sbml(self, sbml_model):
        """ Add this reaction to a SBML model.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`libsbml.Reaction`: SBML reaction

        Raises:
            :obj:`ValueError`: if the reaction has a backward rate law which cannot be exported to SBML
        """
        # create SBML reaction in SBML document
        sbml_rxn = call_libsbml(sbml_model.createReaction)

        # id, name
        call_libsbml(sbml_rxn.setIdAttribute, self.gen_sbml_id())
        call_libsbml(sbml_rxn.setName, self.name)

        # reversibility
        if self.reversible and self.rate_laws.get_one(direction=RateLawDirection.backward) \
                and not are_terms_equivalent(self.submodel.framework, onto['WC:dynamic_flux_balance_analysis']):
            raise ValueError('Reversible reactions with backward rate laws must be split before export to SBML')
        call_libsbml(sbml_rxn.setReversible, self.reversible)

        # dFBA flux bounds
        if are_terms_equivalent(self.submodel.framework, onto['WC:dynamic_flux_balance_analysis']):
            sbml_plugin = call_libsbml(sbml_rxn.getPlugin, 'fbc')
            flux_bounds_units = FluxBounds.Meta.attributes['units'].choices[0]

            if self.flux_bounds:
                bounds = [('Lower', self.flux_bounds.min, -1.), ('Upper', self.flux_bounds.max, 1.)]
            else:
                bounds = [('Lower', float('nan'), -1.), ('Upper', float('nan'), 1.)]

            for bound, value, sense in bounds:
                if isnan(value) or value is None:
                    value = sense * float('inf')

                param_id = "__Reaction__Flux{}Bound__{}".format(bound, self.gen_sbml_id())
                param = LibSbmlInterface.create_parameter(sbml_model, param_id, value,
                                                          flux_bounds_units)
                call_libsbml(getattr(sbml_plugin, 'set' + bound + 'FluxBound'), param_id)

        # identifiers, comments
        LibSbmlInterface.set_commments(self, sbml_rxn)

        # forward rate law
        rl = self.rate_laws.get_one(direction=RateLawDirection.forward)
        if rl:
            rl.export_to_sbml(sbml_model)  # because law = libsbml.KineticLaw(3, 2); reaction.setKineticLaw(law); doesn't work

        # return reaction
        return sbml_rxn

    def export_relations_to_sbml(self, sbml_model, sbml_rxn):
        """ Add relationships to/from object to SBML reaction.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            sbml_rxn (:obj:`libsbml.Reaction`): SBML reaction
        """
        annots = {}

        # participants
        for participant in self.participants:
            if participant.coefficient < 0:
                sbml_part = call_libsbml(sbml_rxn.createReactant)
                call_libsbml(sbml_part.setStoichiometry, -participant.coefficient)
            elif 0 < participant.coefficient:
                sbml_part = call_libsbml(sbml_rxn.createProduct)
                call_libsbml(sbml_part.setStoichiometry, participant.coefficient)
            call_libsbml(sbml_part.setSpecies, participant.species.gen_sbml_id())
            call_libsbml(sbml_part.setConstant, True)

        # rate_units
        annots['rate_units'] = 'rate_units'

        # identifiers
        annots['identifiers'] = 'identifiers'

        # forward rate law
        rl = self.rate_laws.get_one(direction=RateLawDirection.forward)
        if rl:
            modifiers = set(rl.expression.species).difference(set(part.species for part in self.participants))
            for modifier in modifiers:
                call_libsbml(sbml_rxn.addModifier, call_libsbml(sbml_model.getSpecies, modifier.gen_sbml_id()))

        # backward rate law
        rl = self.rate_laws.get_one(direction=RateLawDirection.backward)
        if rl:
            annots['rate_laws.backward.id'] = (('rate_laws', {'direction': rl.direction}), 'id')
            annots['rate_laws.backward.name'] = (('rate_laws', {'direction': rl.direction}), 'name')
            annots['rate_laws.backward.type'] = (('rate_laws', {'direction': rl.direction}), 'type')
            annots['rate_laws.backward.expression'] = (('rate_laws', {'direction': rl.direction}), 'expression')
            annots['rate_laws.backward.units'] = (('rate_laws', {'direction': rl.direction}), 'units')
            annots['rate_laws.backward.identifiers'] = (('rate_laws', {'direction': rl.direction}), 'identifiers')
            annots['rate_laws.backward.comments'] = (('rate_laws', {'direction': rl.direction}), 'comments')

        # annotations
        LibSbmlInterface.set_annotations(self, annots, sbml_rxn)

    def import_from_sbml(self, sbml_rxn):
        """ Load from SBML reaction

        Args:
            sbml (:obj:`libsbml.Reaction`): SBML reaction
        """
        # id, name
        self.id = self.parse_sbml_id(call_libsbml(sbml_rxn.getIdAttribute))
        self.name = call_libsbml(sbml_rxn.getName)

        # reversibility
        self.reversible = call_libsbml(sbml_rxn.getReversible)

        # dFBA flux bounds
        sbml_doc = call_libsbml(sbml_rxn.getSBMLDocument)
        if LibSbmlInterface.call_libsbml(sbml_doc.isSetPackageRequired, 'fbc'):
            sbml_model = call_libsbml(sbml_rxn.getModel)
            sbml_plugin = call_libsbml(sbml_rxn.getPlugin, 'fbc')
            flux_bounds = FluxBounds()
            has_flux_bounds = False
            for bound, attr_name, sense in [('Lower', 'min', -1.), ('Upper', 'max', 1.)]:
                param_id = call_libsbml(getattr(sbml_plugin, 'get' + bound + 'FluxBound'))
                param = call_libsbml(sbml_model.getParameter, param_id)
                if param:
                    _, _, val, flux_bounds.units = LibSbmlInterface.parse_parameter(param)
                    if isinf(val):
                        val = float('nan')
                    else:
                        has_flux_bounds = True
                    setattr(flux_bounds, attr_name, val)
            if has_flux_bounds:
                for rxn in self.model.reactions:
                    if rxn.flux_bounds and rxn.flux_bounds.serialize() == flux_bounds.serialize():
                        flux_bounds = rxn.flux_bounds
                        break
                self.flux_bounds = flux_bounds

        # comments
        LibSbmlInterface.get_commments(self, sbml_rxn)

        # forward rate law
        if call_libsbml(sbml_rxn.isSetKineticLaw):
            rl = self.rate_laws.get_or_create(direction=RateLawDirection.forward)
            rl.model = self.model
            rl.import_from_sbml(call_libsbml(sbml_rxn.getKineticLaw))

    def import_relations_from_sbml(self, sbml_rxn, objs):
        """ Load relationships from SBML reaction

        Args:
            sbml (:obj:`libsbml.Reaction`): SBML reaction
            objs (:obj:`dict`): dictionary that maps WC-Lang types to dictionaries that
                map the ids of WC-Lang objects to WC-Lang objects
        """
        parsed_annots = LibSbmlInterface.parse_annotations(sbml_rxn)
        annots = {}

        # submodel
        self.submodel = self.model.submodels[0]

        # participants
        for num_func, get_func, sense in [(sbml_rxn.getNumReactants, sbml_rxn.getReactant, -1),
                                          (sbml_rxn.getNumProducts, sbml_rxn.getProduct, 1)]:
            for i_part in range(call_libsbml(num_func, returns_int=True)):
                sbml_part = call_libsbml(get_func, i_part)
                species = objs[Species][Species.parse_sbml_id(call_libsbml(sbml_part.getSpecies))]
                coeff = sense * call_libsbml(sbml_part.getStoichiometry)
                part = species.species_coefficients.get_or_create(coefficient=coeff)
                self.participants.append(part)

        # rate units
        annots['rate_units'] = 'rate_units'

        # identifiers
        annots['identifiers'] = 'identifiers'

        # forward rate law
        rl = self.rate_laws.get_one(direction=RateLawDirection.forward)
        if rl:
            rl.import_relations_from_sbml(call_libsbml(sbml_rxn.getKineticLaw), objs)

        # backward rate law
        if 'rate_laws.backward.id' in parsed_annots:
            rl = self.rate_laws.get_or_create(direction=RateLawDirection.backward)
            rl.model = self.model
            annots['rate_laws.backward.id'] = (('rate_laws', {'direction': rl.direction}), 'id')
            annots['rate_laws.backward.name'] = (('rate_laws', {'direction': rl.direction}), 'name')
            annots['rate_laws.backward.type'] = (('rate_laws', {'direction': rl.direction}), 'type')
            annots['rate_laws.backward.expression'] = (('rate_laws', {'direction': RateLawDirection.backward}), 'expression')
            annots['rate_laws.backward.units'] = (('rate_laws', {'direction': rl.direction}), 'units')
            annots['rate_laws.backward.comments'] = (('rate_laws', {'direction': rl.direction}), 'comments')
            annots['rate_laws.backward.identifiers'] = (('rate_laws', {'direction': rl.direction}), 'identifiers')

        # get annotations
        LibSbmlInterface.get_annotations(self, annots, sbml_rxn, objs)


class SpeciesCoefficient(obj_model.Model, SbmlModelMixin):
    """ A tuple of a species and a coefficient

    Attributes:
        species (:obj:`Species`): species
        coefficient (:obj:`float`): coefficient

    Related attributes:

        * reaction (:obj:`Reaction`): reaction
        * observables (:obj:`Observable`): observables
    """
    species = ManyToOneAttribute(Species, related_name='species_coefficients')
    coefficient = FloatAttribute(nan=False)

    class Meta(obj_model.Model.Meta):
        unique_together = (('species', 'coefficient'),)
        attribute_order = ('species', 'coefficient')
        frozen_columns = 1
        tabular_orientation = TabularOrientation.cell
        ordering = ('species', 'coefficient')
        children = {
            'submodel': ('species',),
            'core_model': ('species',),
        }
        child_attrs = {
            'sbml': ('species', 'coefficient'),
            'wc_sim': ('species', 'coefficient'),
        }

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
                species_id = Species._gen_id(match.group(5), compartment.get_primary_attribute())
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


class RateLawExpression(obj_model.Model, Expression, SbmlModelMixin):
    """ Rate law expression

    Attributes:
        expression (:obj:`str`): mathematical expression of the rate law
        _parsed_expression (:obj:`ParsedExpression`): an analyzed `expression`; not an `obj_model.Model`
        species (:obj:`list` of :obj:`Species`): species whose dynamic concentrations are used in the rate law
        parameters (:obj:`list` of :obj:`Parameter`): parameters whose values are used in the rate law
        compartments (:obj:`list` of :obj:`Compartment`): Compartments used by this stop condition expression

    Related attributes:

        * rate_law (:obj:`RateLaw`): the `RateLaw` which uses this `RateLawExpression`
    """
    expression = LongStringAttribute(primary=True, unique=True, default='')
    parameters = ManyToManyAttribute('Parameter', related_name='rate_law_expressions')
    species = ManyToManyAttribute(Species, related_name='rate_law_expressions')
    observables = ManyToManyAttribute(Observable, related_name='rate_law_expressions')
    functions = ManyToManyAttribute(Function, related_name='rate_law_expressions')
    compartments = ManyToManyAttribute(Compartment, related_name='rate_law_expressions')

    class Meta(obj_model.Model.Meta, Expression.Meta):
        attribute_order = ('expression', 'species', 'parameters')
        tabular_orientation = TabularOrientation.cell
        ordering = ('expression',)
        expression_term_models = ('Parameter', 'Species', 'Observable', 'Function', 'Compartment')
        expression_unit_registry = unit_registry
        children = {
            'submodel': ('parameters', 'species', 'observables', 'functions', 'compartments'),
            'core_model': ('parameters', 'species', 'observables', 'functions', 'compartments'),
        }
        child_attrs = {
            'sbml': ('expression', 'parameters', 'species', 'observables', 'functions', 'compartments'),
            'wc_sim': ('expression', 'parameters', 'species', 'observables', 'functions', 'compartments'),
        }

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

    def merge_attrs(self, other, other_objs_in_self, self_objs_in_other):
        """ Merge attributes of two objects

        Args:
            other (:obj:`obj_model.Model`): other model
            other_objs_in_self (:obj:`dict`): dictionary that maps instances of objects in another model to objects
                in a model
            self_objs_in_other (:obj:`dict`): dictionary that maps instances of objects in a model to objects
                in another model
        """
        super(RateLawExpression, self).merge_attrs(other, other_objs_in_self, self_objs_in_other)
        Expression.merge_attrs(self, other, other_objs_in_self, self_objs_in_other)


class RateLaw(obj_model.Model, SbmlModelMixin):
    """ Rate law

    Attributes:
        id (:obj:`str`): identifier equal to `{reaction.id}-{direction.name}`
        name (:obj:`str`): name
        model (:obj:`Model`): model
        reaction (:obj:`Reaction`): reaction
        direction (:obj:`RateLawDirection`): direction
        type (:obj:`pronto.term.Term`): type
        expression (:obj:`RateLawExpression`): expression
        units (:obj:`unit_registry.Unit`): units
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = StringAttribute(primary=True, unique=True)
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='rate_laws')
    reaction = ManyToOneAttribute(Reaction, related_name='rate_laws')
    direction = EnumAttribute(RateLawDirection, default=RateLawDirection.forward)
    type = OntologyAttribute(onto,
                             namespace='WC',
                             terms=onto['WC:rate_law'].rchildren(),
                             default=None, none=True)
    expression = ExpressionManyToOneAttribute(RateLawExpression, min_related=1, min_related_rev=1, related_name='rate_laws')
    units = UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('s^-1'),),
                          default=unit_registry.parse_units('s^-1'))
    identifiers = IdentifierManyToManyAttribute(related_name='rate_laws')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='rate_laws')
    conclusions = ManyToManyAttribute('Conclusion', related_name='rate_laws')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='rate_laws')

    class Meta(obj_model.Model.Meta, ExpressionExpressionTermMeta):
        attribute_order = ('id', 'name', 'reaction', 'direction', 'type',
                           'expression', 'units',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        # unique_together = (('reaction', 'direction'), )
        expression_term_model = RateLawExpression
        expression_term_units = 'units'
        children = {
            'submodel': ('expression', 'identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('expression', 'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'reaction', 'direction', 'type', 'expression', 'units', 'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'reaction', 'direction', 'expression', 'units'),
        }

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
            try:
                calc = self.expression._parsed_expression.test_eval(with_units=True)
            except ParsedExpressionError as error:
                errors.append(InvalidAttribute(self.Meta.attributes['units'], [str(error)]))
            else:
                if hasattr(calc, 'units'):
                    calc_units = calc.units or unit_registry.parse_units('dimensionless')
                else:
                    calc_units = unit_registry.parse_units('dimensionless')

                exp_units = self.units
                if not are_units_equivalent(exp_units, calc_units, check_same_magnitude=True):
                    errors.append(InvalidAttribute(self.Meta.attributes['units'],
                                                   ['Units of "{}" should be "{}" not "{}"'.format(
                                                    self.expression.expression, str(exp_units), str(calc_units))]))
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

    def export_to_sbml(self, sbml_model):
        """ Add this rate law to a SBML reaction.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`libsbml.KineticLaw`: SBML kinetic law
        """
        if self.direction != RateLawDirection.forward:
            return

        sbml_rxn = sbml_model.getReaction(self.reaction.gen_sbml_id())
        if not sbml_rxn or sbml_rxn.getKineticLaw():
            return
        sbml = call_libsbml(sbml_rxn.createKineticLaw)

        # id, name
        call_libsbml(sbml.setIdAttribute, self.gen_sbml_id())
        call_libsbml(sbml.setName, self.name)

        # comments
        LibSbmlInterface.set_commments(self, sbml)

        # return SBML kinetic law
        return sbml

    def export_relations_to_sbml(self, sbml_model, sbml):
        """ Add relationships to/from object to SBML kinetic law.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            sbml (:obj:`libsbml.KineticLaw`): SBML kinetic law
        """
        if self.direction != RateLawDirection.forward:
            return

        sbml_rxn = call_libsbml(sbml_model.getReaction, self.reaction.gen_sbml_id())
        sbml = call_libsbml(sbml_rxn.getKineticLaw)

        # expression
        LibSbmlInterface.set_math(sbml.setMath, self.expression, units_transform='({}) * 1 mole')

        # type, units, identifiers
        annots = ['type', 'units', 'identifiers']
        LibSbmlInterface.set_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml)

    def import_from_sbml(self, sbml):
        """ Load from SBML kinetic law

        Args:
            sbml (:obj:`libsbml.KineticLaw`): SBML kinetic law
        """
        annots = []

        # id, name
        self.id = self.parse_sbml_id(call_libsbml(sbml.getIdAttribute))
        self.name = call_libsbml(sbml.getName)

        # type
        annots.append('type')

        # units
        annots.append('units')

        # comments
        LibSbmlInterface.get_commments(self, sbml)

        # annotations
        LibSbmlInterface.get_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml)

    def import_relations_from_sbml(self, sbml, objs):
        """ Load relationships from SBML kinetic law

        Args:
            sbml (:obj:`libsbml.KineticLaw`): SBML kinetic law
            objs (:obj:`dict`): dictionary that maps WC-Lang types to dictionaries that
                map the ids of WC-Lang objects to WC-Lang objects
        """
        # expression
        self.expression = LibSbmlInterface.get_math(sbml.getMath, self.Meta.expression_term_model,
                                                    objs, units_transform=self._import_relations_from_sbml_units_transform)

        # identifiers
        LibSbmlInterface.get_annotations(self, LibSbmlInterface.gen_nested_attr_paths(['identifiers']), sbml, objs)

    @staticmethod
    def _import_relations_from_sbml_units_transform(formula):
        formula = re.sub(r'^(.*?) \* 1 mole$', r'\1', formula)
        if formula[0] == '(' and formula[-1] == ')':
            formula = formula[1:-1]
        return formula


class DfbaObjSpecies(obj_model.Model, SbmlModelMixin):
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
        units (:obj:`unit_registry.Unit`): units of the value
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
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
    units = UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('M s^-1'), unit_registry.parse_units('mol gDCW^-1 s^-1')),
                          default=unit_registry.parse_units('M s^-1'))
    identifiers = IdentifierManyToManyAttribute(related_name='dfba_obj_species',
                                                verbose_related_name='dFBA objective species')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='dfba_obj_species')
    conclusions = ManyToManyAttribute('Conclusion', related_name='dfba_obj_species',
                                      verbose_related_name='dFBA objective species')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='dfba_obj_species',
                                     verbose_related_name='dFBA objective species')

    class Meta(obj_model.Model.Meta):
        # unique_together = (('dfba_obj_reaction', 'species'), )
        attribute_order = ('id', 'name', 'dfba_obj_reaction',
                           'species', 'value', 'units',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        verbose_name = 'dFBA objective species'
        verbose_name_plural = 'dFBA objective species'
        merge = obj_model.ModelMerge.append
        children = {
            'submodel': ('dfba_obj_reaction', 'species', 'identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('dfba_obj_reaction', 'species', 'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'dfba_obj_reaction', 'species', 'value', 'units', 'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'dfba_obj_reaction', 'species', 'value', 'units'),
        }

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
                ((self.units == unit_registry.parse_units('M s^-1') and
                    self.dfba_obj_reaction.cell_size_units != unit_registry.parse_units('l')) or
                 (self.units == unit_registry.parse_units('mol gDCW^-1 s^-1') and
                    self.dfba_obj_reaction.cell_size_units != unit_registry.parse_units('gDCW'))):
            errors.append(InvalidAttribute(
                self.Meta.attributes['units'],
                ['Units {} are not consistent with cell size units {}'.format(
                    str(self.units), str(self.dfba_obj_reaction.cell_size_units))]))

        if errors:
            return InvalidObject(self, errors)
        return None


class DfbaObjReaction(obj_model.Model, SbmlModelMixin):
    """ A pseudo-reaction used to represent the interface between metabolism and other
    cell processes.

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        submodel (:obj:`Submodel`): submodel that uses this reaction
        units (:obj:`unit_registry.Unit`): rate units
        cell_size_units (:obj:`unit_registry.Unit`): cell size units
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:

        * dfba_obj_expression (:obj:`DfbaObjectiveExpression`): dFBA objectie expression
        * dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): the components of this dFBA objective reaction
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='dfba_obj_reactions', verbose_related_name='dFBA objective reactions')
    submodel = ManyToOneAttribute(Submodel, related_name='dfba_obj_reactions', verbose_related_name='dFBA objective reactions')
    units = UnitAttribute(unit_registry,
                          choices=(unit_registry.parse_units('s^-1'),),
                          default=unit_registry.parse_units('s^-1'))
    cell_size_units = UnitAttribute(unit_registry,
                                    choices=(unit_registry.parse_units('l'), unit_registry.parse_units('gDCW')),
                                    default=unit_registry.parse_units('l'))
    identifiers = IdentifierManyToManyAttribute(related_name='dfba_obj_reactions',
                                                verbose_related_name='dFBA objective reactions')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='dfba_obj_reactions')
    conclusions = ManyToManyAttribute('Conclusion', related_name='dfba_obj_reactions',
                                      verbose_related_name='dFBA objective reactions')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='dfba_obj_reactions', verbose_related_name='dFBA objective reactions')

    class Meta(obj_model.Model.Meta, ExpressionDynamicTermMeta):
        attribute_order = ('id', 'name', 'submodel', 'units', 'cell_size_units',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        indexed_attrs_tuples = (('id',), )
        verbose_name = 'dFBA objective reaction'
        expression_term_units = 'units'
        merge = obj_model.ModelMerge.append
        children = {
            'submodel': ('dfba_obj_species', 'identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('dfba_obj_species', 'identifiers', 'evidence', 'conclusions', 'references'),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'submodel', 'units', 'cell_size_units', 'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'submodel', 'units', 'cell_size_units'),
        }

    def export_to_sbml(self, sbml_model):
        """ Add a dFBA objective reaction to a SBML model.

        DfbaObjReactions are added to the SBML model because they can be used in a dFBA submodel's
        objective function. In fact the default objective function is the submodel's dFBA objective reaction.
        Since SBML does not define DfbaObjReaction as a separate class, DfbaObjReactions are added
        to the SBML model as SBML reactions.
        CheckModel ensures that wc_lang DfbaObjReactions and Reactions have distinct ids.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`libsbml.Reaction`: SBML reaction
        """
        # create SBML reaction in SBML document
        sbml_rxn = call_libsbml(sbml_model.createReaction)
        call_libsbml(sbml_rxn.setReversible, False)

        sbml_plugin = call_libsbml(sbml_rxn.getPlugin, 'fbc')
        for bound, value in [('Lower', -float('inf')), ('Upper', float('inf'))]:
            param_id = "__DfbaObjReaction__Flux{}Bound__{}".format(self.id, bound)
            param = LibSbmlInterface.create_parameter(sbml_model, param_id, value, self.units)
            call_libsbml(getattr(sbml_plugin, 'set' + bound + 'FluxBound'), param_id)

        # id, name
        call_libsbml(sbml_rxn.setIdAttribute, self.gen_sbml_id())
        call_libsbml(sbml_rxn.setName, self.name)

        # comments
        LibSbmlInterface.set_commments(self, sbml_rxn)

        # return SBML reaction
        return sbml_rxn

    def export_relations_to_sbml(self, sbml_model, sbml_rxn):
        """ Add relationships to/from object to SBML reaction.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            sbml_rxn (:obj:`libsbml.Reaction`): SBML reaction
        """
        # participants
        for dfba_obj_species in self.dfba_obj_species:
            if dfba_obj_species.value < 0:
                sbml_part = call_libsbml(sbml_rxn.createReactant)
                coeff = -dfba_obj_species.value
            else:
                sbml_part = call_libsbml(sbml_rxn.createProduct)
                coeff = dfba_obj_species.value
            id = dfba_obj_species.species.gen_sbml_id()
            call_libsbml(sbml_part.setIdAttribute, dfba_obj_species.gen_sbml_id())
            # call_libsbml(sbml_part.setName, dfba_obj_species.name) # because libSBML has a bug in SpeciesReference.setName
            call_libsbml(sbml_part.setSpecies, id)
            call_libsbml(sbml_part.setConstant, True)
            call_libsbml(sbml_part.setStoichiometry, coeff)
            LibSbmlInterface.set_annotations(dfba_obj_species, LibSbmlInterface.gen_nested_attr_paths([
                                             'name', 'units', 'identifiers']), sbml_part)
            LibSbmlInterface.set_commments(dfba_obj_species, sbml_part)

        # units, identifiers
        annots = ['units', 'cell_size_units', 'identifiers']
        LibSbmlInterface.set_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml_rxn)

    def import_from_sbml(self, sbml_rxn):
        """ Load from SBML reaction

        Args:
            sbml (:obj:`libsbml.Reaction`): SBML reaction
        """
        annots = []

        # id, name
        self.id = self.parse_sbml_id(call_libsbml(sbml_rxn.getIdAttribute))
        self.name = call_libsbml(sbml_rxn.getName)

        # units
        annots.extend(['units', 'cell_size_units'])

        # comments
        LibSbmlInterface.get_commments(self, sbml_rxn)

        # annotations
        LibSbmlInterface.get_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml_rxn)

    def import_relations_from_sbml(self, sbml_rxn, objs):
        """ Load relationships from SBML reaction

        Args:
            sbml (:obj:`libsbml.Reaction`): SBML reaction
            objs (:obj:`dict`): dictionary that maps WC-Lang types to dictionaries that
                map the ids of WC-Lang objects to WC-Lang objects
        """
        # submodel
        self.submodel = self.model.submodels[0]

        # participants
        for num_func, get_func, sense in [(sbml_rxn.getNumReactants, sbml_rxn.getReactant, -1),
                                          (sbml_rxn.getNumProducts, sbml_rxn.getProduct, 1)]:
            for i_part in range(call_libsbml(num_func, returns_int=True)):
                sbml_part = call_libsbml(get_func, i_part)
                species = objs[Species][Species.parse_sbml_id(call_libsbml(sbml_part.getSpecies))]
                value = sense * call_libsbml(sbml_part.getStoichiometry)

                dfba_obj_species = self.dfba_obj_species.create()
                dfba_obj_species.id = DfbaObjSpecies.parse_sbml_id(call_libsbml(sbml_part.getIdAttribute))
                # dfba_obj_species.name = call_libsbml(sbml_part.getName) # because libSBML has a bug in SpeciesReference.setName
                dfba_obj_species.model = self.model
                dfba_obj_species.species = species
                dfba_obj_species.value = value
                LibSbmlInterface.get_annotations(dfba_obj_species, LibSbmlInterface.gen_nested_attr_paths([
                                                 'name', 'units', 'identifiers']), sbml_part, objs)
                LibSbmlInterface.get_commments(dfba_obj_species, sbml_part)

        # identifiers
        LibSbmlInterface.get_annotations(self, LibSbmlInterface.gen_nested_attr_paths(['identifiers']), sbml_rxn, objs)


class Parameter(obj_model.Model, SbmlModelMixin):
    """ Parameter

    Attributes:
        id (:obj:`str`): unique identifier per model/submodel
        name (:obj:`str`): name
        model (:obj:`Model`): model
        type (:obj:`pronto.term.Term`): parameter type
        value (:obj:`float`): value
        std (:obj:`float`): standard error of the value
        units (:obj:`unit_registry.Unit`): units of the value and standard error
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:

        * density_compartment (:obj:`Compartment`): compartments whose density is represented by the parameter
        * observable_expressions (:obj:`list` of :obj:`ObservableExpression`): observable expressions
        * function_expressions (:obj:`list` of :obj:`FunctionExpression`): function expressions
        * rate_law_expressions (:obj:`list` of :obj:`RateLawExpression`): rate law expressions
        * stop_condition_expressions (:obj:`list` of :obj:`StopConditionExpression`): stop condition expressions
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='parameters')
    type = OntologyAttribute(onto,
                             namespace='WC',
                             terms=onto['WC:parameter'].rchildren(),
                             default=None, none=True)
    value = FloatAttribute()
    std = FloatAttribute(min=0, verbose_name='Standard error')
    units = UnitAttribute(unit_registry)
    identifiers = IdentifierManyToManyAttribute(related_name='parameters')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='parameters')
    conclusions = ManyToManyAttribute('Conclusion', related_name='parameters')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='parameters')

    class Meta(obj_model.Model.Meta, ExpressionStaticTermMeta):
        attribute_order = ('id', 'name', 'type',
                           'value', 'std', 'units',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references')
        children = {
            'submodel': ('identifiers', 'evidence', 'conclusions', 'references'),
            'core_model': ('identifiers', 'evidence', 'conclusions', 'references'),
        }
        expression_term_value = 'value'
        expression_term_units = 'units'
        child_attrs = {
            'sbml': ('id', 'name', 'model', 'type', 'value', 'std', 'units', 'identifiers', 'comments'),
            'wc_sim': ('id', 'model', 'type', 'value', 'std', 'units'),
        }

    def export_to_sbml(self, sbml_model):
        """ Add this parameter to a SBML model.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model

        Returns:
            :obj:`libsbml.Parameter`: SBML parameter
        """
        sbml_id = self.gen_sbml_id()
        sbml = LibSbmlInterface.create_parameter(sbml_model, sbml_id, self.value, self.units,
                                                 name=self.name)

        LibSbmlInterface.set_commments(self, sbml)

        return sbml

    def export_relations_to_sbml(self, sbml_model, sbml):
        """ Add relationships to/from object to SBML parameter.

        Args:
            sbml_model (:obj:`libsbml.Model`): SBML model
            sbml (:obj:`libsbml.Parameter`): SBML parameter
        """
        annots = ['type', 'std', 'identifiers']
        LibSbmlInterface.set_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml)

    def import_from_sbml(self, sbml):
        """ Load from SBML parameter

        Args:
            sbml (:obj:`libsbml.Parameter`): SBML parameter
        """
        self.id = self.parse_sbml_id(call_libsbml(sbml.getIdAttribute))
        self.name = call_libsbml(sbml.getName)
        self.value = call_libsbml(sbml.getValue)
        self.units = LibSbmlInterface.get_unit(sbml.getUnits)
        LibSbmlInterface.get_commments(self, sbml)

    def import_relations_from_sbml(self, sbml, objs):
        """ Load relationships from SBML parameter

        Args:
            sbml (:obj:`libsbml.Parameter`): SBML parameter
            objs (:obj:`dict`): dictionary that maps WC-Lang types to dictionaries that
                map the ids of WC-Lang objects to WC-Lang objects
        """
        annots = ['type', 'std', 'identifiers']
        LibSbmlInterface.get_annotations(self, LibSbmlInterface.gen_nested_attr_paths(annots), sbml, objs)


class ObservationGenotype(obj_model.Model, SbmlModelMixin):
    """ Genotype of an observation

    Attributes:
        taxon (:obj:`str`): taxon in which the observation was observed
        genetic_variant (:obj:`str`): genetic variant in which the observation was observed

    Related attributes:

        * genotype_observations (:obj:`list` of :obj:`Observation`): observations
    """
    taxon = StringAttribute()
    variant = StringAttribute()

    class Meta(obj_model.Model.Meta):
        tabular_orientation = TabularOrientation.multiple_cells
        unique_together = (('taxon', 'variant', ), )
        attribute_order = ('taxon', 'variant')
        children = {
            'submodel': (),
            'core_model': (),
        }
        child_attrs = {
            'sbml': (),
            'wc_sim': (),
        }

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return '__'.join([self.taxon, self.variant])


class ObservationEnv(obj_model.Model, SbmlModelMixin):
    """ Environment of an observation

    Attributes:
        temp (:obj:`float`): temperature at which the observation was observed
        temp_units (:obj:`unit_registry.Unit`): temperature units
        ph (:obj:`float`): pH at which the observation was observed
        ph_units (:obj:`unit_registry.Unit`): pH units
        growth_media (:obj:`str`): growth media at which the observation was observed
        condition (:obj:`str`): experimental conditions (e.g. control)

    Related attributes:

        * env_observations (:obj:`list` of :obj:`Observation`): observations
    """
    temp = FloatAttribute(nan=True, verbose_name='Temperature')
    temp_units = UnitAttribute(unit_registry,
                               choices=(unit_registry.parse_units('celsius'),),
                               none=True,
                               verbose_name='Temperature units')
    ph = FloatAttribute(nan=True, verbose_name='pH')
    ph_units = UnitAttribute(unit_registry,
                             choices=(unit_registry.parse_units('dimensionless'),),
                             none=True,
                             verbose_name='pH units')
    growth_media = LongStringAttribute()
    condition = LongStringAttribute()

    class Meta(obj_model.Model.Meta):
        tabular_orientation = TabularOrientation.multiple_cells
        unique_together = (('temp', 'temp_units', 'ph', 'ph_units', 'growth_media', 'condition'), )
        attribute_order = ('temp', 'temp_units', 'ph', 'ph_units', 'growth_media', 'condition')
        children = {
            'submodel': (),
            'core_model': (),
        }
        child_attrs = {
            'sbml': (),
            'wc_sim': (),
        }

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return '__'.join([str(self.temp), str(self.temp_units),
                          str(self.ph), str(self.ph_units),
                          self.growth_media, self.condition])

    def validate(self):
        """ Determine if the environment is valid

        * temperature units are defined if the temperature is not None
        * pH units are defined if the pH is not None

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """
        invalid_obj = super(ObservationEnv, self).validate()
        if invalid_obj:
            errors = invalid_obj.attributes
        else:
            errors = []

        if self.temp is not None and not isnan(self.temp) \
                and not isinstance(self.temp_units, unit_registry.Unit):
            errors.append(InvalidAttribute(self.Meta.attributes['temp_units'],
                                           ['Temperature units must be defined']))
        if self.ph is not None and not isnan(self.ph) \
                and not isinstance(self.ph_units, unit_registry.Unit):
            errors.append(InvalidAttribute(self.Meta.attributes['ph_units'],
                                           ['pH units must be defined']))

        if errors:
            return InvalidObject(self, errors)
        return None


class Process(obj_model.Model, SbmlModelMixin):
    """ A process of an observation or conclusion 

    Attributes:
        name (:obj:`str`): procedure which produced the conclusion
        version (:obj:`str`): version of procedure which produced the conclusion

    Related attributes:

        * observation_analysis (:obj:`list` of :obj:`Observation`): observation
        * observation_measurement (:obj:`list` of :obj:`Observation`): observation
        * conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
    """
    name = LongStringAttribute()
    version = StringAttribute()

    class Meta(obj_model.Model.Meta):
        tabular_orientation = TabularOrientation.multiple_cells
        unique_together = (('name', 'version'), )
        attribute_order = ('name', 'version')
        children = {
            'submodel': (),
            'core_model': (),
        }
        child_attrs = {
            'sbml': (),
            'wc_sim': (),
        }

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        return '__'.join([self.name, self.version])


class Observation(obj_model.Model):
    """ Observation

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        value (:obj:`str`): value
        std (:obj:`str`): standard error of the value
        units (:obj:`unit_registry.Unit`): units
        type (:obj:`pronto.term.Term`): type        
        genotype (:obj:`ObservationGenotype`): genotype
        env (:obj:`ObservationEnv`): environment        
        experiment_type (:obj:`str`): type of experiment (e.g. RNA-seq)
        experiment_design (:obj:`str`): experimental design
        data_generation_process (:obj:`Process`): process used to measure data (e.g. deep sequencing)
        data_analysis_process (:obj:`Process`): process used to analyze data (e.g. Cufflinks)
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references

    Related attributes:

        * observation_sets (:obj:`list` of :obj:`ObservationSet`): observation sets
        * evidence (:obj:`list` of :obj:`Evidence`): evidence
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='observations')
    value = StringAttribute()
    std = StringAttribute(verbose_name='Standard error')
    units = UnitAttribute(unit_registry, none=True)
    type = OntologyAttribute(onto,
                             namespace='WC',
                             terms=onto['WC:observation'].rchildren(),
                             default=None, none=True)
    genotype = ManyToOneAttribute(ObservationGenotype, related_name='genotype_observations')
    env = ManyToOneAttribute(ObservationEnv, related_name='env_observations', verbose_name='Environment')
    experiment_type = LongStringAttribute()
    experiment_design = LongStringAttribute()
    data_generation_process = ManyToOneAttribute(Process, related_name='observation_data_generation')
    data_analysis_process = ManyToOneAttribute(Process, related_name='observation_data_analysis')
    identifiers = IdentifierManyToManyAttribute(related_name='observations')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='observations')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'value', 'std', 'units',
                           'type', 'genotype', 'env',
                           'experiment_type', 'experiment_design', 'data_generation_process', 'data_analysis_process',
                           'identifiers', 'comments', 'references')
        children = {
            'submodel': ('genotype', 'env', 'data_generation_process', 'data_analysis_process', 'identifiers', 'references'),
            'core_model': ('genotype', 'env', 'data_generation_process', 'data_analysis_process', 'identifiers', 'references'),
        }
        child_attrs = {
            'sbml': (),
            'wc_sim': (),
        }


class ObservationSet(obj_model.Model, SbmlModelMixin):
    """ Set of co-observed observations

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        observations (:obj:`list` of :obj:`Observation`): observations
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='observation_sets')
    observations = ManyToManyAttribute(Observation, related_name='observation_sets')
    identifiers = IdentifierManyToManyAttribute(related_name='observation_sets')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='observation_sets')

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'observations',
                           'identifiers', 'comments', 'references')
        children = {
            'submodel': ('observations', 'identifiers', 'references'),
            'core_model': ('observations', 'identifiers', 'references'),
        }
        child_attrs = {
            'sbml': (),
            'wc_sim': (),
        }


class Evidence(obj_model.Model):
    """ Observation that supports/disputes an conclusion

    Attributes:
        observation (:obj:`Observation`): observation which supports the conclusion
        type (:obj:`pronto.Term`): how the observation supports the conclusion (e.g. supporting, inconclusive, disputing)
        strength (:obj:`float): how much the observation supports the conclusion
        quality (:obj:`float`): the reliability of the observation

    Related attributes:

        * submodels (:obj:`list` of :obj:`Submodel`): submodels
        * dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        * compartments (:obj:`list` of :obj:`Compartment`): compartments
        * species_types (:obj:`list` of :obj:`SpeciesType`): species types
        * species (:obj:`list` of :obj:`Species`): species
        * distribution_init_concentrations (:obj:`list` of :obj:`DistributionInitConcentration`): initial concentrations
        * observables (:obj:`list` of :obj:`Observable`): observables
        * functions (:obj:`list` of :obj:`Function`): functions
        * stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
        * reactions (:obj:`list` of :obj:`Reaction`): reactions
        * rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        * dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
        * dfba_obj_reactions (:obj:`list` of :obj:`DfbaObjReaction`): dFBA objective reactions
        * parameters (:obj:`list` of :obj:`Parameter`): parameters
        * conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        * changes (:obj:`list` of :obj:`Change`): changes
    """
    observation = ManyToOneAttribute(Observation, related_name='evidence')
    type = OntologyAttribute(onto,
                             namespace='WC',
                             terms=onto['WC:evidence'].rchildren(),
                             default=None, none=True)
    strength = FloatAttribute()
    quality = FloatAttribute()

    class Meta(obj_model.Model.Meta):
        tabular_orientation = TabularOrientation.cell
        attribute_order = ('observation', 'type', 'strength', 'quality')
        children = {
            'submodel': ('observation',),
            'core_model': ('observation',),
        }
        child_attrs = {
            'sbml': (),
            'wc_sim': (),
        }

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: string representation
        """
        args = []

        if are_terms_equivalent(self.type, onto['WC:supporting_evidence']):
            args.append('+')
        elif are_terms_equivalent(self.type, onto['WC:disputing_evidence']):
            args.append('-')
        else:
            args.append('~')

        if self.strength is not None and not isnan(self.strength):
            args.append('s=' + str(self.strength))

        if self.quality is not None and not isnan(self.quality):
            args.append('q=' + str(self.quality))

        return '{}({})'.format(self.observation.serialize(), ', '.join(args))

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`DfbaObjectiveExpression`: cleaned value
            :obj:`InvalidAttribute`: cleaning error
        """
        match = re.match(r'^(.*?)\(([\+\-~])(,([qs])=(.*?))?(,([qs])=(.*?))?\)$', value.replace(' ', ''))
        strength = ''
        quality = ''
        if match:
            observation = match.group(1)
            type = match.group(2)
            if type == '+':
                type = 'supporting_evidence'
            elif type == '-':
                type = 'disputing_evidence'
            else:
                type = 'inconclusive_evidence'

            if match.group(4) and match.group(7) and match.group(4) == match.group(7):
                return (None, InvalidAttribute(cls.Meta.attributes['observation'], ['Invalid syntax']))

            if match.group(4) == 's':
                strength = match.group(5)
            elif match.group(4) == 'q':
                quality = match.group(5)

            if match.group(7) == 's':
                strength = match.group(8)
            elif match.group(7) == 'q':
                quality = match.group(8)
        else:
            return (None, InvalidAttribute(cls.Meta.attributes['observation'],
                                           ['Invalid syntax']))

        errors = []

        observation, error = cls.Meta.attributes['observation'].deserialize(observation, objects)
        if error:
            errors.extend(error.messages)

        type, error = cls.Meta.attributes['type'].deserialize(type)
        if error:
            errors.extend(error.messages)  # pragma: no cover # unreachable because of parsing above

        strength, error = cls.Meta.attributes['strength'].deserialize(strength)
        if error:
            errors.extend(error.messages)

        quality, error = cls.Meta.attributes['quality'].deserialize(quality)
        if error:
            errors.extend(error.messages)

        if errors:
            return (None, InvalidAttribute(cls.Meta.attributes['observation'], errors))

        obj = cls(observation=observation, type=type, strength=strength, quality=quality)
        if cls not in objects:
            objects[cls] = {}
        serialized_val = obj.serialize()
        existing_obj = objects[cls].get(serialized_val, None)
        if existing_obj:
            return (existing_obj, None)
        objects[cls][serialized_val] = obj
        return (obj, None)


class Conclusion(obj_model.Model):
    """ Conclusion of one or more observations

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        value (:obj:`str`): value
        std (:obj:`str`): standard error of the value
        units (:obj:`unit_registry.Unit`): units
        type (:obj:`pronto.term.Term`): type
        process (:obj:`Process`): procedure which produced the conclusion
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        evidence (:obj:`list` of :obj:`Evidence`): evidence that supports/refutes the conclusion
            (e.g. individual observations underlying an average)
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
        authors (:obj:`list` of :obj:`Author`): authors
        date (:obj:`datetime.datetime`): date and time when the conclusion was made

    Related attributes:

        * submodels (:obj:`list` of :obj:`Submodel`): submodels
        * compartments (:obj:`list` of :obj:`Compartment`): compartments
        * species_types (:obj:`list` of :obj:`SpeciesType`): species types
        * species (:obj:`list` of :obj:`Species`): species
        * distribution_init_concentrations (:obj:`list` of :obj:`DistributionInitConcentration`):
          distributions of initial concentrations of species at the beginning of each
          cell cycle
        * observables (:obj:`list` of :obj:`Observable`): observables
        * functions (:obj:`list` of :obj:`Function`): functions
        * reactions (:obj:`list` of :obj:`Reaction`): reactions
        * rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        * dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        * dfba_obj_reactions (:obj:`list` of :obj:`DfbaObjReaction`): dFBA objective reactions
        * dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
        * stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
        * parameters (:obj:`list` of :obj:`Parameter`): parameters
        * changes (:obj:`list` of :obj:`Change`): changes
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='conclusions')
    value = StringAttribute()
    std = StringAttribute(verbose_name='Standard error')
    units = UnitAttribute(unit_registry, none=True)
    type = OntologyAttribute(onto,
                             namespace='WC',
                             terms=onto['WC:conclusion'].rchildren(),
                             default=None, none=True)
    process = ManyToOneAttribute(Process, related_name='conclusions')
    identifiers = IdentifierManyToManyAttribute(related_name='conclusions')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='conclusions')
    comments = CommentAttribute()
    references = ManyToManyAttribute('Reference', related_name='conclusions')
    authors = ManyToManyAttribute('Author', related_name='conclusions')
    date = DateTimeAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'value', 'std', 'units',
                           'type', 'process',
                           'identifiers', 'evidence', 'comments', 'references', 'authors', 'date')
        children = {
            'submodel': ('process', 'identifiers', 'evidence', 'references', 'authors'),
            'core_model': ('process', 'identifiers', 'evidence', 'references', 'authors'),
        }
        child_attrs = {
            'sbml': (),
            'wc_sim': (),
        }


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
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        comments (:obj:`str`): comments

    Related attributes:

        * taxon (:obj:`Taxon`): taxon
        * env (:obj:`Environment`): environment
        * submodels (:obj:`list` of :obj:`Submodel`): submodels
        * compartments (:obj:`list` of :obj:`Compartment`): compartments
        * species_types (:obj:`list` of :obj:`SpeciesType`): species types
        * species (:obj:`list` of :obj:`Species`): species
        * distribution_init_concentrations (:obj:`list` of :obj:`DistributionInitConcentration`):
          distributions of initial concentrations of species at the beginning of
          each cell cycle
        * observables (:obj:`list` of :obj:`Observable`): observables
        * functions (:obj:`list` of :obj:`Function`): functions
        * reactions (:obj:`list` of :obj:`Reaction`): reactions
        * rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        * dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        * dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
        * stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
        * parameters (:obj:`list` of :obj:`Parameter`): parameters
        * changes (:obj:`list` of :obj:`Change`): changes
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute(Model, related_name='references')
    title = StringAttribute()
    author = StringAttribute()
    editor = StringAttribute()
    year = PositiveIntegerAttribute()
    type = OntologyAttribute(onto,
                             namespace='WC',
                             terms=onto['WC:reference'].rchildren(),
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
    identifiers = IdentifierManyToManyAttribute(related_name='references')
    comments = CommentAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'title', 'author', 'editor', 'year', 'type', 'publication', 'publisher',
                           'series', 'volume', 'number', 'issue', 'edition', 'chapter', 'pages',
                           'identifiers', 'comments')
        children = {
            'submodel': ('identifiers',),
            'core_model': ('identifiers',),
        }
        child_attrs = {
            'sbml': (),
            'wc_sim': (),
        }


class Author(obj_model.Model, SbmlModelMixin):
    """ An author of a model

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): full name
        model (:obj:`Model`): model
        last_name (:obj:`str`): last name(s)
        first_name (:obj:`str`): first name(s)
        middle_name (:obj:`str`): middle name(s)
        title (:obj:`str`): title
        organization (:obj:`str`): organization
        email (:obj:`str`): email address
        website (:obj:`str`): website
        address (:obj:`str`): physical address
        identifiers (:obj:`list` of :obj:`Identifier`): identifiers
        comments (:obj:`str`): comments

    Related attributes:

        * conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        * changes (:obj:`list` of :obj:`Change`): changes
    """
    id = SlugAttribute()
    name = StringAttribute(min_length=1)
    model = ManyToOneAttribute(Model, related_name='authors')
    last_name = StringAttribute(min_length=1)
    first_name = StringAttribute()
    middle_name = StringAttribute(min_length=1)
    title = LongStringAttribute()
    organization = LongStringAttribute()
    email = EmailAttribute()
    website = UrlAttribute()
    address = LongStringAttribute()
    identifiers = IdentifierManyToManyAttribute(related_name='authors')
    comments = LongStringAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'last_name', 'first_name', 'middle_name',
                           'title', 'organization',
                           'email', 'website', 'address',
                           'identifiers', 'comments')
        frozen_columns = 2
        children = {
            'submodel': ('identifiers',),
            'core_model': ('identifiers',),
        }
        child_attrs = {
            'sbml': ('id', 'name', 'model',
                     'last_name', 'first_name', 'middle_name',
                     'title', 'organization',
                     'email', 'website', 'address',
                     'identifiers', 'comments'),
            'wc_sim': (),
        }

    def get_identifier(self, namespace):
        """ Get the author's id in a namespace (e.g., `github.user`, `orcid`)

        Args:
            namespace (:obj:`str`): namespace of the identifier to retrieve

        Returns:
            :obj:`str`: user's id in :obj:`namespace`

        Raises:
            :obj:`ValueError`: if the author has multiple ids in :obj:`namespace`
        """
        identifiers = self.identifiers.get(namespace=namespace)
        if len(identifiers) == 1:
            return identifiers[0].id
        if len(identifiers) > 1:
            raise ValueError('Author {} has multiple {} ids'.format(self.id, namespace))
        return None


class Change(obj_model.Model, SbmlModelMixin):
    """ A change to a model

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): full name
        model (:obj:`Model`): model
        type (:obj:`pronto.term.Term`): type
        target (:obj:`str`): target
        target_submodel (:obj:`Submodel`): target submodel
        target_type (:obj:`pronto.term.Term`): target type
        reason (:obj:`str`): reason
        reason_type (:obj:`pronto.term.Term`): type of reason
        intention (:obj:`str`): intention
        intention_type (:obj:`pronto.term.Term`): type of intention
        identifiers (::obj:`list` of :obj:`Identifier`): identifiers
        conclusions (:obj:`list` of :obj:`Conclusion`): conclusions
        comments (:obj:`str`): comments
        references (:obj:`list` of :obj:`Reference`): references
        authors (:obj:`list` of :obj:`Author`): authors
        date (:obj:`datetime.datetime`): date
    """
    id = SlugAttribute()
    name = StringAttribute(min_length=1)
    model = ManyToOneAttribute(Model, related_name='changes')

    type = OntologyAttribute(onto, namespace='WC',
                             terms=onto['WC:change'].rchildren())
    target = LongStringAttribute()
    target_submodel = ManyToOneAttribute(Submodel, related_name='changes')
    target_type = OntologyAttribute(onto, namespace='WC',
                                    terms=onto['WC:target_provenance'].rchildren())
    reason = LongStringAttribute()
    reason_type = OntologyAttribute(onto, namespace='WC',
                                    terms=onto['WC:reason_provenance'].rchildren())
    intention = LongStringAttribute()
    intention_type = OntologyAttribute(onto, namespace='WC',
                                       terms=onto['WC:intention_provenance'].rchildren())
    identifiers = IdentifierManyToManyAttribute(related_name='changes')
    evidence = EvidenceManyToManyAttribute('Evidence', related_name='changes')
    conclusions = ManyToManyAttribute('Conclusion', related_name='changes')
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='changes')
    authors = ManyToManyAttribute('Author', related_name='changes')
    date = DateTimeAttribute()

    class Meta(obj_model.Model.Meta):
        attribute_order = ('id', 'name',
                           'type', 'target', 'target_submodel', 'target_type', 'reason', 'reason_type', 'intention', 'intention_type',
                           'identifiers', 'evidence', 'conclusions', 'comments', 'references',
                           'authors', 'date')
        children = {
            'submodel': ('identifiers', 'evidence', 'conclusions', 'references', 'authors'),
            'core_model': ('identifiers', 'evidence', 'conclusions', 'references', 'authors'),
        }
        child_attrs = {
            'sbml': (),
            'wc_sim': (),
        }


class Identifier(obj_model.Model, SbmlModelMixin):
    """ Reference to an entry in a namespace

    Attributes:
        namespace (:obj:`str`): namespace name
        id (:obj:`str`): id of entry in namespace

    Related attributes:

        * model (:obj:`Model`): model
        * taxon (:obj:`Taxon`): taxon
        * env (:obj:`Environment`): environment
        * submodels (:obj:`list` of :obj:`Submodel`): submodels
        * compartments (:obj:`list` of :obj:`Compartment`): compartments
        * species_types (:obj:`list` of :obj:`SpeciesType`): species types
        * species (:obj:`list` of :obj:`Species`): species
        * distribution_init_concentrations (:obj:`list` of :obj:`DistributionInitConcentration`): distributions
          of initial concentrations of species at the beginning of each cell cycle
        * observables (:obj:`list` of :obj:`Observable`): observables
        * functions (:obj:`list` of :obj:`Function`): functions
        * reactions (:obj:`list` of :obj:`Reaction`): reactions
        * rate_laws (:obj:`list` of :obj:`RateLaw`): rate laws
        * dfba_objs (:obj:`list` of :obj:`DfbaObjective`): dFBA objectives
        * dfba_obj_reactions (:obj:`list` of :obj:`DfbaObjReaction`): dFBA objective reactions
        * dfba_obj_species (:obj:`list` of :obj:`DfbaObjSpecies`): dFBA objective species
        * stop_conditions (:obj:`list` of :obj:`StopCondition`): stop conditions
        * parameters (:obj:`list` of :obj:`Parameter`): parameters
        * references (:obj:`list` of :obj:`Reference`): references
        * authors (:obj:`list` of :obj:`Author`): authors
        * changes (:obj:`list` of :obj:`Change`): changes
    """

    namespace = StringAttribute(min_length=1)
    id = StringAttribute(min_length=1)

    class Meta(obj_model.Model.Meta):
        unique_together = (('namespace', 'id', ), )
        tabular_orientation = TabularOrientation.cell
        attribute_order = ('namespace', 'id')
        frozen_columns = 2
        ordering = ('namespace', 'id', )
        child_attrs = {
            'sbml': ('namespace', 'id'),
            'wc_sim': (),
        }

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return '{}: {}'.format(self.namespace, self.id)

    @classmethod
    def deserialize(cls, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`Identifier`: cleaned value
            :obj:`InvalidAttribute`: cleaning error
        """
        if ': ' not in value:
            return (None, InvalidAttribute(cls.Meta.attributes['id'], ['Invalid format']))

        namespace, _, id = value.strip().partition(': ')
        identifier = cls(namespace=namespace.strip(), id=id.strip())

        if Identifier not in objects:
            objects[Identifier] = {}

        serialized_val = identifier.serialize()
        if serialized_val in objects[Identifier]:
            identifier = objects[Identifier][serialized_val]
        else:
            objects[Identifier][serialized_val] = identifier

        return (identifier, None)


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


class WcLangWarning(UserWarning):
    """ WC-Lang warning """
    pass  # pragma: no cover

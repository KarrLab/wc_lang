""" Data model to represent bio-chemical models.

This module defines classes that represent the schema of a bio-chemical model:

* :obj:`Taxon`
* :obj:`Submodel`
* :obj:`Compartment`
* :obj:`SpeciesType`
* :obj:`Species`
* :obj:`Concentration`
* :obj:`Reaction`
* :obj:`ReactionParticipant`
* :obj:`RateLaw`
* :obj:`RateLawEquation`
* :obj:`Parameter`
* :obj:`Reference`
* :obj:`CrossReference`

These are all instances of `BaseModel`, an alias for `wc_utils.schema.core.Model`.
A bio-chemical model may contain a list of instances of each of these classes, interlinked
by object references. For example, A :obj:`Reaction` will reference its constituent
:obj:`ReactionParticipant` instances, and the :obj:`RateLaw` that describes the reaction's rate.

This module also defines numerous classes that serve as attributes of these classes.

Many classes contain the methods `serialize()` and `deserialize()`, which invert each other.
`serialize()` converts a python object instance into a string representation, whereas
`deserialize()` parses an object's string representation -- as would be stored in a file or spreadsheet
representation of a bio-chemical model -- into a python object instance.
`deserialize()` returns an error when the string representatin cannot be parsed into the
python object.


:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 2016, Karr Lab
:License: MIT
"""

from enum import Enum, EnumMeta
from itertools import chain
from math import ceil, floor, exp, log, log10, isnan
from natsort import natsorted, ns
from six import with_metaclass
from wc_utils.schema.core import (Model as BaseModel,
                                  BooleanAttribute, EnumAttribute, FloatAttribute, IntegerAttribute, PositiveIntegerAttribute,
                                  RegexAttribute, SlugAttribute, StringAttribute, LongStringAttribute, UrlAttribute,
                                  OneToOneAttribute, ManyToOneAttribute, ManyToManyAttribute,
                                  InvalidModel, InvalidObject, InvalidAttribute,
                                  TabularOrientation)
from wc_utils.util.enumerate import CaseInsensitiveEnum, CaseInsensitiveEnumMeta
import re
import sys


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


class TaxonRank(with_metaclass(TaxonRankMeta, Enum)):
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


class SubmodelAlgorithm(CaseInsensitiveEnum):
    """ Submodel algorithms """
    dfba = 1
    ode = 2
    ssa = 3


class SpeciesTypeType(CaseInsensitiveEnum):
    """ Types of species types """
    metabolite = 1
    protein = 2
    rna = 3
    pseudo_species = 4


class RateLawDirection(CaseInsensitiveEnum):
    """ Rate law directions """
    backward = -1
    forward = 1


class ReferenceType(CaseInsensitiveEnum):
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


class OneToOneSpeciesAttribute(OneToOneAttribute):
    """ Species attribute """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(OneToOneSpeciesAttribute, self).__init__('Species',
                                                       related_name=related_name, none=False, related_none=True,
                                                       verbose_name=verbose_name, verbose_related_name=verbose_related_name, help=help)

    def serialize(self, value):
        """ Serialize related object

        Args:
            value (:obj:`Model`): Python representation

        Returns:
            :obj:`str`: simple Python representation
        """
        return value.serialize()

    def deserialize(self, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of `set` of `ReactionParticipant`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        return Species.deserialize(self, value, objects)


class ReactionParticipantsAttribute(ManyToManyAttribute):
    """ Reaction participants """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(ReactionParticipantsAttribute, self).__init__('ReactionParticipant', related_name=related_name,
                                                            verbose_name=verbose_name, verbose_related_name=verbose_related_name, help=help)

    def serialize(self, participants):
        """ Serialize related object

        Args:
            participants (:obj:`set` of `ReactionParticipant`): Python representation of reaction participants

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
                lhs.append(part.serialize(global_comp))
            elif part.coefficient > 0:
                rhs.append(part.serialize(global_comp))

        if global_comp:
            return '[{}]: {} ==> {}'.format(global_comp.get_primary_attribute(), ' + '.join(lhs), ' + '.join(rhs))
        else:
            return '{} ==> {}'.format(' + '.join(lhs), ' + '.join(rhs))

    def deserialize(self, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of `set` of `ReactionParticipant`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        errors = []

        id = '[a-z][a-z0-9_]*'
        stoch = '\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\)'
        gbl_part = '({} )*({})'.format(stoch, id)
        lcl_part = '({} )*({}\[{}\])'.format(stoch, id, id)
        gbl_side = '{}( \+ {})*'.format(gbl_part, gbl_part)
        lcl_side = '{}( \+ {})*'.format(lcl_part, lcl_part)
        gbl_pattern = '^\[({})\]: ({}) ==> ({})$'.format(id, gbl_side, gbl_side)
        lcl_pattern = '^({}) ==> ({})$'.format(lcl_side, lcl_side)

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

        parts = set()

        for part in re.findall('(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*)(\[([a-z][a-z0-9_]*)\])*', lhs, flags=re.I):
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

            coefficient = float(part[1] or 1.)

            if part_errors:
                errors += part_errors
            else:
                spec_primary_attribute = '{}[{}]'.format(
                    species_type.get_primary_attribute(), compartment.get_primary_attribute())
                species, error = Species.deserialize(self, spec_primary_attribute, objects)
                if error:
                    raise ValueError('Invalid species "{}"'.format(spec_primary_attribute))

                if coefficient != 0:
                    rxn_part = ReactionParticipant(species=species, coefficient=-1 * coefficient)
                    if ReactionParticipant not in objects:
                        objects[ReactionParticipant] = {}
                    objects[ReactionParticipant][rxn_part.serialize()] = rxn_part
                    parts.add(rxn_part)

        for part in re.findall('(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*)(\[([a-z][a-z0-9_]*)\])*', rhs, flags=re.I):
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

            coefficient = float(part[1] or 1.)

            if part_errors:
                errors += part_errors
            else:
                spec_primary_attribute = '{}[{}]'.format(
                    species_type.get_primary_attribute(), compartment.get_primary_attribute())
                species, error = Species.deserialize(self, spec_primary_attribute, objects)
                if error:
                    raise ValueError('Invalid species "{}"'.format(spec_primary_attribute))

                if coefficient != 0:
                    rxn_part = ReactionParticipant(species=species, coefficient=coefficient)
                    if ReactionParticipant not in objects:
                        objects[ReactionParticipant] = {}
                    objects[ReactionParticipant][rxn_part.serialize()] = rxn_part
                    parts.add(rxn_part)

        if errors:
            return (None, InvalidAttribute(self, errors))
        return (parts, None)


class RateLawEquationAttribute(OneToOneAttribute):
    """ Rate law equation """

    def __init__(self, related_name='', verbose_name='', verbose_related_name='', help=''):
        """
        Args:
            related_name (:obj:`str`, optional): name of related attribute on `related_class`
            verbose_name (:obj:`str`, optional): verbose name
            verbose_related_name (:obj:`str`, optional): verbose related name
            help (:obj:`str`, optional): help message
        """
        super(RateLawEquationAttribute, self).__init__('RateLawEquation',
                                                       related_name=related_name, none=False, related_none=False,
                                                       verbose_name=verbose_name, verbose_related_name=verbose_related_name, help=help)

    def serialize(self, value):
        """ Serialize related object

        Args:
            participants (:obj:`Model`): Python representation

        Returns:
            :obj:`str`: simple Python representation
        """
        return value.serialize()

    def deserialize(self, value, objects):
        """ Deserialize value

        Args:
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of `object`, `InvalidAttribute` or `None`: tuple of cleaned value and cleaning error
        """
        return RateLawEquation.deserialize(self, value, objects)


class Model(BaseModel):
    """ Model

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        version (:obj:`str`): version number
        wc_lang_version (:obj:`str`): wc_lang version number
        comments (:obj:`str`): comments

        cross_references (:obj:`set` of `CrossReference`): cross references
        taxon (:obj:`Taxon`): taxon
        submodels (:obj:`set` of `Submodel`): submodels
        compartments (:obj:`set` of `Compartment`): compartments
        species_types (:obj:`set` of `SpeciesType`): species types
        parameters (:obj:`set` of `Parameter`): parameters
        references (:obj:`set` of `Reference`): references
    """
    id = SlugAttribute()
    name = StringAttribute()
    version = RegexAttribute(min_length=1, pattern='^[0-9]+\.[0-9+]\.[0-9]+', flags=re.I)
    wc_lang_version = RegexAttribute(min_length=1, pattern='^[0-9]+\.[0-9+]\.[0-9]+', flags=re.I,
                                     verbose_name='wc_lang version')
    comments = LongStringAttribute()

    class Meta(BaseModel.Meta):
        attribute_order = ('id', 'name', 'version', 'wc_lang_version', 'comments')
        tabular_orientation = TabularOrientation.column

    def get_compartments(self):
        """ Get all compartments

        Returns:
            :obj:`set` of `Compartment`: compartments
        """
        return set(self.compartments)

    def get_species_types(self):
        """ Get all species types

        Returns:
            :obj:`set` of `SpeciesType`: species types
        """
        return set(self.species_types)

    def get_submodels(self):
        """ Get all submodels

        Returns:
            :obj:`set` of `Submodel`: submodels
        """
        return set(self.submodels)

    def get_species(self):
        """ Get all species from submodels

        Returns:
            :obj:`set` of `Species`: species
        """
        species = set()

        for submodel in self.submodels:
            species.update(submodel.get_species())

        for concentation in self.get_concentrations():
            species.add(concentation.species)

        return species

    def get_concentrations(self):
        """ Get all concentrations from species types

        Returns:
            :obj:`set` of `Concentration`: concentations
        """
        concentrations = set()
        for species_type in self.species_types:
            for species in species_type.species:
                if species.concentration:
                    concentrations.add(species.concentration)
        return concentrations

    def get_reactions(self):
        """ Get all reactions from submodels

        Returns:
            :obj:`set` of `Reaction`: reactions
        """
        reactions = set()
        for submodel in self.submodels:
            reactions.update(submodel.reactions)
        return reactions

    def get_rate_laws(self):
        """ Get all rate laws from reactions

        Returns:
            :obj:`set` of `RateLaw`: rate laws
        """
        rate_laws = set()
        for reaction in self.get_reactions():
            rate_laws.update(reaction.rate_laws)
        return rate_laws

    def get_parameters(self):
        """ Get all parameters from model and submodels

        Returns:
            :obj:`set` of `Parameter`: parameters
        """
        parameters = set()
        parameters.update(self.parameters)
        for submodel in self.submodels:
            parameters.update(submodel.parameters)
        return parameters

    def get_references(self):
        """ Get all references from model and children

        Returns:
            :obj:`set` of `Reference`: references
        """
        refs = set()

        refs.update(self.references)

        if self.taxon:
            refs.update(self.taxon.references)

        for compartment in self.compartments:
            refs.update(compartment.references)

        for species_type in self.species_types:
            refs.update(species_type.references)

        for concentration in self.get_concentrations():
            refs.update(concentration.references)

        for submodel in self.submodels:
            refs.update(submodel.references)

        for reaction in self.get_reactions():
            refs.update(reaction.references)

        for rate_law in self.get_rate_laws():
            refs.update(rate_law.references)

        for parameter in self.get_parameters():
            refs.update(parameter.references)

        return refs

    def get_component(self, type, id):
        """ Find model component of `type` with `id`

        Args:
            type (:obj:`str`) type of component to find
            id (:obj:`str`): id of component to find

        Returns:
            :obj:`BaseModel`: component with `id`, or `None` if there is no component with `id`=`id`
        """
        types = ['compartment', 'species_type', 'submodel', 'reaction', 'parameter', 'reference']
        if type not in types:
            raise ValueError('Type must be one of {}'.format(', '.join(types)))

        components = getattr(self, 'get_{}s'.format(type))()
        return next((c for c in components if c.id == id), None)


class Taxon(BaseModel):
    """ Biological taxon (e.g. family, genus, species, strain, etc.)

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        rank (:obj:`TaxonRank`): rank
        comments (:obj:`str`): comments
        references (:obj:`set` of `Reference`): references

        cross_references (:obj:`set` of `CrossReference`): cross references
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = OneToOneAttribute('Model', related_name='taxon')
    rank = EnumAttribute(TaxonRank, default=TaxonRank.species)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='taxa')

    class Meta(BaseModel.Meta):
        attribute_order = ('id', 'name',
                           'model',
                           'rank',
                           'comments', 'references')
        tabular_orientation = TabularOrientation.column


class Submodel(BaseModel):
    """ Submodel

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        algorithm (:obj:`SubmodelAlgorithm`): algorithm
        comments (:obj:`str`): comments
        references (:obj:`set` of `Reference`): references

        cross_references (:obj:`set` of `CrossReference`): cross references
        reactions (:obj:`set` of `Reaction`): reactions
        parameters (:obj:`set` of `Parameter`): parameters
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute('Model', related_name='submodels')
    algorithm = EnumAttribute(SubmodelAlgorithm, default=SubmodelAlgorithm.ssa)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='submodels')

    class Meta(BaseModel.Meta):
        attribute_order = ('id', 'name',
                           'model',
                           'algorithm',
                           'comments', 'references')

    def get_species(self):
        """ Get species in reactions

        Returns:
            :obj:`set` of `Species`: species in reactions
        """
        species = set()

        for rxn in self.reactions:
            species.update(rxn.get_species())

        return species


class Compartment(BaseModel):
    """ Compartment

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        initial_volume (:obj:`float`): initial volume(L)
        comments (:obj:`str`): comments
        references (:obj:`set` of `Reference`): references

        cross_references (:obj:`set` of `CrossReference`): cross references
        concentrations (:obj:`set` of `Concentration`): concentrations
        reaction_participants (:obj:`set` of `ReactionParticipant`): reaction participants
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute('Model', related_name='compartments')
    initial_volume = FloatAttribute(min=0)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='compartments')

    class Meta(BaseModel.Meta):
        attribute_order = ('id', 'name',
                           'model',
                           'initial_volume',
                           'comments', 'references')


class SpeciesType(BaseModel):
    """ Species type

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        model (:obj:`Model`): model
        structure (:obj:`str`): structure(InChI for metabolites; sequence for DNA, RNA, proteins)
        empirical_formula (:obj:`str`): empirical formula
        molecular_weight (:obj:`float`): molecular weight
        charge (:obj:`int`): charge
        type (:obj:`SpeciesTypeType`): type
        comments (:obj:`str`): comments
        references (:obj:`set` of `Reference`): references

        cross_references (:obj:`set` of `CrossReference`): cross references
        concentrations (:obj:`set` of `Concentration`): concentrations
        reaction_participants (:obj:`set` of `ReactionParticipant`): reaction participants
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute('Model', related_name='species_types')
    structure = LongStringAttribute()
    empirical_formula = RegexAttribute(pattern='^([A-Z][a-z]?\d*)*$')
    molecular_weight = FloatAttribute(min=0)
    charge = IntegerAttribute()
    type = EnumAttribute(SpeciesTypeType, default=SpeciesTypeType.metabolite)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='species_types')

    class Meta(BaseModel.Meta):
        verbose_name = 'Species type'
        attribute_order = ('id', 'name',
                           'model',
                           'structure', 'empirical_formula', 'molecular_weight', 'charge', 'type',
                           'comments', 'references')

    # todo: move to compiled model
    def is_carbon_containing(self):
        """ Returns `True` is species contains at least one carbon atom.

        Returns:
            :obj:`bool`: `True` is species contains at least one carbon atom.
        """
        return re.match('C[1-9]', self.empirical_formula) is not None


class Species(BaseModel):
    """ Species (tuple of species type, compartment)

    Attributes:
        species_type (:obj:`SpeciesType`): species type
        compartment (:obj:`Compartment`): compartment

        concentration (:obj:`Concentration`): concentration
        reaction_participants (:obj:`set` of `Reaction`): participations in reactions
        rate_law_equations (:obj:`RateLawEquation`): rate law equations
    """
    species_type = ManyToOneAttribute('SpeciesType', related_name='species', none=False)
    compartment = ManyToOneAttribute('Compartment', related_name='species', none=False)

    class Meta(BaseModel.Meta):
        attribute_order = ('species_type', 'compartment')
        frozen_columns = 1
        tabular_orientation = TabularOrientation.inline
        ordering = ('species_type', 'compartment')

    def serialize(self):
        """ Get value of primary attribute

        Returns:
            :obj:`str`: value of primary attribute
        """
        return '{}[{}]'.format(
            self.species_type.get_primary_attribute(),
            self.compartment.get_primary_attribute())

    @classmethod
    def deserialize(cls, attribute, value, objects):
        """ Deserialize value

        Args:
            attribute (:obj:`Attribute`): attribute
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of `object`, `InvalidAttribute` or `None`: tuple of cleaned value and cleaning error
        """
        if cls in objects and value in objects[cls]:
            return (objects[cls][value], None)

        match = re.match('^([a-z][a-z0-9_]*)\[([a-z][a-z0-9_]*)\]$', value, flags=re.I)
        if match:
            errors = []

            if match.group(1) in objects[SpeciesType]:
                species_type = objects[SpeciesType][match.group(1)]
            else:
                errors.append('Species type "{}" is not defined'.format(match.group(1)))

            if match.group(2) in objects[Compartment]:
                compartment = objects[Compartment][match.group(2)]
            else:
                errors.append('Compartment "{}" is not defined'.format(match.group(2)))

            if errors:
                return (None, InvalidAttribute(attribute, errors))
            else:
                obj = cls(species_type=species_type, compartment=compartment)
                if cls not in objects:
                    objects[cls] = {}
                objects[cls][obj.serialize()] = obj
                return (obj, None)

        return (None, InvalidAttribute(attribute, ['Invalid species']))


class Concentration(BaseModel):
    """ Species concentration

    Attributes:
        species (:obj:`Species`): species
        value (:obj:`float`): value (M)
        comments (:obj:`str`): comments
        references (:obj:`set` of `Reference`): references
    """
    species = OneToOneSpeciesAttribute(related_name='concentration')
    value = FloatAttribute(min=0)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='concentrations')

    class Meta(BaseModel.Meta):
        attribute_order = ('species',
                           'value',
                           'comments', 'references')
        frozen_columns = 1
        ordering = ('species',)

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return self.species.serialize()


class Reaction(BaseModel):
    """ Reaction

    Attributes:
        id (:obj:`str`): unique identifier
        name (:obj:`str`): name
        submodel (:obj:`Submodel`): submodel that reaction belongs to
        participants (:obj:`set` of `ReactionParticipant`): participants
        reversible (:obj:`bool`): indicates if reaction is thermodynamically reversible
        comments (:obj:`str`): comments
        references (:obj:`set` of `Reference`): references

        cross_references (:obj:`set` of `CrossReference`): cross references
        rate_laws (:obj:`set` of `RateLaw`): rate laws
    """
    id = SlugAttribute()
    name = StringAttribute()
    submodel = ManyToOneAttribute('Submodel', related_name='reactions')
    participants = ReactionParticipantsAttribute(related_name='reactions')
    reversible = BooleanAttribute()
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='reactions')

    class Meta(BaseModel.Meta):
        attribute_order = ('id', 'name',
                           'submodel',
                           'participants', 'reversible',
                           'comments', 'references')

    def get_species(self):
        """ Get species

        Returns:
            :obj:`set`: set of `Species`
        """
        species = set()

        for part in self.participants:
            species.add(part.species)

        for rate_law in self.rate_laws:
            species.update(rate_law.equation.modifiers)

        return species


class ReactionParticipant(BaseModel):
    """ Species in a reaction

    Attributes:
        species (:obj:`Species`): species
        coefficient (:obj:`float`): coefficient

        reaction (:obj:`Reaction`): reaction
    """
    species = ManyToOneAttribute('Species', related_name='reaction_participants')
    coefficient = FloatAttribute(nan=False)

    class Meta(BaseModel.Meta):
        attribute_order = ('species', 'coefficient')
        frozen_columns = 1
        tabular_orientation = TabularOrientation.inline
        ordering = ('species',)

    def serialize(self, compartment=None):
        """ Serialize related object

        Args:
            compartment (:obj:`Compartment`, optional): global compartment

        Returns:
            :obj:`str`: simple Python representation
        """
        return self._serialize(self.species, self.coefficient, compartment=compartment)

    @staticmethod
    def _serialize(species, coefficient, compartment=None):
        """ Serialize values

        Args:
            species (:obj:`Species`): species
            coefficient (:obj:`float`): coefficient
            compartment (:obj:`Compartment`, optional): global compartment

        Returns:
            :obj:`str`: simple Python representation
        """
        coefficient = float(coefficient)

        if abs(coefficient) == 1:
            coefficient_str = ''
        elif coefficient % 1 == 0 and abs(coefficient) < 1000:
            coefficient_str = '({:.0f}) '.format(abs(coefficient))
        else:
            coefficient_str = '({:e}) '.format(abs(coefficient))

        if compartment:
            return '{}{}'.format(coefficient_str, species.species_type.get_primary_attribute())
        else:
            return '{}{}'.format(coefficient_str, species.serialize())

    @classmethod
    def deserialize(cls, attribute, value, objects, compartment=None):
        """ Deserialize value

        Args:
            attribute (:obj:`Attribute`): attribute
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model
            compartment (:obj:`Compartment`, optional): compartment

        Returns:
            :obj:`tuple` of `set` of `ReactionParticipant`, `InvalidAttribute` or `None`: tuple of cleaned value
                and cleaning error
        """
        errors = []

        if compartment:
            pattern = '^(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*)$'
        else:
            pattern = '^(\(((\d*\.?\d+|\d+\.)(e[\-\+]?\d+)?)\) )*([a-z][a-z0-9_]*\[[a-z][a-z0-9_]*\])$'

        match = re.match(pattern, value, flags=re.I)
        if match:
            errors = []

            coefficient = float(match.group(2) or 1.)

            if compartment:
                species_id = '{}[{}]'.format(match.group(5), compartment.get_primary_attribute())
            else:
                species_id = match.group(5)

            species, error = Species.deserialize(attribute, species_id, objects)
            if error:
                return error

            serial_val = cls._serialize(species, coefficient)
            if cls in objects and serial_val in objects[cls]:
                return (objects[cls][serial_val], None)

            obj = cls(species=species, coefficient=coefficient)
            if cls not in objects:
                objects[cls] = {}
            objects[cls][obj.serialize()] = obj
            return (obj, None)

        else:
            attr = cls.Meta.attributes['species']
            return (None, InvalidAttribute(attr, ['Invalid reaction participant']))


class RateLaw(BaseModel):
    """ Rate law

    Attributes:
        reaction (:obj:`Reaction`): reaction
        direction (:obj:`RateLawDirection`): direction
        equation (:obj:`RateLawEquation`): equation
        k_cat (:obj:`float`): v_max (reactions enz^-1 s^-1)
        k_m (:obj:`float`): k_m (M)
        comments (:obj:`str`): comments
        references (:obj:`set` of `Reference`): references
    """

    reaction = ManyToOneAttribute('Reaction', related_name='rate_laws')
    direction = EnumAttribute(RateLawDirection, default=RateLawDirection.forward)
    equation = RateLawEquationAttribute(related_name='rate_law')
    k_cat = FloatAttribute(min=0, nan=True)
    k_m = FloatAttribute(min=0, nan=True)
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='rate_laws')

    class Meta(BaseModel.Meta):
        attribute_order = ('reaction', 'direction',
                           'equation', 'k_cat', 'k_m',
                           'comments', 'references')
        unique_together = (('reaction', 'direction'), )
        ordering = ('reaction', 'direction',)

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return '{}.{}'.format(self.reaction.serialize(), self.direction.name)


class RateLawEquation(BaseModel):
    """ Rate law equation

    Attributes:
        expression (:obj:`str`): expression
        modifiers (:obj:`set` of `Species`): modifiers

        rate_law (:obj:`RateLaw`): rate law
    """
    expression = LongStringAttribute()
    modifiers = ManyToManyAttribute('Species', related_name='rate_law_equations')

    class Meta(BaseModel.Meta):
        """
        Attributes:
            valid_functions (:obj:`tuple` of `str`): tuple of names of functions that can be used in `law`
        """
        attribute_order = ('expression', 'modifiers')
        tabular_orientation = TabularOrientation.inline
        valid_functions = (ceil, floor, exp, pow, log, log10, min, max)
        ordering = ('rate_law',)

    def serialize(self):
        """ Generate string representation

        Returns:
            :obj:`str`: value of primary attribute
        """
        return self.expression

    @classmethod
    def deserialize(cls, attribute, value, objects):
        """ Deserialize value

        Args:
            attribute (:obj:`Attribute`): attribute
            value (:obj:`str`): String representation
            objects (:obj:`dict`): dictionary of objects, grouped by model

        Returns:
            :obj:`tuple` of `object`, `InvalidAttribute` or `None`: tuple of cleaned value and cleaning error
        """
        modifiers = set()
        errors = []
        pattern = '(^|[^a-z0-9_])({}\[{}\])([^a-z0-9_]|$)'.format(SpeciesType.id.pattern[1:-1],
                                                                  Compartment.id.pattern[1:-1])

        for match in re.findall(pattern, value, flags=re.I):
            species, error = Species.deserialize(attribute, match[1], objects)
            if error:
                errors += error.messages
            else:
                modifiers.add(species)

        if errors:
            attr = cls.Meta.attributes['expression']
            return (None, InvalidAttribute(attribute, errors))

        # return value
        obj = cls(expression=value, modifiers=modifiers)
        if cls not in objects:
            objects[cls] = {}
        objects[cls][obj.serialize()] = obj
        return (obj, None)

    def validate(self):
        """ Determine if the object is valid

        Returns:
            :obj:`InvalidObject` or None: `None` if the object is valid,
                otherwise return a list of errors as an instance of `InvalidObject`
        """

        pattern = '(^|[^a-z0-9_])({}\[{}\])([^a-z0-9_]|$)'.format(SpeciesType.id.pattern[1:-1],
                                                                  Compartment.id.pattern[1:-1])

        """ check that all named entities are defined """
        modifier_ids = set((x.serialize() for x in self.modifiers))
        species_ids = set([x[1] for x in re.findall(pattern, self.expression, flags=re.I)])

        if modifier_ids != species_ids:
            errors = []

            for id in modifier_ids.difference(species_ids):
                errors.append('Extraneous modifier "{}"'.format(id))

            for id in species_ids.difference(modifier_ids):
                errors.append('Undefined modifier "{}"'.format(id))

            attr = self.__class__.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, errors)
            return InvalidObject(self, [attr_err])

        """ Check that rate law evaluates """
        # setup name space and format expression for python evaluation
        local_ns = {f.__name__: f for f in self.__class__.Meta.valid_functions}

        if self.rate_law:
            if not isnan(self.rate_law.k_cat):
                local_ns['k_cat'] = 1.

            if not isnan(self.rate_law.k_m):
                local_ns['k_m'] = 1.

        local_ns['s'] = []
        py_expression = self.expression
        for i_species, species_id in enumerate(set([x[1] for x in re.findall(pattern, self.expression, flags=re.I)])):
            local_ns['s'].append(1.)
            py_expression = py_expression.replace(species_id, 's[{}]'.format(i_species))

        # try evaluating law
        try:
            eval(py_expression, {'__builtins__': None}, local_ns)
        except Exception as error:
            msg = 'Invalid function: {}\n  {}'.format(self.expression, str(error).replace('\n', '\n  '))
            attr = self.__class__.Meta.attributes['expression']
            attr_err = InvalidAttribute(attr, [msg])
            return InvalidObject(self, [attr_err])

        """ return `None` to indicate valid object """
        return None


class Parameter(BaseModel):
    """ Parameter

    Attributes:
        id (:obj:`str`): unique identifier per model/submodel
        name (:obj:`str`): name
        model (:obj:`Model`): model
        submodel (:obj:`Submodel`): submodel
        value (:obj:`float`): value
        units (:obj:`str`): units of value
        comments (:obj:`str`): comments
        references (:obj:`set` of `Reference`): references
    """
    id = SlugAttribute(unique=False)
    name = StringAttribute()
    model = ManyToOneAttribute('Model', related_name='parameters')
    submodels = ManyToManyAttribute('Submodel', related_name='parameters')
    value = FloatAttribute(min=0)
    units = StringAttribute()
    comments = LongStringAttribute()
    references = ManyToManyAttribute('Reference', related_name='parameters')

    class Meta(BaseModel.Meta):
        attribute_order = ('id', 'name',
                           'model', 'submodels',
                           'value', 'units',
                           'comments', 'references')

    @classmethod
    def validate_unique(cls, objects):
        """ Validate attribute uniqueness. Make sure that parameter ids are unique
        within each model/submodel.

        Args:
            objects (:obj:`set` of `Parameter`): set of objects

        Returns:
            :obj:`InvalidModel`: list of invalid attributes and their errors
        """

        # validate attributes uniqueness
        error = super(Parameter, cls).validate_unique(objects)
        if error:
            errors = error.attributes
        else:
            errors = []

        # check for id uniqueness within each model/submodel
        model_ids = set()
        model_rep_ids = set()
        submodel_ids = {}
        submodel_rep_ids = {}
        for obj in objects:
            if obj.model:
                if obj.id in model_ids:
                    model_rep_ids.add(obj.id)
                else:
                    model_ids.add(obj.id)

            for submodel in obj.submodels:
                if submodel not in submodel_ids:
                    submodel_ids[submodel] = set()

                if obj.id in submodel_ids[submodel]:
                    if submodel not in submodel_rep_ids:
                        submodel_rep_ids[submodel] = set()
                    submodel_rep_ids[submodel].add(obj.id)
                else:
                    submodel_ids[submodel].add(obj.id)

        if model_rep_ids or submodel_rep_ids:
            msg = 'IDs must be unique within each model/submodel. The following IDs were repeated:'

            for id in model_rep_ids:
                msg += '\n- {}'.format(id)

            for submodel, rep_ids in submodel_rep_ids.items():
                msg += '\n- {}:'.format(submodel.id)
                for id in rep_ids:
                    msg += '\n  - {}'.format(id)

            errors.append(InvalidAttribute(cls.Meta.attributes['id'], [msg]))

        # return
        if errors:
            return InvalidModel(cls, errors)
        return None


class Reference(BaseModel):
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

        cross_references (:obj:`set` of `CrossReference`): cross references
        taxa (:obj:`set` of `Taxon`): taxa
        submodels (:obj:`set` of `Submodel`): submodels
        compartments (:obj:`set` of `Compartment`): compartments
        species_types (:obj:`set` of `SpeciesType`): species types
        concentrations (:obj:`set` of `Concentration`): concentrations
        reactions (:obj:`set` of `Reaction`): reactions
        rate_laws (:obj:`set` of `RateLaw`): rate laws
        parameters (:obj:`set` of `Parameter`): parameters
    """
    id = SlugAttribute()
    name = StringAttribute()
    model = ManyToOneAttribute('Model', related_name='references')
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

    class Meta(BaseModel.Meta):
        attribute_order = ('id', 'name', 'model',
                           'title', 'author', 'editor', 'year', 'type', 'publication', 'publisher',
                           'series', 'volume', 'number', 'issue', 'edition', 'chapter', 'pages',
                           'comments')


class CrossReference(BaseModel):
    """ Cross reference to a database

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
    url = UrlAttribute()
    model = ManyToOneAttribute('Model', related_name='cross_references')
    taxon = ManyToOneAttribute('Taxon', related_name='cross_references')
    submodel = ManyToOneAttribute('Submodel', related_name='cross_references')
    compartment = ManyToOneAttribute('Compartment', related_name='cross_references')
    species_type = ManyToOneAttribute('SpeciesType', related_name='cross_references')
    reaction = ManyToOneAttribute('Reaction', related_name='cross_references')
    reference = ManyToOneAttribute('Reference', related_name='cross_references')

    class Meta(BaseModel.Meta):
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

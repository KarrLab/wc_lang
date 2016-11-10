""" Data model to represent models. 

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2017-07-25
:Copyright: 2016, Karr Lab
:License: MIT

.. todo:: implement2 __eq__, __ne__ operators
"""

from itertools import chain
from operator import attrgetter
import math


class Model(object):
    """ Represents a model containing:

    * compartments
    * species
    * submodels
    * reactions
    * parameters
    * references

    Attributes:
        compartments (:obj:`list`): list of compartments
        species (:obj:`list`): list of species
        submodels (:obj:`list`): list of submodels
        reactions (:obj:`list`): list of reactions
        parameters (:obj:`list`): list of parameters
        references (:obj:`list`): list of references
    """

    def __init__(self, compartments=None, species=None, submodels=None, reactions=None, parameters=None, references=None):
        """ Construct a model

        Args:
            compartments (:obj:`list`, optional): list of compartments
            species (:obj:`list`), optional: list of species
            submodels (:obj:`list`, optional): list of submodels
            reactions (:obj:`list`), optional: list of reactions
            parameters (:obj:`list`, optional): list of parameters
            references (:obj:`list`, optional): list of references
        """

        self.compartments = compartments or []
        self.species = species or []
        self.submodels = submodels or []
        self.reactions = reactions or []
        self.parameters = parameters or []
        self.references = references or []

    def get_component_by_id(self, id, component_type=''):
        """ Find model component with id.

        Args:
            id (:obj:`str`): id of component to find
            component_type (:obj:`str`, optional): type of component to search for; if empty search over all components

        Returns:
            :obj:`object`: component with id, or `None` if there is no component with the id
        """

        # components to search over
        if component_type in ['compartments', 'species', 'submodels', 'reactions', 'parameters', 'references']:
            components = getattr(self, component_type)
        elif not component_type:
            components = chain(self.compartments, self.species, self.submodels,
                               self.reactions, self.parameters, self.references)
        else:
            raise Exception('Invalid component type "{}"'.format(component_type))

        # find component
        return next((component for component in components if component.id == id), None)

    def summary(self):
        '''Provide a string summary of the contents of a model.'''
        counts = []
        for attr in 'submodels compartments species reactions parameters references'.split():
            counts.append("{}: {}".format(attr, len(getattr(self, attr))))
        return "Model contains:\n{}".format('\n'.join(counts))


class Submodel(object):
    """ Represents a submodel.

    Attributes:
        id (:obj:`str`): unique id
        name (:obj:`str`): name
        algorithm (:obj:`str`): algorithm
        species (:obj:`list`): list of species in submodel
        reactions (:obj:`list`): list of reactions in submodel
        parameters (:obj:`list`): list of parameers in submodel
    """

    def __init__(self, id, name='', algorithm='SSA', species=None, reactions=None, parameters=None):
        """ Construct a submodel.

        Args:
            id (:obj:`str`): unique id
            name (:obj:`str`, optional): name
            algorithm (:obj:`str`, optional): algorithm
            species (:obj:`list`, optional): list of species in submodel
            reactions (:obj:`list`, optional): list of reactions in submodel
            parameters (:obj:`list`, optional): list of parameers in submodel
        """

        self.id = id
        self.name = name
        self.algorithm = algorithm
        self.species = species or []
        self.reactions = reactions or []
        self.parameters = parameters or []


class Compartment(object):
    """ Represents a compartment. 

    Attributes:
        id (:obj:`str`): unique id
        name (:obj:`str`): name
        initial_volume (:obj:`float`): initial volume (L)
        comments (:obj:`str`): comments
    """

    def __init__(self, id, name='', initial_volume=float('nan'), comments=''):
        """ Construct a compartment.

        Args:
            id (:obj:`str`): unique id
            name (:obj:`str`, optional): name
            initial_volume (:obj:`float`, optional): initial volume (L)
            comments (:obj:`str`, optional): comments
        """
        self.id = id
        self.name = name
        self.initial_volume = initial_volume
        self.comments = comments


class Species(object):
    """ Represents a species

    Attributes:
        id                 (:obj:`str`):   id
        name               (:obj:`str`):   name
        structure          (:obj:`str`):   structure (InChI for metabolites; sequence for DNA, RNA, proteins)
        empirical_formula  (:obj:`str`):   empirical formula
        molecular_weight   (:obj:`float`): molecular weight
        charge             (:obj:`int`):   charge
        type               (:obj:`str`):   type (metabolite, RNA, protein)
        concentrations     (:obj:`list`):  list of comparments in which the species has non-zero concentrations
        cross_refs         (:obj:`list`):  list of cross references
        comments           (:obj:`str`):   comments
    """

    def __init__(self, id, name='', structure='', empirical_formula='', molecular_weight=float('nan'),
                 charge=float('nan'), type='', concentrations=None, cross_refs=None, comments=''):
        """ Construct a species.

        Args:
            id                 (:obj:`str`):             id
            name               (:obj:`str`, optional):   name
            structure          (:obj:`str`, optional):   structure (InChI for metabolites; sequence for DNA, RNA, proteins)
            empirical_formula  (:obj:`str`, optional):   empirical formula
            molecular_weight   (:obj:`float`, optional): molecular weight
            charge             (:obj:`int`, optional):   charge
            type               (:obj:`str`, optional):   type (metabolite, RNA, protein)
            concentrations     (:obj:`list`, optional):  list of comparments in which the species has non-zero concentrations
            cross_refs         (:obj:`list`, optional):  list of cross references
            comments           (:obj:`str`, optional):   comments
        """

        self.id = id
        self.name = name
        self.structure = structure
        self.empirical_formula = empirical_formula
        self.molecular_weight = molecular_weight
        self.charge = charge
        self.type = type
        self.concentrations = concentrations or []
        self.cross_refs = cross_refs or []
        self.comments = comments

    def is_carbon_containing(self):
        """ Returns `True` is species contains at least one carbon atom.

        Returns:
            :obj:`bool`: `True` is species contains at least one carbon atom.
        """
        return self.empirical_formula and self.empirical_formula.upper().find('C') != -1


class Reaction(object):
    """ Represents a reaction. 

    Attributes:
        id            (:obj:`str`): id
        name          (:obj:`str`): name
        submodel      (:obj:`wc_lang.core.Submodel`): submodel that reaction belongs to
        reversible    (:obj:`bool`): indicates if reaction is thermodynamically reversible
        participants  (:obj:`list`): list of reaction participants
        enzyme        (:obj:`wc_lang.core.SpeciesCompartment`): enzyme
        rate_law      (:obj:`wc_lang.core.RateLaw`): rate law
        vmax          (:obj:`float`): vmax
        km            (:obj:`float`): km
        cross_refs    (:obj:`list`): list of cross references to external databases
        comments      (:obj:`str`): comments
    """

    def __init__(self, id, name='', submodel=None, reversible=False, participants=None,
                 enzyme=None, rate_law=None, vmax=float('nan'), km=float('nan'), cross_refs=None, comments=''):
        """ Construct a reaction.

        Args:
            id            (:obj:`str`): id
            name          (:obj:`str`, optional): name
            submodel      (:obj:`wc_lang.core.Submodel`, optional): submodel that reaction belongs to
            reversible    (:obj:`bool`, optional): indicates if reaction is thermodynamically reversible
            participants  (:obj:`list`, optional): list of reaction participants
            enzyme        (:obj:`wc_lang.core.SpeciesCompartment`, optional): enzyme
            rate_law      (:obj:`wc_lang.core.RateLaw`, optional): rate law
            vmax          (:obj:`float`, optional): vmax
            km            (:obj:`float`, optional): km
            cross_refs    (:obj:`list`, optional): list of cross references to external databases
            comments      (:obj:`str`, optional): comments
        """

        self.id = id
        self.name = name
        self.submodel = submodel
        self.reversible = reversible
        self.participants = participants or []
        self.enzyme = enzyme
        self.rate_law = rate_law
        self.vmax = vmax
        self.km = km
        self.cross_refs = cross_refs or []
        self.comments = comments

    def __str__(self):
        """ Generate string representation of reaction stoichometry.

        Returns:
            :obj:`str`: string representation of reaction stoichometry
        """
        global_comp = self.participants[0].compartment
        for part in self.participants:
            if part.compartment != global_comp:
                global_comp = None
                break

        lhs = []
        rhs = []
        for part in self.participants:
            if part.coefficient < 0:
                part_str = ''
                if part.coefficient != -1:
                    if math.ceil(part.coefficient) == part.coefficient:
                        part_str += '({:d}) '.format(-part.coefficient)
                    else:
                        part_str += '({:e}) '.format(-part.coefficient)
                part_str += part.species.id
                if global_comp is None:
                    part_str += '[{}]'.format(part.compartment.id)
                lhs.append(part_str)
            else:
                part_str = ''
                if part.coefficient != 1:
                    if math.ceil(part.coefficient) == part.coefficient:
                        part_str += '({:d}) '.format(part.coefficient)
                    else:
                        part_str += '({:e}) '.format(part.coefficient)
                part_str += part.species.id
                if global_comp is None:
                    part_str += '[{}]'.format(part.compartment.id)
                rhs.append(part_str)

        stoich_str = ''
        if global_comp is not None:
            stoich_str += '[{}]: '.format(global_comp.id)
        stoich_str += '{0} {1}==> {2}'.format(' + '.join(lhs), '<' if self.reversible else '', ' + '.join(rhs))

        return stoich_str


class Parameter(object):
    """ Represents a model parameter.

    Attributes:
        id        (:obj:`str`): id
        name      (:obj:`str`): name
        submodel  (:obj:`wc_lang.core.Submodel`): submodel
        value     (:obj:`float`): value
        units     (:obj:`str`): units of value
        comments  (:obj:`str`): comments    
    """

    def __init__(self, id, name='', submodel=None, value=float('nan'), units='', comments=''):
        """ Construct a parameter.

        Args:
            id        (:obj:`str`): id
            name      (:obj:`str`, optional): name
            submodel  (:obj:`wc_lang.core.Submodel`, optional): submodel
            value     (:obj:`float`, optional): value
            units     (:obj:`str`, optional): units of value
            comments  (:obj:`str`, optional): comments
        """
        self.id = id
        self.name = name
        self.submodel = submodel
        self.value = value
        self.units = units
        self.comments = comments


class Reference(object):
    """ Represents a reference.

    Attributes:
        id          (:obj:`str`): id
        name        (:obj:`str`): name
        cross_refs  (:obj:`list`): list of cross references to external sources
        comments    (:obj:`str`): comments
    """

    def __init__(self, id, name='', cross_refs=None, comments=''):
        """ Construct a reference.

        Args:
            id          (:obj:`str`): id
            name        (:obj:`str`, optional): name
            cross_refs  (:obj:`list`, optional): list of cross references to external sources
            comments    (:obj:`str`, optional): comments
        """
        self.id = id
        self.name = name
        self.cross_refs = cross_refs or []
        self.comments = comments


class Concentration(object):
    """ Represents a concentration in a compartment.

    Attributes:
        compartment  (:obj:`wc_lang.core.Compartment`): compartment in which concentration occurs
        value        (:obj:`float`): concentration in compartment
    """

    def __init__(self, compartment, value=float('nan')):
        """ Construct a concentration

        Args:
            compartment  (:obj:`wc_lang.core.Compartment`): compartment in which concentration occurs
            value        (:obj:`float`, optional): concentration in compartment
        """
        self.compartment = compartment
        self.value = value


class SpeciesCompartment(object):
    """ Represents a participant (species/compartment tuple) in a submodel.

    Attributes:
        species      (:obj:`wc_lang.core.Species`): species
        compartment  (:obj:`wc_lang.core.Compartment`): compartment
    """

    def __init__(self, species, compartment):
        """ Construct a species/compartment tuple.

        Args:
            species      (:obj:`wc_lang.core.Species`): species
            compartment  (:obj:`wc_lang.core.Compartment`): compartment
        """
        self.species = species
        self.compartment = compartment

    @property
    def id(self):
        return '{0}[{1}]'.format(self.species.id, self.compartment.id)

    @property
    def name(self):
        return '{0} ({1})'.format(self.species.name, self.compartment.name)


class ReactionParticipant(object):
    """ Represents a participant (species/compartment tuple) in a reaction.

    Attributes:
        species      (:obj:`wc_lang.core.Species`): species
        compartment  (:obj:`wc_lang.core.Compartment`): compartment
        coefficient  (:obj:`float`): coefficient of species/compartment in reaction
    """

    def __init__(self, species, compartment, coefficient=float('nan')):
        """ Construct a reaction participant (tuple of species, compartment, coefficient).

        Args:
            species      (:obj:`wc_lang.core.Species`): species
            compartment  (:obj:`wc_lang.core.Compartment`): compartment
            coefficient  (:obj:`float`, optional): coefficient of species/compartment in reaction
        """
        self.species = species
        self.compartment = compartment
        self.coefficient = coefficient

    @property
    def id(self):
        return '{0}[{1}]'.format(self.species.id, self.compartment.id)

    @property
    def name(self):
        return '{0} ({1})'.format(self.species.name, self.compartment.name)

    def __str__(self):
        return '({:f}) {:s}[{:s}]'.format(self.coefficient, self.species.id, self.compartment.id)


class RateLaw(object):
    """ Represents a rate law.

    Attributes:
        law  (:obj:`str`):rate law
    """

    def __init__(self, law=''):
        """ Construct a rate law.

        Args:
            law  (:obj:`str`):rate law
        """
        self.law = law

    def get_modifiers(self, model):
        """ Get modifiers (species/compartment tuples) of rate law.

        Args:
            model  (:obj:`wc_lang.core.Model`): model

        Returns:
            :obj:`list`: list of modifiers
        """

        #.. todo:: improve parsing: fully parse rate law and check that all named entities are defined in model
        modifiers = []
        for spec in model.species:
            for comp in model.compartments:
                id = '{0}[{1}]'.format(spec.id, comp.id)
                if self.law.find(id) != -1:
                    modifiers.append(id)

        # return list of modifiers
        return modifiers


class CrossReference(object):
    """ Represents a cross reference to an external database.

    Attributes:
        source  (:obj:`str`): ID of external source e.g. DOI, URI, URL, database name
        id      (:obj:`str`): ID of reference object in external source
    """

    def __init__(self, source, id):
        """ Construct a cross reference.

        Args:
            source  (:obj:`str`): ID of external source e.g. DOI, URI, URL, database name
            id      (:obj:`str`): ID of reference object in external source
        """

        self.source = source
        self.id = id

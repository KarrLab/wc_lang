""" Tests of input/output.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-11-10
:Copyright: 0016, Karr Lab
:License: MIT
"""

from wc_lang.core import (Model, Taxon, TaxonRank, Submodel, Reaction, SpeciesType, SpeciesTypeType, Species, Compartment,
                          ReactionParticipant, Parameter, Reference, ReferenceType, CrossReference,
                          RateLaw, RateLawEquation, SubmodelAlgorithm, Concentration)
from wc_lang.io import ExcelIo
import os
import tempfile
import unittest


class TestIo(unittest.TestCase):

    def setUp(self):
        self.model = mdl = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1b')

        mdl.taxon = Taxon(id='taxon', name='test taxon', rank=TaxonRank['species'])

        self.comp_0 = comp_0 = mdl.compartments.create(id='comp_0', name='compartment 0', initial_volume=1.25)
        self.comp_1 = comp_1 = mdl.compartments.create(id='comp_1', name='compartment 1', initial_volume=2.5)
        self.compartments = compartments = [comp_0, comp_1]

        self.species_types = species_types = []
        self.species = species = []
        self.concentrations = concentrations = []
        for i in range(8):
            spec_type = mdl.species_types.create(
                id='spec_type_{}'.format(i),
                name='species type {}'.format(i),
                type=SpeciesTypeType['metabolite'],
                structure='C' * (i + 1),
                empirical_formula='C' + str(i + 1),
                molecular_weight=12 * (i + 1),
                charge=i + 1)
            species_types.append(spec_type)

            if i != 3:
                spec = Species(species_type=spec_type, compartment=comp_0)
            else:
                spec = Species(species_type=spec_type, compartment=comp_1)
            species.append(spec)

            conc = Concentration(species=spec, value=3 * i)
            concentrations.append(conc)

        self.submdl_0 = submdl_0 = mdl.submodels.create(
            id='submodel_0', name='submodel 0', algorithm=SubmodelAlgorithm['ssa'])
        self.submdl_1 = submdl_1 = mdl.submodels.create(
            id='submodel_1', name='submodel 1', algorithm=SubmodelAlgorithm['ssa'])
        self.submdl_2 = submdl_2 = mdl.submodels.create(
            id='submodel_2', name='submodel 2', algorithm=SubmodelAlgorithm['dfba'])
        self.submodels = submodels = [submdl_0, submdl_1, submdl_2]

        self.rxn_0 = rxn_0 = submdl_0.reactions.create(id='rxn_0', name='reaction 0')
        rxn_0.participants.create(species=species[0], coefficient=-2)
        rxn_0.participants.create(species=species[1], coefficient=-3)
        rxn_0.participants.create(species=species[2], coefficient=1)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[5].serialize()),
            modifiers=set(species[5:6]))
        rate_law_0 = rxn_0.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.rxn_1 = rxn_1 = submdl_1.reactions.create(id='rxn_1', name='reaction 1')
        rxn_1.participants.create(species=species[0], coefficient=-2)
        rxn_1.participants.create(species=species[1], coefficient=-3)
        rxn_1.participants.create(species=species[3], coefficient=2)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[6].serialize()),
            modifiers=set(species[6:7]))
        rate_law_1 = rxn_1.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.rxn_2 = rxn_2 = submdl_2.reactions.create(id='rxn_2', name='reaction 2')
        rxn_2.participants.create(species=species[0], coefficient=-2)
        rxn_2.participants.create(species=species[1], coefficient=-3)
        rxn_2.participants.create(species=species[4], coefficient=1)
        equation = RateLawEquation(
            expression='k_cat * {0} / (k_m + {0})'.format(species[7].serialize()),
            modifiers=set(species[7:8]))
        rate_law_2 = rxn_2.rate_laws.create(equation=equation, k_cat=2, k_m=1)

        self.reactions = [rxn_0, rxn_1, rxn_2]
        self.rate_laws = [rate_law_0, rate_law_1, rate_law_2]

        self.parameters = parameters = []
        self.references = references = []
        self.cross_references = cross_references = []
        for i in range(3):
            param = mdl.parameters.create(
                id='param_{}'.format(i), name='parameter {}'.format(i),
                value=i * 4, units='dimensionless')
            param.submodels = submodels[i:i + 1]
            parameters.append(param)

            ref = param.references.create(
                id='ref_{}'.format(i), name='reference {}'.format(i),
                type=ReferenceType['misc'])
            references.append(ref)

            x_ref = ref.cross_references.create(database='x', id='y' * (i + 1),
                                                url='http://x.com/{}'.format('y' * (i + 1)))
            cross_references.append(x_ref)

        _, self.filename = tempfile.mkstemp(suffix='.xlsx')

    def tearDown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)

    def test_write_read(self):
        ExcelIo.write(self.filename, self.model)
        model = ExcelIo.read(self.filename)

        self.assertEqual(model, self.model)
        self.assertEqual(self.model.difference(model), '')

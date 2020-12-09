""" Test creation of implicit exchange reactions for dFBA submodels.

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-29
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from wc_lang import Model, ChemicalStructure
from wc_lang.transform import CreateImplicitDfbaExchangeReactionsTransform
from wc_utils.util.chem import EmpiricalFormula
from wc_onto import onto
import unittest
import wc_lang.config.core


class CreateImplicitDfbaExchangeReactionsTransformTestCase(unittest.TestCase):
    def test(self):
        config = wc_lang.config.core.get_config()['wc_lang']['dfba']

        model = Model()
        submodel = model.submodels.create(id='submdl', name='submodel', framework=onto['WC:dynamic_flux_balance_analysis'])

        comps = [
            model.compartments.create(id='c', name='cytosol'),
            model.compartments.create(id='d', name='dna'),
            model.compartments.create(id='e', name='extracellular space'),
        ]
        sts = [
            model.species_types.create(id='st_1', name='species type 1',
                                       structure=ChemicalStructure(empirical_formula=EmpiricalFormula('C'))),
            model.species_types.create(id='st_2', name='species type 2',
                                       structure=ChemicalStructure(empirical_formula=EmpiricalFormula('H'))),
            model.species_types.create(id='st_3', name='species type 3',
                                       structure=ChemicalStructure(empirical_formula=EmpiricalFormula('O'))),
        ]
        specs = []
        for st in sts:
            spec_comps = []
            specs.append(spec_comps)
            for comp in comps:
                spec_comps.append(model.species.create(species_type=st, compartment=comp))

        rxn = model.reactions.create(submodel=submodel)
        rxn.participants.create(species=specs[0][0])
        rxn.participants.create(species=specs[0][1])
        rxn.participants.create(species=specs[0][2])

        rxn = model.reactions.create(submodel=submodel)
        rxn.participants.create(species=specs[1][0])
        rxn.participants.create(species=specs[1][1])
        rxn.participants.create(species=specs[1][2])

        rxn = model.reactions.create(submodel=submodel)
        rxn.participants.create(species=specs[2][0])
        rxn.participants.create(species=specs[2][1])

        transform = CreateImplicitDfbaExchangeReactionsTransform()
        transform.run(model)

        self.assertEqual(len(model.reactions), 5)
        self.assertEqual(len(submodel.reactions), 5)
        self.assertNotEqual(model.reactions.get_one(id='__dfba_ex_submdl_st_1_e'), None)
        self.assertNotEqual(model.reactions.get_one(id='__dfba_ex_submdl_st_2_e'), None)
        self.assertEqual(model.reactions.get_one(id='__dfba_ex_submdl_st_3_e'), None)

        rxn = model.reactions.get_one(id='__dfba_ex_submdl_st_1_e')
        self.assertEqual(rxn.name, 'dFBA exchange (submodel, species type 1, extracellular space)')
        self.assertEqual(len(rxn.participants), 1)
        self.assertEqual(rxn.participants[0].species, specs[0][2])
        self.assertEqual(rxn.participants[0].coefficient, 1.)
        self.assertEqual(rxn.reversible, True)
        self.assertEqual(rxn.flux_bounds.min, -config['flux_bounds']['ex_carbon'])
        self.assertEqual(rxn.flux_bounds.max, config['flux_bounds']['ex_carbon'])

        rxn = model.reactions.get_one(id='__dfba_ex_submdl_st_2_e')
        self.assertEqual(rxn.flux_bounds.min, -config['flux_bounds']['ex_no_carbon'])
        self.assertEqual(rxn.flux_bounds.max, config['flux_bounds']['ex_no_carbon'])

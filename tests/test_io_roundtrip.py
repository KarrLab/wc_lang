""" Tests of io round tripping

:Author: Arthur Goldberg, Arthur.Goldberg@mssm.edu
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-18
:Copyright: 2018, Karr Lab
:License: MIT
"""

import obj_model
import os
import pytest
import re
import shutil
import tempfile
import unittest

from wc_lang.core import (ReactionRateUnit, Model, SpeciesCoefficient, Expression, Species, Observable, Function,
                          DistributionInitConcentration, ConcentrationUnit, RateLaw, RateLawDirection, RateLawExpression,
                          Parameter)
from wc_lang.io import Reader, Writer
from wc_utils.util.chem import EmpiricalFormula


class RoundTripTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_write_error(self):
        model = Model(id='test_model', version='0.0.0')
        comp = model.compartments.create(id='compartment_1')
        comp.init_density = model.parameters.create(id='density_compartment_1', value=1100, units='g l^-1')
        species_type_1 = model.species_types.create(
            id='species_type_1',
            empirical_formula=EmpiricalFormula('CHO'),
            charge=1)
        species_type_2 = model.species_types.create(
            id='species_type_2',
            empirical_formula=EmpiricalFormula('C2H2O2'),
            charge=2)
        species_1 = comp.species.create(species_type=species_type_1,
                                        id='', model=model)
        species_2 = comp.species.create(species_type=species_type_2,
                                        id='', model=model)
        submdl = model.submodels.create(id='submodel_1')

        # create a DistributionInitConcentration so that Species are provided to ExpressionAttribute.deserialize()
        species_1.distribution_init_concentration = DistributionInitConcentration(
            model=model,
            mean=1, units=ConcentrationUnit.M)
        species_1.distribution_init_concentration.id = species_1.distribution_init_concentration.gen_id()
        objects = {Species: {}}
        objects[Species][species_1.id] = species_1
        observable_1 = Expression.make_obj(model, Observable, 'observable_1', species_1.id, objects)

        rxn_species_coeffs = [
            species_1.species_coefficients.create(coefficient=-2.),
            species_2.species_coefficients.create(coefficient=1.),
        ]
        rxn = submdl.reactions.create(id='reaction_1', model=model)
        rxn.participants.extend(rxn_species_coeffs)
        rl = rxn.rate_laws.create(direction=RateLawDirection.forward,
                                  units=ReactionRateUnit['s^-1'],
                                  model=model)
        rl.id = rl.gen_id()
        param_1 = model.parameters.create(id='param_1', value=1., units='s^-1')
        rl.expression, error = RateLawExpression.deserialize('param_1', {Parameter: {'param_1': param_1}})
        self.assertEqual(error, None)

        errors = obj_model.Validator().run(model, get_related=True)
        self.assertNotEqual(errors, None)

        filename = os.path.join(self.tmp_dir, 'model.xlsx')
        with pytest.warns(obj_model.io.IoWarning, match='objects are not valid'):
            Writer().run(filename, model, set_repo_metadata_from_path=False)
        with self.assertRaisesRegex(ValueError, re.escape('contains error(s)')):
            Reader().run(filename)

    def test_successful_roundtrip(self):
        """
        Args:
            species_coefficient_creation_method (:obj:`str`): name of method to use to get or create
                an instance of :obj:`wc_lang.SpeciesCoefficient`
        """
        model = Model(id='test_model', version='0.0.0')
        comp = model.compartments.create(id='compartment_1')
        comp.init_density = model.parameters.create(id='density_compartment_1', value=1100, units='g l^-1')
        species_type_1 = model.species_types.create(
            id='species_type_1',
            empirical_formula=EmpiricalFormula('CHO'),
            charge=1)
        species_type_2 = model.species_types.create(
            id='species_type_2',
            empirical_formula=EmpiricalFormula('C3H3O3'),
            charge=3)
        species_1 = comp.species.create(species_type=species_type_1,
                                        model=model)
        species_2 = comp.species.create(species_type=species_type_2,
                                        model=model)
        species_1.id = species_1.gen_id()
        species_2.id = species_2.gen_id()
        submdl = model.submodels.create(id='submodel_1')

        # create a DistributionInitConcentration so that Species are provided to ExpressionAttribute.deserialize()
        species_1.distribution_init_concentration = DistributionInitConcentration(
            model=model,
            mean=1, units=ConcentrationUnit.M)
        species_1.distribution_init_concentration.id = species_1.distribution_init_concentration.gen_id()
        objects = {Species: {}}
        objects[Species][species_1.id] = species_1
        observable_1 = Expression.make_obj(model, Observable, 'observable_1', species_1.id, objects)
        objects = {Observable: {'observable_1': observable_1}}
        observable_2 = Expression.make_obj(model, Observable, 'observable_2', 'obs_1', objects)

        param_1 = model.parameters.create(id='param_1', value=1., units='dimensionless')
        param_2 = model.parameters.create(id='param_2', value=1., units='s^-1')
        objects = {
            Parameter: {
                'param_1': param_1,
                'param_2': param_2,
            },
        }
        func_1 = Expression.make_obj(model, Function, 'func_1', 'param_1', objects)

        rxn_species_coeffs = [
            species_1.species_coefficients.get_or_create(coefficient=-3.),
            species_2.species_coefficients.get_or_create(coefficient=1.),
        ]
        rxn = submdl.reactions.create(id='reaction_1', model=model)
        rxn.participants.extend(rxn_species_coeffs)
        rl = rxn.rate_laws.create(direction=RateLawDirection.forward,
                                  units=ReactionRateUnit['s^-1'],
                                  model=model)
        rl.id = rl.gen_id()
        rl.expression, error = RateLawExpression.deserialize('param_2', objects)
        self.assertEqual(error, None)

        errors = obj_model.Validator().run(model, get_related=True)
        self.assertEqual(errors, None, str(errors))

        filename = os.path.join(self.tmp_dir, 'model.xlsx')
        Writer().run(filename, model, set_repo_metadata_from_path=False)
        model_2 = Reader().run(filename)[Model][0]

        self.assertEqual(model.difference(model_2), '')
        self.assertTrue(model_2.is_equal(model))
        self.assertTrue(model.is_equal(model_2))

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

from wc_lang.core import (Model, SpeciesCoefficient, Expression, Species, Observable, Function,
                          Concentration, ConcentrationUnit, RateLaw, RateLawDirection, RateLawExpression,
                          Parameter)
from wc_lang.io import Reader, Writer


class RoundTripTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_write_error(self):
        model = Model(id='test_model', version='0.0.0')
        comp = model.compartments.create(id='compartment_1')
        species_type_1 = model.species_types.create(
            id='species_type_1',
            empirical_formula='CHO',
            charge=1)
        species_type_2 = model.species_types.create(
            id='species_type_2',
            empirical_formula='C2H2O2',
            charge=2)
        species_1 = comp.species.create(species_type=species_type_1,
                                        id='', model=model)
        species_2 = comp.species.create(species_type=species_type_2,
                                        id='', model=model)
        submdl = model.submodels.create(id='submodel_1')

        # create a Concentration so that Species are provided to ExpressionAttribute.deserialize()
        species_1.concentration = Concentration(
            id=Concentration.gen_id(species_1.id),
            model=model,
            mean=1, units=ConcentrationUnit.M)
        objects = {Species: {}}
        objects[Species][species_1.id] = species_1
        observable_1 = Expression.make_obj(model, Observable, 'observable_1', species_1.id, objects)

        rxn_species_coeffs = [
            species_1.species_coefficients.create(coefficient=-2.),
            species_2.species_coefficients.create(coefficient=1.),
        ]
        rxn = submdl.reactions.create(id='reaction_1', model=model)
        rxn.rate_laws.create(id=RateLaw.gen_id(rxn.id, RateLawDirection.forward.name),
                             direction=RateLawDirection.forward, model=model,
                             expression=RateLawExpression(expression='1.'))
        rxn.participants.extend(rxn_species_coeffs)

        errors = obj_model.Validator().run(model, get_related=True)
        self.assertNotEqual(errors, None)

        filename = os.path.join(self.tmp_dir, 'model.xlsx')
        with pytest.warns(obj_model.io.IoWarning, match='objects are not valid'):
            Writer().run(model, filename, set_repo_metadata_from_path=False)
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
        species_type_1 = model.species_types.create(
            id='species_type_1',
            empirical_formula='CHO',
            charge=1)
        species_type_2 = model.species_types.create(
            id='species_type_2',
            empirical_formula='C3H3O3',
            charge=3)
        species_1 = comp.species.create(species_type=species_type_1,
                                        id=Species.gen_id(species_type_1.id, comp.id),
                                        model=model)
        species_2 = comp.species.create(species_type=species_type_2,
                                        id=Species.gen_id(species_type_2.id, comp.id),
                                        model=model)
        submdl = model.submodels.create(id='submodel_1')

        # create a Concentration so that Species are provided to ExpressionAttribute.deserialize()
        species_1.concentration = Concentration(
            id=Concentration.gen_id(species_1.id),
            model=model,
            mean=1, units=ConcentrationUnit.M)
        objects = {Species: {}}
        objects[Species][species_1.id] = species_1
        observable_1 = Expression.make_obj(model, Observable, 'observable_1', species_1.id, objects)
        objects = {Observable: {'observable_1': observable_1}}
        observable_2 = Expression.make_obj(model, Observable, 'observable_2', 'obs_1', objects)

        param = model.parameters.create(id='param_1', value=1.)
        objects = {Parameter: {'param_1': param}}
        fun_1 = Expression.make_obj(model, Function, 'fun_1', 'param_1', objects)

        rxn_species_coeffs = [
            species_1.species_coefficients.get_or_create(coefficient=-3.),
            species_2.species_coefficients.get_or_create(coefficient=1.),
        ]
        rxn = submdl.reactions.create(id='reaction_1', model=model)
        rxn.participants.extend(rxn_species_coeffs)
        rxn.rate_laws.create(id=RateLaw.gen_id(rxn.id, RateLawDirection.forward.name),
                             direction=RateLawDirection.forward, model=model,
                             expression=RateLawExpression(expression='1.'))

        errors = obj_model.Validator().run(model, get_related=True)
        self.assertEqual(errors, None, str(errors))

        filename = os.path.join(self.tmp_dir, 'model.xlsx')
        Writer().run(model, filename, set_repo_metadata_from_path=False)
        model_2 = Reader().run(filename)

        self.assertEqual(model.difference(model_2), '')
        self.assertTrue(model_2.is_equal(model))
        self.assertTrue(model.is_equal(model_2))

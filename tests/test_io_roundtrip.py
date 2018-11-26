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

from wc_lang.core import (Model, SpeciesCoefficient, ExpressionMethods, Species, Observable, Function,
                          Concentration, ConcentrationUnit, RateLaw, RateLawDirection, RateLawEquation,
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
        species_type = model.species_types.create(id='species_type_1')
        species = comp.species.create(species_type=species_type)
        species.id = ''
        species.model = model
        submdl = model.submodels.create(id='submodel_1')

        # create a Concentration so that Species are provided to ExpressionAttribute.deserialize()
        species.concentration = Concentration(
            id=Concentration.gen_id(species.id),
            model=model,
            value=1, units=ConcentrationUnit.M)
        objects = {Species: {}}
        objects[Species][species.id] = species
        observable_1 = ExpressionMethods.make_obj(model, Observable, 'observable_1', species.id, objects)

        coefficient = 1
        rxn_species_coeff = species.species_coefficients.create(coefficient=coefficient)
        rxn = submdl.reactions.create(id='reaction_1', model=model)
        rxn.rate_laws.create(id=RateLaw.gen_id(rxn.id, RateLawDirection.forward.name),
                             direction=RateLawDirection.forward, model=model,
                             equation=RateLawEquation(expression='1.'))
        rxn.participants.append(rxn_species_coeff)

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
        species_type = model.species_types.create(id='species_type_1')
        species = comp.species.create(species_type=species_type)
        species.id = species.gen_id(species_type.id, comp.id)
        species.model = model
        submdl = model.submodels.create(id='submodel_1')

        # create a Concentration so that Species are provided to ExpressionAttribute.deserialize()
        species.concentration = Concentration(
            id=Concentration.gen_id(species.id),
            model=model,
            value=1, units=ConcentrationUnit.M)
        objects = {Species: {}}
        objects[Species][species.id] = species
        observable_1 = ExpressionMethods.make_obj(model, Observable, 'observable_1', species.id, objects)
        objects = {Observable: {'observable_1': observable_1}}
        ExpressionMethods.make_obj(model, Observable, 'observable_2', 'obs_1', objects)

        param = model.parameters.create(id='param_1')
        objects = {Parameter: {'param_1': param}}
        ExpressionMethods.make_obj(model, Function, 'fun_1', 'param_1', objects)

        coefficient = 1

        rxn_species_coeff = species.species_coefficients.get_or_create(coefficient=coefficient)
        rxn = submdl.reactions.create(id='reaction_1', model=model)
        rxn.participants.append(rxn_species_coeff)
        rxn.rate_laws.create(id=RateLaw.gen_id(rxn.id, RateLawDirection.forward.name),
                             direction=RateLawDirection.forward, model=model,
                             equation=RateLawEquation(expression='1.'))

        errors = obj_model.Validator().run(model, get_related=True)
        self.assertEqual(errors, None, str(errors))

        filename = os.path.join(self.tmp_dir, 'model.xlsx')
        Writer().run(model, filename, set_repo_metadata_from_path=False)
        model_2 = Reader().run(filename)

        self.assertEqual(model.difference(model_2), '')
        self.assertTrue(model_2.is_equal(model))
        self.assertTrue(model.is_equal(model_2))

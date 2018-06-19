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
import shutil
import tempfile
import unittest

from wc_lang.core import Model, SpeciesCoefficient
from wc_lang.io import Reader, Writer


class RoundTripTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    #@unittest.expectedFailure
    def test_create(self):
        model = Model(id='test_model', version='0.0.0')
        comp = model.compartments.create(id='compartment_1')
        species_type = model.species_types.create(id='species_type_1')
        species = comp.species.create(species_type=species_type)
        submdl = model.submodels.create(id='submodel_1', compartment=comp)

        coefficient = 1

        obs_species_coeff = species.species_coefficients.create(coefficient=coefficient)
        obs = model.observables.create(id='observable_1')
        obs.species.append(obs_species_coeff)

        rxn_species_coeff = species.species_coefficients.create(coefficient=coefficient)
        rxn = submdl.reactions.create(id='reaction_1')
        rxn.participants.append(rxn_species_coeff)

        self.assertNotEqual(obj_model.Validator().run(model, get_related=True), None)

        filename = os.path.join(self.tmp_dir, 'model.xlsx')
        with pytest.warns(obj_model.io.IoWarning, match='objects are not valid'):
            Writer().run(model, filename, set_repo_metadata_from_path=False)
        model_2 = Reader().run(filename)

        self.assertNotEqual(model.difference(model_2), '')
        self.assertFalse(model_2.is_equal(model))
        self.assertFalse(model.is_equal(model_2))

    def test_get_or_create(self):
        """
        Args:
            species_coefficient_creation_method (:obj:`str`): name of method to use to get or create 
                an instance of :obj:`wc_lang.SpeciesCoefficient`
        """
        model = Model(id='test_model', version='0.0.0')
        comp = model.compartments.create(id='compartment_1')
        species_type = model.species_types.create(id='species_type_1')
        species = comp.species.create(species_type=species_type)
        submdl = model.submodels.create(id='submodel_1', compartment=comp)

        coefficient = 1

        obs_species_coeff = species.species_coefficients.get_or_create(coefficient=coefficient)
        obs = model.observables.create(id='observable_1')
        obs.species.append(obs_species_coeff)

        rxn_species_coeff = species.species_coefficients.get_or_create(coefficient=coefficient)
        rxn = submdl.reactions.create(id='reaction_1')
        rxn.participants.append(rxn_species_coeff)

        self.assertEqual(obj_model.Validator().run(model, get_related=True), None)

        filename = os.path.join(self.tmp_dir, 'model.xlsx')
        Writer().run(model, filename, set_repo_metadata_from_path=False)
        model_2 = Reader().run(filename)

        self.assertEqual(model.difference(model_2), '')
        self.assertTrue(model_2.is_equal(model))
        self.assertTrue(model.is_equal(model_2))

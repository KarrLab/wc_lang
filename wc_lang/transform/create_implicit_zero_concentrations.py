""" Create the implicit zero concentrations in a model

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-11-28
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from .core import Transform
from wc_lang.core import Concentration, ConcentrationUnit


class CreateImplicitZeroConcentrationsTransform(Transform):
    """ Create the implicit zero concentrations in a model
    """

    class Meta(object):
        id = 'CreateImplicitZeroConcentrations'
        label = 'Create the implicit zero concentrations in a model'

    def run(self, model):
        """ Transform model

        Args:
            model (:obj:`Model`): model

        Returns:
            :obj:`Model`: same model, but transformed
        """
        for species in model.get_species():
            if species.concentration is None:
                species.concentration = Concentration(
                    id=Concentration.gen_id(species.id),
                    model=model,
                    species=species,
                    mean=0.0, std=0.0, units=ConcentrationUnit.molecules)

        return model

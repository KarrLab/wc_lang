""" Test round-trip `wc_lang` model migration

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2019-07-22
:Copyright: 2019, Karr Lab
:License: MIT
"""
from obj_model.migrate import MigrationWrapper


class ReverseParameterTransformation(MigrationWrapper):
    """ These transformations and the ones in ./reverse.py invert each other
    """

    def prepare_existing_models(self, migrator, existing_models):
        pass

    def modify_migrated_models(self, migrator, migrated_models):
        # decrement the value of the 'carbonExchangeRate' Parameter
        for migrated_model in migrated_models:
            if isinstance(migrated_model, migrator.migrated_defs['Parameter']) and\
                migrated_model.id == 'carbonExchangeRate':
                migrated_model.value -= 1

transformations = ReverseParameterTransformation()
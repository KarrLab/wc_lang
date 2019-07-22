""" Test round-trip `wc_lang` model migration

:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2019-07-22
:Copyright: 2019, Karr Lab
:License: MIT
"""
from obj_model.migrate import MigrationWrapper


class ForwardParameterTransformation(MigrationWrapper):
    """ These transformations and the ones in ./reverse.py invert each other
    """

    def prepare_existing_models(self, migrator, existing_models):
        # increment the value of the 'carbonExchangeRate' Parameter
        for existing_model in existing_models:
            if isinstance(existing_model, migrator.existing_defs['Parameter']) and\
                existing_model.id == 'carbonExchangeRate':
                existing_model.value += 1

    def modify_migrated_models(self, migrator, migrated_models):
        pass


transformations = ForwardParameterTransformation()
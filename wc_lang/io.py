""" Reading and writing models to/from files.

* Comma separated values (.csv)
* Excel (.xlsx)
* Tab separated values (.tsv)

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-12-05
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang.core import (Model, Taxon, Submodel, Compartment, SpeciesType,
                          Concentration, Reaction, RateLaw, Parameter, Reference,
                          CrossReference)
from wc_lang.util import get_models
from wc_utils.schema import io


class Writer(object):
    """ Write model to file(s) """

    def run(self, path, model=None):
        """ Write model to file(s)

        Args:
            path (:obj:`str`): path to file(s)
            model (:obj:`Model`): model
        """
        models = [
            Model, Taxon,
            Submodel, Compartment, SpeciesType, Concentration,
            Reaction, RateLaw, Parameter,
            Reference, CrossReference,
        ]

        kwargs = {
            'language': 'wc_lang',
            'creator': '{}.{}'.format(self.__class__.__module__, self.__class__.__name__),
        }
        if model:
            objects = set((model,))
            kwargs['title'] = model.id
            kwargs['description'] = model.name
            kwargs['version'] = model.version
        else:
            objects = set()

        io.Writer().run(path, objects, models, **kwargs)


class Reader(object):
    """ Read model from file(s) """

    def run(self, path):
        """ Read model from file(s)

        Args:
            path (:obj:`str`): path to file(s)

        Returns:
            :obj:`Model`: model
        """
        objects = io.Reader().run(path, get_models(inline=False))

        if not objects[Model]:
            return None

        if len(objects[Model]) > 1:
            raise ValueError('Model file "{}" should only define one model'.format(path))

        return objects[Model].pop()


def create_template(path):
    """ Create file with model template, including row and column headings

    Args:
        path (:obj:`str`): path to file(s)
    """
    Writer().run(path)

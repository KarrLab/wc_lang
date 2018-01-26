""" Reading and writing models to/from files.

Supported file types:

* Comma separated values (.csv)
* Excel (.xlsx)
* Tab separated values (.tsv)

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-12-05
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang.core import (Model, Taxon, Submodel, ObjectiveFunction, Compartment, SpeciesType,
                          Concentration, Reaction, RateLaw, BiomassComponent,
                          BiomassReaction, Parameter, Reference, DatabaseReference)
from wc_lang.util import get_models
from obj_model import io


class Writer(object):
    """ Write model to file(s) """

    model_order = [
        Model, Taxon, Submodel, ObjectiveFunction, Compartment, SpeciesType,
        Concentration, Reaction, RateLaw, BiomassComponent,
        BiomassReaction, Parameter, Reference, DatabaseReference
    ]

    def run(self, path, model=None):
        """ Write model to file(s)

        Args:
            path (:obj:`str`): path to file(s)
            model (:obj:`Model`): model
        """

        kwargs = {
            'language': 'wc_lang',
            'creator': '{}.{}'.format(self.__class__.__module__, self.__class__.__name__),
        }
        if model:
            objects = [model]
            kwargs['title'] = model.id
            kwargs['description'] = model.name
            kwargs['version'] = model.version
        else:
            objects = []

        io.Writer().run(path, objects, self.model_order, **kwargs)


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


def convert(source, destination):
    """ Convert among Excel (.xlsx), comma separated (.csv), and tab separated model formats (.tsv)

    Read a model from the `source` files(s) and write it to the `destination` files(s). A path to a
    delimiter separated set of models must be represented by a Unix glob pattern (with a \*) that
    matches all delimiter separated files.

    Args:
        source (:obj:`str`): path to source file(s)
        destination (:obj:`str`): path to save converted file
    """
    io.convert(source, destination, models=Writer.model_order)


def create_template(path):
    """ Create file with model template, including row and column headings

    Args:
        path (:obj:`str`): path to file(s)
    """
    Writer().run(path)

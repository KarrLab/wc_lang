""" Classes for reading and writing models to/from files.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-12-05
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang.core import (Model, Taxon, Submodel, Compartment, SpeciesType,
                          Concentration, Reaction, RateLaw, Parameter, Reference,
                          CrossReference)
from wc_utils.schema.io import ExcelIo as BaseExcelIo


class ExcelIo(object):
    """ Read/write model to/from Excel workbooks """

    @classmethod
    def write(cls, filename, model):
        """ Write model to file

        Args:
            filename (:obj:`str`): path to file
            model (:obj:`Model`): model
        """
        BaseExcelIo.write(filename, set((model,)), [
            Model, Taxon, Submodel, Compartment, SpeciesType,
            Concentration, Reaction, RateLaw, Parameter, Reference,
            CrossReference],
            title=model.id, description=model.name, version=model.version,
            language='wc_lang', creator=cls.__module__ + '.' + cls.__name__)

    @classmethod
    def read(cls, filename):
        """ Read model from file

        Args:
            filename (:obj:`str`): path to file

        Returns:
            :obj:`Model`: model
        """
        objects = BaseExcelIo.read(filename, [
            Model, Taxon, Submodel, Compartment, SpeciesType,
            Concentration, Reaction, RateLaw, Parameter, Reference,
            CrossReference, ])

        if not objects[Model]:
            return None

        if len(objects[Model]) > 1:
            raise ValueError('Model file "{}" should only define one model'.format(filename))

        return objects[Model].pop()

    @classmethod
    def create_template(cls, filename):
        """ Create file with template (i.e. row, column headings)

        Args:
            filename (:obj:`str`): path to file
        """
        BaseExcelIo.create_template(filename, [
            Model, Taxon, Submodel, Compartment, SpeciesType,
            Concentration, Reaction, RateLaw, Parameter, Reference,
            CrossReference],
            language='wc_lang', creator=cls.__module__ + '.' + cls.__name__)

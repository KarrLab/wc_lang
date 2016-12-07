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
        """ Write model to Excel workbook

        Args:
            filename (:obj:`str`): path to Excel file
            model (:obj:`Model`): model
            keywords (:obj:`str`, optional): keywords
        """
        BaseExcelIo.write(filename, set((model,)), [
            Model, Taxon, Submodel, Compartment, SpeciesType,
            Concentration, Reaction, RateLaw, Parameter, Reference,
            CrossReference],
            title=model.id, description=model.name, version=model.version,
            language='wc_lang', creator=cls.__module__ + '.' + cls.__name__)

    @classmethod
    def read(cls, filename):
        """ Read model from Excel workbook

        Args:
            filename (:obj:`str`): path to Excel file

        Returns:
            :obj:`Model`: model
        """
        objects = BaseExcelIo.read(filename, set((
            Model, Taxon, Submodel, Compartment, SpeciesType,
            Concentration, Reaction, RateLaw, Parameter, Reference,
            CrossReference, )))
        return objects[Model].pop()

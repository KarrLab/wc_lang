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

from wc_lang import core
from wc_lang import util
import obj_model
import os
import wc_lang


class Writer(object):
    """ Write model to file(s) """

    model_order = [
        core.Model, core.Taxon, core.Submodel, core.Compartment, core.SpeciesType,
        core.Concentration, core.Observable, core.Function, core.Reaction, core.RateLaw, core.BiomassComponent,
        core.BiomassReaction, core.Parameter, core.StopCondition, core.Reference, core.DatabaseReference
    ]

    def run(self, model, path):
        """ Write model to file(s)

        Args:            
            model (:obj:`core.Model`): model
            path (:obj:`str`): path to file(s)
        """
        _, ext = os.path.splitext(path)
        obj_model.io.get_writer(ext)().run(path, [model], models=self.model_order, 
            language='wc_lang',
            creator='{}.{}'.format(self.__class__.__module__, self.__class__.__name__),
            title=model.id,
            description=model.name,
            version=model.version)


class Reader(object):
    """ Read model from file(s) """

    def run(self, path, strict=True):
        """ Read model from file(s)

        Args:
            path (:obj:`str`): path to file(s)
            strict (:obj:`str`, optional): if :obj:`True`, validate that the the model file(s) strictly follow the
                :obj:`obj_model` serialization format:

                * The worksheets are in the expected order
                * There are no missing worksheets
                * There are no extra worksheets
                * The columns are in the expected order
                * There are no missing columns
                * There are no extra columns

        Returns:
            :obj:`core.Model`: model

        Raises:
            :obj:`ValueError`: if :obj:`path` defines multiple models
        """
        _, ext = os.path.splitext(path)
        reader = obj_model.io.get_reader(ext)()
        
        kwargs = {}
        if isinstance(reader, obj_model.io.WorkbookReader) and not strict:
            kwargs['ignore_missing_sheets'] = True
            kwargs['ignore_extra_sheets'] = True
            kwargs['ignore_sheet_order'] = True
            kwargs['ignore_missing_attributes'] = True
            kwargs['ignore_extra_attributes'] = True
            kwargs['ignore_attribute_order'] = True
        objects = reader.run(path, models=Writer.model_order, **kwargs)

        if not objects[core.Model]:
            return None

        if len(objects[core.Model]) > 1:
            raise ValueError('Model file "{}" should only define one model'.format(path))

        return objects[core.Model].pop()


def convert(source, destination, strict=True):
    """ Convert among Excel (.xlsx), comma separated (.csv), and tab separated (.tsv) file formats

    Read a model from the `source` files(s) and write it to the `destination` files(s). A path to a
    delimiter separated set of models must be represented by a Unix glob pattern (with a \*) that
    matches all delimiter separated files.

    Args:
        source (:obj:`str`): path to source file(s)
        destination (:obj:`str`): path to save converted file
        strict (:obj:`str`, optional): if :obj:`True`, validate that the the model file(s) strictly follow the
                :obj:`obj_model` serialization format:

                * The worksheets are in the expected order
                * There are no missing worksheets
                * There are no extra worksheets
                * The columns are in the expected order
                * There are no missing columns
                * There are no extra columns
    """
    kwargs = {}
    if not strict:
        kwargs['ignore_missing_sheets'] = True
        kwargs['ignore_extra_sheets'] = True
        kwargs['ignore_sheet_order'] = True
        kwargs['ignore_missing_attributes'] = True
        kwargs['ignore_extra_attributes'] = True
        kwargs['ignore_attribute_order'] = True
    obj_model.io.convert(source, destination, models=Writer.model_order, **kwargs)


def create_template(path):
    """ Create file with model template, including row and column headings

    Args:
        path (:obj:`str`): path to file(s)
    """
    model = core.Model(id='template', name='Template', version=wc_lang.__version__)
    Writer().run(model, path)

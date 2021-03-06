""" Reading and writing models to/from files.

Supported file types:

* Comma separated values (.csv)
* Excel (.xlsx)
* Tab separated values (.tsv)

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2016-12-05
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang import core
from wc_utils.util.string import indent_forest
from wc_utils.util.list import get_count_limited_class
import obj_tables
import obj_tables.io
import obj_tables.utils
import os
import wc_lang
import wc_lang.config.core
from wc_utils.util import git


def get_root_model(models):
    """ Get a root model class from an iterator over `obj_tables.Model` classes

    The root model is obtained by name, rather than from `core.Model`, because `models` may have
    been obtained from `obj_tables` migration, which loads other versions of `wc_lang.core` in
    private imports.

    Args:
        models (:obj:`iterator`): subclasses of :obj:`obj_tables.Model` being read or written

    Returns:
        :obj:`type`: a subclass of `obj_tables.Model` whose name is `root_model_name`
    """
    # the singleton, root obj_tables.Model in wc_lang.core
    ROOT_MODEL_NAME = 'Model'
    return get_count_limited_class(models, ROOT_MODEL_NAME)


class Writer(obj_tables.io.Writer):
    """ Write model to file(s) """

    MODELS = (
        core.Model, core.Taxon, core.Environment,
        core.Submodel, core.Compartment, core.SpeciesType, core.Species,
        core.DistributionInitConcentration, core.Observable, core.Function,
        core.Reaction, core.RateLaw,
        core.DfbaObjective, core.DfbaObjReaction, core.DfbaObjSpecies,
        core.Parameter, core.StopCondition,
        core.Observation, core.ObservationSet, core.Conclusion, core.Reference, core.Author, core.Change,
    )

    def run(self, path, model, models=None, get_related=True, include_all_attributes=False, validate=None,
            title=None, description=None, keywords=None, version=None, language=None, creator=None,
            write_schema=False, write_toc=True, extra_entries=0, data_repo_metadata=False, schema_package=None,
            protected=True):
        """ Write a list of model classes to an Excel file, with one worksheet for each model, or to
            a set of .csv or .tsv files, with one file for each model.

        Args:
            path (:obj:`str`): path to write file(s)
            model (:obj:`type`): a `wc_lang` Model that's a subclass of `obj_tables.Model`
            models (:obj:`list` of :obj:`Model`, optional): models in the order that they should
                appear as worksheets; all models which are not in `models` will
                follow in alphabetical order
            get_related (:obj:`bool`, optional): if :obj:`True`, write object and all related objects
            include_all_attributes (:obj:`bool`, optional): if :obj:`True`, export all attributes including those
                not explictly included in `Model.Meta.attribute_order`
            validate (:obj:`bool`, optional): if :obj:`True`, validate the data
            title (:obj:`str`, optional): title
            description (:obj:`str`, optional): description
            keywords (:obj:`str`, optional): keywords
            version (:obj:`str`, optional): version
            language (:obj:`str`, optional): language
            creator (:obj:`str`, optional): creator
            write_schema (:obj:`bool`, optional): if :obj:`True`, include additional worksheet with schema
            write_toc (:obj:`bool`, optional): if :obj:`True`, include additional worksheet with table of contents
            extra_entries (:obj:`int`, optional): additional entries to display
            data_repo_metadata (:obj:`bool`, optional): if :obj:`True`, try to write metadata information
                about the file's Git repo; the repo must be current with origin, except for the file
            schema_package (:obj:`str`, optional): the package which defines the `obj_tables` schema
                used by the file; if not :obj:`None`, try to write metadata information about the
                the schema's Git repository: the repo must be current with origin
            protected (:obj:`bool`, optional): if :obj:`True`, protect the worksheet
        """
        if models is None:
            models = self.MODELS
        root_model = get_root_model(models)

        if issubclass(self.get_writer(path), obj_tables.io.WorkbookWriter):
            self.validate_implicit_relationships(root_model)
            self.validate_implicit_relationships_are_set(model, root_model)

        config = wc_lang.config.core.get_config()['wc_lang']['io']
        if validate is None:
            validate = config['validate']

        # default metadata for exported file
        if title is None:
            title = model.id
        if description is None:
            description = model.name
        if version is None:
            version = model.version
        if language is None:
            language = 'wc_lang'
        if creator is None:
            creator = '{}.{}'.format(self.__class__.__module__, self.__class__.__name__)

        super(Writer, self).run(path, model, schema_name='wc_lang', models=models, get_related=get_related,
                                include_all_attributes=include_all_attributes, validate=validate,
                                title=title, description=description, version=version, language=language,
                                creator=creator,
                                write_schema=write_schema, write_toc=write_toc,
                                extra_entries=extra_entries,
                                data_repo_metadata=data_repo_metadata, schema_package=schema_package,
                                protected=protected)

    @classmethod
    def validate_implicit_relationships(cls, root_model):
        """ Check that relationships to :obj:`root_model` do not need to be explicitly exported because
            they can be inferred by :obj:`Reader.run`. This is necessary to enable the relationships
            to :obj:`root_model` to not be exported in workbooks, and instead added by :obj:`Reader.run`.

        Raises:
            :obj:`Exception`: if there are relationships from :obj:`core.Model` or one-to-many or
                many-to-many relationships to :obj:`root_model`
        """
        for attr in root_model.Meta.attributes.values():
            if isinstance(attr, obj_tables.RelatedAttribute) and \
                    attr.related_class.__name__ != 'Identifier':
                raise Exception('Relationships from `Model` not supported')

        for attr in root_model.Meta.related_attributes.values():
            if not isinstance(attr, (obj_tables.OneToOneAttribute, obj_tables.ManyToOneAttribute)):
                raise Exception('Only one-to-one and many-to-one relationships are supported to `{}`'.format(
                    root_model.__name__))

    def validate_implicit_relationships_are_set(self, model, root_model):
        """ Check that there is only one instance of :obj:`root_model` and that each relationship to
        :obj:`root_model` is set. This is necessary to enable the relationships to :obj:`root_model`
        to not be exported in workbooks, and instead added by :obj:`Reader.run`.

        Args:
            model (:obj:`obj_tables.Model`): the root model instance
            root_model (:obj:`type`): the type of `model`, a subclass of :obj:`obj_tables.Model`

        Raises:
            :obj:`ValueError`: if there are multiple instances of `root_model` in the object graph
        """
        for obj in model.get_related():
            for attr in obj.Meta.attributes.values():
                if isinstance(attr, obj_tables.RelatedAttribute) and \
                        attr.related_class == root_model:
                    if getattr(obj, attr.name) != model:
                        raise ValueError('{}.{} must be set to the instance of `Model`'.format(obj.__class__.__name__, attr.name))


class Reader(obj_tables.io.Reader):
    """ Read model from file(s) """

    MODELS = Writer.MODELS

    def run(self, path, models=None,
            ignore_missing_models=None, ignore_extra_models=None, ignore_sheet_order=None,
            include_all_attributes=False, ignore_missing_attributes=None, ignore_extra_attributes=None,
            ignore_attribute_order=None, validate=None):
        """ Read a list of model objects from file(s) and, optionally, validate them

        Args:
            path (:obj:`str`): path to file(s)
            models (:obj:`types.TypeType` or :obj:`list` of :obj:`types.TypeType`, optional): type
                of object to read or list of types of objects to read
            ignore_missing_models (:obj:`bool`, optional): if :obj:`False`, report an error if a worksheet/
                file is missing for one or more models
            ignore_extra_models (:obj:`bool`, optional): if :obj:`True` and all `models` are found, ignore
                other worksheets or files
            ignore_sheet_order (:obj:`bool`, optional): if :obj:`True`, do not require the sheets to be provided
                in the canonical order
            include_all_attributes (:obj:`bool`, optional): if :obj:`True`, export all attributes including those
                not explictly included in `Model.Meta.attribute_order`
            ignore_missing_attributes (:obj:`bool`, optional): if :obj:`False`, report an error if a
                worksheet/file doesn't contain all of attributes in a model in `models`
            ignore_extra_attributes (:obj:`bool`, optional): if :obj:`True`, do not report errors if
                attributes in the data are not in the model
            ignore_attribute_order (:obj:`bool`): if :obj:`True`, do not require the attributes to be provided
                in the canonical order
            validate (:obj:`bool`, optional): if :obj:`True`, validate the data

        Returns:
            :obj:`dict`: model objects grouped by `obj_tables.Model` class

        Raises:
            :obj:`ValueError`: if the file defines zero or multiple models or the model defined in the file(s) is
                invalid
        """
        if models is None:
            models = self.MODELS

        root_model = get_root_model(models)
        if issubclass(self.get_reader(path), obj_tables.io.WorkbookReader):
            Writer.validate_implicit_relationships(root_model)

        config = wc_lang.config.core.get_config()['wc_lang']['io']
        if ignore_missing_models is None:
            ignore_missing_models = not config['strict']
        if ignore_extra_models is None:
            ignore_extra_models = not config['strict']
        if ignore_sheet_order is None:
            ignore_sheet_order = not config['strict']
        if ignore_missing_attributes is None:
            ignore_missing_attributes = not config['strict']
        if ignore_extra_attributes is None:
            ignore_extra_attributes = not config['strict']
        if ignore_attribute_order is None:
            ignore_attribute_order = not config['strict']

        objects = super(Reader, self).run(path, schema_name='wc_lang', models=models,
                                          ignore_missing_models=ignore_missing_models,
                                          ignore_extra_models=ignore_extra_models,
                                          ignore_sheet_order=ignore_sheet_order,
                                          include_all_attributes=include_all_attributes,
                                          ignore_missing_attributes=ignore_missing_attributes,
                                          ignore_extra_attributes=ignore_extra_attributes,
                                          ignore_attribute_order=ignore_attribute_order,
                                          group_objects_by_model=True,
                                          validate=False)

        # check that file only has 1 wc_lang Model instance
        root_model = get_root_model(objects)
        if len(objects[root_model]) != 1:
            raise ValueError('"{}" should define one model'.format(path))
        model = objects[root_model][0]

        # add implicit relationships to `Model`
        if issubclass(self.get_reader(path), obj_tables.io.WorkbookReader):
            for cls, cls_objects in objects.items():
                for attr in cls.Meta.attributes.values():
                    if isinstance(attr, obj_tables.RelatedAttribute) and \
                            attr.related_class == root_model:
                        for cls_obj in cls_objects:
                            setattr(cls_obj, attr.name, model)

        # validate
        config = wc_lang.config.core.get_config()['wc_lang']['io']
        if (validate is not None and validate) or (validate is None and config['validate']):
            objs = []
            for cls_objs in objects.values():
                objs.extend(cls_objs)

            errors = obj_tables.Validator().validate(objs)
            if errors:
                raise ValueError(
                    indent_forest(['The model cannot be loaded because it fails to validate:', [errors]]))

        # return model
        return objects


def convert(source, destination, protected=True):
    """ Convert among Excel (.xlsx), comma separated (.csv), and tab separated (.tsv) file formats

    Read a model from the `source` files(s) and write it to the `destination` files(s). A path to a
    delimiter separated set of models must be represented by a Unix glob pattern (with a \\*) that
    matches all delimiter separated files.

    Args:
        source (:obj:`str`): path to source file(s)
        destination (:obj:`str`): path to save converted file
        protected (:obj:`bool`, optional): if :obj:`True`, protect the worksheet
    """
    model = Reader().run(source)[core.Model][0]
    Writer().run(destination, model, data_repo_metadata=False, protected=protected)


def create_template(path, write_schema=False, write_toc=True,
                    extra_entries=10, data_repo_metadata=True,
                    protected=True):
    """ Create file with model template, including row and column headings

    Args:
        path (:obj:`str`): path to file(s)
        write_schema (:obj:`bool`, optional): if :obj:`True`, include additional worksheet with schema
        write_toc (:obj:`bool`, optional): if :obj:`True`, include additional worksheet with table of contents
        extra_entries (:obj:`int`, optional): additional entries to display
        data_repo_metadata (:obj:`bool`, optional): if :obj:`True`, try to write metadata information
            about the file's Git repo
        protected (:obj:`bool`, optional): if :obj:`True`, protect the worksheet
    """
    model = core.Model(id='template', name='Template', version=wc_lang.__version__)
    Writer().run(path, model,
                 write_schema=write_schema, write_toc=write_toc,
                 extra_entries=extra_entries,
                 data_repo_metadata=data_repo_metadata,
                 protected=protected)

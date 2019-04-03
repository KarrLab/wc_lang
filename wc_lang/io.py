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
from wc_utils.util.string import indent_forest
import obj_model.utils
import obj_model
import os
import wc_lang
import wc_lang.config.core


class Writer(obj_model.io.Writer):
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
            extra_entries=0, set_repo_metadata_from_path=True):
        """ Write a list of model classes to an Excel file, with one worksheet for each model, or to
            a set of .csv or .tsv files, with one file for each model.

        Args:
            path (:obj:`str`): path to write file(s)
            model (:obj:`core.Model`): model
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
            extra_entries (:obj:`int`, optional): additional entries to display
            set_repo_metadata_from_path (:obj:`bool`, optional): if :obj:`True`, set the Git repository metadata (URL,
                branch, revision) for the model from the parent directory of :obj:`core_path`
        """
        if issubclass(self.get_writer(path), obj_model.io.WorkbookWriter):
            self.validate_implicit_relationships()
            self.validate_implicit_relationships_are_set(model)

        if models is None:
            models = self.MODELS

        config = wc_lang.config.core.get_config()['wc_lang']['io']
        if validate is None:
            validate = config['validate']

        # default meta data for exported model
        if set_repo_metadata_from_path:
            obj_model.utils.set_git_repo_metadata_from_path(model, path)

        # default meta data for exported file
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

        super(Writer, self).run(path, model, models=models, get_related=get_related,
                                include_all_attributes=include_all_attributes, validate=validate,
                                title=title, description=description, version=version, language=language, creator=creator,
                                extra_entries=extra_entries)

    @classmethod
    def validate_implicit_relationships(cls):
        """ Check that relationships to :obj:`core.Model` do not need to be explicitly exported because they can be inferred
        by :obj:`Reader.run`. This is necessary to enable the relationships to :obj:`core.Model` to not be exported in workbooks, and
        instead added by :obj:`Reader.run`.

        Raises:
            :obj:`Exception`: if there are relationships from :obj:`core.Model` or one-to-many or many-to-many relationships
                to :obj:`core.Model`
        """
        for attr in core.Model.Meta.attributes.values():
            if isinstance(attr, obj_model.RelatedAttribute) and \
                    attr.related_class != core.Identifier:
                raise Exception('Relationships from `Model` not supported')

        for attr in core.Model.Meta.related_attributes.values():
            if not isinstance(attr, (obj_model.OneToOneAttribute, obj_model.ManyToOneAttribute)):
                raise Exception('Only one-to-one and many-to-one relationships are supported to `Model`')

    def validate_implicit_relationships_are_set(self, model):
        """ Check that there is only one instance of :obj:`core.Model` and that each relationship to :obj:`core.Model`
        is set. This is necessary to enable the relationships to :obj:`core.Model` to not be exported in workbooks, and
        instead added by :obj:`Reader.run`.

        Args:
            model (:obj:`core.Model`): model

        Raises:
            :obj:`ValueError`: if there are multiple instances of :obj:`core.Model` in the object graph
        """
        for obj in model.get_related():
            for attr in obj.Meta.attributes.values():
                if isinstance(attr, obj_model.RelatedAttribute) and \
                        attr.related_class == core.Model:
                    if getattr(obj, attr.name) != model:
                        raise ValueError('{}.{} must be set to the instance of `Model`'.format(obj.__class__.__name__, attr.name))


class Reader(obj_model.io.Reader):
    """ Read model from file(s) """

    MODELS = Writer.MODELS

    def run(self, path, models=None,
            ignore_missing_sheets=None, ignore_extra_sheets=None, ignore_sheet_order=None,
            include_all_attributes=False, ignore_missing_attributes=None, ignore_extra_attributes=None, ignore_attribute_order=None,
            group_objects_by_model=True, validate=None):
        """ Read a list of model objects from file(s) and, optionally, validate them

        Args:
            path (:obj:`str`): path to file(s)
            models (:obj:`types.TypeType` or :obj:`list` of :obj:`types.TypeType`, optional): type
                of object to read or list of types of objects to read
            ignore_missing_sheets (:obj:`bool`, optional): if :obj:`False`, report an error if a worksheet/
                file is missing for one or more models
            ignore_extra_sheets (:obj:`bool`, optional): if :obj:`True` and all `models` are found, ignore
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
            group_objects_by_model (:obj:`bool`, optional): if :obj:`True`, group decoded objects by their
                types
            validate (:obj:`bool`, optional): if :obj:`True`, validate the data

        Returns:
            :obj:`dict`: model objects grouped by `obj_model.Model` class

        Raises:
            :obj:`ValueError`: if the file defines zero or multiple models or the model defined in the file(s) is
                invalid
        """
        if issubclass(self.get_reader(path), obj_model.io.WorkbookReader):
            Writer.validate_implicit_relationships()

        if models is None:
            models = self.MODELS

        config = wc_lang.config.core.get_config()['wc_lang']['io']
        if ignore_missing_sheets is None:
            ignore_missing_sheets = not config['strict']
        if ignore_extra_sheets is None:
            ignore_extra_sheets = not config['strict']
        if ignore_sheet_order is None:
            ignore_sheet_order = not config['strict']
        if ignore_missing_attributes is None:
            ignore_missing_attributes = not config['strict']
        if ignore_extra_attributes is None:
            ignore_extra_attributes = not config['strict']
        if ignore_attribute_order is None:
            ignore_attribute_order = not config['strict']

        objects = super(Reader, self).run(path, models=models,
                                          ignore_missing_sheets=ignore_missing_sheets,
                                          ignore_extra_sheets=ignore_extra_sheets,
                                          ignore_sheet_order=ignore_sheet_order,
                                          include_all_attributes=include_all_attributes,
                                          ignore_missing_attributes=ignore_missing_attributes,
                                          ignore_extra_attributes=ignore_extra_attributes,
                                          ignore_attribute_order=ignore_attribute_order,
                                          group_objects_by_model=group_objects_by_model,
                                          validate=False)

        # check that file only has 1 model
        if len(objects[core.Model]) != 1:
            raise ValueError('"{}" should define one model'.format(path))
        model = objects[core.Model][0]

        # add implicit relationships to `Model`
        if issubclass(self.get_reader(path), obj_model.io.WorkbookReader):
            for cls, cls_objects in objects.items():
                for attr in cls.Meta.attributes.values():
                    if isinstance(attr, obj_model.RelatedAttribute) and \
                            attr.related_class == core.Model:
                        for cls_obj in cls_objects:
                            setattr(cls_obj, attr.name, model)

        # validate
        config = wc_lang.config.core.get_config()['wc_lang']['io']
        if (validate is not None and validate) or (validate is None and config['validate']):
            objs = []
            for cls_objs in objects.values():
                objs.extend(cls_objs)

            errors = obj_model.Validator().validate(objs)
            if errors:
                raise ValueError(
                    indent_forest(['The model cannot be loaded because it fails to validate:', [errors]]))

        # return model
        return objects


def convert(source, destination):
    """ Convert among Excel (.xlsx), comma separated (.csv), and tab separated (.tsv) file formats

    Read a model from the `source` files(s) and write it to the `destination` files(s). A path to a
    delimiter separated set of models must be represented by a Unix glob pattern (with a \\*) that
    matches all delimiter separated files.

    Args:
        source (:obj:`str`): path to source file(s)
        destination (:obj:`str`): path to save converted file
    """
    model = Reader().run(source)[core.Model][0]
    Writer().run(destination, model, set_repo_metadata_from_path=False)


def create_template(path, extra_entries=10, set_repo_metadata_from_path=True):
    """ Create file with model template, including row and column headings

    Args:
        path (:obj:`str`): path to file(s)
        extra_entries (:obj:`int`, optional): additional entries to display
        set_repo_metadata_from_path (:obj:`bool`, optional): if :obj:`True`, set the Git repository metadata (URL,
            branch, revision) for the model from the parent directory of :obj:`core_path`
    """
    model = core.Model(id='template', name='Template', version=wc_lang.__version__)
    Writer().run(path, model, extra_entries=extra_entries, set_repo_metadata_from_path=set_repo_metadata_from_path)

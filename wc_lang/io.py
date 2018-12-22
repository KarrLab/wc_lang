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
from wc_utils.util.string import indent_forest
import obj_model
import os
import wc_lang
import wc_lang.config.core


class Writer(object):
    """ Write model to file(s) """

    model_order = [
        core.Model, core.Taxon,
        core.Submodel, core.Compartment, core.SpeciesType, core.Species,
        core.DistributionInitConcentration, core.Observable, core.Function,
        core.Reaction, core.RateLaw,
        core.DfbaObjective, core.DfbaNetReaction, core.DfbaNetSpecies,
        core.Parameter, core.StopCondition,
        core.Evidence, core.Reference,
    ]

    def run(self, model, path, set_repo_metadata_from_path=True):
        """ Write model to file(s)

        Args:
            model (:obj:`core.Model`): model
            path (:obj:`str`): path to file(s)
            set_repo_metadata_from_path (:obj:`bool`, optional): if :obj:`True`, set the Git repository metadata (URL,
                branch, revision) for the model from the parent directory of :obj:`core_path`
        """
        self.validate_implicit_relationships()

        # check that there is only 1 :obj:`Model`and that each relationship to :obj:`Model` is set. This is necessary to
        # enable the relationships to :obj:`Model` to be implicit in the Excel output and added by :obj:`Reader.run`
        for obj in model.get_related():
            for attr in obj.Meta.attributes.values():
                if isinstance(attr, obj_model.RelatedAttribute) and \
                        attr.related_class == core.Model:
                    if getattr(obj, attr.name) != model:
                        raise ValueError('{}.{} must be set to the instance of `Model`'.format(obj.__class__.__name__, attr.name))

        # set Git repository metadata from the parent directories of :obj:`core_path`
        if set_repo_metadata_from_path:
            util.set_git_repo_metadata_from_path(model, path)

        # write objects
        _, ext = os.path.splitext(path)
        writer = obj_model.io.get_writer(ext)()

        kwargs = {}
        if isinstance(writer, obj_model.io.WorkbookWriter):
            kwargs['include_all_attributes'] = False

        writer.run(path, [model], models=self.model_order,
                   language='wc_lang',
                   creator='{}.{}'.format(self.__class__.__module__, self.__class__.__name__),
                   title=model.id,
                   description=model.name,
                   version=model.version,
                   **kwargs)

    @classmethod
    def validate_implicit_relationships(cls):
        """ Check that relationships to :obj:`core.Model` do not need to be explicitly written to
        workbooks because they can be inferred by :obj:`Reader.run`
        """
        for attr in core.Model.Meta.attributes.values():
            if isinstance(attr, obj_model.RelatedAttribute) and \
                    attr.related_class != core.DatabaseReference:
                raise Exception('Relationships from `Model` not supported')

        for attr in core.Model.Meta.related_attributes.values():
            if not isinstance(attr, (obj_model.OneToOneAttribute, obj_model.ManyToOneAttribute)):
                raise Exception('Only one-to-one and many-to-one relationships are supported to `Model`')


class Reader(object):
    """ Read model from file(s) """

    def run(self, path):
        """ Read model from file(s)

        Args:
            path (:obj:`str`): path to file(s)

        Returns:
            :obj:`core.Model`: model

        Raises:
            :obj:`ValueError`: if :obj:`path` defines multiple models
        """
        config = wc_lang.config.core.get_config()

        Writer.validate_implicit_relationships()

        # read objects from file
        _, ext = os.path.splitext(path)
        reader = obj_model.io.get_reader(ext)()

        kwargs = {}
        if isinstance(reader, obj_model.io.WorkbookReader):
            kwargs['include_all_attributes'] = False
            if not config['wc_lang']['io']['strict']:
                kwargs['ignore_missing_sheets'] = True
                kwargs['ignore_extra_sheets'] = True
                kwargs['ignore_sheet_order'] = True
                kwargs['ignore_missing_attributes'] = True
                kwargs['ignore_extra_attributes'] = True
                kwargs['ignore_attribute_order'] = True
        objects = reader.run(path, models=Writer.model_order, validate=False, **kwargs)

        # check that file only has 0 or 1 models
        if not objects[core.Model]:
            for cls, cls_objects in objects.items():
                if cls_objects:
                    raise ValueError('"{}" cannot contain instances of `{}` without an instance of `Model`'.format(
                        path, cls.__name__))
            return None

        elif len(objects[core.Model]) > 1:
            raise ValueError('"{}" should define one model'.format(path))

        else:
            model = objects[core.Model].pop()

        # add implicit relationships to `Model`
        for cls, cls_objects in objects.items():
            for attr in cls.Meta.attributes.values():
                if isinstance(attr, obj_model.RelatedAttribute) and \
                        attr.related_class == core.Model:
                    for cls_obj in cls_objects:
                        setattr(cls_obj, attr.name, model)

        # validate
        objs = []
        for cls_objs in objects.values():
            objs.extend(cls_objs)

        errors = obj_model.Validator().validate(objs)
        if errors:
            raise ValueError(
                indent_forest(['The model cannot be loaded because it fails to validate:', [errors]]))

        # return model
        return model


def convert(source, destination):
    """ Convert among Excel (.xlsx), comma separated (.csv), and tab separated (.tsv) file formats

    Read a model from the `source` files(s) and write it to the `destination` files(s). A path to a
    delimiter separated set of models must be represented by a Unix glob pattern (with a \\*) that
    matches all delimiter separated files.

    Args:
        source (:obj:`str`): path to source file(s)
        destination (:obj:`str`): path to save converted file
    """
    model = Reader().run(source)
    Writer().run(model, destination, set_repo_metadata_from_path=False)


def create_template(path, set_repo_metadata_from_path=True):
    """ Create file with model template, including row and column headings

    Args:
        path (:obj:`str`): path to file(s)
        set_repo_metadata_from_path (:obj:`bool`, optional): if :obj:`True`, set the Git repository metadata (URL,
            branch, revision) for the model from the parent directory of :obj:`core_path`
    """
    model = core.Model(id='template', name='Template', version=wc_lang.__version__)
    Writer().run(model, path, set_repo_metadata_from_path=set_repo_metadata_from_path)

""" Command line programs for manipulating model definitions

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2016-12-07
:Copyright: 2016, Karr Lab
:License: MIT
"""

from obj_tables.migrate import schema_repo_migration_controllers
from wc_lang import transform
from wc_lang.core import Model
from wc_lang.io import Writer, Reader, convert, create_template
from wc_utils.workbook.io import read as read_workbook
import cement
import obj_tables
import os
import sys
import wc_lang
import wc_lang.sbml.io
import wc_utils.workbook


class BaseController(cement.Controller):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = "Command line utilities for managing whole-cell model definitions"
        help = "Command line utilities for managing whole-cell model definitions"
        arguments = [
            (['-v', '--version'], dict(action='version', version=wc_lang.__version__)),
        ]

    @cement.ex(hide=True)
    def _default(self):
        self._parser.print_help()


class CutSubmodelsController(cement.Controller):
    """ Cut submodels into separate models """
    class Meta:
        label = 'cut-submodels'
        description = 'Cut submodels into separate models'
        help = 'Cut submodels into separate models'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['in_file'], dict(type=str, help='Path to model to cut')),
            (['out_dir'], dict(type=str, help='Directory to save cut submodels')),
            (['--unprotected'], dict(action='store_true', default=False,
                                     help='If set, do not protect the outputted workbook')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs

        # read model
        model = Reader().run(args.in_file)[Model][0]

        # split submodels into separate models
        core, submodels = model.submodels.gen_models()

        # create output directory, if it doesn't exist
        if not os.path.isdir(args.out_dir):
            os.makedirs(args.out_dir)

        # save separated submodels to file
        Writer().run(os.path.join(args.out_dir, 'core.xlsx'), core, data_repo_metadata=False,
                     protected=(not args.unprotected))
        for submodel in submodels:
            Writer().run(os.path.join(args.out_dir, '{}.xlsx'.format(
                submodel.submodels[0].id)), submodel, data_repo_metadata=False,
                protected=(not args.unprotected))


class MergeModelsController(cement.Controller):
    """ Merge models """
    class Meta:
        label = 'merge-models'
        description = 'Merge multiple models'
        help = 'Merge multiple models'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['-p', '--primary'], dict(dest='primary_path', type=str, help='Path to base for merged model')),
            (['-s', '--secondary'], dict(dest='secondary_paths', type=str, nargs='*', help='Path to models to merge into primary model')),
            (['-o', '--out'], dict(dest='out_path', type=str, help='Path to save merged model')),
            (['--unprotected'], dict(action='store_true', default=False,
                                     help='If set, do not protect the outputted workbook')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs

        primary_model = Reader().run(args.primary_path)[Model][0]

        for secondary_path in args.secondary_paths:
            secondary_model = Reader().run(secondary_path)[Model][0]
            primary_model.merge(secondary_model)

        Writer().run(args.out_path, primary_model, data_repo_metadata=False,
                     protected=(not args.unprotected))


class ValidateController(cement.Controller):
    """ Validate model and display errors """

    class Meta:
        label = 'validate'
        description = 'Validate model and display errors'
        help = 'Validate model and display errors'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path'], dict(type=str, help='Path to model')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        try:
            Reader().run(args.path)  # reader already does validation
            print('Model is valid')
        except ValueError as exception:
            raise SystemExit('Model is invalid: ' + str(exception))


class DifferenceController(cement.Controller):
    """ Display difference between two model definitions """

    class Meta:
        label = 'difference'
        description = 'Get difference between two model definitions'
        help = 'Get difference between two model definitions'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path_1'], dict(type=str, help='Path to first model')),
            (['path_2'], dict(type=str, help='Path to second model')),
            (['--compare-files'], dict(dest='compare_files', default=False, action='store_true',
                                       help='If true, compare models; otherwise compare files directly')),
            (['--compare-metadata-in-files'], dict(dest='compare_metadata_in_files', default=False, action='store_true',
                                                   help='If true, compare metadata (tables of contents, header rows)')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs

        if args.compare_files:
            model1 = read_workbook(args.path_1)
            model2 = read_workbook(args.path_2)
            if not args.compare_metadata_in_files:
                self.remove_metadata(model1)
                self.remove_metadata(model2)

            diff = model1.difference(model2)

        else:
            model1 = Reader().run(args.path_1)[Model][0]
            model2 = Reader().run(args.path_2)[Model][0]
            diff = model1.difference(model2)

        if diff:
            print(diff)
        else:
            print('Models are identical')

    @staticmethod
    def remove_metadata(model):
        """ Remove metadata from model

        Args:
            model (:obj:`wc_utils.workbook.Workbook`): model
        """
        model.pop('!' + obj_tables.TOC_SHEET_NAME, None)
        for sheet in model.values():
            for row in list(sheet):
                if row and isinstance(row[0], str) and row[0].startswith('!!'):
                    sheet.remove(row)
                else:
                    break


transform_list = ''
for trans in transform.get_transforms().values():
    transform_list += '\n  {}: {}'.format(trans.Meta.id, trans.Meta.label)


class TransformController(cement.Controller):
    """ Apply one, or more, transforms to a model and save the result """

    class Meta:
        label = 'transform'
        description = 'Apply one, or more, transforms to a model and save the result'
        help = 'Apply one, or more, transforms to a model and save the result'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['source'], dict(type=str, help='Path to model')),
            (['dest'], dict(type=str, help='Path to save transformed model')),
            (['--transform'], dict(dest='transforms', action='append',
                                   help='Model transform:' + transform_list)),
            (['--unprotected'], dict(action='store_true', default=False,
                                     help='If set, do not protect the outputted workbook')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs

        if not args.transforms:
            raise SystemExit('Please select at least one transform')

        # read model
        model = Reader().run(args.source)[Model][0]

        # apply transforms
        transforms = transform.get_transforms()
        for id in args.transforms:
            cls = transforms[id]
            instance = cls()
            instance.run(model)

        # write model
        Writer().run(args.dest, model, data_repo_metadata=False,
                     protected=(not args.unprotected))


class NormalizeController(cement.Controller):
    """ Normalize model """

    class Meta:
        label = 'normalize'
        description = 'Normalize model'
        help = 'Normalize model'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['source'], dict(type=str, help='Path to model')),
            (['--dest'], dict(default='', type=str, help='Path to save normalized model')),
            (['--unprotected'], dict(action='store_true', default=False,
                                     help='If set, do not protect the outputted workbook')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        model = Reader().run(args.source)[Model][0]
        if args.dest:
            Writer().run(args.dest, model, data_repo_metadata=False,
                         protected=(not args.unprotected))
        else:
            Writer().run(args.source, model, data_repo_metadata=False,
                         protected=(not args.unprotected))


class ConvertController(cement.Controller):
    """ Convert model among Excel (.xlsx), comma separated (.csv), JavaScript Object Notation (.json),
    tab separated (.tsv), and Yet Another Markup Language (.yaml, .yml) formats """

    class Meta:
        label = 'convert'
        description = 'Convert model among .csv, .json, .tsv, .xlsx, .yaml, and .yml formats'
        help = 'Convert model among .csv, .json, .tsv, .xlsx, .yaml, and .yml formats'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['source'], dict(type=str, help='Path to model')),
            (['dest'], dict(type=str, help='Path to save model in converted format')),
            (['--unprotected'], dict(action='store_true', default=False,
                                     help='If set, do not protect the outputted workbook')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        convert(args.source, args.dest, protected=(not args.unprotected))


class CreateTemplateController(cement.Controller):
    """ Create file with model template (i.e. create file with row and column labels) """

    class Meta:
        label = 'create-template'
        description = 'Create file with model template: blank file(s) with row and column labels'
        help = 'Create file with model template: blank file(s) with row and column labels'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path'], dict(type=str, help='Path to save model template')),
            (['--ignore-repo-metadata'], dict(dest='data_repo_metadata', default=True, action='store_false',
                                              help=('If set, do not set the Git repository metadata for the knowledge base from '
                                                    'the parent directory of `path`'))),
            (['--unprotected'], dict(action='store_true', default=False,
                                     help='If set, do not protect the outputted workbook')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        create_template(args.path, data_repo_metadata=args.data_repo_metadata, protected=(not args.unprotected))


class UpdateVersionMetadataController(cement.Controller):
    """ Update the version metadata (repository URL, branch, revision; wc_lang version) of a model """

    class Meta:
        label = 'update-version-metadata'
        description = 'Update the version metadata (repository URL, branch, revision; wc_lang version) of a model'
        help = 'Update the version metadata (repository URL, branch, revision; wc_lang version) of a model'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path'], dict(type=str, help='Path to model')),
            (['--ignore-repo-metadata'], dict(dest='data_repo_metadata', default=True, action='store_false',
                                              help=('If set, do not set the Git repository metadata for the knowledge base from '
                                                    'the parent directory of `path-core`'))),
            (['--unprotected'], dict(action='store_true', default=False,
                                     help='If set, do not protect the outputted workbook')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        model = Reader().run(args.path)[Model][0]
        model.wc_lang_version = wc_lang.__version__
        Writer().run(args.path, model, data_repo_metadata=args.data_repo_metadata,
                     protected=(not args.unprotected))


class ExportController(cement.Controller):
    """ Export a model to SBML """

    class Meta:
        label = 'export'
        description = 'Export a model to SBML'
        help = 'Export a model to SBML'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['in_path'], dict(type=str, help='Path to model to export')),
            (['out_dir'], dict(type=str, help='Directory to save exported model')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        model = Reader().run(args.in_path)[Model][0]
        wc_lang.sbml.io.SbmlWriter().run(model, args.out_dir)


class ImportController(cement.Controller):
    """ Import a model from SBML """

    class Meta:
        label = 'import'
        description = 'Import a model from SBML'
        help = 'Import a model from SBML'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['in_dir'], dict(type=str, help='Directory with model to import')),
            (['out_path'], dict(type=str, help='Path to save model')),
            (['--unprotected'], dict(action='store_true', default=False,
                                     help='If set, do not protect the outputted workbook')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        model = wc_lang.sbml.io.SbmlReader().run(args.in_dir)
        Writer().run(args.out_path, model, data_repo_metadata=False,
                     protected=(not args.unprotected))


class App(cement.App):
    """ Command line application """
    class Meta:
        label = 'wc-lang'
        base_controller = 'base'
        handlers = [
            BaseController,
            CutSubmodelsController,
            MergeModelsController,
            ValidateController,
            DifferenceController,
            TransformController,
            NormalizeController,
            ConvertController,
            CreateTemplateController,
            UpdateVersionMetadataController,
            ExportController,
            ImportController,
        ] + schema_repo_migration_controllers


def main():
    with App() as app:
        app.run()

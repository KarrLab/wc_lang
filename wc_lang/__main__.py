""" Command line programs for manipulating model definitions

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-12-07
:Copyright: 2016, Karr Lab
:License: MIT
"""

from wc_lang import transform
from wc_lang.io import Writer, Reader, convert, create_template
from wc_utils.workbook.io import read as read_workbook
import cement
import sys
import wc_lang


class BaseController(cement.Controller):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = "Command line utilities for managing whole-cell model definitions"
        arguments = [
            (['-v', '--version'], dict(action='version', version=wc_lang.__version__)),
        ]

    @cement.ex(hide=True)
    def _default(self):
        self._parser.print_help()


class ValidateController(cement.Controller):
    """ Validate model definition and display errors """

    class Meta:
        label = 'validate'
        description = 'Validate model definition and display errors'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path'], dict(type=str, help='Path to model definition')),
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
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path_1'], dict(type=str, help='Path to first model definition')),
            (['path_2'], dict(type=str, help='Path to second model definition')),
            (['--compare-files'], dict(dest='compare_files', default=False, action='store_true',
                                       help='If true, compare models; otherwise compare files directly')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs

        if args.compare_files:
            model1 = read_workbook(args.path_1)
            model2 = read_workbook(args.path_2)
            diff = model1.difference(model2)

        else:
            model1 = Reader().run(args.path_1)
            model2 = Reader().run(args.path_2)
            diff = model1.difference(model2)

        if diff:
            print(diff)
        else:
            print('Models are identical')


transform_list = ''
for trans in transform.get_transforms().values():
    transform_list += '\n  {}: {}'.format(trans.Meta.id, trans.Meta.label)


class TransformController(cement.Controller):
    """ Apply one, or more, transforms to a model and save the result """

    class Meta:
        label = 'transform'
        description = 'Apply one, or more, transforms to a model and save the result'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['source'], dict(type=str, help='Path to model definition')),
            (['dest'], dict(type=str, help='Path to save transformed model definition')),
            (['--transform'], dict(dest='transforms', action='append',
                                   help='Model transform:' + transform_list)),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs

        if not args.transforms:
            raise SystemExit('Please select at least one transform')

        # read model
        model = Reader().run(args.source)

        # apply transforms
        transforms = transform.get_transforms()
        for id in args.transforms:
            cls = transforms[id]
            instance = cls()
            instance.run(model)

        # write model
        Writer().run(model, args.dest, set_repo_metadata_from_path=False)


class NormalizeController(cement.Controller):
    """ Normalize model definition """

    class Meta:
        label = 'normalize'
        description = 'Normalize model definition'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['source'], dict(type=str, help='Path to model definition')),
            (['--dest'], dict(default='', type=str, help='Path to save normalized model definition')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        model = Reader().run(args.source)
        if args.dest:
            Writer().run(model, args.dest, set_repo_metadata_from_path=False)
        else:
            Writer().run(model, args.source, set_repo_metadata_from_path=False)


class ConvertController(cement.Controller):
    """ Convert model definition among Excel (.xlsx), comma separated (.csv), JavaScript Object Notation (.json),
    tab separated (.tsv), and Yet Another Markup Language (.yaml, .yml) formats """

    class Meta:
        label = 'convert'
        description = 'Convert model definition among .csv, .json, .tsv, .xlsx, .yaml, and .yml formats'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['source'], dict(type=str, help='Path to model definition')),
            (['dest'], dict(type=str, help='Path to save model in converted format')),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        convert(args.source, args.dest)


class CreateTemplateController(cement.Controller):
    """ Create file with model template (i.e. create file with row and column labels) """

    class Meta:
        label = 'create-template'
        description = 'Create file with model template: blank file(s) with row and column labels'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path'], dict(type=str, help='Path to save model template')),
            (['--ignore-repo-metadata'], dict(dest='set_repo_metadata_from_path', default=True, action='store_false',
                                              help=('If set, do not set the Git repository metadata for the knowledge base from '
                                                    'the parent directory of `path`'))),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        create_template(args.path, set_repo_metadata_from_path=args.set_repo_metadata_from_path)


class UpdateVersionMetadataController(cement.Controller):
    """ Update the version metadata (repository URL, branch, revision; wc_lang version) of a model """

    class Meta:
        label = 'update-version-metadata'
        description = 'Update the version metadata (repository URL, branch, revision; wc_lang version) of a model'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path'], dict(type=str, help='Path to model')),
            (['--ignore-repo-metadata'], dict(dest='set_repo_metadata_from_path', default=True, action='store_false',
                                              help=('If set, do not set the Git repository metadata for the knowledge base from '
                                                    'the parent directory of `path-core`'))),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        model = Reader().run(args.path)
        model.wc_lang_version = wc_lang.__version__
        Writer().run(model, args.path, set_repo_metadata_from_path=args.set_repo_metadata_from_path)


class App(cement.App):
    """ Command line application """
    class Meta:
        label = 'wc_lang'
        base_controller = 'base'
        handlers = [
            BaseController,
            ValidateController,
            DifferenceController,
            TransformController,
            NormalizeController,
            ConvertController,
            CreateTemplateController,
            UpdateVersionMetadataController,
        ]


def main():
    with App() as app:
        app.run()

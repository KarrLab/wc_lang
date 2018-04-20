""" Command line programs for manipulating model definitions

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-12-07
:Copyright: 2016, Karr Lab
:License: MIT
"""

from cement.core.foundation import CementApp
from cement.core.controller import CementBaseController, expose
from wc_lang import transform
from wc_lang.io import Writer, Reader, convert, create_template
from wc_utils.workbook.io import read as read_workbook
import wc_lang


class BaseController(CementBaseController):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = "Command line utilities for managing whole-cell model definitions"

    @expose(help='Get version')
    def get_version(self):
        """ Get version """
        print(wc_lang.__version__)


class ValidateController(CementBaseController):
    """ Validate model definition and display errors """

    class Meta:
        label = 'validate'
        description = 'Validate model definition and display errors'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path'], dict(type=str, help='Path to model definition')),
            (['--sloppy'], dict(dest='strict', default=True, action='store_false',
                                help='If set, do not validate the format of the model file(s)')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        try:
            Reader().run(args.path, strict=args.strict)  # reader already does validation
            print('Model is valid')
        except ValueError as exception:
            raise ValueError('Model is invalid: ' + str(exception))


class DifferenceController(CementBaseController):
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
            (['--sloppy'], dict(dest='strict', default=True, action='store_false',
                                help='If set, do not validate the format of the model file(s)')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs

        if args.compare_files:
            model1 = read_workbook(args.path_1)
            model2 = read_workbook(args.path_2)
            diff = model1.difference(model2)

        else:
            model1 = Reader().run(args.path_1, strict=args.strict)
            model2 = Reader().run(args.path_2, strict=args.strict)
            diff = model1.difference(model2)

        if diff:
            print(diff)
        else:
            print('Models are identical')


transform_list = ''
for trans in transform.get_transforms().values():
    transform_list += '\n  {}: {}'.format(trans.Meta.id, trans.Meta.label)


class TransformController(CementBaseController):
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
            (['--sloppy'], dict(dest='strict', default=True, action='store_false',
                                help='If set, do not validate the format of the model file(s)')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs

        if not args.transforms:
            raise ValueError('Please select at least one transform')

        # read model
        model = Reader().run(args.source, strict=args.strict)

        # apply transforms
        transforms = transform.get_transforms()
        for id in args.transforms:
            cls = transforms[id]
            instance = cls()
            instance.run(model)

        # write model
        Writer().run(model, args.dest)


class NormalizeController(CementBaseController):
    """ Normalize model definition """

    class Meta:
        label = 'normalize'
        description = 'Normalize model definition'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['source'], dict(type=str, help='Path to model definition')),
            (['--dest'], dict(default='', type=str, help='Path to save normalized model definition')),
            (['--sloppy'], dict(dest='strict', default=True, action='store_false',
                                help='If set, do not validate the format of the model file(s)')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        model = Reader().run(args.source, strict=args.strict)
        if args.dest:
            Writer().run(model, args.dest)
        else:
            Writer().run(model, args.source)


class ConvertController(CementBaseController):
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
            (['--sloppy'], dict(dest='strict', default=True, action='store_false',
                                help='If set, do not validate the format of the model file(s)')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        convert(args.source, args.dest, strict=args.strict)


class CreateTemplateController(CementBaseController):
    """ Create file with model template (i.e. create file with row and column labels) """

    class Meta:
        label = 'create-template'
        description = 'Create file with model template: blank file(s) with row and column labels'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path'], dict(type=str, help='Path to save model template')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        create_template(args.path)


class UpdateWcLangVersionController(CementBaseController):
    """ Update wc_lang_version of model """

    class Meta:
        label = 'update-wc-lang-version'
        description = 'Update wc_lang_version of model'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['path'], dict(type=str, help='Path to model')),
            (['--sloppy'], dict(dest='strict', default=True, action='store_false',
                                help='If set, do not validate the format of the model file(s)')),
        ]

    @expose(hide=True)
    def default(self):
        args = self.app.pargs
        model = Reader().run(args.path, strict=args.strict)
        model.wc_lang_version = wc_lang.__version__
        Writer().run(model, args.path)


class App(CementApp):
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
            UpdateWcLangVersionController,
        ]


def main():
    with App() as app:
        app.run()

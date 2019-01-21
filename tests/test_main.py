""" Tests of command line program

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2016-12-07
:Copyright: 2016, Karr Lab
:License: MIT
"""

from capturer import CaptureOutput
from obj_model import Validator
from os import path
from shutil import rmtree
from tempfile import mkdtemp
from wc_lang import __main__
from wc_lang import Model, Parameter
from wc_lang.io import Writer, Reader
from wc_utils.util.units import unit_registry
import datetime
import mock
import unittest
import wc_lang


class TestCli(unittest.TestCase):

    def setUp(self):
        self.tempdir = mkdtemp()

    def tearDown(self):
        rmtree(self.tempdir)

    def test_get_version(self):
        with CaptureOutput() as capturer:
            with __main__.App(argv=['-v']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(capturer.get_text(), wc_lang.__version__)

        with CaptureOutput() as capturer:
            with __main__.App(argv=['--version']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(capturer.get_text(), wc_lang.__version__)

    def test_merge_models(self):
        in_paths = [
            path.join(self.tempdir, 'in-0.xlsx'),
            path.join(self.tempdir, 'in-1.xlsx'),
        ]
        out_path = path.join(self.tempdir, 'out.xlsx')

        # write models
        model_0 = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1')
        model_0.species_types.create(id='a')
        model_1 = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1')
        model_1.species_types.create(id='b')
        Writer().run(in_paths[0], model_0, set_repo_metadata_from_path=False)
        Writer().run(in_paths[1], model_1, set_repo_metadata_from_path=False)

        # merge models
        with __main__.App(argv=['merge', '-p', in_paths[0], '-s', in_paths[1], '-o', out_path]) as app:
            app.run()

        # read merged model
        merged_model = Reader().run(out_path)[Model][0]

        # verify merged model
        self.assertEqual(len(merged_model.species_types), 2)
        self.assertNotEqual(merged_model.species_types.get_one(id='a'), None)
        self.assertNotEqual(merged_model.species_types.get_one(id='b'), None)

    def test_validate(self):
        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1')
        self.assertEqual(Validator().run(model, get_related=True), None)
        filename = path.join(self.tempdir, 'model.xlsx')
        Writer().run(filename, model, set_repo_metadata_from_path=False)

        with CaptureOutput() as capturer:
            with __main__.App(argv=['validate', filename]) as app:
                app.run()
            self.assertEqual(capturer.get_text(), 'Model is valid')

    def test_validate_exception(self):
        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1')
        model.parameters.append(Parameter(id='param_1', value=1., units=unit_registry.parse_units('dimensionless')))
        model.parameters.append(Parameter(id='param_1', value=1., units=unit_registry.parse_units('dimensionless')))

        self.assertNotEqual(Validator().run(model, get_related=True), None)
        filename = path.join(self.tempdir, 'model.xlsx')
        Writer().run(filename, model, set_repo_metadata_from_path=False)

        with self.assertRaisesRegex(SystemExit, '^Model is invalid: '):
            with __main__.App(argv=['validate', filename]) as app:
                app.run()

    def test_difference(self):
        now = datetime.datetime.now().replace(microsecond=0)

        model1 = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0', created=now, updated=now)
        filename1 = path.join(self.tempdir, 'model1.xlsx')
        Writer().run(filename1, model1, set_repo_metadata_from_path=False)

        model2 = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0', created=now, updated=now)
        filename2 = path.join(self.tempdir, 'model2.xlsx')
        Writer().run(filename2, model2, set_repo_metadata_from_path=False)

        model3 = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1', created=now, updated=now)
        filename3 = path.join(self.tempdir, 'model3.xlsx')
        Writer().run(filename3, model3, set_repo_metadata_from_path=False)

        with CaptureOutput() as capturer:
            with __main__.App(argv=['difference', filename1, filename2]) as app:
                app.run()
            self.assertEqual(capturer.get_text(), 'Models are identical')

        with CaptureOutput() as capturer:
            with __main__.App(argv=['difference', filename1, filename2, '--compare-files']) as app:
                app.run()
            self.assertEqual(capturer.get_text(), 'Models are identical')

        with CaptureOutput() as capturer:
            with __main__.App(argv=['difference', filename1, filename3]) as app:
                app.run()
            diff = ('Objects (Model: "model", Model: "model") have different attribute values:\n  '
                    '`wc_lang_version` are not equal:\n    0.0.0 != 0.0.1')
            self.assertEqual(capturer.get_text(), diff)

        with CaptureOutput() as capturer:
            with __main__.App(argv=['difference', filename1, filename3, '--compare-files']) as app:
                app.run()
            diff = 'Sheet Model:\n  Row 7:\n    Cell B: 0.0.0 != 0.0.1'
            self.assertEqual(capturer.get_text(), diff)

    def test_transform(self):
        source = path.join(self.tempdir, 'source.xlsx')
        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0')
        Writer().run(source, model, set_repo_metadata_from_path=False)

        dest = path.join(self.tempdir, 'dest.xlsx')
        with __main__.App(argv=['transform', source, dest, '--transform', 'MergeAlgorithmicallyLikeSubmodels']) as app:
            app.run()

        self.assertTrue(path.isfile(dest))

    def test_transform_exception(self):
        source = path.join(self.tempdir, 'source.xlsx')
        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0')
        Writer().run(source, model, set_repo_metadata_from_path=False)

        dest = path.join(self.tempdir, 'dest.xlsx')
        with self.assertRaisesRegex(SystemExit, 'Please select at least one transform'):
            with __main__.App(argv=['transform', source, dest]) as app:
                app.run()

    def test_normalize(self):
        filename_xls_1 = path.join(self.tempdir, 'model-1.xlsx')
        filename_xls_2 = path.join(self.tempdir, 'model-2.xlsx')

        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0')
        Writer().run(filename_xls_1, model, set_repo_metadata_from_path=False)

        # with same destination
        with __main__.App(argv=['normalize', filename_xls_1]) as app:
            app.run()

        model2 = Reader().run(filename_xls_1)[Model][0]
        self.assertTrue(model2.is_equal(model))

        # with different destination
        with __main__.App(argv=['normalize', filename_xls_1, '--dest', filename_xls_2]) as app:
            app.run()

        model2 = Reader().run(filename_xls_2)[Model][0]
        self.assertTrue(model2.is_equal(model))

    def test_convert(self):
        filename_xls = path.join(self.tempdir, 'model.xlsx')
        filename_csv = path.join(self.tempdir, 'model-*.csv')

        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0')
        Writer().run(filename_xls, model, set_repo_metadata_from_path=False)

        with __main__.App(argv=['convert', filename_xls, filename_csv]) as app:
            app.run()

        self.assertTrue(path.isfile(path.join(self.tempdir, 'model-Model.csv')))

    def test_create_template(self):
        filename = path.join(self.tempdir, 'template.xlsx')

        with __main__.App(argv=['create-template', filename, '--ignore-repo-metadata']) as app:
            app.run()

        self.assertTrue(path.isfile(filename))

    def test_update_version_metadata(self):
        filename = path.join(self.tempdir, 'model.xlsx')

        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0')
        self.assertNotEqual(model.wc_lang_version, wc_lang.__version__)
        Writer().run(filename, model, set_repo_metadata_from_path=False)

        with __main__.App(argv=['update-version-metadata', filename, '--ignore-repo-metadata']) as app:
            app.run()

        model = Reader().run(filename)[Model][0]
        self.assertEqual(model.wc_lang_version, wc_lang.__version__)

    def test_migrate(self):
        model = Model(id='model', name='test model', version='0.0.1', wc_lang_version='0.0.1')
        filename = path.join(self.tempdir, 'model.xlsx')
        Writer().run(filename, model, set_repo_metadata_from_path=False)

        version = '0.0.2'
        filename_2 = path.join(self.tempdir, 'model-2.xlsx')
        with __main__.App(argv=['migrate', filename, version,
                                '--ignore-repo-metadata', '--out-path', filename_2]) as app:
            app.run()
        model_2 = Reader().run(filename_2)[Model][0]
        self.assertEqual(model_2.wc_lang_version, version)

        version = '0.0.3'
        with __main__.App(argv=['migrate', filename, version,
                                '--ignore-repo-metadata']) as app:
            app.run()
        model_2 = Reader().run(filename)[Model][0]
        self.assertEqual(model_2.wc_lang_version, version)

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['wc_lang', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: wc_lang')

        with mock.patch('sys.argv', ['wc_lang']):
            with CaptureOutput() as capturer:
                __main__.main()
                self.assertRegex(capturer.get_text(), 'usage: wc_lang')

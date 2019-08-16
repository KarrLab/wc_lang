""" Tests of command line program

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
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
import git
import mock
import os
import re
import unittest
import wc_lang


class TestCli(unittest.TestCase):

    def setUp(self):
        self.tempdir = mkdtemp()

    def tearDown(self):
        rmtree(self.tempdir)

    def test_get_version(self):
        with CaptureOutput(relay=False) as capturer:
            with __main__.App(argv=['-v']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(capturer.get_text(), wc_lang.__version__)

        with CaptureOutput(relay=False) as capturer:
            with __main__.App(argv=['--version']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(capturer.get_text(), wc_lang.__version__)

    def test_cut_submodels(self):
        timestamp = datetime.datetime(2018, 1, 1, 12, 0, 0)
        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1',
                      created=timestamp, updated=timestamp)
        model.submodels.create(id='submodel_1')
        model.submodels.create(id='submodel_2')
        model.references.create(id='ref_0')
        model.references.create(id='ref_1', submodels=model.submodels[0:1])
        model.references.create(id='ref_2', submodels=model.submodels[1:2])
        model.references.create(id='ref_3', submodels=model.submodels[0:2])

        in_path = path.join(self.tempdir, 'in.xlsx')
        Writer().run(in_path, model, data_repo_metadata=False)

        out_path = path.join(self.tempdir, 'out')
        with __main__.App(argv=['cut-submodels', in_path, out_path]) as app:
            app.run()

        model_0 = Reader().run(path.join(out_path, 'core.xlsx'))[Model][0]
        model_1 = Reader().run(path.join(out_path, 'submodel_1.xlsx'))[Model][0]
        model_2 = Reader().run(path.join(out_path, 'submodel_2.xlsx'))[Model][0]

        exp_model_0 = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1',
                            created=timestamp, updated=timestamp)
        exp_model_0.references.create(id='ref_0')
        self.assertTrue(model_0.is_equal(exp_model_0))

        exp_model_1 = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1',
                            created=timestamp, updated=timestamp)
        exp_model_1.submodels.create(id='submodel_1')
        exp_model_1.references.create(id='ref_1', submodels=exp_model_1.submodels)
        exp_model_1.references.create(id='ref_3', submodels=exp_model_1.submodels)
        self.assertTrue(model_1.is_equal(exp_model_1))

        exp_model_2 = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1',
                            created=timestamp, updated=timestamp)
        exp_model_2.submodels.create(id='submodel_2')
        exp_model_2.references.create(id='ref_2', submodels=exp_model_2.submodels)
        exp_model_2.references.create(id='ref_3', submodels=exp_model_2.submodels)
        self.assertTrue(model_2.is_equal(exp_model_2))

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
        Writer().run(in_paths[0], model_0, data_repo_metadata=False)
        Writer().run(in_paths[1], model_1, data_repo_metadata=False)

        # merge models
        with __main__.App(argv=['merge-models', '-p', in_paths[0], '-s', in_paths[1], '-o', out_path]) as app:
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
        Writer().run(filename, model, data_repo_metadata=False)

        with CaptureOutput(relay=False) as capturer:
            with __main__.App(argv=['validate', filename]) as app:
                app.run()
            self.assertEqual(capturer.get_text(), 'Model is valid')

    def test_validate_exception(self):
        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1')
        model.parameters.append(Parameter(id='param_1', value=1., units=unit_registry.parse_units('dimensionless')))
        model.parameters.append(Parameter(id='param_1', value=1., units=unit_registry.parse_units('dimensionless')))

        self.assertNotEqual(Validator().run(model, get_related=True), None)
        filename = path.join(self.tempdir, 'model.xlsx')
        Writer().run(filename, model, data_repo_metadata=False)

        with self.assertRaisesRegex(SystemExit, '^Model is invalid: '):
            with __main__.App(argv=['validate', filename]) as app:
                app.run()

    def test_difference(self):
        now = datetime.datetime.now().replace(microsecond=0)

        model1 = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0', created=now, updated=now)
        filename1 = path.join(self.tempdir, 'model1.xlsx')
        Writer().run(filename1, model1, data_repo_metadata=False)

        model2 = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0', created=now, updated=now)
        filename2 = path.join(self.tempdir, 'model2.xlsx')
        Writer().run(filename2, model2, data_repo_metadata=False)

        model3 = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.1', created=now, updated=now)
        filename3 = path.join(self.tempdir, 'model3.xlsx')
        Writer().run(filename3, model3, data_repo_metadata=False)

        with CaptureOutput(relay=False) as capturer:
            with __main__.App(argv=['difference', filename1, filename2]) as app:
                app.run()
            self.assertEqual(capturer.get_text(), 'Models are identical')

        with CaptureOutput(relay=False) as capturer:
            with __main__.App(argv=['difference', filename1, filename2, '--compare-files']) as app:
                app.run()
            self.assertEqual(capturer.get_text(), 'Models are identical')

        with CaptureOutput(relay=False) as capturer:
            with __main__.App(argv=['difference', filename1, filename3]) as app:
                app.run()
            diff = ('Objects (Model: "model", Model: "model") have different attribute values:\n  '
                    '`wc_lang_version` are not equal:\n    0.0.0 != 0.0.1')
            self.assertEqual(capturer.get_text(), diff)

        with CaptureOutput(relay=False) as capturer:
            with __main__.App(argv=['difference', filename1, filename3, '--compare-files']) as app:
                app.run()
            diff = 'Sheet Model:\n  Row 7:\n    Cell B: 0.0.0 != 0.0.1'
            self.assertEqual(capturer.get_text(), diff)

    def test_transform(self):
        source = path.join(self.tempdir, 'source.xlsx')
        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0')
        Writer().run(source, model, data_repo_metadata=False)

        dest = path.join(self.tempdir, 'dest.xlsx')
        with __main__.App(argv=['transform', source, dest, '--transform', 'MergeAlgorithmicallyLikeSubmodels']) as app:
            app.run()

        self.assertTrue(path.isfile(dest))

    def test_transform_exception(self):
        source = path.join(self.tempdir, 'source.xlsx')
        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0')
        Writer().run(source, model, data_repo_metadata=False)

        dest = path.join(self.tempdir, 'dest.xlsx')
        with self.assertRaisesRegex(SystemExit, 'Please select at least one transform'):
            with __main__.App(argv=['transform', source, dest]) as app:
                app.run()

    def test_normalize(self):
        filename_xls_1 = path.join(self.tempdir, 'model-1.xlsx')
        filename_xls_2 = path.join(self.tempdir, 'model-2.xlsx')

        model = Model(id='model', name='test model', version='0.0.1a', wc_lang_version='0.0.0')
        Writer().run(filename_xls_1, model, data_repo_metadata=False)

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
        Writer().run(filename_xls, model, data_repo_metadata=False)

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
        Writer().run(filename, model, data_repo_metadata=False)

        with __main__.App(argv=['update-version-metadata', filename, '--ignore-repo-metadata']) as app:
            app.run()

        model = Reader().run(filename)[Model][0]
        self.assertEqual(model.wc_lang_version, wc_lang.__version__)

    def test_export_import(self):
        filename = path.join(path.dirname(__file__), 'fixtures', 'sbml-io.xlsx')
        filename_transformed = path.join(path.dirname(__file__), 'fixtures', 'sbml-io-transformed.xlsx')
        sbml_dir = self.tempdir
        filename_2 = path.join(self.tempdir, 'export.xlsx')

        with __main__.App(argv=['export', filename, sbml_dir]) as app:
            app.run()

        with __main__.App(argv=['import', sbml_dir, filename_2]) as app:
            app.run()

        model = Reader().run(filename_transformed)[Model][0]
        model_2 = Reader().run(filename_2)[Model][0]
        self.assertTrue(model_2.is_equal(model))

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['wc-lang', '--help']):
            with CaptureOutput(relay=False):
                with self.assertRaises(SystemExit) as context:
                    __main__.main()
                    self.assertRegex(context.Exception, 'usage: wc-lang')

        with mock.patch('sys.argv', ['wc-lang']):
            with CaptureOutput(relay=False) as capturer:
                __main__.main()
                self.assertRegex(capturer.get_text(), 'usage: wc-lang')

    def test_migration_handlers(self):
        with CaptureOutput(relay=False) as capturer:
            with __main__.App(argv=['make-changes-template']) as app:
                app.run()
                m = re.search(r"Created and added template schema changes file: '(.+)'",
                    capturer.get_text())
                filename = m.group(1)
                self.assertTrue(os.path.isfile(filename))
                # remove schema changes file from the git repo
                repo = git.Repo('.')
                repo.index.remove([filename])
                os.remove(filename)
                dir = os.path.dirname(filename)
                try:
                    os.rmdir(dir)
                except OSError as ex:
                    pass

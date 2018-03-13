""" Tests of the model generator

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_lang import model_gen
import unittest
import wc_kb


class TestModelGenerator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.knowledge_base = kb = wc_kb.KnowledgeBase()

    def test___init__(self):
        gen = model_gen.ModelGenerator(self.knowledge_base)
        self.assertEqual(gen.knowledge_base, self.knowledge_base)

    def test_init_model_generator(self):
        gen = model_gen.ModelGenerator(
            self.knowledge_base,
            component_generators=[model_gen.InitModelGenerator],
            options={'InitModelGenerator': {'id': 'test', 'version': 'test version'}})
        model = gen.run()

        self.assertEqual(model.id, 'test')
        self.assertEqual(model.version, 'test version')

    def test_run(self):
        gen = model_gen.ModelGenerator(
            self.knowledge_base,
            component_generators=[model_gen.InitModelGenerator],
            options={
                'InitModelGenerator': {
                    'id': 'test',
                },
            },
        )
        model = gen.run()

        self.assertEqual(model.id, 'test')
        self.assertEqual(model.version, None)
        self.assertEqual(model.species_types, [])

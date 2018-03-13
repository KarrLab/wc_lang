""" Base classes for generating wc_lang-formatted models from a knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import abc
import six
import wc_kb
import wc_lang.core


class ModelComponentGenerator(six.with_metaclass(abc.ABCMeta, object)):
    """ Generator of a component of a model

    Attributes:
        knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
        model (:obj:`wc_lang.core.Model`): model
    """

    def __init__(self, knowledge_base, model):
        """
        Args:
            knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
            model (:obj:`wc_lang.core.Model`): model
        """
        self.knowledge_base = knowledge_base
        self.model = model

    @abc.abstractmethod
    def run(self):
        """ Run the generator """
        pass  # pragma: no cover


class InitModelGenerator(ModelComponentGenerator):
    """ Initialize a model

    Attributes:
        id (:obj:`str`): model identifier
        version (:obj:`str`): model version
    """

    def __init__(self, knowledge_base, model, id=None, version=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
            model (:obj:`wc_lang.core.Model`): model
            id (:obj:`str`, optional): model identifier
            version (:obj:`str`, optional): model version
        """
        super(InitModelGenerator, self).__init__(knowledge_base, model)
        self.id = id
        self.version = version

    def run(self):
        """ Run the generator """
        self.model.id = self.id
        self.model.version = self.version


class ModelGenerator(object):
    """ Machinery for generating a model (an instance of :obj:`wc_lang.core.Model`) from a knowledge base
    (an instance of :obj:`wc_kb.KnowledgeBase`)

    Attributes:
        knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
        component_generators (:obj:`list` of :obj:`ModelComponentGenerator`): model component generators
        options (:obj:`dict`, optional): dictionary of options whose keys are methods and values are
            optional arguments to the methods
    """

    DEFAULT_MODEL_COMPONENT_GENARATORS = ()

    def __init__(self, knowledge_base, component_generators=None, options=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
            component_generators (:obj:`tuple` of :obj:`ModelComponentGenerator`): model component generators
            options (:obj:`dict`, optional): dictionary of options whose keys are method names and values are
                optional arguments to the methods
        """
        self.knowledge_base = knowledge_base
        self.component_generators = component_generators or self.DEFAULT_MODEL_COMPONENT_GENARATORS
        self.options = options or {}

    def run(self):
        """ Generate a model from a knowledge base

        Returns:
            :obj:`wc_lang.core.Model`: model
        """
        # create model
        model = wc_lang.core.Model()

        # run component generators
        for gen_cls in self.component_generators:
            gen = gen_cls(self.knowledge_base, model, **self.options.get(gen_cls.__name__, {}))
            gen.run()

        # return model
        return model

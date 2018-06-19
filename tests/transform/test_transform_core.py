""" Tests of model transforms.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_lang.transform import get_transforms, MergeAlgorithmicallyLikeSubmodelsTransform
import unittest


class TransformCoreTestCase(unittest.TestCase):
    """ Test model transforms """

    def test_get_transforms(self):
        transforms = get_transforms()
        self.assertEqual(transforms[MergeAlgorithmicallyLikeSubmodelsTransform.Meta.id],
                         MergeAlgorithmicallyLikeSubmodelsTransform)

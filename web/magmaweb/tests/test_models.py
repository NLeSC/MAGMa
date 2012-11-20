import unittest
import magmaweb.models as mm


class TestReactionSequence(unittest.TestCase):
    def setUp(self):
        self.rs = mm.ReactionSequence()

    def test_set(self):
        reactions = ['PARENT', 'CHILD']
        r = self.rs.process_bind_param(reactions, 'sqlite')
        self.assertEqual(r, 'PARENT\nCHILD')

    def test_set_none(self):
        reactions = None
        r = self.rs.process_bind_param(reactions, 'sqlite')
        self.assertIsNone(r)

    def test_get(self):
        reactions = 'PARENT\nCHILD'
        r = self.rs.process_result_value(reactions, 'sqlite')
        self.assertEqual(r, ['PARENT', 'CHILD'])

    def test_get_empty(self):
        reactions = ''
        r = self.rs.process_result_value(reactions, 'sqlite')
        self.assertEqual(r, [])

    def test_get_none(self):
        reactions = None
        r = self.rs.process_result_value(reactions, 'sqlite')
        self.assertEqual(r, None)

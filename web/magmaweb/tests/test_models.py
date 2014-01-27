import unittest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from magmaweb.models import ReactionSequence, fill_molecules_reactionsequence, Metabolite, Reaction


class TestReactionSequence(unittest.TestCase):
    reactions = {
       'reactantof': [{
          'esterase': {'nr': 123, 'nrp': 45}
       }],
       'productof': [{
          'theogallin': {'nr': 678, 'nrp': 90}
       }]
    }
    reactions_as_json = u'{"reactantof": [{"esterase": {"nr": 123, "nrp": 45}}], "productof": [{"theogallin": {"nr": 678, "nrp": 90}}]}'

    def setUp(self):
        self.rs = ReactionSequence()

    def test_set(self):
        r = self.rs.process_bind_param(self.reactions, 'sqlite')
        self.assertEqual(r, self.reactions_as_json)

    def test_set_none(self):
        reactions = None
        r = self.rs.process_bind_param(reactions, 'sqlite')
        self.assertIsNone(r)

    def test_get(self):
        r = self.rs.process_result_value(self.reactions_as_json, 'sqlite')
        self.assertEqual(r, self.reactions)

    def test_get_empty(self):
        reactions = u''
        r = self.rs.process_result_value(reactions, 'sqlite')
        self.assertEqual(r, {})

    def test_get_none(self):
        reactions = None
        r = self.rs.process_result_value(reactions, 'sqlite')
        self.assertEqual(r, None)

    def test_get_badjson2empty(self):
        reactions = u'PARENT'
        r = self.rs.process_result_value(reactions, 'sqlite')
        self.assertEqual(r, {})


class TestReactionSequenceFiller(unittest.TestCase):
    def setUp(self):
        # default use a in memory db
        url = 'sqlite://'
        engine = create_engine(url)
        self.Session = sessionmaker(bind=engine)
        from magmaweb.models import Base
        Base.metadata.create_all(engine)  # @UndefinedVariable

    def getReactionSequence(self, metid):
        return self.Session().query(Metabolite).get(metid).reactionsequence

    def test_noreactions(self):
        session = self.Session()
        session.add(Metabolite(
            metid=1, nhits=1
        ))
        session.flush()

        fill_molecules_reactionsequence(session)

        self.assertDictEqual(self.getReactionSequence(1), {})

    def test_singlereaction(self):
        """
        1 -1> 2
        """
        session = self.Session()
        session.add(Metabolite(metid=1, nhits=1))
        session.add(Metabolite(metid=2, nhits=1))
        session.add(Reaction(reactid=1, name='esterase',
                             reactant=1, product=2))
        session.flush()

        fill_molecules_reactionsequence(session)

        expected1 = {
            u'reactantof': {
                u'esterase': {
                    u'nr': 1,
                    u'nrp': 1
                }
            }
        }
        self.assertDictEqual(self.getReactionSequence(1), expected1)
        expected2 = {
            u'productof': {
                u'esterase': {
                    u'nr': 1,
                    u'nrp': 1
                }
            }
        }
        self.assertDictEqual(self.getReactionSequence(2), expected2)

    def test_reactionson2molecules(self):
        """
        1 -1> 2 -2> 1
        """
        session = self.Session()
        session.add(Metabolite(metid=1, nhits=1))
        session.add(Metabolite(metid=2, nhits=1))
        session.add(Reaction(reactid=1, name='esterase',
                             reactant=1, product=2))
        session.add(Reaction(reactid=2, name='theogallin',
                             reactant=2, product=1))
        session.flush()

        fill_molecules_reactionsequence(session)

        expected1 = {
            u'reactantof': {
                u'esterase': {
                    u'nr': 1,
                    u'nrp': 1
                }
            },
            u'productof': {
                u'theogallin': {
                    u'nr': 1,
                    u'nrp': 1
                }
            }
        }
        self.assertDictEqual(self.getReactionSequence(1), expected1)
        expected2 = {
            u'reactantof': {
                u'theogallin': {
                    u'nr': 1,
                    u'nrp': 1
                }
            },
            u'productof': {
                u'esterase': {
                    u'nr': 1,
                    u'nrp': 1
                }
            }
        }
        self.assertDictEqual(self.getReactionSequence(2), expected2)

    def test_reactionson4molecules(self):
        """
        1 -1> 2 -3> 4
        1 -2> 3 -4> 4
        """
        session = self.Session()
        session.add(Metabolite(metid=1, nhits=1))
        session.add(Metabolite(metid=2, nhits=1))
        session.add(Metabolite(metid=3, nhits=0))
        session.add(Metabolite(metid=4, nhits=1))
        session.add(Reaction(reactid=1, name='esterase',
                             reactant=1, product=2))
        session.add(Reaction(reactid=2, name='theogallin',
                             reactant=1, product=3))
        session.add(Reaction(reactid=3, name='dehydrox',
                             reactant=2, product=4))
        session.add(Reaction(reactid=4, name='reduc',
                             reactant=3, product=4))
        session.flush()

        fill_molecules_reactionsequence(session)

        expected1 = {
            u'reactantof': {u'esterase': {u'nr': 1, u'nrp': 1},
            u'theogallin': {u'nr': 1, u'nrp': 0}}
        }
        self.assertDictEqual(self.getReactionSequence(1), expected1)
        expected2 = {
            u'productof': {u'esterase': {u'nr': 1, u'nrp': 1}},
            u'reactantof': {u'dehydrox': {u'nr': 1, u'nrp': 1}}
        }
        self.assertDictEqual(self.getReactionSequence(2), expected2)
        expected3 = {
            u'productof': {u'theogallin': {u'nr': 1, u'nrp': 1}},
            u'reactantof': {u'reduc': {u'nr': 1, u'nrp': 1}}
        }
        self.assertDictEqual(self.getReactionSequence(3), expected3)
        expected4 = {
            u'productof': {u'dehydrox': {u'nr': 1, u'nrp': 1},
            u'reduc': {u'nr': 1, u'nrp': 0}}
        }
        self.assertDictEqual(self.getReactionSequence(4), expected4)

    def test_reactionson1reaction2products(self):
        """
        1 -1> 2
        1 -1> 3
        """
        session = self.Session()
        session.add(Metabolite(metid=1, nhits=1))
        session.add(Metabolite(metid=2, nhits=1))
        session.add(Metabolite(metid=3, nhits=1))
        session.add(Reaction(reactid=1, name='esterase',
                             reactant=1, product=2))
        session.add(Reaction(reactid=2, name='esterase',
                             reactant=1, product=3))
        session.flush()

        fill_molecules_reactionsequence(session)

        expected1 = {
            u'reactantof': {u'esterase': {u'nr': 2, u'nrp': 2}},
        }
        self.assertDictEqual(self.getReactionSequence(1), expected1)
        expected2 = {
            u'productof': {u'esterase': {u'nr': 1, u'nrp': 1}},
        }
        self.assertDictEqual(self.getReactionSequence(2), expected2)
        expected3 = {
            u'productof': {u'esterase': {u'nr': 1, u'nrp': 1}},
        }
        self.assertDictEqual(self.getReactionSequence(3), expected3)

    def test_reactionson1reaction2reactants(self):
        """
        1 -1> 2
        3 -1> 2
        """
        session = self.Session()
        session.add(Metabolite(metid=1, nhits=1))
        session.add(Metabolite(metid=2, nhits=1))
        session.add(Metabolite(metid=3, nhits=1))
        session.add(Reaction(reactid=1, name='esterase',
                             reactant=1, product=2))
        session.add(Reaction(reactid=2, name='esterase',
                             reactant=3, product=2))
        session.flush()

        fill_molecules_reactionsequence(session)

        expected1 = {
            u'reactantof': {u'esterase': {u'nr': 1, u'nrp': 1}},
        }
        self.assertDictEqual(self.getReactionSequence(1), expected1)
        expected2 = {
            u'productof': {u'esterase': {u'nr': 2, u'nrp': 2}},
        }
        self.assertDictEqual(self.getReactionSequence(2), expected2)
        expected3 = {
            u'reactantof': {u'esterase': {u'nr': 1, u'nrp': 1}},
        }
        self.assertDictEqual(self.getReactionSequence(3), expected3)

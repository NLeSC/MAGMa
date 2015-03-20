import unittest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from magmaweb.models import ReactionSequence, Molecule, Reaction
from magmaweb.models import fill_molecules_reactions


class TestReactionSequence(unittest.TestCase):
    reactions = {
       'products': [{
          'esterase': {'nr': 123, 'nrp': 45}
       }],
       'reactants': [{
          'theogallin': {'nr': 678, 'nrp': 90}
       }]
    }
    reactions_json = '{"reactants": [{"theogallin": {"nr": 678, "nrp": 90}}], '
    reactions_json += '"products": [{"esterase": {"nr": 123, "nrp": 45}}]}'

    def setUp(self):
        self.rs = ReactionSequence()

    def test_set(self):
        r = self.rs.process_bind_param(self.reactions, 'sqlite')
        self.assertEqual(r, self.reactions_json)

    def test_set_none(self):
        reactions = None
        r = self.rs.process_bind_param(reactions, 'sqlite')
        self.assertIsNone(r)

    def test_get(self):
        r = self.rs.process_result_value(self.reactions_json, 'sqlite')
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


class TestFillMoleculeReactions(unittest.TestCase):
    def setUp(self):
        # default use a in memory db
        url = 'sqlite://'
        engine = create_engine(url)
        self.Session = sessionmaker(bind=engine)
        from magmaweb.models import Base
        Base.metadata.create_all(engine)  # @UndefinedVariable

    def getReactionSequence(self, molid):
        return self.Session().query(Molecule).get(molid).reactionsequence

    def test_noreactions(self):
        session = self.Session()
        session.add(Molecule(
            molid=1, nhits=1
        ))
        session.flush()

        fill_molecules_reactions(session)

        self.assertDictEqual(self.getReactionSequence(1), {})

    def test_singlereaction(self):
        """
        1 -1> 2
        """
        session = self.Session()
        session.add(Molecule(molid=1, nhits=1))
        session.add(Molecule(molid=2, nhits=1))
        session.add(Reaction(reactid=1, name='esterase',
                             reactant=1, product=2))
        session.flush()

        fill_molecules_reactions(session)

        expected1 = {
            u'products': {
                u'esterase': {
                    u'nr': 1,
                    u'nrp': 1
                }
            }
        }
        self.assertDictEqual(self.getReactionSequence(1), expected1)
        expected2 = {
            u'reactants': {
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
        session.add(Molecule(molid=1, nhits=1))
        session.add(Molecule(molid=2, nhits=1))
        session.add(Reaction(reactid=1, name='esterase',
                             reactant=1, product=2))
        session.add(Reaction(reactid=2, name='theogallin',
                             reactant=2, product=1))
        session.flush()

        fill_molecules_reactions(session)

        expected1 = {
            u'products': {
                u'esterase': {
                    u'nr': 1,
                    u'nrp': 1
                }
            },
            u'reactants': {
                u'theogallin': {
                    u'nr': 1,
                    u'nrp': 1
                }
            }
        }
        self.assertDictEqual(self.getReactionSequence(1), expected1)
        expected2 = {
            u'products': {
                u'theogallin': {
                    u'nr': 1,
                    u'nrp': 1
                }
            },
            u'reactants': {
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
        session.add(Molecule(molid=1, nhits=1))
        session.add(Molecule(molid=2, nhits=1))
        session.add(Molecule(molid=3, nhits=0))
        session.add(Molecule(molid=4, nhits=1))
        session.add(Reaction(reactid=1, name='esterase',
                             reactant=1, product=2))
        session.add(Reaction(reactid=2, name='theogallin',
                             reactant=1, product=3))
        session.add(Reaction(reactid=3, name='dehydrox',
                             reactant=2, product=4))
        session.add(Reaction(reactid=4, name='reduc',
                             reactant=3, product=4))
        session.flush()

        fill_molecules_reactions(session)

        expected1 = {
            u'products': {u'esterase': {u'nr': 1, u'nrp': 1},
                          u'theogallin': {u'nr': 1}}
        }
        self.assertDictEqual(self.getReactionSequence(1), expected1)
        expected2 = {
            u'reactants': {u'esterase': {u'nr': 1, u'nrp': 1}},
            u'products': {u'dehydrox': {u'nr': 1, u'nrp': 1}}
        }
        self.assertDictEqual(self.getReactionSequence(2), expected2)
        expected3 = {
            u'reactants': {u'theogallin': {u'nr': 1, u'nrp': 1}},
            u'products': {u'reduc': {u'nr': 1, u'nrp': 1}}
        }
        self.assertDictEqual(self.getReactionSequence(3), expected3)
        expected4 = {
            u'reactants': {u'dehydrox': {u'nr': 1, u'nrp': 1},
                           u'reduc': {u'nr': 1}}
        }
        self.assertDictEqual(self.getReactionSequence(4), expected4)

    def test_reactionson1reaction2reactants(self):
        """
        1 -1> 2
        1 -1> 3
        """
        session = self.Session()
        session.add(Molecule(molid=1, nhits=1))
        session.add(Molecule(molid=2, nhits=1))
        session.add(Molecule(molid=3, nhits=1))
        session.add(Reaction(reactid=1, name='esterase',
                             reactant=1, product=2))
        session.add(Reaction(reactid=2, name='esterase',
                             reactant=1, product=3))
        session.flush()

        fill_molecules_reactions(session)

        expected1 = {
            u'products': {u'esterase': {u'nr': 2, u'nrp': 2}},
        }
        self.assertDictEqual(self.getReactionSequence(1), expected1)
        expected2 = {
            u'reactants': {u'esterase': {u'nr': 1, u'nrp': 1}},
        }
        self.assertDictEqual(self.getReactionSequence(2), expected2)
        expected3 = {
            u'reactants': {u'esterase': {u'nr': 1, u'nrp': 1}},
        }
        self.assertDictEqual(self.getReactionSequence(3), expected3)

    def test_reactionson1reaction2products(self):
        """
        1 -1> 2
        3 -1> 2
        """
        session = self.Session()
        session.add(Molecule(molid=1, nhits=1))
        session.add(Molecule(molid=2, nhits=1))
        session.add(Molecule(molid=3, nhits=1))
        session.add(Reaction(reactid=1, name='esterase',
                             reactant=1, product=2))
        session.add(Reaction(reactid=2, name='esterase',
                             reactant=3, product=2))
        session.flush()

        fill_molecules_reactions(session)

        expected1 = {
            u'products': {u'esterase': {u'nr': 1, u'nrp': 1}},
        }
        self.assertDictEqual(self.getReactionSequence(1), expected1)
        expected2 = {
            u'reactants': {u'esterase': {u'nr': 2, u'nrp': 2}},
        }
        self.assertDictEqual(self.getReactionSequence(2), expected2)
        expected3 = {
            u'products': {u'esterase': {u'nr': 1, u'nrp': 1}},
        }
        self.assertDictEqual(self.getReactionSequence(3), expected3)

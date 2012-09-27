import unittest
from pyramid import testing


class Test_jsonhtml_renderer_factory(unittest.TestCase):
    def setUp(self):
        self.config = testing.setUp()

    def tearDown(self):
        testing.tearDown()

    def _callFUT(self, name):
        from magmaweb import jsonhtml_renderer_factory
        return jsonhtml_renderer_factory(name)

    def test_body(self):
        renderer = self._callFUT(None)
        result = renderer({'a':1}, {})
        self.assertEqual(result, '{"a": 1}')

    def test_content_type(self):
        request = testing.DummyRequest()
        renderer = self._callFUT(None)
        renderer({'a':1}, {'request':request})
        self.assertEqual(request.response.content_type, 'text/html')

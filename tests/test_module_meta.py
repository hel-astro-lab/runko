import unittest

class runko_module_meta(unittest.TestCase):
    def test_runko_is_importable(self):
        import importlib.util
        runko_spec = importlib.util.find_spec("runko")
        found = runko_spec is not None
        self.assertTrue(found)

if __name__ == "__main__":
    unittest.main()

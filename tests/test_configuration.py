import unittest
import runko

class configuration(unittest.TestCase):
    def setUp(self):
        from pathlib import Path
        source_dir = Path(__file__).resolve().parent
        self.conf_A = runko.Configuration(str(source_dir) + "/configuration_A.ini")

    def test_all_sections_are_read(self):
        self.assertEqual(self.conf_A.int_var, 42)
        self.assertEqual(self.conf_A.string_var, "foo bar")
        self.assertEqual(self.conf_A.float_var, 4.2)
        self.assertEqual(self.conf_A.string_array_var, ["foo", "bar"])
        self.assertTrue(self.conf_A.true_var)
        self.assertFalse(self.conf_A.false_var)

    def test_configuration_can_be_edited_dynamically(self):
        import copy

        conf = copy.deepcopy(self.conf_A)

        self.assertIsNone(conf.foo)
        conf.foo = "bar"
        self.assertEqual(conf.foo, "bar")

    def test_special_attributes_raise_attribute_error(self):
        with self.assertRaises(AttributeError):
            self.conf_A.__foobar__

    def test_empty_conf(self):
        conf = runko.Configuration(None)

        self.assertIsNone(conf.foo)
        conf.foo = "bar"
        self.assertEqual(conf.foo, "bar")

if __name__ == "__main__":
    unittest.main()

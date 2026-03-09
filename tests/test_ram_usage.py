import unittest


class TestGetRssKB(unittest.TestCase):
    """Tests for the cross-platform RAM measurement helper."""

    def test_returns_positive_integer(self):
        from runko.ram_usage import get_rss_kB
        rss = get_rss_kB()
        self.assertIsInstance(rss, int)
        self.assertGreater(rss, 0)

    def test_returns_reasonable_value(self):
        """RSS of a Python process should be at least a few MB."""
        from runko.ram_usage import get_rss_kB
        rss = get_rss_kB()
        # Python interpreter alone is typically 10-50 MB;
        # sanity-check that we get at least 1 MB (1024 kB).
        self.assertGreater(rss, 1024)

    def test_allocate_increases_rss(self):
        """Allocating a large block should increase reported RSS."""
        from runko.ram_usage import get_rss_kB
        before = get_rss_kB()
        # Allocate ~100 MB and touch it so it is resident
        block = bytearray(100 * 1024 * 1024)
        after = get_rss_kB()
        self.assertGreater(after, before)
        del block


if __name__ == "__main__":
    unittest.main()

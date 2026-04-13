import unittest
import numpy as np

from runko.hilbert import Hilbert2D, Hilbert3D

CONFIGS_2D = [(1, 1), (2, 2), (3, 3), (2, 1), (1, 2), (3, 2), (2, 3), (4, 2)]
CONFIGS_3D = [(1, 1, 1), (2, 2, 2), (2, 1, 1), (1, 2, 1), (1, 1, 2),
              (3, 2, 1), (1, 2, 3), (2, 2, 1)]


class TestHilbert2D(unittest.TestCase):

    def test_roundtrip(self):
        for m0, m1 in CONFIGS_2D:
            with self.subTest(m0=m0, m1=m1):
                gen = Hilbert2D(m0, m1)
                for x in range(1 << m0):
                    for y in range(1 << m1):
                        h = gen.hindex(x, y)
                        self.assertEqual(gen.inv(h), (x, y))

    def test_reverse_roundtrip(self):
        for m0, m1 in CONFIGS_2D:
            with self.subTest(m0=m0, m1=m1):
                gen = Hilbert2D(m0, m1)
                n_total = 1 << (m0 + m1)
                for h in range(n_total):
                    self.assertEqual(gen.hindex(*gen.inv(h)), h)

    def test_permutation(self):
        for m0, m1 in CONFIGS_2D:
            with self.subTest(m0=m0, m1=m1):
                gen = Hilbert2D(m0, m1)
                n_total = 1 << (m0 + m1)
                indices = set()
                for x in range(1 << m0):
                    for y in range(1 << m1):
                        indices.add(gen.hindex(x, y))
                self.assertEqual(indices, set(range(n_total)))


class TestHilbert3D(unittest.TestCase):

    def test_roundtrip(self):
        for m0, m1, m2 in CONFIGS_3D:
            with self.subTest(m0=m0, m1=m1, m2=m2):
                gen = Hilbert3D(m0, m1, m2)
                for x in range(1 << m0):
                    for y in range(1 << m1):
                        for z in range(1 << m2):
                            h = gen.hindex(x, y, z)
                            self.assertEqual(gen.inv(h), (x, y, z))

    def test_reverse_roundtrip(self):
        for m0, m1, m2 in CONFIGS_3D:
            with self.subTest(m0=m0, m1=m1, m2=m2):
                gen = Hilbert3D(m0, m1, m2)
                n_total = 1 << (m0 + m1 + m2)
                for h in range(n_total):
                    self.assertEqual(gen.hindex(*gen.inv(h)), h)

    def test_permutation(self):
        for m0, m1, m2 in CONFIGS_3D:
            with self.subTest(m0=m0, m1=m1, m2=m2):
                gen = Hilbert3D(m0, m1, m2)
                n_total = 1 << (m0 + m1 + m2)
                indices = set()
                for x in range(1 << m0):
                    for y in range(1 << m1):
                        for z in range(1 << m2):
                            indices.add(gen.hindex(x, y, z))
                self.assertEqual(indices, set(range(n_total)))

    def test_vec_matches_scalar(self):
        for m0, m1, m2 in CONFIGS_3D:
            with self.subTest(m0=m0, m1=m1, m2=m2):
                gen = Hilbert3D(m0, m1, m2)
                xs = np.arange(1 << m0, dtype=np.int64)
                ys = np.arange(1 << m1, dtype=np.int64)
                zs = np.arange(1 << m2, dtype=np.int64)
                xx, yy, zz = np.meshgrid(xs, ys, zs, indexing='ij')
                x_flat = xx.ravel()
                y_flat = yy.ravel()
                z_flat = zz.ravel()

                vec_result = gen.hindex_vec(x_flat, y_flat, z_flat)

                scalar_result = np.array([
                    gen.hindex(int(x_flat[i]), int(y_flat[i]), int(z_flat[i]))
                    for i in range(len(x_flat))
                ], dtype=np.int64)

                np.testing.assert_array_equal(vec_result, scalar_result)


if __name__ == "__main__":
    unittest.main()

import unittest

import sys
sys.path.append('python')

import plasmatools as pl





class Basics(unittest.TestCase):

    def setUp(self):
        self.s1 = pl.Sheet()

    def test_resize(self):
        self.s1.resize(10,5)
        self.assertEqual(self.s1.Ni, 10)
        self.assertEqual(self.s1.Nj, 5)

    def test_container(self):
        Ni = 20
        Nj = 10
        self.s1.resize(Ni, Nj)

        self.assertEqual( len(self.s1.iGrid), Ni )
        self.assertEqual( len(self.s1.jGrid), Nj )
        self.assertEqual( len(self.s1.values), Ni*Nj )


class Arithmetics(unittest.TestCase):

    Ni = 10
    Nj = 5

    def setUp(self):
        self.s1 = pl.Sheet()
        self.s1.resize(self.Ni, self.Nj)

        self.s2 = pl.Sheet()
        self.s2.resize(self.Ni, self.Nj)

        for i in range(self.Ni):
            for j in range(self.Nj):
                self.s1.load_value(i, j, 1.0)
                self.s2.load_value(i, j, 2.0)


    # require that every element in sheet is equal to the val
    # this also unpacks the blocks
    def assertAll(self, s1, val):
        for i in range(self.Ni):
            for j in range(self.Nj):
                block = s1.get_block(i,j)
                self.assertEqual( block[0], val )


    def test_iadd1(self):
        self.s1 += self.s2
        self.assertAll(self.s1, 3.0)


    def test_iadd2(self):
        self.s2 += self.s1
        self.assertAll(self.s2, 3.0)


    def test_add(self):
        s3 = self.s1 + self.s2
        self.assertAll(s3, 3.0)


    def test_sub(self):
        s3 = self.s2 - self.s1
        self.assertAll(s3, 1.0)


    def test_mul(self):
        s3 = self.s2 * 2.0
        self.assertAll(s3, 4.0)





if __name__ == '__main__':
    unittest.main()

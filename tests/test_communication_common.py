import unittest
from pyrunko.tools import comm_mode as m


class communication_common(unittest.TestCase):

    def test_communication_modes_have_legacy_values(self):
        modes = [m.emf_E, m.emf_B, m.emf_J, m.pic_particle, m.pic_particle_extra]

        # According to scheduler mpi solver case.
        self.assertEqual(m.emf_J.value, 0)
        self.assertEqual(m.emf_E.value, 1)
        self.assertEqual(m.emf_B.value, 2)
        self.assertEqual(m.pic_particle.value, 3)
        self.assertEqual(m.pic_particle_extra.value, 4)

    def test_communication_modes_are_comparable(self):
        modes = [m.emf_E, m.emf_B, m.emf_J, m.pic_particle, m.pic_particle_extra]

        for i in range(len(modes)):
            self.assertEqual(modes[i], modes[i])
            for j in range(i + 1, len(modes)):
                self.assertNotEqual(modes[i], modes[j])
                self.assertNotEqual(modes[j], modes[i])


if __name__ == "__main__":
    unittest.main()

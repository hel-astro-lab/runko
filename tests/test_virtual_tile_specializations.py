import unittest

import runko

class virtual_tile_specializations(unittest.TestCase):

    def test_emf_virtual_tile(self):
        self.assertEqual(runko.emf.threeD.VirtualTile,
                         runko.emf.threeD.Tile.virtual_tile_specialization())

        self.assertEqual(runko.emf.threeD.Tile,
                         runko.emf.threeD.Tile.canonical_type())

        self.assertEqual(runko.emf.threeD.Tile,
                         runko.emf.threeD.VirtualTile.canonical_type())

    def test_pic_virtual_tile(self):
        # self.assertEqual(runko.pic.threeD.VirtualTile,
        #                  runko.pic.threeD.Tile.virtual_tile_specialization())

        self.assertEqual(runko.pic.threeD.Tile,
                         runko.pic.threeD.Tile.canonical_type())

        # self.assertEqual(runko.pic.threeD.Tile,
        #                  runko.pic.threeD.VirtualTile.canonical_type())



if __name__ == "__main__":
    unittest.main()

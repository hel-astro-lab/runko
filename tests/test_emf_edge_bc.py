import unittest
import numpy as np
import runko


def make_tile(tile_idx=(0, 0, 0), Nx=4, Ny=4, Nz=4, NxMesh=10, NyMesh=11, NzMesh=13):
    config = runko.Configuration(None)
    config.Nx = Nx
    config.Ny = Ny
    config.Nz = Nz
    config.NxMesh = NxMesh
    config.NyMesh = NyMesh
    config.NzMesh = NzMesh
    config.xmin = 0
    config.ymin = 0
    config.zmin = 0
    config.cfl = 1
    config.field_propagator = "FDTD2"
    return config, runko.emf.threeD.Tile(tile_idx, config)


nonzero_E = lambda x, y, z: (1.0, 2.0, 3.0)
nonzero_B = lambda x, y, z: (4.0, 5.0, 6.0)
nonzero_J = lambda x, y, z: (0.1, 0.2, 0.3)
zero = lambda x, y, z: (0, 0, 0)


class emf_edge_bc(unittest.TestCase):

    def test_apply_edge_bc_E_only(self):
        """apply_edge_bc with emf_E sets E only; B and J unchanged."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, nonzero_B, nonzero_J)

        _, (Bx0, By0, Bz0), (Jx0, Jy0, Jz0) = tile.get_EBJ()

        bc = runko.emf.threeD.edge_bc(direction=0, side=0, position=5.0)
        tile.apply_edge_bc(bc, runko.tools.comm_mode.emf_E)

        _, (Bx, By, Bz), (Jx, Jy, Jz) = tile.get_EBJ()

        self.assertTrue(np.allclose(Bx, Bx0))
        self.assertTrue(np.allclose(By, By0))
        self.assertTrue(np.allclose(Bz, Bz0))
        self.assertTrue(np.allclose(Jx, Jx0))
        self.assertTrue(np.allclose(Jy, Jy0))
        self.assertTrue(np.allclose(Jz, Jz0))

    def test_apply_edge_bc_B_only(self):
        """apply_edge_bc with emf_B sets B only; E and J unchanged."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, nonzero_B, nonzero_J)

        (Ex0, Ey0, Ez0), _, (Jx0, Jy0, Jz0) = tile.get_EBJ()

        bc = runko.emf.threeD.edge_bc(direction=0, side=0, position=5.0,
                                       Bx=10.0, By=20.0, Bz=30.0)
        tile.apply_edge_bc(bc, runko.tools.comm_mode.emf_B)

        (Ex, Ey, Ez), _, (Jx, Jy, Jz) = tile.get_EBJ()

        self.assertTrue(np.allclose(Ex, Ex0))
        self.assertTrue(np.allclose(Ey, Ey0))
        self.assertTrue(np.allclose(Ez, Ez0))
        self.assertTrue(np.allclose(Jx, Jx0))
        self.assertTrue(np.allclose(Jy, Jy0))
        self.assertTrue(np.allclose(Jz, Jz0))

    def test_apply_edge_bc_J_only(self):
        """apply_edge_bc with emf_J sets J only; E and B unchanged."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, nonzero_B, nonzero_J)

        (Ex0, Ey0, Ez0), (Bx0, By0, Bz0), _ = tile.get_EBJ()

        bc = runko.emf.threeD.edge_bc(direction=0, side=0, position=5.0,
                                       Jx=0.5, Jy=0.6, Jz=0.7)
        tile.apply_edge_bc(bc, runko.tools.comm_mode.emf_J)

        (Ex, Ey, Ez), (Bx, By, Bz), _ = tile.get_EBJ()

        self.assertTrue(np.allclose(Ex, Ex0))
        self.assertTrue(np.allclose(Ey, Ey0))
        self.assertTrue(np.allclose(Ez, Ez0))
        self.assertTrue(np.allclose(Bx, Bx0))
        self.assertTrue(np.allclose(By, By0))
        self.assertTrue(np.allclose(Bz, Bz0))

    def test_left_x_edge_sets_correct_cells(self):
        """Left x-edge: cells at and behind position are set, rest unchanged."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, zero, zero)

        bc = runko.emf.threeD.edge_bc(direction=0, side=0, position=5.0,
                                       Ex=99.0, Ey=99.0, Ez=99.0)
        tile.apply_edge_bc(bc, runko.tools.comm_mode.emf_E)

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        # tile (0,0,0) has mins=(0,0,0), maxs=(10,11,13)
        # position=5.0, side=0 (left): width = floor(5.0 - 0) + 1 = 6
        # edge region in full lattice: [0, halo+6) = [0, 9) in lattice coords
        # non-halo region starts at halo=3, so numpy x-indices 0..5 are set
        iw = 6
        self.assertTrue(np.all(Ex[:iw, :, :] == 99.0))
        self.assertTrue(np.all(Ey[:iw, :, :] == 99.0))
        self.assertTrue(np.all(Ez[:iw, :, :] == 99.0))

        # cells in front are untouched
        self.assertTrue(np.allclose(Ex[iw:, :, :], 1.0))
        self.assertTrue(np.allclose(Ey[iw:, :, :], 2.0))
        self.assertTrue(np.allclose(Ez[iw:, :, :], 3.0))

    def test_right_x_edge_sets_correct_cells(self):
        """Right x-edge: cells at and beyond position are set, rest unchanged."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, zero, zero)

        # NxMesh=10, tile maxs[0]=10
        # position=7.0, side=1 (right): width = 10 - floor(7.0 - 0) = 3
        # edge region: [halo + 10 - 3, halo + 10 + halo) = [10, 16) in lattice
        # non-halo: numpy x-indices 7..9 are set
        bc = runko.emf.threeD.edge_bc(direction=0, side=1, position=7.0,
                                       Ex=88.0, Ey=88.0, Ez=88.0)
        tile.apply_edge_bc(bc, runko.tools.comm_mode.emf_E)

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        iw = 7
        self.assertTrue(np.all(Ex[iw:, :, :] == 88.0))
        self.assertTrue(np.all(Ey[iw:, :, :] == 88.0))
        self.assertTrue(np.all(Ez[iw:, :, :] == 88.0))

        # cells before position untouched
        self.assertTrue(np.allclose(Ex[:iw, :, :], 1.0))
        self.assertTrue(np.allclose(Ey[:iw, :, :], 2.0))
        self.assertTrue(np.allclose(Ez[:iw, :, :], 3.0))

    def test_y_direction_edge(self):
        """Y-direction edge sets correct y-cells."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, zero, zero)

        # NyMesh=11, tile maxs[1]=11
        # direction=1, side=0, position=4.0: width = floor(4-0)+1 = 5
        bc = runko.emf.threeD.edge_bc(direction=1, side=0, position=4.0,
                                       Ex=77.0, Ey=77.0, Ez=77.0)
        tile.apply_edge_bc(bc, runko.tools.comm_mode.emf_E)

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        iw = 5
        self.assertTrue(np.all(Ex[:, :iw, :] == 77.0))
        self.assertTrue(np.allclose(Ex[:, iw:, :], 1.0))

    def test_z_direction_edge(self):
        """Z-direction edge sets correct z-cells."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, zero, zero)

        # NzMesh=13, tile maxs[2]=13
        # direction=2, side=1, position=10.0: width = 13 - floor(10-0) = 3
        bc = runko.emf.threeD.edge_bc(direction=2, side=1, position=10.0,
                                       Ex=66.0, Ey=66.0, Ez=66.0)
        tile.apply_edge_bc(bc, runko.tools.comm_mode.emf_E)

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        iw = 10
        self.assertTrue(np.all(Ex[:, :, iw:] == 66.0))
        self.assertTrue(np.allclose(Ex[:, :, :iw], 1.0))

    def test_component_mask_Ey_Ez_only(self):
        """E_components=0b110 sets only Ey/Ez, leaves Ex untouched."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, zero, zero)

        (Ex0, _, _), _, _ = tile.get_EBJ()

        bc = runko.emf.threeD.edge_bc(direction=0, side=0, position=5.0,
                                       E_components=0b110)
        tile.apply_edge_bc(bc, runko.tools.comm_mode.emf_E)

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        iw = 6
        # Ex unchanged everywhere
        self.assertTrue(np.allclose(Ex, Ex0))

        # Ey and Ez zeroed in edge region
        self.assertTrue(np.all(Ey[:iw, :, :] == 0))
        self.assertTrue(np.all(Ez[:iw, :, :] == 0))

        # Ey and Ez untouched in front
        self.assertTrue(np.allclose(Ey[iw:, :, :], 2.0))
        self.assertTrue(np.allclose(Ez[iw:, :, :], 3.0))

    def test_multi_tile_fully_behind_edge(self):
        """Tile fully behind edge position gets entirely set."""
        # Tile (0,0,0): mins=(0,0,0), maxs=(10,11,13)
        # position=15.0 >= tile_max=10 → width = Nd = 10 (entire tile)
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, zero, zero)

        bc = runko.emf.threeD.edge_bc(direction=0, side=0, position=15.0,
                                       Ex=55.0, Ey=55.0, Ez=55.0)
        tile.apply_edge_bc(bc, runko.tools.comm_mode.emf_E)

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        self.assertTrue(np.all(Ex == 55.0))
        self.assertTrue(np.all(Ey == 55.0))
        self.assertTrue(np.all(Ez == 55.0))

    def test_multi_tile_fully_in_front_unaffected(self):
        """Tile fully in front of edge position is not affected."""
        # Tile (1,0,0): mins=(10,0,0), maxs=(20,11,13)
        # position=5.0, side=0 → position <= tile_min=10 → skip
        _, tile = make_tile(tile_idx=(1, 0, 0))
        tile.set_EBJ(nonzero_E, zero, zero)

        (Ex0, Ey0, Ez0), _, _ = tile.get_EBJ()

        bc = runko.emf.threeD.edge_bc(direction=0, side=0, position=5.0,
                                       Ex=55.0, Ey=55.0, Ez=55.0)
        tile.apply_edge_bc(bc, runko.tools.comm_mode.emf_E)

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        self.assertTrue(np.allclose(Ex, Ex0))
        self.assertTrue(np.allclose(Ey, Ey0))
        self.assertTrue(np.allclose(Ez, Ez0))

    def test_registered_bcs_apply(self):
        """apply_edge_bcs applies all registered BCs."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, zero, zero)

        bc = runko.emf.threeD.edge_bc(direction=0, side=0, position=5.0,
                                       Ex=42.0, Ey=42.0, Ez=42.0)
        tile.register_edge_bc(bc)
        tile.apply_edge_bcs(runko.tools.comm_mode.emf_E)

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        iw = 6
        self.assertTrue(np.all(Ex[:iw, :, :] == 42.0))

    def test_no_registered_bcs_is_noop(self):
        """apply_edge_bcs with no registered BCs is a no-op."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, zero, zero)

        (Ex0, Ey0, Ez0), _, _ = tile.get_EBJ()

        tile.apply_edge_bcs(runko.tools.comm_mode.emf_E)

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        self.assertTrue(np.allclose(Ex, Ex0))
        self.assertTrue(np.allclose(Ey, Ey0))
        self.assertTrue(np.allclose(Ez, Ez0))

    def test_stacking_later_overwrites_earlier(self):
        """Two BCs on same region: later overwrites earlier in overlap."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, zero, zero)

        # First BC: set left 6 cells to 11.0
        bc1 = runko.emf.threeD.edge_bc(direction=0, side=0, position=5.0,
                                        Ex=11.0, Ey=11.0, Ez=11.0)
        # Second BC: set left 4 cells to 22.0 (overwrites part of bc1)
        bc2 = runko.emf.threeD.edge_bc(direction=0, side=0, position=3.0,
                                        Ex=22.0, Ey=22.0, Ez=22.0)

        tile.apply_edge_bc(bc1, runko.tools.comm_mode.emf_E)
        tile.apply_edge_bc(bc2, runko.tools.comm_mode.emf_E)

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        # First 4 cells overwritten by bc2
        self.assertTrue(np.all(Ex[:4, :, :] == 22.0))
        # Cells 4..5 still from bc1
        self.assertTrue(np.all(Ex[4:6, :, :] == 11.0))
        # Rest untouched
        self.assertTrue(np.allclose(Ex[6:, :, :], 1.0))

    def test_stacking_xy_corner(self):
        """x-direction + y-direction BCs compose into a corner."""
        _, tile = make_tile()
        tile.set_EBJ(nonzero_E, zero, zero)

        bc_x = runko.emf.threeD.edge_bc(direction=0, side=0, position=3.0,
                                         Ex=11.0, Ey=11.0, Ez=11.0)
        bc_y = runko.emf.threeD.edge_bc(direction=1, side=0, position=3.0,
                                         Ex=22.0, Ey=22.0, Ez=22.0)

        tile.apply_edge_bc(bc_x, runko.tools.comm_mode.emf_E)
        tile.apply_edge_bc(bc_y, runko.tools.comm_mode.emf_E)

        (Ex, _, _), _, _ = tile.get_EBJ()

        # Corner (x<4, y<4): overwritten by bc_y (applied second)
        self.assertTrue(np.all(Ex[:4, :4, :] == 22.0))

        # x-only strip (x<4, y>=4): set by bc_x only
        self.assertTrue(np.all(Ex[:4, 4:, :] == 11.0))

        # y-only strip (x>=4, y<4): set by bc_y only
        self.assertTrue(np.all(Ex[4:, :4, :] == 22.0))

        # Outside both: original value
        self.assertTrue(np.allclose(Ex[4:, 4:, :], 1.0))

    def test_conducting_bc_behind_wall_equivalent(self):
        """Edge BC reproduces the old reflector_wall_field_bc behavior."""
        _, tile = make_tile()

        nonzero = lambda x, y, z: (1.0, 2.0, 3.0)
        tile.set_EBJ(nonzero, zero, zero)

        (Ex0, _, _), _, _ = tile.get_EBJ()

        # Replicate: zero Ey, Ez behind wall at x=5.0 using edge_bc
        bc = runko.emf.threeD.edge_bc(direction=0, side=0, position=5.0,
                                       E_components=0b110)
        tile.apply_edge_bc(bc, runko.tools.comm_mode.emf_E)

        (Ex, Ey, Ez), _, _ = tile.get_EBJ()

        iw = 6  # first non-halo index in front of wall
        self.assertTrue(np.all(Ey[:iw, :, :] == 0))
        self.assertTrue(np.all(Ez[:iw, :, :] == 0))
        self.assertTrue(np.all(Ey[iw:, :, :] != 0))
        self.assertTrue(np.all(Ez[iw:, :, :] != 0))
        self.assertTrue(np.allclose(Ex, Ex0))


if __name__ == "__main__":
    unittest.main()

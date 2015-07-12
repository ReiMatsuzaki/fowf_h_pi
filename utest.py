import unittest
import fowf
import fowf_h_pi

class TestGrid(unittest.TestCase):
    def setUp(self):
        self.grid = fowf.Grid(0.1, 100)

    def test_accessor(self):
        self.assertEqual(100, self.grid.num())
        self.assertEqual(0.1, self.grid.h())
        self.assertEqual(0.1, self.grid.xs_i(0))
        self.assertAlmostEqual(0.3, self.grid.xs_i(2))

    def test_copy(self):
        g2 = self.grid
        self.assertEqual(100, g2.num())


class TestHydrogenPI(unittest.TestCase):
    def setUp(self):
        self.h_pi = fowf.HydrogenPI(0.43)

    def test_accessor(self):
        self.assertEqual(      0,    self.h_pi.l0())
        self.assertEqual(      1,    self.h_pi.l1())
        self.assertAlmostEqual(0.43, self.h_pi.ene())


class TestDiscreteFunc(unittest.TestCase):
    def setUp(self):
        g = fowf.Grid(0.1, 100)
        self.func = fowf.discretize_func([1,2,3,4], g)
        
    def test_accessor(self):
        self.assertAlmostEqual(0.1, self.func.xs_i(0))
        self.assertAlmostEqual(0.3, self.func.xs_i(2))
        self.assertAlmostEqual(1.0, self.func.ys_i(0))
        self.assertAlmostEqual(2.0, self.func.ys_i(1))

class TestMain(unittest.TestCase):
    def setUp(self):
        grid = fowf.Grid(0.1, 1000)
        h_pi = fowf.HydrogenPI(0.5)
        driv_sys = fowf.DrivSys(grid, h_pi) 
        self.wave_func = driv_sys.solve()
        self.res_data= fowf_h_pi.calculate_fowf(grid, h_pi)
    def test_print(self):
        print self.wave_func.xs_i(0)
        print self.wave_func.ys_i(0)
    def test_np(self):
        res = fowf_h_pi.DiscreteFunc_to_ndarray(self.wave_func)
        self.assertAlmostEqual(res[1][1], self.res_data[1][1])
        print self.wave_func.ys_i(0)
    def test_print2(self):
        grid = fowf.Grid(0.1, 1000)
        h_pi = fowf.HydrogenPI(0.5)
        wf   = fowf_h_pi.calculate_fowf(grid, h_pi)
        print "/////////////"
        print wf[1]




if __name__ == '__main__':
    unittest.main()
    


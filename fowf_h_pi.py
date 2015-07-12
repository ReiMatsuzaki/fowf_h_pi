import fowf
import numpy as np

def DiscreteFunc_to_ndarray(f):

    num = f.size()
    def x_yr_yi(i):
        x = f.xs_i(i)
        y = f.ys_i(i)
        return [x, y.real, y.imag]

    x_yr_yi_list = [[0.0, 0.0, 0.0]] + [ x_yr_yi(i) for i in range(num) ]
    res = np.array(x_yr_yi_list)
    return res

def calculate_fowf(g, h_pi):
    driv_sys = fowf.DrivSys(g, h_pi)
    func = driv_sys.solve()
    return DiscreteFunc_to_ndarray(func)



# Copyright of the Board of Trustees of Columbia University in the City of New York

'''
Testing for spiral-fMRI repository
Author: Marina Manso Jimeno
Last Updated: 05/11/2020
'''

import unittest
import numpy as np
import scipy.io as sio
import math
import matplotlib.pyplot as plt

from Acq.make_slr_rf import make_slr_rf
from Acq.trapwave2 import trapwave2
from Acq.make_crusher import make_crusher
from Acq.makeSystemlength import makeSystemlength
from Acq.interpolation import rf_interpolate, grad_interpolate
from Acq.archimedian import archimedian_spiral
from Acq.make_balanced import make_balanced

from pypulseq.opts import Opts

system = {'gamma': 42576000, 'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s',
                 'grad_raster_time': 10e-6, 'rf_raster_time': 10e-6}

class testMakeSLRPulse(unittest.TestCase):

    def test_flip_angle_range(self):
        with self.assertRaises(ValueError): make_slr_rf(flip_angle=-10,slice_thickness=10,time_bw_product=1.5,duration=1e-3,ncycles=0,system=system)
        with self.assertRaises(ValueError): make_slr_rf(flip_angle=500, slice_thickness=10, time_bw_product=1.5,
                                                  duration=1e-3, ncycles=0, system=system)

    def test_flip_angle(self):
        flip_in = 90
        rf, _ = make_slr_rf(flip_angle=flip_in, slice_thickness=10, time_bw_product=1.5,
                    duration=1e-3, ncycles=0, system=system)
        # TODO: Is there a way to calculate the flip angle from the RF waveform?

    def test_tbw_product_val(self):
        with self.assertRaises(ValueError): make_slr_rf(flip_angle=90, slice_thickness=10, time_bw_product=1,
                                                  duration=1e-3, ncycles=0, system=system)

    def test_tbw_product(self):
        rf, _ = make_slr_rf(flip_angle=90, slice_thickness=10, time_bw_product=8,
                            duration=1e-3, ncycles=0, system=system)
        # TODO: Is there a way to calculate the bandwidth from the RF waveform?
        '''dur1 = rf.shape[1] * 10e-6
        t = rf * 10e-6
        BW = 0.2 # 1 / 1e-3
        tbw = dur1 * BW
        freq = np.fft.fftfreq(t.shape[-1])
        plt.plot(np.flipud(freq), np.fft.fft(rf[0]))
        plt.show()'''

    def test_duration_val(self):
        with self.assertRaises(ValueError): make_slr_rf(flip_angle=90, slice_thickness=10, time_bw_product=8,
                                                  duration=0, ncycles=0, system=system)

    def test_duration(self):
        dur_in = 1
        rf, _ = make_slr_rf(flip_angle=90, slice_thickness=10, time_bw_product=8,
                            duration=dur_in, ncycles=0, system=system)
        dur_out = rf.shape[0] * system['rf_raster_time']
        self.assertEqual(dur_in, np.round(dur_out,3))

    def test_slice_thickness_val(self):
        with self.assertRaises(ValueError): make_slr_rf(flip_angle=90, slice_thickness=-10, time_bw_product=8,
                                                  duration=1e-3, ncycles=0, system=system)

    def test_slice_thickness(self):
        sl_thick_in = 10e-3
        tbwp = 8
        dur = 1e-3
        rf, gss = make_slr_rf(flip_angle=90, slice_thickness=sl_thick_in, time_bw_product=tbwp,
                            duration=dur, ncycles=0, system=system)

        gplat = gss.max() # Hz/m
        sl_thick_out = tbwp / dur / gplat
        self.assertEqual(sl_thick_in, np.round(sl_thick_out,3))

    def test_makeslrpulse_matlabvspython(self):
        rf_mat = sio.loadmat('testing_assets/rf_slr.mat')['rf'] * system['gamma'] / 10000  # from G to Hz
        gss_mat = sio.loadmat('testing_assets/gss_slr.mat')['gex'] * system['gamma'] / 100 # from G/cm to Hz/m
        systemmat = {'gamma': 42576000, 'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s',
                     'grad_raster_time': 4e-6, 'rf_raster_time': 4e-6}
        rf_pyth, gss_pyth = make_slr_rf(flip_angle=10, slice_thickness=0.144, time_bw_product=8, duration=1e-3, ncycles=0, system=systemmat)
        self.assertAlmostEqual(rf_mat.all(),rf_pyth.all())
        self.assertEqual(gss_mat.all(),gss_pyth.all())

class testTrap2Wave(unittest.TestCase):
    def test_area_val(self):
        with self.assertRaises(ValueError): trapwave2(area=0,mxg=system['max_grad'],mxs=system['max_slew'],rasterTime=system['grad_raster_time'])

    def test_area(self):
        area_in = 10
        gtrap = trapwave2(area=area_in,mxg=system['max_grad'],mxs=system['max_slew'],rasterTime=system['grad_raster_time'])

        base1 = gtrap.shape[1]
        base2 = len(np.where(gtrap == gtrap.max())[1])
        height = gtrap.max()
        area_out = (base1 + base2) * system['grad_raster_time'] / 2 * height # Hz/m*sec
        self.assertEqual(area_in, round(area_out))

    def test_maxGrad(self):
        maxGrad = 50
        gtrap = trapwave2(area=1, mxg=maxGrad, mxs=system['max_slew'],
                          rasterTime=system['grad_raster_time'])
        self.assertLess(gtrap.max(), maxGrad)

    def test_maxSlew(self):
        maxSlew = 150
        gtrap = trapwave2(area=1, mxg=system['max_grad'], mxs=maxSlew,
                          rasterTime=system['grad_raster_time'])
        slew = np.diff(gtrap) / system['grad_raster_time'] * 1e-3
        self.assertLess(slew.max(), maxSlew)

    def test_rasterTime_val(self):
        with self.assertRaises(ValueError): trapwave2(area=1, mxg=system['max_grad'], mxs=system['max_slew'],
                          rasterTime=0)

class testMakeCrusher(unittest.TestCase):
    def test_ncycles_val(self):
        with self.assertRaises(ValueError): make_crusher(ncycles=0,opslthick=1e-3,gzarea=0,max_grad=system['max_grad'],max_slew=system['max_slew'],sys_lims=system)

    def test_opsthick_val(self):
        with self.assertRaises(ValueError): make_crusher(ncycles=1, opslthick=0, gzarea=0,
                                                         max_grad=system['max_grad'], max_slew=system['max_slew'],                                               sys_lims=system)
    def test_missing_gamma(self):
        sys_lims = {'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s',
                 'grad_raster_time': 10e-6, 'rf_raster_time': 10e-6}
        with self.assertRaises(ValueError): make_crusher(ncycles=1, opslthick=1e-3, gzarea=0,
                                                         max_grad=system['max_grad'], max_slew=system['max_slew'],
                                                sys_lims=sys_lims)
    def test_makecrusher_matlabvspython(self):
        gcrush_mat = sio.loadmat('testing_assets/gcrusher.mat')['gspoil1'] * system['gamma'] / 100 # G/cm to Hz/m
        systemmat = {'gamma': 42576000, 'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s',
                     'grad_raster_time': 4e-6, 'rf_raster_time': 4e-6}
        gcrush_pyth = make_crusher(ncycles=1, opslthick=0.3333e-2, gzarea=0, max_grad=32, max_slew=0.9*130/math.sqrt(2),sys_lims=systemmat)
        self.assertEqual(gcrush_mat.all(), gcrush_pyth.all())

class testMakeSystemLength(unittest.TestCase):
    def test_rater_val(self):
        with self.assertRaises(ValueError): makeSystemlength(g=np.ones((50,1)),raster=0)
        with self.assertRaises(ValueError): makeSystemlength(g=np.ones((50, 1)), raster=-1)

    def test_g_shape_out(self):
        g1 = makeSystemlength(g=np.ones((1,50)),raster=4e-6)
        self.assertEqual(g1.shape[1], 1)
        g2 = makeSystemlength(g=np.ones((50, 1)), raster=4e-6)
        self.assertEqual(g2.shape[1], 1)
        g4 = makeSystemlength(g=np.ones((50,)), raster=4e-6)
        self.assertEqual(g4.shape[1], 1)

    def test_g_len_out(self):
        gGE = makeSystemlength(g=np.ones((43, 1)), raster=4e-6)
        self.assertEqual(len(gGE) % 4, 0)
        gSiemens = makeSystemlength(g=np.ones((43, 1)), raster=1e-6)
        self.assertEqual(len(gSiemens) % 2, 0)

    def test_remove_negzero(self):
        g = makeSystemlength(g= - np.zeros((4,1)), raster=4e-6)
        self.assertNotEqual(g.all(), (-np.zeros((4,1))).all)

class testInterpolation(unittest.TestCase):
    def test_rfinterp_rasterval(self):
        with self.assertRaises(ValueError): rf_interpolate(rf_signal=np.ones((20,1)),rf_raster_time=0)
        with self.assertRaises(ValueError): rf_interpolate(rf_signal=np.ones((20, 1)), rf_raster_time=-1)

    def test_rfinterp_signalout_len(self):
        rf_in = np.ones((20,1))
        rf_largeRaster = rf_interpolate(rf_signal=np.ones((20,1)),rf_raster_time=10e-6)
        self.assertLess(len(rf_largeRaster),len(rf_in))
        rf_smallRaster = rf_interpolate(rf_signal=np.ones((20,1)))
        self.assertGreater(len(rf_smallRaster), len(rf_in))

    def test_gradinterp_rasterval(self):
        with self.assertRaises(ValueError): grad_interpolate(grad_signal=np.ones((20,1)),grad_raster_time=0)
        with self.assertRaises(ValueError): grad_interpolate(grad_signal=np.ones((20, 1)), grad_raster_time=-1)

    def test_gradinterp_signalout_len(self):
        grad_in = np.ones((20,1))
        grad_largeRaster = grad_interpolate(grad_signal=grad_in)
        self.assertLess(len(grad_largeRaster),len(grad_in))
        grad_smallRaster = grad_interpolate(grad_signal=grad_in,grad_raster_time=2e-6)
        self.assertGreater(len(grad_smallRaster),len(grad_in))

class testArchimedianSpiral(unittest.TestCase):
    def test_rmax(self):
        rmax_in = 150.0
        k, _ = archimedian_spiral(smax=system['max_slew'],gmax=system['max_grad'],T=system['grad_raster_time'],N=3,FOV=0.24,rmax=rmax_in)
        self.assertLess(np.real(k.max()), rmax_in)

    def test_max_grad(self):
        gmax_in = system['max_grad']
        _, g = archimedian_spiral(smax=system['max_slew'], gmax=gmax_in, T=system['grad_raster_time'], N=3,
                                  FOV=0.24, rmax=150)
        gmax_out = np.real(g.max()) / system['gamma'] * 1e3 # mT/m
        self.assertLess(gmax_out, gmax_in)

    def test_max_slew(self):
        smax_in = system['max_slew']
        T_in = system['grad_raster_time']
        _, g = archimedian_spiral(smax=smax_in, gmax=system['max_grad'], T=T_in, N=3,
                                  FOV=0.24, rmax=150)
        s_out = np.diff(g) / system['grad_raster_time']
        smax_out = np.real(s_out.max()) / system['gamma']  # mT/m
        self.assertLess(smax_out, smax_in)

    def test_N_val(self):
        with self.assertRaises(ValueError) : archimedian_spiral(smax=system['max_slew'], gmax=system['max_grad'], T=system['grad_raster_time'], N=0,
                                  FOV=0.24, rmax=150)

    def test_rmax_val_zero(self):
        with self.assertRaises(ValueError): archimedian_spiral(smax=system['max_slew'], gmax=system['max_grad'], T=system['grad_raster_time'], N=3,
                                  FOV=0.24, rmax=0)

    def test_rmax_val_neg(self):
        kneg, gneg = archimedian_spiral(smax=system['max_slew'], gmax=system['max_grad'], T=system['grad_raster_time'], N=3,
                                  FOV=0.24, rmax=-150)
        kpos, gpos = archimedian_spiral(smax=system['max_slew'], gmax=system['max_grad'], T=system['grad_raster_time'],
                                        N=3,FOV=0.24, rmax=150)
        self.assertEqual(kpos.all(), kneg.all())
        self.assertEqual(gpos.all(), gneg.all())

    def test_FOV_val(self):
        with self.assertRaises(ValueError): archimedian_spiral(smax=system['max_slew'], gmax=system['max_grad'], T=system['grad_raster_time'],
                                        N=3, FOV=0, rmax=150)
        with self.assertRaises(ValueError): archimedian_spiral(smax=system['max_slew'], gmax=system['max_grad'], T=system['grad_raster_time'],
                                        N=3, FOV=-0.10, rmax=150)

    def test_archspiral_matlabvspython(self):
        g_mat = sio.loadmat('testing_assets/arch_spiral.mat')['g'] * system['gamma'] / 100  # from G/cm to Hz/m
        _, g_pyth = archimedian_spiral(smax=0.99 * 130, gmax=0.99 * 32, T=4e-6, N=3, FOV=0.24, rmax=150)
        self.assertEqual(g_mat.all(),g_pyth.all())

class testMakeBalanced(unittest.TestCase):
    def test_g_shape(self):
        with self.assertRaises(ValueError): make_balanced(np.ones((4,500)), system['max_grad'], system['max_slew'], system)
        g1 = make_balanced(np.ones((1,500)),system['max_grad'], system['max_slew'], system)
        g2 = make_balanced(np.ones((500,)), system['max_grad'], system['max_slew'], system)
        g3 = make_balanced(np.ones((500,1)), system['max_grad'], system['max_slew'], system)
        self.assertEqual(g1.all(), g2.all(), g3.all())

    def test_missingelement_missing(self):
        sys1 = {'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s',
                 'grad_raster_time': 10e-6, 'rf_raster_time': 10e-6}
        sys2 = {'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s', 'rf_raster_time': 10e-6}
        with self.assertRaises(ValueError): make_balanced(np.ones((50,1)),system['max_grad'], system['max_slew'], sys1)
        with self.assertRaises(ValueError): make_balanced(np.ones((50, 1)),system['max_grad'], system['max_slew'], sys2)

    def test_g_allzeros(self):
        with self.assertRaises(ValueError): make_balanced(g=np.zeros((50,1)), max_grad=system['max_grad'], max_slew=system['max_slew'], system=system)

    def test_makebalanced_matlabvspython(self):
        gbal_mat = sio.loadmat('testing_assets/arch_spiral_gxbalanced.mat')['gx'] * system['gamma'] / 100  # from G/cm to Hz/m
        _, g_pyth = archimedian_spiral(smax=0.99 * 130, gmax=0.99 * 32, T=4e-6, N=3, FOV=0.24, rmax=150)
        gx = np.real(np.vstack((0,0,g_pyth.reshape(len(g_pyth),1))))
        systemmat = {'gamma': 42576000, 'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s',
                 'grad_raster_time': 4e-6}
        gbal_pyth = make_balanced(g=gx, max_grad=32, max_slew=130/math.sqrt(2), system=systemmat)
        self.assertEqual(gbal_mat.all(),gbal_pyth.all())
if __name__ == "__main__":
    unittest.main()
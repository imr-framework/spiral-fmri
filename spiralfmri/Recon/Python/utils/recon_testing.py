'''
Testing for spiral-fMRI repository (Image recon)
Author: Marina Manso Jimeno
Last Updated: 06/10/2020
'''
import configparser
import scipy.io as sio
import unittest
import numpy as np

from spiralfmri.Recon.Python.utils.recon_opts import data_check_file_inputs, data_check_getparams_input, rfphs_compensate, spiral_alignment
from spiralfmri.Recon.Python.recon_SoS import recon_SoS

path_settings = '../../../settings.ini'
config = configparser.ConfigParser()
config.read(path_settings)
config_data = config['RECON TEST']

dat = np.load(config_data['path_rawdat_file'])
ktraj = np.load(config_data['path_ktraj_file'])
dat_4recon = np.load(config_data['path_rawdat_4recon'])
ktraj_4recon = np.load(config_data['path_ktraj_4recon'])

class TestReconOpts(unittest.TestCase):
    def test_input_dims(self):
        dat1 = dat
        ktraj1 = ktraj[:100, :, :, :]
        with self.assertRaises(ValueError): data_check_file_inputs(raw_dat=dat1, ktraj=ktraj1)

        dat2 = dat
        ktraj2 = ktraj[:, :50, :, :]
        with self.assertRaises(ValueError): data_check_file_inputs(raw_dat=dat2, ktraj=ktraj2)

        dat3 = dat
        ktraj3 = ktraj[:, :, :18, :]
        with self.assertRaises(ValueError): data_check_file_inputs(raw_dat=dat3, ktraj=ktraj3)

    def test_input_getparams_test(self):
        params = {'acq_params': {'nz': 54}, 'fmri': {'nframes': 21}}
        dat1 = dat[:, :50, :, :]
        with self.assertRaises(ValueError): data_check_getparams_input(raw_dat=dat1, params=params)

        dat2 = dat[:, :, :19, :]
        with self.assertRaises(ValueError): data_check_getparams_input(raw_dat=dat2, params=params)

    def test_rfphs_compensate(self):
        dat1 = dat
        angle = 20
        corr_dat = rfphs_compensate(raw_dat=dat1, rfphs_angle=angle)
        self.assertEqual(dat.shape, corr_dat.shape)
        self.assertEqual(corr_dat.dtype, 'complex128')

    def test_spiral_alignment(self):
        dat1 = dat
        nleafs = 4
        dat_aligned = spiral_alignment(dat1, nleafs)
        self.assertEqual(dat_aligned.shape[0], dat1.shape[0])
        self.assertEqual(dat_aligned.shape[1], nleafs)
        self.assertEqual(dat_aligned.shape[2],dat1.shape[1])
        self.assertEqual(dat_aligned.shape[3], dat1.shape[2] - nleafs + 1)
        self.assertEqual(dat_aligned.shape[-1], dat1.shape[-1])

class TestReconSoS(unittest.TestCase):
    def test_reconsos_output(self):
        dat1 = np.expand_dims(dat_4recon[:, :, 0, 0, :], axis=2)
        ktraj1 = ktraj_4recon
        N = 72
        recon = recon_SoS(dat1, ktraj1, N)
        self.assertEqual(recon.shape[0],recon.shape[1],N)
        self.assertEqual(recon.shape[-1], dat1.shape[2])
        self.assertGreaterEqual(recon.min(), 0)


if __name__ == "__main__":
    unittest.main()
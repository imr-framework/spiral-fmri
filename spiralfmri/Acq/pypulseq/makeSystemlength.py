import numpy as np

def makeSystemlength(g, raster):
    '''
    Ensures length is divisible by 4 for GE and by 2 for Siemens by padding with zeros if necessary

    Parameters
    ----------
    g : numpy.ndarray
        gradient/RF waveform
    raster : float
        gradient/RF raster time

    Returns
    -------
    g : numpy.ndarray
        Waveform with the proper length
    '''
    if raster <= 0:
        raise ValueError('This is not a valid raster time')

    if raster == 4e-6:
        sampling = int(raster * 1e6) # GE
    else:
        sampling = 2 # Siemens

    if len(g.shape) == 1:
        g = g.reshape(len(g),1)

    if g.shape[0] < g.shape[1]:
        g = g.reshape(-1,1)

    if g.shape[0] % sampling:
        nrows, ncols = g.shape
        g = np.vstack((g, np.zeros((sampling - (g.shape[0] % sampling), ncols))))
        # g = np.concatenate((g, np.zeros((nrows, sampling-(g.shape[1] % sampling)))), axis=1)

    g[np.where(g == -0)] = 0

    return g
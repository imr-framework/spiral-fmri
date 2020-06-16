'''
Helper functions for spiral-fmri reconstruction
Author : Marina Manso Jimeno
Last updated : 06/08/2020
'''
import math
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import ImageGrid

def disp_timepoint(time_point):
    '''
    Displays in a grid the images corresponding to one time point

    Parameters
    ----------
    time_point : numpy.ndarray
        Array of images corresponding to one time-point. Dimensions (matrix size(N x N) x Nslices)
    '''
    N = time_point.shape[0]
    nslices = time_point.shape[-1]
    sq = math.sqrt(nslices)
    if sq == round(sq):
        ncols = sq
        nrows = sq
    else:
        nrows = round(sq)
        left = nslices - nrows ** 2
        ncols = nrows + math.ceil(left / nrows)
        extra = nrows * ncols - nslices

    im_list = np.moveaxis(time_point, -1, 0)
    if extra:
        im_list = np.vstack((im_list, np.zeros((extra, N, N))))


    fig = plt.figure()
    grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(nrows, ncols),  # creates 2x2 grid of axes
                     axes_pad=0,  # pad between axes in inch.
                     )

    for ax, im in zip(grid, [im for im in im_list]):
        # Iterating over the grid returns the Axes.
        ax.imshow(im, cmap='gray', vmin=im_list.min(), vmax=im_list.max())
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)


    plt.show()


"""
The code here takes a scene image and convolves with the NIRISS WFSS "PSF"
image to produce a simulated dispersed scene.
"""
import numpy
import astropy.io.fits as fits
import scipy.signal as signal
from numpy.fft import rfftn, irfftn, fftshift, ifftshift


def wfss_scene(scene_image, filtername, grismname, x0, y0, path='./'):
    """
    Convolve a scene image with the WFSS PSF and return dispersed image over
    the 2322x2322 pixel POM image area.

    Parameters
    ----------

    scene_image:  A numpy 2-d image (float) of an imaging scene to disperse. 
                  must be 2322x2322 pixels or larger

    filtername:   A string, one of the WFSS blocking filter names

    grimsname:    A string, the NIRISS GR150 grism name, either 'GR150R' or
                  'GR150C' 

    x0:           An integer value, the lower left corner x pixel value for 
                  the POM image read-out area

    y0:           An integer value, the lower left corner y pixel value for 
                  the POM image read-out area

    path:         An optional string value, the path to the WFSS PSF images

    Returns
    -------

    outimage:     A numpy 2-d image (float) of the dispersed scene; size 
                  2322x2322 pixels, or None if there is an issue

    The scene image is multiplied by the spot mask before the convolution in 
    the area that corresponds to the output pixels.
    """
    grisms = ['GR150R', 'GR150C']
    filters = ['F090W', 'F115W', 'F140M', 'F150W', 'F158M', 'F200W']
    if not grismname.upper() in grisms:
        print('Error: bad grism name %s passed to wfss_scene' % (grismname))
        return None
    if not filtername.upper() in filters:
        print('Error: bad grism name %s passed to wfss_scene' % (filtername))
        return None
    imshape = scene_image.shape
    if (x0 < 0) or (y0 < 0) or (x0+2322 > imshape[1]) or (y0+2322 > imshape[0]):
        print('Error in wfss_scene: bad image offset values' + \
              '(%d, %d) passed to the routine.' % (x0, y0))
        return None
    try:
        if path[-1] != '/':
            path = path+'/'
        spotmask = fits.getdata(path+'occulting_spots_mask.fits')
    except:
        print('Error: the occulting spot mask was not found.')
        spotmask = numpy.zeros((2048, 2048), dtype=numpy.float32)+1.
    psfname = filtername+'_'+grismname+'_psfimage.fits'
    psfname = psfname.lower()
    try:
        psfimage = fits.getdata(path+psfname)
    except:
        print('Error: PSF image %s not found in directory %s.' % (psfname,
                                                                  path))
        return None
    field_image = numpy.copy(scene_image[y0:y0+2322, x0:x0+2322])
    y1 = y0+137
    x1 = x0+137
    field_image[y1:y1+2048, x1:x1+2048] = \
        field_image[y1:y1+2048,x1:x1+2048]*spotmask
    newimage = signal.fftconvolve(field_image, psfimage, mode='same')
    
    # scale by the grism total throughput of 0.8
    newimage = newimage*0.8
    return newimage

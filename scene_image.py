"""
Code to make a scene image from input files for stars and galaxies.

Given a pointing and suitable catalogue files, one for stars and one for 
galaxies as in with Mirage, generate a scene image in the standard orientation.

Actually two scene images are made, one with the stars as true point sources 
and one where the first image is convolved with the PSF function.

Routines in this file

make_star_image:  Make a star scene image from a Mirage star list file, with 
                  the PSF assumed to be a single pixel (for later convolution 
                  with the real PSF)

make_galaxy_image:  Make a galaxy scene image from a Mirage galaxy list file, 
                    based on the Sersic parameters 

do_convolve:  Convolve the ideal scene image with a suitable NIRISS PSF image 
              to make the proper scene image

rotate_image:  Use the scipy.ndimage.rotate function to rotate a scene image

generate_image: Make an ideal star image from a list of stellar positions and 
                total signal values

get_pixel:   Calculate the pixel position of a given (RA, Dec) sky position 
             with respect to a given image reference position (RA0, Dec0), 
             using pysiaf.

rel_pos:   Calculate the pixel position of a given (RA, Dec) sky position 
           with respect to an given image reference position (RA0, Dec0) 
           using direct geometry.

"""
import sys
import math
import numpy
import pysiaf
import astropy.io.fits as fits
import scipy.signal as signal
import scipy.ndimage as ndimage
from astropy.modeling.models import Sersic2D

def make_star_image(star_file_name, position, filter1, path='./',
                    simple=False):
    """
    Make the star scene image from an input file of positions/brightnesses, 
    each star on a single pixel.

    Parameters
    ----------

    star_file_name:  a string variable giving the filename for the list of 
                     stars to use, with (RA, Dec) sky positions and the 
                     NIRISS magnitudes

    position:        a two-element float tuple with the (RA, Dec) values in 
                     decimal degrees for the centre of the field

    filter1:         a string variable giving the NIRISS filter name for which 
                     to make the scene image

    path:            an optional string variable pointing to where the PSF 
                     images are stored, defaults to the current directory

    simple:          an optional Boolean variable, if True use the simple 
                     sky->pixel calculation, if False, the default, use 
                     pysiaf

    Returns
    -------

    scene_image:    a numpy float array of the brightness values (ADU/s) for 
                    the stars in the scene.  Stars are true point sources in 
                    this image, occupying a single pixel.  The image is 3631 
                    pixels square with pixel (1816, 1816) being the reference 
                    position.  The image is undistorted.  Each star is put at 
                    the nearest pixel position.  The image is in the standard
                    orientation.  The NIRISS detector field is then pixels 
                    [792:2840, 792:2840] by assumption.  The larger POM field 
                    is assumed to cover pixels [655:2976, 655:2976].

                    None is returned if there is an issue.

    star_list:      A tuple of the [ra, dec, signal] values found within the 
                    scene, each element being a numpy array; if no stars are 
                    found in the field then the values are None.
    """
    blank = [None, None, None]
    mag0 = [1.243877e+11, 1.041117e+11, 3.256208e+10, 6.172448e+10,
            2.877868e+10, 4.245261e+10, 2.472333e+10, 1.626810e+10,
            2.930977e+09, 1.902929e+09, 9.948760e+09, 1.703311e+09]
    aboff = [0.48790, 0.74678, 1.07748, 1.17275, 1.27628, 1.65571,
             2.25786, 2.77045, 2.91546, 3.14410, 3.19330, 3.37689]
    fnames = ['F090W', 'F115W', 'F140M', 'F150W', 'F158M', 'F200W',
              'F277W', 'F356W', 'F380M', 'F430M', 'F444W', 'F480M']
    mag0 = numpy.asarray(mag0)
    aboff = numpy.asarray(aboff)
    aboff = numpy.power(10., aboff*0.4)
    mag0 = mag0/1.612
    if not filter1.upper() in fnames:
        print('Filter name not recognized.')
        return None, blank
    for loop in range(len(fnames)):
        if filter1 == fnames[loop]:
            findex = loop
    if findex > 5:
        print('Filter %d is not a WFSS blocking filter.')
        return None, blank
    try:
        infile = open(star_file_name, 'r')
        lines = infile.readlines()
        infile.close()
        abmag_flag = False
        for line in lines:
            if ('pixel' in line) and ('#' in line):
                print('Error: star positions must be in (RA, Dec) form.')
                return None, blank
            if ('abmag' in line) and ('#' in line):
                abmag_flag = True
            if 'x_or_RA' in line:
                racol =-1
                deccol = -1
                magcol = -1
                target = 'niriss_' + filter1.lower() + '_magnitude'
                line = line.strip('\n')
                values = line.split()
                for loop in range(len(values)):
                    if 'x_or_RA' in values[loop]:
                        racol = loop
                    if 'y_or_Dec' in values[loop]:
                        deccol = loop
                    if target in values[loop]:
                        magcol = loop
                if (racol < 0) or (deccol < 0) or (magcol < 0):
                    print('Unable to parse columns in star file.')
                    return None, blank
                try:
                    ravalues = numpy.loadtxt(star_file_name, usecols=[racol,],
                                             comments=['#', 'index'])
                    decvalues = numpy.loadtxt(star_file_name, usecols=[deccol,],
                                              comments=['#', 'index'])
                    magvalues = numpy.loadtxt(star_file_name, usecols=[magcol,],
                                              comments=['#', 'index'])
                    break
                except:
                    print('Unable to read the position/magnitude values.')
                    return None, blank
        lines = 0
        signal = magvalues*0.
        for loop in range(len(magvalues)):
            signal[loop] = mag0[findex]/(10.**(magvalues[loop]*0.4))
            if abmag_flag:
                signal[loop] = signal[loop]*aboff[findex]
        star_list = [ravalues, decvalues, signal]
        scene_image, new_star_list = generate_image(star_list, position,
                                                    simple=simple)
        return scene_image, new_star_list
    except:
        print('An error occurred generating the scene image.')
        return None, blank

def make_galaxy_image(galaxy_file_name, position, filter1, path='./', 
                      simple=False):
    """
    Make the galaxy scene image from an input parameter file, in standard 
    orientation (not convolved with a PSF)

    Parameters
    ----------

    galaxy_file_name:  a string variable giving the name for the file of 
                       galaxy parameters to use, with (RA, Dec) sky positions, 
                       shape parameters, and the NIRISS magnitudes

    postion:         a two-element real tuple with the (RA, Dec) values in 
                     decimal degrees for the centre of the field

    filter1:         a string variable giving the NIRISS filter name for which 
                     to make the scene image

    path:            an optional string variable pointing to where the PSF 
                     images are stored, defaults to the current directory

    simple:          an optional Boolean variable, if True use the simple 
                     sky->pixel calculation, if False, the default, use 
                     pysiaf

    Returns
    -------

    scene_image:    a numpy float array of the brightness values (ADU/s) for 
                    the galaxies in the scene.  The image is 3631 pixels square 
                    with pixel (1816, 1816) being the reference position.  The 
                    image is undistorted.  Each galaxy is put at the nearest 
                    pixel position.  The NIRISS detector field is at pixels  
                    [792:2840, 792:2840] by assumption.  The larger POM field 
                    is assumed to cover pixels [655:2976, 655:2976].

                    None is returned if there is an issue.

    """
    instrument = 'NIRISS'
    aperture = 'NIS_CEN'
    mag0 = [1.243877e+11, 1.041117e+11, 3.256208e+10, 6.172448e+10,
            2.877868e+10, 4.245261e+10, 2.472333e+10, 1.626810e+10,
            2.930977e+09, 1.902929e+09, 9.948760e+09, 1.703311e+09]
    aboff = [0.48790, 0.74678, 1.07748, 1.17275, 1.27628, 1.65571,
             2.25786, 2.77045, 2.91546, 3.14410, 3.19330, 3.37689]
    fnames = ['F090W', 'F115W', 'F140M', 'F150W', 'F158M', 'F200W',
              'F277W', 'F356W', 'F380M', 'F430M', 'F444W', 'F480M']
    mag0 = numpy.asarray(mag0)
    aboff = numpy.asarray(aboff)
    aboff = numpy.power(10., aboff*0.4)
    mag0 = mag0/1.612
    if not filter1.upper() in fnames:
        print('Filter name (%s) not recognized.' % (filter1))
        return None
    for loop in range(len(fnames)):
        if filter1 == fnames[loop]:
            findex = loop
            break
    if findex > 5:
        print('Filter %d is not a WFSS blocking filter.')
        return None
    try:
        infile = open(galaxy_file_name, 'r')
        lines = infile.readlines()
        infile.close()
        abmag_flag = False
        for line in lines:
            if ('pixel' in line) and ('#' in line):
                print('Error: galaxy positions must be in (RA, Dec) form.')
                return None
            if ('abmag' in line) and ('#' in line):
                abmag_flag = True
            if 'x_or_RA' in line:
                racol =-1
                deccol = -1
                pacol = -1
                radcol = -1
                indexcol = -1
                ellipcol = -1
                magcol = -1
                target = 'niriss_' + filter1.lower() + '_magnitude'
                line = line.strip('\n')
                values = line.split()
                for loop in range(len(values)):
                    if 'x_or_RA' in values[loop]:
                        racol = loop
                    if 'y_or_Dec' in values[loop]:
                        deccol = loop
                    if target in values[loop]:
                        magcol = loop
                    if 'ellipticity' in values[loop]:
                        ellipcol = loop
                    if 'sersic_index' in values[loop]:
                        indexcol = loop
                    if 'pos_angle' in values[loop]:
                        pacol = loop
                    if 'radius' in values[loop]:
                        radcol = loop
                cols = numpy.asarray([racol, deccol, radcol, ellipcol, pacol,
                                     indexcol, magcol], dtype=numpy.int16)
                if numpy.min(cols) < 0:
                    print('Unable to parse columns in galaxies file.')
                    return None
        try:
            ravalues = numpy.loadtxt(galaxy_file_name, usecols=[racol,],
                                     comments=['#', 'index'])
            decvalues = numpy.loadtxt(galaxy_file_name, usecols=[deccol,],
                                      comments=['#', 'index'])
            magvalues = numpy.loadtxt(galaxy_file_name, usecols=[magcol,],
                                      comments=['#', 'index'])
            radvalues = numpy.loadtxt(galaxy_file_name, usecols=[radcol,],
                                      comments=['#', 'index'])
            ellipvalues = numpy.loadtxt(galaxy_file_name, usecols=[ellipcol,],
                                        comments=['#', 'index'])
            pavalues = numpy.loadtxt(galaxy_file_name, usecols=[pacol,],
                                     comments=['#', 'index'])
            indexvalues = numpy.loadtxt(galaxy_file_name, usecols=[indexcol,],
                                        comments=['#', 'index'])
            signal = mag0[findex]/numpy.power(10., magvalues*0.4)
            signal = signal/1.612
            if abmag_flag:
                signal = signal*aboff[findex]
        except:
            print('Unable to read the position/magnitude/shape values.')
            return None
        galaxy_list = [ravalues, decvalues, magvalues, signal, 
                       radvalues, ellipvalues, pavalues, indexvalues]
        gimage = generate_galaxy_image(galaxy_list, position, simple=simple)
        return gimage
    except:
#        print('Some error occured trying to read the galaxies parameter file.')
        return None

def generate_galaxy_image(galaxy_list, position, rotation=0., simple=False):
    """
    Do the work of making a galaxies scene image.

    Parameters
    ----------

    galaxy_list:   A list of numpy 1-d float arrays holding the catalogue 
                   values for the scene (RA, DEC, mag, signal, radii, 
                   ellipticity, angle, sersic index)

    position:      A two-element list giving the image center (RA, Dec) in 
                   degrees 

    rotation:      An optional float value, the rotation angle in degrees E 
                   of N

    simple:        A boolean value, if True use the simple projection to get 
                   pixel positions, if False use pysiaf.  The latter is the 
                   default.

    Returns
    -------

    galaxy_image:   A 3631x3631 numpy float array, the scene image.
    """
    instrument = 'NIRISS'
    aperture = 'NIS_CEN'
    x, y = numpy.meshgrid(numpy.arange(301), numpy.arange(301))
    xcen = 150
    ycen = 150
    amplitude = 1.
    twopi = 2.*3.14159265358979
    ravalues = galaxy_list[0]
    decvalues = galaxy_list[1]
    signal = galaxy_list[3]
    radvalues = galaxy_list[4]/0.0656
    ellipvalues = galaxy_list[5]
    pavalues = galaxy_list[6]*twopi/360.0
    indvalues = galaxy_list[7]
    galaxy_image = numpy.zeros((3631, 3631), dtype=numpy.float32)
    for loop in range(len(indvalues)):
        mod = Sersic2D(amplitude=1., r_eff=radvalues[loop],
                       n=indvalues[loop], x_0=xcen, y_0=ycen,
                       ellip=ellipvalues[loop],
                       theta=pavalues[loop]-rotation)
        image = mod(x,y)
        tsig = numpy.sum(image)
        if tsig > 0.:
            image = image/numpy.sum(image)
        if simple:
            xpix, ypix = relpos(ravalues[loop], decvalues[loop],
                                position[0], position[1], rotation,
                                0.0656)
            xpix = xpix + 1023.5
            ypix = ypix + 1023.5
        else:
            xpix, ypix = get_pixel(ravalues[loop], decvalues[loop],
                                   position[0], position[1], rotation,
                                   instrument, aperture)
        nxpix = int(xpix)+792
        nypix = int(ypix)+792
        if (nxpix >= 0) and (nxpix < 3631) and (nypix >= 0) and (nypix < 3631):
            x0 = nxpix - 150
            y0 = nypix - 150
            x1 = nxpix + 151
            y1 = nypix + 151
            xmin = 0
            xmax = 301
            ymin = 0
            ymax = 301
            if x0 < 0:
                xmin = -x0
                x0 = 0
            if y0 < 0:
                ymin = -y0
                y0 = 0
            if y1 > 3631:
                ymax = 301 - (y1-3631)
                y1 = 3631
            if x1 > 3631:
                xmax = 301 - (x1-3631)
                x1 = 3631
            try:
                galaxy_image[y0:y1, x0:x1] = galaxy_image[y0:y1, x0:x1] + \
                    image[ymin:ymax, xmin:xmax]*signal[loop]
            except:
                pass
    return galaxy_image


def do_convolve(scene_image, filter1, path):
    """
    Convolve a scene image with a PSF file for the filter name given.

    Parameters
    ----------

    scene_image:   a two-dimensional numpy float array, the input image to be 
                   convovled with a PSF image

    filter1:       a string giving the NIRISS filter name to use

    path:          a string giving the path to the imaging PSF images

    Returns
    -------

    convolved_image:  a two-dimensional numpy float array of the convolved 
                      image, of the same dimensions are scene_image, or 
                      None is there is an issue.
    """
    psfname = 'niriss_NIS_x1024_y1024_' + filter1.lower() + '_predicted_0_0p00_0p00.fits'
    if path[-1] != '/':
        path = path+'/'
    try:
        psf_image = fits.getdata(path+psfname)
    except:
        print('Failed to read PSF image %d.' % (path+psfname))
        return None
    convolved_image = signal.fftconvolve(scene_image, psf_image, mode='same')
    return convolved_image

def rotate_image(scene_image, angle):
    """
    Rotate the input scene_image by the angle given in degrees.

    Rotations are assumed to be counter-clockwise from the y axis

    Parameters
    ----------

    scene_image:   a numpy two-dimensional array of float values

    angle:         a float value, the angle of rotation in degrees

    Returns:

    rotated_image:  a numpy two-dimensional array of float values, of the 
                    same dimensions as the scene_image array, containing the 
                    rotated version of the image

    The rotation is done using the scipy ndimage package.
    """
    term = angle/360.
    offset = math.floor(term)
    rotangle = angle - offset*360.
    if rotangle == 0.:
        return scene_image
    zmin = numpy.min(scene_image)
    sh1 = scene_image.shape
    rotated_image = ndimage.rotate(scene_image, -rotangle, cval=zmin, order=5)
    sh2 = rotated_image.shape
    if sh2 != sh1:
        xmin = (sh2[1]-sh1[1]) // 2
        xmax = xmin+sh1[1]
        ymin = (sh2[0]-sh1[0]) // 2
        ymax = ymin+sh1[0]
        rotated_image = numpy.copy(rotated_image[ymin:ymax, xmin:xmax])
# ndimage.rotate tends to produce negative artifacts in the image, screen 
# these out to the minimum in the original image
    inds = numpy.where(rotated_image < zmin)
    rotated_image[inds] = zmin
    return rotated_image
# note this needs to be tested for accuracy.  Sharp PSFs may cause issues in
# this code....

def generate_image(star_list, position, rotation=0., simple=False):
    """
    Do the work of making a star scene image.  Each star is one pixel in size.

    Parameters
    ----------

    star_list:     A list of numpy 1-d float arrays holding the catalogue 
                   values for the scene (RA, DEC, mag, signal)

    position:      A two-element list giving the image center (RA, Dec) in 
                   degrees 

    rotation:      An optional float value, the rotation angle in degrees E 
                   of N

    simple:        A boolean value, if True use the simple projection to get 
                   pixel positions, if False use pysiaf.  The latter is the 
                   default.

    Returns
    -------

    scene_image:   A 3631x3631 numpy float array, the scene image.

    new_star_list:   A new list, same structure as star_list, containts the 
                     values for the stars within the field
    """
    nout = 0
    scene_image = numpy.zeros((3631, 3631), dtype=numpy.float32)
    ravalues = star_list[0]
    decvalues = star_list[1]
    signal = star_list[2]
    newravalues = []
    newdecvalues = []
    newsignal = []
    instrument = 'NIRISS'
    aperture = 'NIS_CEN'
    for loop in range(len(ravalues)):
        if simple:
            xpix, ypix = relpos(ravalues[loop], decvalues[loop],
                                position[0], position[1], rotation,
                                0.0656)
            xpix = xpix + 1023.5
            ypix = ypix + 1023.5
        else:
            xpix, ypix = get_pixel(ravalues[loop], decvalues[loop],
                                   position[0], position[1], rotation,
                                   instrument, aperture)
        nxpix = int(xpix)+792
        nypix = int(ypix)+792
        if (nxpix >= 0) and (nxpix < 3631) and (nypix >= 0) and (nypix < 3631):
            scene_image[nypix, nxpix] = scene_image[nypix, nxpix]+signal[loop]
            nout = nout + 1
            newravalues.append(ravalues[loop])
            newdecvalues.append(decvalues[loop])
            newsignal.append(signal[loop])
    new_star_list = [numpy.asarray(newravalues, dtype=numpy.float32),
                     numpy.asarray(newdecvalues, dtype=numpy.float32),
                     numpy.asarray(newsignal, dtype=numpy.float32)]
    return scene_image, new_star_list

def get_pixel(ratarget, dectarget, ra0, dec0, rotation, instrument, aperture):
    """
    Calculate the pixel position of a target for a given aperture and sky 
    postion plus orientation.

    Parameters
    ----------

    ratarget:   a float value, the target RA in decimal degrees

    dectarget:  a float value, the target Dec in decimal degrees

    ra0:        a float value, the field pointing RA in decimal degrees

    dec0:       a float value, the field pointing Dec in decimal degrees

    rotation:   a float value, the field rotation in decimal degrees E of N

    instrument: a string variable giving the instrument name (e.g. 'NIRISS')

    aperture:   a string variable giving the instrument aperture name (e.g. 
                'NIS_CEN')

    Returns
    -------

    xpixel:    the object x pixel position, a float value

    ypixel:    the object y pixel position, a float value

    """
    dtor = 3.14159265358979/180.
    siaf_instance = pysiaf.Siaf(instrument)
    siaf = siaf_instance[aperture]
    v2_arcsec = siaf.V2Ref
    v3_arcsec = siaf.V3Ref
    v2 = v2_arcsec*dtor/3600.
    v3 = v3_arcsec*dtor/3600.
    ra_ref = ra0*dtor
    dec_ref = dec0*dtor
    pa_v3 = rotation*dtor
    mat1 = numpy.array([[numpy.cos(ra_ref) * numpy.cos(dec_ref),
        -numpy.sin(ra_ref) * numpy.cos(pa_v3) + numpy.cos(ra_ref) * numpy.sin(dec_ref) * numpy.sin(pa_v3),
        -numpy.sin(ra_ref) * numpy.sin(pa_v3) - numpy.cos(ra_ref) * numpy.sin(dec_ref) * numpy.cos(pa_v3)],
        [numpy.sin(ra_ref) * numpy.cos(dec_ref),
         numpy.cos(ra_ref) * numpy.cos(pa_v3) + numpy.sin(ra_ref) * numpy.sin(dec_ref) * numpy.sin(pa_v3),
         numpy.cos(ra_ref) * numpy.sin(pa_v3) - numpy.sin(ra_ref) * numpy.sin(dec_ref) * numpy.cos(pa_v3)],
        [numpy.sin(dec_ref), -numpy.cos(dec_ref) * numpy.sin(pa_v3),
         numpy.cos(dec_ref) * numpy.cos(pa_v3)]])

    X = -(mat1[2, 0] * numpy.cos(v2) + mat1[2, 1] * numpy.sin(v2)) * numpy.sin(v3) + mat1[2, 2] * numpy.cos(v3)
    Y = (mat1[0, 0] *  mat1[1, 2] - mat1[1, 0] * mat1[0, 2]) * numpy.cos(v2) + \
      (mat1[0, 1] * mat1[1, 2] - mat1[1, 1] * mat1[0, 2]) * numpy.sin(v2)
    local_roll = numpy.rad2deg(numpy.arctan2(Y, X))
    if local_roll < 0:
        local_roll = local_roll+360.
    attitude_matrix = pysiaf.utils.rotations.attitude(v2_arcsec, v3_arcsec, ra0, dec0, local_roll)
    loc_v2, loc_v3 = pysiaf.utils.rotations.getv2v3(attitude_matrix, ratarget, dectarget)
    xpixel, ypixel = siaf.tel_to_sci(loc_v2, loc_v3)
    return xpixel, ypixel

def relpos(ra1, dec1, ra0, dec0, rotation, pixelsize):
    """
    Calculate the offset from position (RA0, Dec0) to (RA1, Dec1) in pixels.

    Parameters:

    ra1:    A float value, the RA of position 1 in decimal degrees

    dec1:   A float value, the Dec of position 1 in decimal degrees

    ra0:    A float value, the RA of the reference position in decimal degrees

    dec0:    A float value, the Dec of the reference position in decimal 
             degrees

    rotation:   A float value, the rotation angle in degrees E of N

    pixelsize:  A float value, the "pixel" size in arc-seconds

    Returns
    -------

    delxpix:  A float value, the x offset in pixels

    delypix:  A float value, the y offset in pixels
    """
    dtor = math.radians(1.)
    if (ra1 == ra0) and (dec1 == dec0):
        return 0., 0.
    angle=math.atan2(math.sin((ra1-ra0)*dtor),
                     math.cos(dec0*dtor)* \
                     math.tan(dec1*dtor)- \
                     math.sin(dec0*dtor)*math.cos((ra1-ra0)*dtor))/dtor
    if angle > 360.:
        angle=angle-360.
    if angle < 0.:
        angle=angle+360.
    angle = angle - rotation
    arcdist=math.sin(dec0*dtor)*math.sin(dec1*dtor)+math.cos(dec0*dtor)*math.cos(dec1*dtor)*math.cos(dtor*(ra1-ra0))
    if abs(arcdist) > 1.:
        arcdist = 1.
    arcdist=math.acos(arcdist)
    arcdist=arcdist/dtor
    if arcdist < 0.:
        arcdist=arcdist+180.
    arcdist = arcdist*(3600./pixelsize)
    delxpix = -arcdist*math.sin(angle*dtor)
    delypix = arcdist*math.cos(angle*dtor)
    return delxpix, delypix

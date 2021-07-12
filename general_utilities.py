"""
This file contains general utilities for the plot code and image display code.

get_subimage:   extract a subimage from a (2-D) image based on zoom parameters

save_data_set_values:    save plot (x, y) values to an ascii file

save_png_figure:   save a plot to a PNG file

save_ps_figure:    save a plot to a postscipt file

hybrid_transform:   apply an IRAF style hybrid logarithmic transformation to 
                    data values

inverse_hybrid_transform:   apply the inverse of the hubrid logarithmic 
                            transformation to a data value

hybrid_labels:   generate plot labels for the hybrid logarithmic case

parse_data_input_text:   parse a set of input lines to extract data values

slope_calculation:  carry out a standard least squares data fit to (x, y) 
                    values

list_polynomial_fitpars:   write data fit parameters to a file

line_range:   find the next comment line in a line list

put_value:   assign a value to Tkinter text field

separator_line:   make Tkinter separator line object

put_yes_no:  make a Tkinter two choice radio button object

add_yes_no_field:   utility routine to add a yes/no radio button and return the 
                    associated Tkinter int variable (1 for yes, 0 for no)

put_message:   append a string to a scrolled text box or text box

range_check:   given an angle in degrees, return the equivalent angle in the 
               primary range 0 to 360 degrees

"""

import math
import numpy
import matplotlib
import tkinter as Tk
import tkinter.filedialog
from astropy.io import fits as fits

def get_subimage(image, zoom):
    """
    Extract a subimage according to the zoom parameters.

    Parameters
    ----------

    image:   A numpy 2-D array of values, by assumption; if not the program
             returns this value unaltered

    zoom:    A three element list with the zoom parameters (zoom value,
             lower left corner pixel x value, lower left corner pixel y
             value); all values are integers, zoom parameter being > 1
             and the corner values >= 0 to affect the image

    Returns
    -------

    subimage:  A numpy 2-D array of values, a sub-set of original image, or
               the image itself if the zoom[0] value is 1 or if any of the
               inputs do not match what is expected
    """
    try:
        sh1 = image.shape
    except AttributeError:
        return image
    if (zoom[0] <= 1) or (zoom[1] < 0) or (zoom[2] < 0) or \
       (zoom[1] >= sh1[1]) or (zoom[2] >= sh1[0]) or (len(sh1) != 2):
        return image
    else:
        npixelx = sh1[1] // zoom[0]
        npixely = sh1[0] // zoom[0]
        if npixely == 0:
            npixelx = 1
        if npixely == 0:
            npixely = 1
        subimage = numpy.copy(image[
            zoom[2]:zoom[2]+npixely, zoom[1]:zoom[1]+npixelx])
        return subimage


def save_data_set_values(xvalues, yvalues, labelstring=None):
    """
    Save plot values to an ascii output file.

    Parameters
    ----------

        xvalues :  A numpy array of float or integer values, the x plot values

        yvalues :  A numpy array of float or integer values, the y plot values
                   (the same length as xvalues, by assumption)

        labelstring :  An optional string to write to the top of the file.
                       If the value is None, nothing is written.

    Returns
    -------

        None

    """
    outfilename = tkinter.filedialog.asksaveasfilename()
    outfile = open(outfilename, 'w')
    if labelstring is not None:
        print(labelstring, file=outfile)
    for loop in range(len(xvalues)):
        print('%20.6g %20.6g' % (xvalues[loop], yvalues[loop]), file=outfile)
    outfile.close()


def save_png_figure(fig):
    """
    Save the current plot as a PNG file.

    Parameters
    ----------
        fig :   A matplotlib Figure instance, the figure to be
                written out

    Returns
    -------

        None

    """
    outfile = tkinter.filedialog.asksaveasfilename(filetypes=[('PNG',
                                                               '.png')])
    if isinstance(outfile, type('string')):
        s1 = outfile.split('.')
    if 'png' not in s1[-1]:
        outfile = outfile + '.png'
    fig.savefig(outfile, format="PNG")


def save_ps_figure(fig):
    """
    Save the current plot as a Postscript file.

    Parameters
    ----------
        fig :  A matplotlib Figure instance, the figure to be
               written out

    Returns
    -------

        None

    """
    outfile = tkinter.filedialog.asksaveasfilename(filetypes=[('PS',
                                                               '.ps')])
    if isinstance(outfile, type('string')):
        s1 = outfile.split('.')
        if 'ps' not in s1[-1]:
            outfile = outfile + '.ps'
        fig.savefig(outfile, format="PS")


def hybrid_transform(datavalues):
    """
    Apply an IRAF-style hybrid log transformation.

    This routine transforms a set of input data values for the hybrid
    log plot:

    Numbers with absolute value > 10 have the sign times the base
    10 log of the absolute value i.e. 100 -> 2.0, -100 -> -2.0.

    Numbers with absolute value <= 10 are scaled down by a factor of
    10 to lie in the range from +1.0 to -1.0 i.e 8 -> 0.8 and -6 -> -0.6.

    Parameter
    ---------
        datavalues :   A numpy array of numbers (real or integer), or a
                       single float value

    Returns
    -------
        newdatavalues :   A numpy array of the same length as datavalues,
                          with the transformed numbers, or a single float
                          value with the transformed mumber

    """
    try:
        newdatavalues = datavalues.astype(numpy.float32)
        inds = numpy.where(numpy.abs(datavalues) < 10.)
        newdatavalues[inds] = datavalues[inds]/10.
        inds = numpy.where(datavalues >= 10.)
        newdatavalues[inds] = numpy.log10(newdatavalues[inds])
        inds = numpy.where(datavalues <= -10.)
        newdatavalues[inds] = -1.0*numpy.log10(numpy.abs(newdatavalues[inds]))
    except Exception:
        if abs(datavalues) >= 10.:
            newdatavalues = math.log10(abs(datavalues))
            if datavalues < 0.:
                newdatavalues = -1.*newdatavalues
        else:
            newdatavalues = datavalues/10.
    return newdatavalues


def inverse_hybrid_transform(value):
    """
    Transform back from the IRAF-style hybrid log values.

    This takes the hybrid log value and transforms it back to the
    actual value.  That value is returned.  Unlike the hybrid_transform
    function, this works on single values not a numpy array.  That is because
    one is not going to have a data array in the hybrid transformation form.

    Parameters
    ----------
        value :  A real number to be transformed back from the hybrid
                 log scaling

    Returns
    -------
        newvalue :   The associated value that maps to the given hybrid
                     log value.

    """
    if value < 0.:
        workvalue = -1.0*value
        sign = -1.0
    else:
        workvalue = value
        sign = +1.0
    if workvalue < 1.0:
        newvalue = 10.*workvalue
    else:
        newvalue = 10.**workvalue
    newvalue = sign * newvalue
    return newvalue


def hybrid_labels(datavalues):
    """
    Create labels for the hybrid log case.

    This routine makes hybrid log labels for the range of values in a
    given data set: it assumes that the range is large enough that
    integer ticks between  10 and +10 are suitable.  There is no reason
    to use the hybridlog option if the range is smaller than this, or
    indeed if the range is less than a couple of orders of magnitude.

    Parameters
    ----------
        datavalues :     A one-dimensional numpy array of float data values.

    Returns
    -------
        rangeout :  A list of the positions for the tick marks (floats).

        ticksout :  A list of the tick labels (strings).

    Both the return values are used with set_xticks/set_xticklabels
    (main matplotlib) or xticks (pyplot) or the y axis equivalents.
    """
    ticklabels = []
    datamax = numpy.max(datavalues)
    datamin = numpy.min(datavalues)
    baserange = numpy.arange(-10, 12, 2)/10.0
    subvalues = numpy.arange(2, 11)
    subvalues = numpy.log10(subvalues)
    n1 = int(datamin)
    if datamin < 0.:
        n1 = n1 - 1
    n2 = int(datamax)
    if n1 < -1:
        for loop in range(-1, n1-1, -1):
            for n1 in range(len(subvalues)-1, -1, -1):
                baserange = numpy.insert(baserange, 0, loop-subvalues[n1])
    if n2 > 1:
        for loop in range(1, n2+1):
            for n1 in range(len(subvalues)):
                baserange = numpy.append(baserange, loop+subvalues[n1])
    ticklabels = []
    tickzero = None
    for loop in range(len(baserange)):
        if (tickzero is None) and (datamin < baserange[loop]) and \
           (datamax > baserange[loop]):
            tickzero = baserange[loop]
    for loop in range(len(baserange)):
        if baserange[loop] < datamin:
            baserange[loop] = tickzero
            ticklabels.append('')
        elif baserange[loop] > datamax:
            baserange[loop] = tickzero
            ticklabels.append('')
        else:
            if abs(baserange[loop]) < 1.1:
                ticklabels.append('%.0f' % (baserange[loop]*10.0))
            else:
                if (baserange[loop] > 1.1) and (math.floor(
                        baserange[loop]) == baserange[loop]):
                    ticklabels.append(r'$10^%d$' % (baserange[loop]))
                elif (baserange[loop] < -1.1) and (math.floor(
                        baserange[loop]) == baserange[loop]):
                    ticklabels.append(r'$-10^%d$' % (abs(baserange[loop])))
                else:
                    ticklabels.append('')
    rangeout = []
    ticksout = []
    for loop, value in enumerate(baserange):
        if not math.isnan(value):
            rangeout.append(value)
            ticksout.append(ticklabels[loop])
    return rangeout, ticksout


def parse_data_input_text(text):
    """
    This routine parses a list of text lines to extract the numerical
    values for a set, used with the text entry option.

    Parameters
    ----------
    text : a list of strings

    Returns
    -------
    xvalues : a list of float values, the x data values for a set

    dxvalues1 : a list of float values, the lower x error values for a set

    dxvalues2 : a list of float values, the upper x error values for a set

    yvalues : a list of float values, the y data values for a set

    dyvalues1 : a list of float values, the lower y error values for a set

    dyvalues2 : a list of float values, the upper y error values for a set

    errorflag : boolean value, flags whether the uncertaintes are defined

    """
    xvalues = []
    dxvalues1 = []
    dxvalues2 = []
    yvalues = []
    dyvalues1 = []
    dyvalues2 = []
    lines = text.split('\n')
    for line in lines:
        if '#' in line:
            pass
        else:
            values = line.split()
            numbers = []
            errorflag = False
            for loop in range(len(values)):
                try:
                    v1 = float(values[loop])
                    numbers.append(v1)
                except ValueError:
                    pass
            if len(numbers) == 2:
                xvalues.append(numbers[0])
                yvalues.append(numbers[1])
                dxvalues1.append(0.0)
                dxvalues2.append(0.0)
                dyvalues1.append(0.0)
                dyvalues2.append(0.0)
            elif len(numbers) == 4:
                xvalues.append(numbers[0])
                yvalues.append(numbers[2])
                dxvalues1.append(numbers[1])
                dxvalues2.append(numbers[1])
                dyvalues1.append(numbers[3])
                dyvalues2.append(numbers[3])
                errorflag = True
            elif len(numbers) == 6:
                xvalues.append(numbers[0])
                yvalues.append(numbers[3])
                dxvalues1.append(numbers[1])
                dxvalues2.append(numbers[2])
                dyvalues1.append(numbers[4])
                dyvalues2.append(numbers[5])
                errorflag = True
            elif len(numbers) > 2:
                xvalues.append(numbers[0])
                yvalues.append(numbers[1])
                dxvalues1.append(0.0)
                dxvalues2.append(0.0)
                dyvalues1.append(0.0)
                dyvalues2.append(0.0)
    return xvalues, dxvalues1, dxvalues2, yvalues, dyvalues1, \
        dyvalues2, errorflag


def slope_calculation(xdata, ydata, yerrors=None):
    """
    Calculate a standard least-squares linear fit to input data.

    This is for the case where the standard library function fails, as is
    found to be the case in some instances.

    The treatment is as in "Numerical Recipes in C" (second edition) section
    15.2.

    Parameters
    ----------

    xdata:  A one-dimensional numpy float or integer array of the x values

    ydata:  A one-dimensional numpy float or integer array of the y values,
            which must be the same length as the xdata array

    yerrors:  An optional one-dimensional numpy float array of the y value
              uncertainties, which must be the same length as the xdata array
              if defined.  If not defined, the values are set to a constant.
              The error values must be strictly positive, non-zero.

    Returns
    -------

    slope:           A float value, the best fit slope

    intercept:       A float value, the best fit intercept

    slopeerror:      A float value, the uncertainty in the best fit slope
                     (standard deviation estimate)

    intercepterror:  A float value, the uncertainty in the best fit intercept
                     (standard deviation estimate)

    covariance:      A float value, the covariance between the fit parameters

    correlation:     A float value, the correlation coefficient of the fit

    """
    if yerrors is None:
        yerrors = ydata*0. + 1.
    else:
        inds = numpy.where(yerrors > 0.)
        mean1 = numpy.mean(yerrors[inds])
        inds = numpy.where(yerrors <= 0.)
        yerrors[inds] = mean1
    if (xdata.shape != ydata.shape) or (len(xdata.shape) > 1) or \
       (xdata.shape != yerrors.shape):
        return None, None, None, None, None, None
    invvariance = 1./(yerrors*yerrors)
    sum1 = numpy.sum(invvariance)
    sum2 = numpy.sum(xdata*invvariance)
    sum3 = numpy.sum(ydata*invvariance)
    xmean = sum2/sum1
    t1 = (xdata - xmean) / yerrors
    sum4 = numpy.sum(t1*t1)
    slope = numpy.sum(t1*ydata/yerrors)/sum4
    intercept = (sum3 - (slope*sum2))/sum1
    interceptvariance = (1. + (sum2*sum2)/(sum1*sum4))/sum1
    slopevariance = 1./sum4
    covariance = -1.*sum2/(sum1*sum4)
    correlation = covariance/numpy.sqrt(slopevariance*interceptvariance)
    return slope, intercept, math.sqrt(slopevariance), \
        math.sqrt(interceptvariance), covariance, correlation


def list_polynomial_fitpars(fit_type, fit_order, fitpars,
                            filename='fit_values.txt'):
    """
    This code writes out the fit parameters to  an ascii output file

    Parameters
    ----------
    fit_type : integer variable, flags the function used in the fit

    fit_order : integer variable, gives the order of the fit function

    fitpars : numpy float array, the fit parameters

    filename : an optional string variable that gives the file name for the 
               output; default is "fit_values.txt"

    Returns
    -------

    None

    """
    outfile = open(filename, 'a')
    print('Order %d %s polynomial fit:' % (fit_order, fit_type), file=outfile)
    print('Parameter    Value', file=outfile)
    for loop in range(len(fitpars)):
        print('%3d %g' % (loop, fitpars[loop]), file=outfile)
    print(' ', file=outfile)
    outfile.close()


def line_range(lines, ind1, comment_flag='#'):
    """
    Find a range of data lines within a line list.

    Given an input line list and a starting index, subsequent lines are
    examined to see where the next comment line is.  Comment lines are
    assumed to start with the # character by default, or one can set this
    with the comment_flag variable in the call.  Lines that are not comments
    are assumed to be data lines.  The index of the next comment line is
    returned, or the index that gives a range to the end of the line list
    where there is no such comment line after the index specified.

    Parameters
    ----------
        lines :  A list of input lines (assumed to be from the readlines()
                 function)

        ind1 :   A starting index in the list of lines

        comment_flag:   An optional string variable that marks comment lines

    Returns
    -------
        n1 : an integer value for the next comment line (assumed to
             start with '#') in the list of input lines, or the index
             for the length of the line list if no other comment line is
             found
    """
    ncomment = len(comment_flag)
    for n1 in range(ind1+1, len(lines)):
        if comment_flag in lines[n1][0:ncomment]:
            return n1
    return len(lines)


def round_float(plotgui, value, minimum_flag):
    """
    Round a floating point value to the nearest significant figure.

    If the plotgui.matplotlib_rounding flag is False this performs rounding of
    floating point values to convenient values for such things as
    calculating plot ranges.  In some cases the matplotlib rounding is
    not very good so this is provided as a way to get better results.
    When one is using a logarithmic scaling the matplotlib ropunding is
    not always good, for example, and this routine can do better.  In
    other cases this routine does not work as well, such as when the
    total range of values on the plot is small.

    Parmeters
    ---------

        value :         A floating point number to be rounded off

        minimum_flag :  A boolean value value to determine whether the
                        value should be rounded down.  If True values
                        are rounded down as one wishes for the minimum
                        value in a plot, and if False values are rounded
                        up as one wishes for the maximum value in a plot.

    Returns
    -------
       rounded_value :    The rounded off floating point number calculated
                          from `value', or the input value in the case
                          where the self.matplotlib_routine flag is True.

    The code differs from use of the floor and ceil functions in that
    it tries to round to the nearest significant figure for a given
    value.  So passing in a value of 1.2e-08, for exmaple, with
    minimum_flag = False returns 2.0e-08, and if the flag is True it
    returns 1.0e-08.

    The code uses the math log10, floor, and ceil functions.

    """
    if value == 0.:
        return value
    if value < 0.:
        sign = -1.
        value1 = -value
    else:
        sign = 1.
        value1 = value
    power = math.log10(value1)
    if power < 0.:
        exp = int(power-1.)
    else:
        exp = int(power)
    shift = 10.**exp
    x = value1/shift
    delx = 1.0
    if x < 1.7:
        x = x*10.
        shift = shift/10.
    elif x < 2.5:
        x = x*5.
        shift = shift/5.
    if (minimum_flag) and sign > 0.:
        x = math.floor(x)
    elif (minimum_flag) and sign < 0.:
        x = math.ceil(x)
    elif (not minimum_flag) and sign > 0.:
        x = math.ceil(x)
    elif (not minimum_flag) and sign < 0.:
        x = math.floor(x)
    rounded_value = x*shift*sign
    # If the rounded value is very close to the input value, offset
    # by one unit in x...not normally part of the routine, but needed
    # for matplotlib plots because of symbols close to the edges of
    # the plot.
    ratio = abs(value/rounded_value)
    if (ratio > 0.97) and (ratio < 1.03):
        if (minimum_flag) and sign > 0.:
            rounded_value = (x-delx)*shift*sign
        elif (minimum_flag) and sign < 0.:
            rounded_value = (x+delx)*shift*sign
        elif (not minimum_flag) and sign > 0.:
            rounded_value = (x+delx)*shift*sign
        elif (not minimum_flag) and sign < 0.:
            rounded_value = (x-delx)*shift*sign
    return rounded_value


def put_value(value, field):
    """
    Place a value in a widgit text field.

    Any current contents of the field are deleted.

    Parameters
    ----------

        value :  the string value to be placed in the text field

        field :  the tkinter text field variable where the string is to
                 be put

    Returns
    -------

        None

    """
    try:
        s1 = field.get()
        field.delete(0, last=len(s1))
        field.insert(0, str(value))
    except Exception:
        pass

def separator_line(parent, w1, h1, pad, flag, packvalue=Tk.TOP,
                   gridvalues=None):
    """
    Create the Tkinter canvas object making a separator line.

    This is a utility routine to make a separator line within a Tkinter
    frame.

    Parameters
    ----------

        parent :  A Tkinter Frame variable, that will contain the line

        w1 :      An integer value for the line canvas width (pixels)

        h1 :      An integer value for the line canvas height (pixels)

        pad :     An integer value for the line padding (pixels)

        flag :    A Boolean value, the line is horizontal if the
                  value is True, vertical otherwise

        packvalue :  An optional Tkinter pack direction (one of Tk.LEFT,
                     Tk.RIGHT, Tk.TOP, or Tk.BOTTOM) with default value
                     of Tk.TOP; if None use the gridvalues

        gridvalues : An optional list of seven grid parameters (row, column, 
                     rowspan, columnspan, sticky_string, padx, pady) to 
                     use if the packvalue is None; the first 4 values are 
                     integers, the next is a string using Tk for tkinter as in
                     'Tk.E+Tk.W' for example, and the last 2 are integers

    Returns
    -------

        linecanvas:  the Tkinter Canvas variable for the line, in case
                     the user wishes to modify it later

    For a vertical line normally the height will be much larger than the
    width, as in

    sl = self.separator_line(frame, 10, 300, 5, False)

    while for a horizontal line normally the width will be much larger
    than the height as in

    sl = self.separator_line(frame, 500, 5, 5, True)

    """
    lincanvas = Tk.Canvas(parent, height=h1, width=w1)
    if packvalue is None:
        lincanvas.grid(row=gridvalues[0], column=gridvalues[1],
                       rowspan=gridvalues[2], columnspan=gridvalues[3],
                       sticky=gridvalues[4], padx=gridvalues[5],
                       pady=gridvalues[7])
    else:
        lincanvas.pack(side=packvalue, fill=Tk.BOTH, expand=Tk.YES)
    if flag:
        lincanvas.create_line(pad, h1/2, w1-pad, h1/2)
    else:
        lincanvas.create_line(w1/2, pad, w1/2, h1-pad)
    return lincanvas

def put_yes_no(root, var, labels, flag):
    """
    Create a Tkinter yes/no radio button.

    This is a utility routine to make a yes/no radio button in Tkinter.
    The required variables are passed to the code and it produces the
    button pair.

    Parameters
    ----------

        root :  The tkinter frame variable to hold the buttons

        var :   The tkinter IntVar that is used to communicate with the
                buttons

        labels : A two element list with strings for the two states of
                 the buttons, first "yes" then "no".

        flag :  A boolean value that causes the code to set the yes field
                (the first of the two) if it is True.  If the value is
                False the no field (the second of the two) is set instead.

    Returns
    -------

        None

    """
    yesfield = Tk.Radiobutton(root, text=labels[0], variable=var, value=1)
    yesfield.grid(row=0, column=0, sticky=Tk.W)
    nofield = Tk.Radiobutton(root, text=labels[1], variable=var, value=0)
    nofield.grid(row=0, column=1, sticky=Tk.W)
    if flag:
        yesfield.select()
    else:
        nofield.select()

def add_yes_no_field(frame, text, flag):
    """
    Add a yes/no radio button pair to a frame

    Parmeters
    ---------

    frame:   a tkinter Frame variable that will hold the radio buttons

    text:    a tuple with text strings for the two buttons

    """
    field1 = Tk.Frame(frame)
    field1.pack(side=Tk.TOP)
    label = Tk.Label(field1, text=text)
    label.pack(side=Tk.LEFT)
    newvar = Tk.IntVar()
    b1 = Tk.Frame(field1)
    put_yes_no(b1, newvar, ['Yes', 'No'], flag)
    b1.pack()
    return newvar

def put_message(textbox, outstr):
    """
    Append a string to a message box.

    Parameters
    ----------

    textbox:   a Tkinter textbox or scrolledtext variable wherein to put the 
               message string

    outstr:    a string variable, that is appended to the text of textbox
    """
    textbox.insert(Tk.END, outstr)
    textbox.see(Tk.END)


def range_check(angle):
    """
    Check that an angle in degrees is in the primary range 0 to 360.

    Parameters
    ----------

    angle:    a float value, an angle in degrees

    Returns
    -------

    newangle:   a float value, the angle in the range 0 to 360 degrees

    If the angle value passed is inf, -inf, of NaN a value of zero is returned.
    """
    if not math.isfinite(angle):
        return 0.
    if (angle >= 0.) and (angle <= 360.):
        return angle
    term = angle/360.
    offset = math.floor(term)
    newangle = angle - offset*360.
    return newangle

def save_fits(image):
    if True:
#    try:
        outfile = tkinter.filedialog.asksaveasfilename(
            filetypes=[('FITS', '*.fits')], title='Output FITS File Name')
        print(outfile)
        if (isinstance(outfile, type('string'))) and (len(outfile) > 0):
            hdu = fits.PrimaryHDU(image)
            hdu.writeto(outfile)
            values = outfile.split('/')
            filename = values[-1]
            print('Have saved the current image to: %s.\n' % (filename))
#    except:
#        pass

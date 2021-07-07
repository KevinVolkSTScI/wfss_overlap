#! /usr/bin/env python
#
"""
Tool to access the possibility of overlap of spectra for a field of view.

This code takes as input either a pair or Mirage source input files (one for 
stars and one for galaxies) or a scene image made from such input files and 
allows the user to investigate the resulting scene at arbitrary rotations 
to determine if spectral order overlap will be an issue for any given source 
in the field.

"""
import sys
import os
import numpy
import tkinter as Tk
import tkinter.ttk
import tkinter.messagebox
import tkinter.filedialog
from tkinter.scrolledtext import ScrolledText
import astropy.io.fits as fits
import fits_image_display
import general_utilities
import scene_image
import wfss_scene
BGCOL = '#F8F8FF'
FILTERNAMES = ['F090W', 'F115W', 'F140M', 'F140W', 'F158M', 'F200W']
GRISMNAMES = ['GR150R', 'GR150C']


class WFSSOverlap(Tk.Frame):
    """
    This class brings up the WFSS Overlap Tool window.

    Parameters
    ----------

    Tk.Frame:   The base class of the object, matching a Tkinter root or
                Toplevel variable

    Returns
    -------

    The class variable is returned, effectively.
    """


    def __init__(self, parent=None, **args):
        self.scene_image = None
        self.convolved_image = None
        self.filtername = 'F090W'
        self.grismname = 'GR150R'
        self.position = None
        self.imagewin = None
        self.last_type = None
        if parent is not None:
            # initialize the window and make the plot area.
            Tk.Frame.__init__(self, parent, args)
            self.root = parent
            self.make_wfss_window()


    def make_wfss_window(self):
        """
        Make the main control window.

        Parameters
        ----------

        No parameters.

        Returns
        -------

        None
        """
        if self.root is not None:
            wfsswindow = self.root
        else:
            wfsswindow = Tk.Toplevel()
            wfsswindow.title("WFSS Overlap Tool")
        wfsswindow.config(bg=BGCOL)
        frame1 = Tk.Frame(wfsswindow)
        frame1.pack(side=Tk.TOP)
        self.message_area = ScrolledText(
            frame1, height=5, width=100, bd=1, relief=Tk.RIDGE, wrap=Tk.NONE)
        self.message_area.pack(side=Tk.TOP)
        filterframe = Tk.Frame(frame1)
        filterframe.pack(side=Tk.TOP)
        self.filtervar = Tk.IntVar()
        for loop in range(len(FILTERNAMES)):
            rb1 = Tk.Radiobutton(filterframe, text=FILTERNAMES[loop],
                                 variable=self.filtervar, value=loop,
                                 command=self.set_filter)
            rb1.pack(side=Tk.LEFT)
        self.filtervar.set(0)
        self.filtername = FILTERNAMES[0]
        grismframe = Tk.Frame(frame1)
        grismframe.pack(side=Tk.TOP)
        self.grismvar = Tk.IntVar()
        for loop in range(len(GRISMNAMES)):
            rb1 = Tk.Radiobutton(grismframe, text=GRISMNAMES[loop],
                                 variable=self.grismvar, value=loop,
                                 command=self.set_grism)
            rb1.pack(side=Tk.LEFT)
        self.grismvar.set(0)
        self.grismname = GRISMNAMES[0]
        fields = Tk.Frame(frame1)
        fields.pack()
        lab1 = Tk.Label(fields, text='Background (ADU/s):')
        lab1.grid(row=0, column=0)
        self.background_entry = Tk.Entry(fields, width=30)
        self.background_entry.grid(row=0, column=1)
        self.background_entry.insert(0, '0.1')
        lab1 = Tk.Label(fields, text='Point Source File Name:')
        lab1.grid(row=1, column=0)
        self.star_entry = Tk.Entry(fields, width=30)
        self.star_entry.grid(row=1, column=1)
        lab1 = Tk.Label(fields, text='Extended Source File Name:')
        lab1.grid(row=2, column=0)
        self.galaxy_entry = Tk.Entry(fields, width=30)
        self.galaxy_entry.grid(row=2, column=1)
        lab1 = Tk.Label(fields, text='File Path:')
        lab1.grid(row=3, column=0)
        self.path_entry = Tk.Entry(fields, width=30)
        self.path_entry.grid(row=3, column=1)
        wd = os.getcwd()
        if len(wd) == 0:
            wd = './'
        if wd[-1] != '/':
            wd = wd+'/'
        self.path_entry.insert(0, wd)
        self.filepath1 = Tk.Button(
            fields, text='Select File',
            command=lambda: self.select_directory(
                self.star_entry, self.path_entry))
        self.filepath1.grid(row=1, column=2)
        self.filepath2 = Tk.Button(
            fields, text='Select File',
            command=lambda: self.select_directory(
                self.galaxy_entry, None))
        self.filepath2.grid(row=2, column=2)
        self.filepath3 = Tk.Button(
            fields, text='Select Directory',
            command=lambda: self.select_directory(None, self.path_entry))
        self.filepath3.grid(row=3, column=2)
        self.filepath4 = Tk.Button(
            fields, text='Select Directory',
            command=lambda: self.select_directory(None, self.psf_path_entry))
        self.filepath4.grid(row=4, column=2)
        lab1 = Tk.Label(fields, text='WFSS PSF Path:')
        lab1.grid(row=4, column=0)
        self.psf_path_entry = Tk.Entry(fields, width=30)
        self.psf_path_entry.grid(row=4, column=1)
        wd = os.getcwd()
        if len(wd) == 0:
            wd = './'
        if wd[-1] != '/':
            wd = wd+'/'
        self.psf_path_entry.insert(0, wd)
        lab1 = Tk.Label(fields, text='Sky Position (degrees):')
        lab1.grid(row=5, column=0)
        self.position_entry = Tk.Entry(fields, width=30)
        self.position_entry.grid(row=5, column=1)
        buttonframe = Tk.Frame(frame1)
        buttonframe.pack()
        b1 = Tk.Button(
            buttonframe, text='Compute Scene', command=self.make_scene)
        b1.pack(side=Tk.LEFT)
        b1 = Tk.Button(
            buttonframe, text='Save Scene Image', command=self.save_scene)
        b1.pack(side=Tk.LEFT)
        b1 = Tk.Button(
            buttonframe, text='Read Scene Image', command=self.load_scene)
        b1.pack(side=Tk.LEFT)
        b1 = Tk.Button(
            buttonframe, text='Display Scene Image', command=self.display_scene)
        b1.pack(side=Tk.LEFT)
        f1 = Tk.Frame(frame1)
        f1.pack()
        typeframe = Tk.Frame(f1)
        typeframe.pack(side=Tk.LEFT)
        lab1 = Tk.Label(typeframe, text='Display type:')
        lab1.pack(side=Tk.LEFT)
        tframe = Tk.Frame(typeframe)
        tframe.pack(side=Tk.LEFT)
        self.typevar = Tk.IntVar()
        rb1 = Tk.Radiobutton(tframe, text='Scene', variable=self.typevar, 
                             value=0, command=lambda: self.redisplay(None))
        rb1.pack(side=Tk.LEFT)
        rb1 = Tk.Radiobutton(tframe, text='POM', variable=self.typevar, 
                             value=1, command=lambda: self.redisplay(None))
        rb1.pack(side=Tk.LEFT)
        rb1 = Tk.Radiobutton(tframe, text='Image Area', variable=self.typevar, 
                             value=2, command=lambda: self.redisplay(None))
        rb1.pack(side=Tk.LEFT)
        self.typevar.set(2)
        showframe = Tk.Frame(f1)
        showframe.pack(side=Tk.LEFT)
        lab1 = Tk.Label(showframe, text='   Image type:')
        lab1.pack(side=Tk.LEFT)
        tframe = Tk.Frame(showframe)
        tframe.pack()
        self.imagevar = Tk.IntVar()
        rb1 = Tk.Radiobutton(tframe, text='Direct', variable=self.imagevar, 
                             value=0, command=lambda: self.redisplay(None))
        rb1.pack(side=Tk.LEFT)
        rb1 = Tk.Radiobutton(tframe, text='Dispersed', variable=self.imagevar, 
                             value=1, command=lambda: self.redisplay(None))
        rb1.pack(side=Tk.LEFT)
        rb1 = Tk.Radiobutton(tframe, text='Unconvolved',
                             variable=self.imagevar, 
                             value=2, command=lambda: self.redisplay(None))
        rb1.pack(side=Tk.LEFT)
        self.imagevar.set(1)
        anglefield = Tk.Frame(frame1)
        anglefield.pack()
        self.anglevar = Tk.DoubleVar()
        self.angle_slider = Tk.Scale(anglefield, orient=Tk.HORIZONTAL,
                                     length = 720, from_=0., to=360.,
                                     resolution=0.01,
                                     label='Rotation Angle (degrees)', 
                                     variable=self.anglevar)
        self.angle_slider.bind("<ButtonRelease-1>", self.apply_angle)
        self.angle_slider.bind("<KeyPress>", self.apply_angle)
        self.angle_slider.pack()
        stepframe = Tk.Frame(frame1)
        stepframe.pack()
 #       lab1 = Tk.Label(stepframe, text='Step Mode')
 #       lab1.pack(side=Tk.LEFT)
        self.step_button = Tk.Button(
            stepframe, text='Step +', command=lambda: self.run_step(True))
        self.step_button.pack(side=Tk.LEFT)
        self.step_button = Tk.Button(
            stepframe, text='Step -', command=lambda: self.run_step(False))
        self.step_button.pack(side=Tk.LEFT)
        lab1 = Tk.Label(stepframe, text='Minimum')
        lab1.pack(side=Tk.LEFT)
        self.min_angle_entry = Tk.Entry(stepframe, width=7)
        self.min_angle_entry.pack(side=Tk.LEFT)
        self.min_angle_entry.insert(0, '0.0')
        lab1 = Tk.Label(stepframe, text='Maximum')
        lab1.pack(side=Tk.LEFT)
        self.max_angle_entry = Tk.Entry(stepframe, width=7)
        self.max_angle_entry.pack(side=Tk.LEFT)
        self.max_angle_entry.insert(0, '360.0')
        lab1 = Tk.Label(stepframe, text='Step')
        lab1.pack(side=Tk.LEFT)
        self.angle_step_entry = Tk.Entry(stepframe, width=7)
        self.angle_step_entry.pack(side=Tk.LEFT)
        self.angle_step_entry.insert(0, '0.5')
#        b1 = Tk.Button(
#            stepframe, text='Movie', command=self.run_movie)
#        b1.pack(side=Tk.LEFT)
        self.step_button.pack(side=Tk.LEFT)
        buttonframe = Tk.Frame(frame1)
        buttonframe.pack()
        b1 = Tk.Button(buttonframe, text='Close', command=self.root.destroy)
        b1.pack(side=Tk.LEFT)


    def set_filter(self):
        """
        Set the filter name according to the radio button selection.

        Parameters
        ----------

        No parameters

        Returns
        -------

        None
        """
        var = self.filtervar.get()
        self.filtername = FILTERNAMES[var]


    def set_grism(self):
        """
        Set the grism name according to the radio button selection.

        Parameters
        ----------

        No parameters

        Returns
        -------

        None
        """
        var = self.grismvar.get()
        self.grismname = GRISMNAMES[var]


    def apply_angle(self, event):
        self.redisplay(None)
        general_utilities.put_message(
                self.message_area,
                'Angle set to: %f degrees.\n' % (self.anglevar.get()))
        
    def redisplay(self, event):
        image_option = self.imagevar.get()
        if self.imagewin is None:
            pass
        else:
            angle = self.anglevar.get()
            rotated_image = scene_image.rotate_image(self.scene_image, angle)
            self.generate_image(rotated_image)
            self.imagewin.displayImage(angle=angle)


    def make_scene(self):
        """
        Create the scene image.

        Parameters
        ----------

        None

        Returns
        -------

        None
        """
        if True:
#        try:
            background = float(self.background_entry.get())
            if background < 0.:
                background = 0.
            general_utilities.put_message(
                self.message_area,
                'Background level set to: %f ADU/s.\n' % (background))
            starname = self.star_entry.get()
            galaxies_image = None
            extname = self.galaxy_entry.get()
            position_str = self.position_entry.get()
            psf_path = self.psf_path_entry.get()
            if len(psf_path) == 0:
                psf_path = './'
            if psf_path[-1] != '/':
                psf_path = psf_path+'/'
            if ',' in position_str:
                position = position_str.split(',')
            else:
                position_str = position_str.lstrip(' ')
                position_str = position_str.rstrip(' ')
                position = position_str.split(' ')
            newpos = []
            for loop in range(len(position)):
                try:
                    p1 = float(position[loop])
                    newpos.append(p1)
                except:
                    pass
            if len(newpos) >= 2:
                position = [newpos[0], newpos[1]]
                position = numpy.asarray(position)
            else:
                general_utilities.put_message(
                    self.message_area, 
                    'Position string not properly defined, will use mean of star positions.\n')
                position = None
            path = self.path_entry.get()
            if len(path) == 0:
                path = './'
            if path[-1] != '/':
                path = path+'/'
            if position is None:
                position = self.mean_position(path+starname)
                outstr = '%15.8f %15.8f' % (position[0], position[1])
                general_utilities.put_value(outstr, self.position_entry)
            stars_image, star_list = scene_image.make_star_image(
                path+starname, position, self.filtername, path=psf_path)
            if not stars_image is None:
                general_utilities.put_message(
                    self.message_area, 
                    'Have made star scene image from file %s.\n' % (starname))
            else:
                general_utilities.put_message(
                    self.message_area, 
                    'Error in making the star scene image from file %s.\n' % (
                        starname))
                return
            try:
                galaxies_image = scene_image.make_galaxy_image(
                    path+extname, position, self.filtername, path=psf_path)
            except:
                galaxies_image = None
            if not galaxies_image is None:
                general_utilities.put_message(
                    self.message_area, 
                    'Have made galaxies scene image from file %s.\n' % (extname))
                stars_image = stars_image+galaxies_image
            stars_image = stars_image+background
            self.scene_image = stars_image
#        except:
#            general_utilities.put_message(
#                self.message_area, 
#                'An error encountered making the scene image.\n')


    def save_scene(self):
        """
        Save the current scene image to a FITS file.

        Parameters
        ----------

        None

        Returns
        -------

        None
        """
        try:
            outfile = tkinter.filedialog.asksaveasfilename(
                filetypes=[('FITS', '*.fits')], title='Output FITS File Name')
            if (isinstance(outfile, type('string'))) and (len(outfile) > 0):
                hdu = fits.PrimaryHDU(self.scene_image)
                hdu.writeto(outfile)
            values = outfile.split('/')
            filename = values[-1]
            general_utilities.put_message(
                self.message_area, 
                'Have saved the current scene image to: %s.\n' % (filename))
        except:
            pass


    def load_scene(self):
        """
        Load a scene image from a FITS file.

        Parameters
        ----------

        None

        Returns
        -------

        None
        """
        try:
            infile = tkinter.filedialog.askopenfilename(
                filetypes=[('FITS', '*.fits')], title='Output FITS File Name')
            if (isinstance(infile, type('string'))) and (len(infile) > 0):
                try:
                    newimage = fits.getdata(infile)
                except:
                    try:
                        newimage = fits.getdata(infile, ext=1)
                    except:
                        newimage = None
            else:
                newimage = None
            dims = newimage.shape
            if (dims[0] == 3631) and (dims[1] == 3631) and (len(dims) == 2):
                self.scene_image = newimage
                general_utilities.put_message(
                    self.message_area, 
                    'Have read the image from: %s.\n' % (infile))
            else:
                general_utilities.put_message(
                    self.message_area, 
                    'Could not read the image from: %s.\n' % (inname))
        except:
            pass


    def display_scene(self):
        """
        Make the output image at the current rotation and display it.

        Parameters
        ----------

        None

        Returns
        -------

        None
        """
        option = False
        if self.imagewin is None:
            win1 = Tk.Toplevel(self.root)
            self.imagewin = fits_image_display.ImageGUI(win1)
            self.imagefilename = ' '
            self.indpi = 300
            self.imagewin.make_image_window()
            option = True
        if self.imagewin is None:
            return
        angle = self.anglevar.get()
        display_option = self.typevar.get()
        image_option = self.imagevar.get()
        total_option = image_option*10+display_option
        work_image = scene_image.rotate_image(self.scene_image, angle)
        if work_image is None:
            return
        self.generate_image(work_image)
        if total_option != self.last_type:
            option = True
            self.last_type = total_option
        self.imagewin.displayImage(getrange=option, angle=angle)


    def generate_image(self, work_image):
        """
        Take a rotated unconvolved image and make the output image requested.

        Parameters
        ----------

        work_image:   A numpy 2-D float array of dimension [3631, 36310]

        Returns
        -------

        new_image:   A numpy 2-D float array, the requested image
        """
        psf_path = self.psf_path_entry.get()
        image_option = self.imagevar.get()
        if image_option == 0:
            work_image = scene_image.do_convolve(
                work_image, self.filtername, psf_path)
            self.extract_image(work_image)
        elif image_option == 1:
            pom_image = numpy.copy(work_image[655:2977, 655:2977])
            dispersed_image = wfss_scene.wfss_scene(
                pom_image, self.filtername, self.grismname, 0, 0,
                path=psf_path)
            big_image = work_image*0.
            big_image[655:2977, 655:2977] = dispersed_image
            self.extract_image(big_image)
        else:
            self.extract_image(work_image)


    def extract_image(self, work_image):
        """
        Extract the image are according the typevar value, and set the 
        image display to this area.

        Parameters
        ----------

        work_image:   a scene image, numpy 2-d float array 3631x3631 pixels

        Returns
        -------

        None
        """
        display_option = self.typevar.get()
        if display_option == 0:
            newimage = work_image
        elif display_option == 1:
            newimage = work_image[655:2977, 655:2977]
        else:
            newimage = work_image[792:2840, 792:2840]
        self.imagewin.image = newimage


    def select_directory(self, fileentry, pathentry):
        """
        Query for a file path and extract the name/directory.

        Parameters
        ----------

        fileentry:  A tkinter entry field where the file name will be 
                    displayed.
        pathentry:  A tkinter entry field where the directory path will be 
                    displayed.

        Returns
        -------

        None
        """
        filename = tkinter.filedialog.askopenfilename()
        if filename is None:
            return
        values = filename.split('/')
        if len(values) == 0:
            dirpath = None
        else:
            dirpath = filename.replace(values[-1], '')
            filename = values[-1]
        if (pathentry is not None) and (dirpath is not None):
            general_utilities.put_value(dirpath, pathentry)
        if fileentry is not None:
            general_utilities.put_value(filename, fileentry)


    def mean_position(self, filename):
        """
        Calculate the mean sky position from a Mirage star or galaxies input 
        catalogue file.  RAs and Decs must be in columns 2 and 3 (1 and 2 in 
        python notation).

        Parameters
        ----------

        filename:    a string variable, the file name to read

        Returns
        -------

        meanra:   a float value, the mean of column 2 of the file, or None if 
                  there is an issue

        meandec:  a float value, the mean of column 3 of the file, or None if 
                  there is an issue
        """
        try:
            ra = numpy.loadtxt(
                filename, comments=('#', 'index'), usecols=(1, ))
            dec = numpy.loadtxt(
                filename, comments=('#', 'index'), usecols=(2, ))
            meanra = numpy.mean(ra)
            meandec = numpy.mean(dec)
            return meanra, meandec
        except:
            return None, None


    def run_step(self, forward):
        """
        Move the angle one step forwards or backwards.

        Parameters
        ----------

        forward:   a Boolean value for whether to offset forwards; if False, 
                   go backwards

        Returns
        -------

        None
        """
        try:
            angmin = float(self.min_angle_entry.get())
            angmin = general_utilities.range_check(angmin)
            angmax = float(self.max_angle_entry.get())
            angmax = general_utilities.range_check(angmax)
            angstep = float(self.angle_step_entry.get())
            nsteps = int((angmax - angmin)/angstep)
            if (angmin >= angmax) or (angstep <= 0.01) or (nsteps < 1):
                raise ValueError
        except:
            str1 = 'Error found in the "step" range parameters.'
            str1 = str1 + '  Please check your inputs.\n'
            general_utilities.put_message(
                self.message_area, str1)
            return
        current_angle = self.anglevar.get()
        if (current_angle < angmin) or (current_angle > angmax):
            self.anglevar.set(angmin)
            current_angle = angmin
        else:
            if forward:
                current_angle = current_angle+angstep
            else:
                current_angle = current_angle-angstep
            if current_angle < angmin:
                current_angle = angmax
            if current_angle > angmax:
                current_angle = angmin
            self.anglevar.set(current_angle)
        self.apply_angle(None)


    def run_movie(self):
        """
        Code intended to run the tool in movie mode; currently not functional
        """
        angmin = float(self.min_angle_entry.get())
        angmin = general_utilities.range_check(angmin)
        angmax = float(self.max_angle_entry.get())
        angmax = general_utilities.range_check(angmax)
        angstep = float(self.angle_step_entry.get())
        nsteps = int((angmax - angmin)/angstep)
        for loop in range(nsteps):
            self.run_step()
            if (loop > 0):
                newang = self.anglevar.get()+angstep
                if newang > angmax:
                    return


if __name__ == "__main__":
    # create the window
    root = Tk.Tk()
    root.title('WFSS Overlap Tool')
    wfsstool = WFSSOverlap(root)
    root.mainloop()

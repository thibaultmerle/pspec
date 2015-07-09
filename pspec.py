#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import struct
import argparse as ap

try:
    import pylab as pl
except ImportError:
    raise ImportError('pylab module is not installed.')

try:
    import pyfits as pf
except ImportError:
    raise ImportError('pyfits module is not installed.')

#import pdb; pdb.set_trace()

#----------------------------------------------------------------------
def load_fits(finp, verb, extn):
    '''Convert fits file in ascii file
       extn: extension number in case of multiple fits
    '''
    try:
        hdulist = pf.open(finp)
    except:
        print "Input filename does not exist or is not in FITS format."
        quit(1)

    head = hdulist[extn].header
    flux = pl.array(hdulist[extn].data).flatten()
    hdulist.close()

    try:
        ctype = head["CTYPE1"]
    except:
        print "CTYPE1 is undefined in "+finp
    try:
        cunit = head["CUNIT1"]
    except:
        print "Warning: CUNIT1 is undefined in " + finp, '...',
        cunit = "Angstrom"
    try:
        crpix = head["CRPIX1"]
    except:
        crpix = 1.
    try:
        crval = head["CRVAL1"]
    except:
        print "Warning: CRVAL1 is undefined in " + finp
        quit(1)
    try:
        cdelt = head["CDELT1"]
    except:
        print "Warning: CDELT1 is undefined in " + finp
        try:
            cdelt = head["CD1_1"] # ESO keyword used if CDELT1 is undefined
        except:
            print "CD1_1 is undefined in "+finp
            quit(1)
    try:
        naxis = head["NAXIS1"]
    except:
        naxis = pl.size(flux)

    if ctype not in ["", None, "WAVELENGTH", "WAVELENGTH [A]", "AWAV", "LINEAR", "ANGSTROM"]:
        print "CTYPE1 = "+ctype+" not yet implemented."
        quit(1)

    wave_start = crval - (crpix - 1)*cdelt
    wave_final = wave_start + (naxis-1)*cdelt
    wave = pl.linspace(wave_start, wave_final, naxis)
    wave_delta = wave[1]-wave[0]

    if abs(wave_delta - cdelt) > 1e-8:
        print "wavelength interval problem."
        quit(1)

    if verb:
        print " Number of points   :", naxis
        print " Initial wavelength :", wave[0]
        print " Final wavelength   :", wave[naxis-1]
        print " Step               :", cdelt, wave_delta

    return pl.asarray([wave, flux]).T

#----------------------------------------------------------------------
def cosmic_rejection(x, y, n):
    '''Try to reject cosmic from a spectrum
    '''
    bla = False
    blabla = False

    if n == 0:
        print " Warning: sigma for rejection = 0. Take 0.01."
        n = 0.1

    msk = abs(y-y.mean()) > n * y.std(ddof=1)

    xrej = x[msk]
    yrej = y[msk]

    if bla:
        print " Rejection of points outside", n, "sigma around the mean."
        print " Number of rejected points:", xrej.size, '/', x.size

    if blabla:
        print " Rejected points:"
        print xrej
        print yrej

    msk = pl.asarray([not i for i in msk])

    return msk, xrej, yrej

#----------------------------------------------------------------------
def broadgauss(x, y, sigma):
    '''Gaussian function for broadening
    '''

    bla = True
    plot = False
    c = 299792458.

    if bla:
        print " sigma = ", round(sigma, 4), " km/s"

    sigma = sigma * 1.0e3/c * pl.mean(x)   # sigma in Å

    if bla:
        print " sigma = ", round(sigma, 3), "  Å "

    xk = x - pl.mean(x)

    g = make_gauss(1, 0, sigma)
    yk = [g(i) for i in xk]

    if bla:
        print " Integral of the gaussian function: ", pl.trapz(yk, xk).__format__('5.3')

    if plot:
        pl.figure(2)
        pl.plot(xk, yk, '+-')
        pl.show()
    #if bla: print" size y:", y.size
    y = pl.convolve(y, yk, mode='same')
    #if bla: print" size y:", y.size

    return y/max(y)

#----------------------------------------------------------------------
def make_gauss(N, mu, sd):
    '''Normal gaussian function called by broadgauss function
    '''
    k = N / (sd * pl.sqrt(2*pl.pi))
    s = -1.0 / (2 * sd * sd)

    def f(x):
        return k * pl.exp(s * (x - mu)*(x - mu))
    return f

#----------------------------------------------------------------------
def linear_norm(x, y, msk, eps=0.003, deps=0.001, nmin=2, nwin=3):
    '''Linear normalization of a slice of a spectra,
       assuming that the slice is centered on the line to normalized.
    '''

    bla = False
    blabla = False

    x = x[msk]
    y = y[msk]

    n = int((len(y)/2.))
    yl = y[:n]
    yr = y[n+1:]

    # Criteria on the left of the central wavelength
    epsl, epsr = eps, eps
    while 1:
        critl = abs(max(yl)-yl) / max(yl)
        idx_yl = pl.where(critl <= epsl)[0]
        idx_yl = idx_yl.astype(int)
        if blabla:
            print " epsl:", epsl
            print " idx_yl, yl:", idx_yl, [y[i] for i in idx_yl]
        if pl.size(idx_yl) >= nmin:
            break
        else:
            epsl = epsl + deps

   # Criteria on the right of the central wavelength
    while 1:
        critr = abs(max(yr)-yr) / max(yr)
        idx_yr = pl.where(critr <= epsr)[0] + n
        idx_yr = idx_yr.astype(int)
        if blabla:
            print " epsr:", epsr
            print "idx_yr, yr:", idx_yr, [y[i] for i in idx_yr]
        if pl.size(idx_yr) >= nmin:
            break
        else:
            epsr = epsr + deps

    idx_y = pl.concatenate([idx_yl, idx_yr])

    if bla:
        print " nmin, nwin =", nmin, nwin
        print " Number of selected left continuum points:  ", idx_yl.size, "/", n
        print " Number of selected right continuum points: ", idx_yr.size, "/", n
        print " Number of selected continuum points:       ", idx_y.size, "/", y.size

    xs = [x[i] for i in idx_y]
    ys = [y[i] for i in idx_y]

    xs, ys = pl.asarray(xs), pl.asarray(ys)
    n_xs = xs.size

    # Mean value around selected points
    for ind, val in enumerate(ys):
        i = idx_y[ind] - nwin
        j = idx_y[ind] + nwin
        if i < 0:
            i = 0
        if j > len(y):
            j = len(y)
        ys[ind] = y[i:j].mean()

    if blabla:
        print "xs, ys", xs, ys

    A = pl.concatenate([xs, pl.ones(n_xs)])
    A = A.reshape((2, n_xs))
    w = pl.linalg.lstsq(A.T, ys)[0]

    # Test if one of the value of w is a nan
    if pl.any(w != w):
        print "Pb with linalg.lstsq. Try to reduce eps or nmin."
        quit(1)

    a, b = w[0], w[1]

    if blabla:
        print "a =", a, "b =", b

    return a, b, xs, ys

#----------------------------------------------------------------------
def get_data(ifname):
    '''Get spectrum for a given input filename
    '''

    bla = False

    i = ifname.rfind('.')

    if i == -1:
        name = ifname
        ext = ''
    else:
        name = ifname[:i]
        ext = ifname[i+1:]

    if bla:
        print name, ext

    if ext in ['dat', 'txt', 'convol', 'spec', 'tmp', 'cvl', '']:
        data = pl.loadtxt(ifname)
    elif ext == 'bin':
        data = []
        ifile = open(ifname, 'rb')
        while True:
            record = ifile.read(16)
            if len(record) != 16:
                break
            record = struct.unpack('dd', record)
            data.append(record)
        data = pl.asarray(data)
    elif ext in ['fits', 'fit']:
        data = load_fits(ifname, False, extn)
    else:
        print "Unknown format of input spectra."
        quit(1)

    # Remove possible nan in the input data
    idx, = pl.where(data[:, 1] == data[:, 1]) # idx is a tuple
    data = data[idx, :]
    idx, = pl.where(data[:, 0] == data[:, 0])
    data = data[idx, :]

    return data

#----------------------------------------------------------------------
def uniform_wave(x, y, xmin, xmax, n=6000):
    '''Set an uniform x scale by interpolation of x
    '''

    eps = 0.001 # Tolerance on the constant step
    bla = False

    if x.size >= n:
        print "Warning: the number of wavelengths is greater than the default number of interpolated points."
        print "       : use -nu option to increase the number of interpolated points."
        n = x.size

    print " Number of points:             ", x.size

    dx = x[:x.size-1] - wave[1:x.size]
    dxmax = abs(max(dx))
    dxmin = abs(min(dx))

    if bla:
        print min(dx), max(dx)

    crit = dxmax - dxmin < eps

    if crit:
        print " Uniformed step of:         ", dxmax.__format__('7.4f')
    else:
        print " Non-uniformed step."
        print " Min and max steps:         ", dxmin.__format__('7.4f'), dxmax.__format__('7.4f')

    xlin = pl.linspace(xmin, xmax, n)
    y = pl.interp(xlin, x, y)
    x = xlin

    print " Number of interpolated points:", x.size.__format__('7g')

    return x, y

#----------------------------------------------------------------------
def select_data(i, data, lbdmin, lbdmax, obslist, lbdunitlist, hshiftlist):
    '''Selection of spectrum in a given wavelength range
    '''

    bla = False
    global p

    # Case where wavelength are not in Angstrom but in nm
    if lbdunitlist[i] in ['a', 'Å', 'A', '', None]:
        pass
    elif lbdunitlist[i] == 'nm':
        data[:, 0] = data[:, 0]*10
    else:
        print "Wavelength unit unknown. Try Å|A|a or nm."
        quit(1)

    if hshiftlist[i] in [None, '']:
        hshift = 0.0
    else:
        hshift = float(hshiftlist[i])

    data[:, 0] = data[:, 0] + hshift

    crit = pl.logical_and(data[:, 0] <= lbdmax, data[:, 0] >= lbdmin)

    if not crit.any():
        # Case of 1 spectrum only without any idea of wavelength unit
        #if len(obslist) == 1:
        data[:, 0] = data[:, 0]*10
        crit = pl.logical_and(data[:, 0] <= lbdmax, data[:, 0] >= lbdmin)
        print "wavelenght of the spectrum are outside the selection range."
        print "Maybe wavelengths are in nm ... convert wavelengths in Å and try again..."
        if not crit.any():
        #      print "Wavelength range outside the spectrum."
        #      quit(1)
        # Skip this spectrum
        #else:
            print "Required (or default) wavelength range outside of this spectrum."
            print " Number of wavelength points:", len(data[:, 0]).__format__('11g')
            print "  Lambda min:                ", min(data[:, 0]/10).__format__('9.3f')
            print "  Lambda max:                ", max(data[:, 0]/10).__format__('9.3f')
            return None, None, p

    if bla:
        print " Number of wavelength points:", len(data[:, 0]).__format__('11g')
        print "  Lambda min:                ", min(data[:, 0]).__format__('9.3f')
        print "  Lambda max:                ", max(data[:, 0]).__format__('9.3f')

    x = pl.compress(crit, data[:, 0])
    y = pl.compress(crit, data[:, 1])

    # Case where there is a third column in input ASCII data considered as the absolute flux
    if abs_flux:
        try:
            y = pl.compress(crit, data[:, 2])
        except:
            pass

    if bla:
        print " Number of selected points:  ", len(x).__format__('11g')
        print "  Selected lambda min:       ", min(x).__format__('9.3f')
        print "  Selected lambda max:       ", max(x).__format__('9.3f')
        if hshift != 0.0:
            print " Including a shift of:       ", hshift.__format__('9.3f')

    #Flag to plot spectra
    p = True

    return x, y, p

#----------------------------------------------------------------------
def select_ll(llfile, lbdmin, lbdmax, lbdrange):
    '''Select lines to show on the plot by vertical lines
    '''

    bla = False

    try:
        llfile = pl.loadtxt(llname, dtype='str', comments='#', delimiter='\n')
    except IOError:
        print "Linelist file does not exist."
        quit(1)
        #return None, None, None

    elt_ll = [line[0:9] for line in llfile]
    lbd_ll = [float(line[9:18]) for line in llfile]

    try:
        lgf_ll = [float(line[18:]) for line in llfile]
    except:
        lgf_ll = None

    elt_ll = pl.asarray(elt_ll)
    lbd_ll = pl.asarray(lbd_ll)
    lgf_ll = pl.asarray(lgf_ll)

    crit = pl.logical_and(lbd_ll <= lbdmax, lbd_ll >= lbdmin)

    elt_ll_set = pl.compress(crit, elt_ll)
    lbd_ll_set = pl.compress(crit, lbd_ll)

    try:
        lgf_ll_set = pl.compress(crit, lgf_ll)
    except:
        lgf_ll_set = None

    nid = lbd_ll_set.size
    dlbd_ll_set = lbd_ll_set[1:] - lbd_ll_set[0:nid-1]
    dlbd_mean = pl.mean(dlbd_ll_set)

    if bla:
        print ""
        print "  Number of line identification in the selected range:", nid
        print "  Number of identification/Å:", (nid/lbdrange).__format__('7.2f')
        print "  Mean interval [Å]:         ", dlbd_mean.__format__('7.2f')

    loggfmin = -1

    while nid > 100:
        crit = lgf_ll_set > loggfmin
        if not crit:
            break

        elt_ll_set = pl.compress(crit, elt_ll_set)
        lbd_ll_set = pl.compress(crit, lbd_ll_set)
        lgf_ll_set = pl.compress(crit, lgf_ll_set)

        nid = lbd_ll_set.size

        dlbd_ll_set = lbd_ll_set[1:] - lbd_ll_set[0:nid-1]
        dlbd_mean = pl.mean(dlbd_ll_set)

        if bla:
            print "Number of line with log gf >", loggfmin, "to display:", nid
            print "Number of identification/Å:", nid/lbdrange
            print "Mean interval [Å]:", dlbd_mean

        loggfmin = loggfmin + 0.2

    return elt_ll_set, lbd_ll_set, dlbd_mean

#----------------------------------------------------------------------
def get_ilist(fname):
    '''Get input list of spectra filenames
    '''

    ilist = []

    for itm in fname:
    # fname is a list of spectra, itm is a name of a spectrum
        try:
            ifile = get_data(itm)
            ilist.append(itm)
        # fname is a list of list of spectra, itm is a name of a list of spectra
        except ValueError:
            ifile = pl.loadtxt(itm, dtype='str', comments='#', delimiter='\n')

            if ifile.size == 1:
                ilist = ilist + list([str(ifile), ])
            else:
                ilist = ilist + list(ifile)

    #ilist = sorted(list(set(ilist)))
    ilist = [itm for ind, itm in enumerate(ilist) if itm not in ilist[:ind]]

    print "ilist:"
    for i in ilist:
        print i

    return ilist

#----------------------------------------------------------------------
def get_ip(ilist):
    '''Get input parameters from the configuration file
    '''

    bla = False

    # Initialize empty lists for the configuration files
    obslist, namelist, lbdunitlist = [], [], []
    symlist, hshiftlist, savelist = [], [], []
    broadlist, grouplist, alplist = [], [], []
    normlist = []

    #If an input file is given (default name: speclist.txt)
    #Then get the obslist, namelist, lbdunitlist, symlist and hshiftlist
    for itm in ilist:
        i = itm.split('|')
        if not os.path.isfile(i[0].strip()):
            print "File", i[0].strip(), "does NOT exist."
            ilist.remove(itm)
            if len(ilist) == 0:
                quit(1)
            continue
        obslist.append(i[0].strip())
        try:
            namelist.append(i[1].strip())
        except:
            namelist.append(i[0].strip())
        try:
            lbdunitlist.append(i[2].strip())
        except:
            lbdunitlist.append(None)
        try:
            symlist.append(i[3].strip())
        except:
            symlist.append(None)
        try:
            hshiftlist.append(i[4].strip())
        except:
            hshiftlist.append(None)
        try:
            savelist.append(i[5].strip())
        except:
            savelist.append(False)
        try:
            broadlist.append(i[6].strip())
        except:
            broadlist.append(False)
        try:
            grouplist.append(i[7].strip())
        except:
            grouplist.append(False)
        try:
            alplist.append(float(i[8].strip()))
        except:
            alplist.append(0)
        try:
            normlist.append(i[9].strip())
        except:
            normlist.append(None)

    if bla:
        print ""
        print "obslist:"
        for i in obslist:
            print i
        print ""
        print "namelist:"
        for i in namelist:
            print i
        print ""
        print "lbdunitlist:"
        print lbdunitlist
        print "symlist:"
        print symlist
        print "hshiftlist:"
        print hshiftlist
        print "savelist:"
        print savelist
        print "broadlist:"
        print broadlist
        print "grouplist:"
        print grouplist
        print "alplist:"
        print alplist
        print "normlist:"
        print normlist

    return obslist, namelist, lbdunitlist, symlist, hshiftlist, savelist, broadlist, grouplist, alplist, normlist

#----------------------------------------------------------------------
if __name__ == "__main__":

    # Internal parameters
    dft_filename = ['speclist.txt', ]
    dft_lbd = 6569.214    # Default central wavelength on Fe I line in the red wing of Halpha
    dft_lbdrange = 20.001 # Default wavelength range
    dft_dlbd_bdn = 3.0    # Default larger of the border for broadening
    dft_sig_rej = 3.0     # Default sigma for rejecting cosmics
    dft_ln = 1.0         # Y value for the name of lines when -l is used

    ll0name = '/home/tmerle/development/pspec/ll/ll_fraunhofer.dat'
    llname = '/home/tmerle/development/pspec/ll/ll_moore.dat'

    description = 'pspec.py is a versatile command linetool to plot one/several observed/theoretical spectra \
                  \n (desirable unit in Å but also manage nm). \
                  \n Generates a plot output file (default: plot.png).'
    epilog = '2014-10-03 T. Merle pspec.py version 1.2'

    parser = ap.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('filename', nargs='*', default=None, help='Filenames of the spectra or lists of the spectra [fits, ascii or binary]')
    parser.add_argument('-w', '--wavelength', type=float, default=dft_lbd, help='Central wavelength in Angström')
    parser.add_argument('-r', '--range', type=float, default=dft_lbdrange, help='Wavelenght range in Angström')
    parser.add_argument('-wmin', '--wavemin', type=float, default=None, help='Minimum wavelength in Angström')
    parser.add_argument('-wmax', '--wavemax', type=float, default=None, help='Maximum wavelength in Angström')
    parser.add_argument('-ymin', '--ymin', type=float, default=None, help='Minimum ordinate value')
    parser.add_argument('-ymax', '--ymax', type=float, default=None, help='Maximum ordinate value')
    parser.add_argument('-n', '--normalize', action='count', default=None, help='Simple normalization [n constant | nn linear]')
    parser.add_argument('-vs', '--vshift', type=float, default=0., help='Vertical shift for normalize spectra')
    parser.add_argument('-b', '--broadening', type=float, default=None, help='standard deviation of a normalized gaussian function [km/s] for convolution (overides broadening value if given in the input file)')
    parser.add_argument('-l', '--linelist', nargs='?', const=llname, help='Display linelist (default: '+llname+')')
    parser.add_argument('-ln', '--linename', type=float, default=dft_ln, help='Ordinate position of line names when -l option is used (default:'+str(dft_ln)+')')
    parser.add_argument('-nu', '--nb_of_uniform_points', type=float, default=None, help='Number of desired uniform wavelength points')
    parser.add_argument('-nop', '--noplot', action='store_true', default=False, help='No interactive plot')
    parser.add_argument('-nog', '--nogrid', action='store_true', default=False, help='Do not display the grid and the selected continuum points')
    parser.add_argument('-nol', '--nolegend', action='store_true', help='Do not display the legend')
    parser.add_argument('-lin', '--legend_inside', action='store_true', help='Display the legend in the plot')
    parser.add_argument('-nocl', '--no_center_line', action='store_true', default=False, help='Do not display a vertical line at the center of the wavelength range')
    parser.add_argument('-baw', '--black_and_white', action='store_true', default=False, help='Plot in black and white')
    parser.add_argument('-wtf', '--writetofile', action='store_true', default=False, help='write to ascii file')
    parser.add_argument('-o', '--output', default=None, help='Output plot name with extension [png | eps | pdf]')
    parser.add_argument('-t', '--title', type=str, default=None, help='Title of the plot')
    parser.add_argument('-s', '--size', nargs='*', type=float, default=[12, 6], help='Size of the graphic [xsize, ysize] (default=[12, 6])')
    parser.add_argument('-a', '--abs', action='store_true', default=False, help='Use data of a third column considered as absolute flux')
    parser.add_argument('-e', '--ext_number', default=0, type=int, help='Extension number in case of FITS input file')

    print ""

    arguments = parser.parse_args()

    # Line command optionnal arguments
    fname = arguments.filename
    lbd = arguments.wavelength
    lbdmin0 = arguments.wavemin
    lbdmax0 = arguments.wavemax
    ymin = arguments.ymin
    ymax = arguments.ymax
    lbdrange = arguments.range
    norm = arguments.normalize
    vshift = arguments.vshift
    broad = arguments.broadening
    nu = arguments.nb_of_uniform_points
    nop = arguments.noplot
    nog = arguments.nogrid
    nol = arguments.nolegend
    lin = arguments.legend_inside
    nocl = arguments.no_center_line
    baw = arguments.black_and_white
    ll = arguments.linelist
    ln = arguments.linename
    wtf = arguments.writetofile
    output = arguments.output
    title = arguments.title
    figsize = arguments.size
    abs_flux = arguments.abs
    extn = arguments.ext_number

    if lbd != dft_lbd and (lbdmin0 or lbdmax0):
        print "Pb: -w and (-wmin and -wmax) are not compatible"
        quit(1)

    # Default list of spectra
    if not fname:
        fname = dft_filename
        print 'Default file name of list of observations:', fname
        print ''

    # Check the existence of the given line command argument's files
    for itm in fname:
        if not os.path.isfile(itm):
            print "File", itm, "does NOT exist."
            fname.remove(itm)
            if len(fname) == 0:
                quit(1)

    # Get the list [of list] of spectra
    ilist = get_ilist(fname)

    # Get input parameter for each spectrum from the ilist
    dumb = get_ip(ilist)
    obslist, namelist, lbdunitlist, symlist, \
    hshiftlist, savelist, broadlist, grouplist, alplist, normlist = dumb

    #Get the wavelength range if only one spectrum
    if len(ilist) == 1 and lbd == dft_lbd and lbdrange == dft_lbdrange and not lbdmin0 and not lbdmax0:
        a = get_data(obslist[0])
        print pl.shape(a)
        lbdmin0 = min(a[:, 0])
        lbdmax0 = max(a[:, 0])
        lbdrange = lbdmax0 - lbdmin0
        lbd = (lbdmin0 + lbdmax0) / 2.
    elif lbdmin0 and lbdmax0:
        if lbdmin0 > lbdmax0:
            print "Pb: -wmin > -wmax"
            quit(1)
        lbd = (lbdmin0 + lbdmax0) / 2.
        lbdrange = lbdmax0 - lbdmin0
    #Or use specified/default values
    else:
        lbdmin0 = lbd - lbdrange/2.
        lbdmax0 = lbd + lbdrange/2.

   # Check the existence of the linelist file
    if ll:
        if ll != llname:
            llname = ll

  #      if lbdrange > 600.:
  #         llname = ll0name
        if not os.path.isfile(ll):
            print "File", ll, "does NOT exist."
            quit(1)

    print ""
    print "Selected/Default central wavelength:", format(lbd, '7.2f'), "Å"
    print "Selected/Default wavelength range:  ", format(lbdrange, '7.2f'), "Å"

    # Figure size option
    if pl.size(figsize) == 2:
        figsize_x, figsize_y = figsize[0], figsize[1]
    elif pl.size(figsize) == 1:
        figsize_x, figsize_y = figsize[0], figsize[0]
    else:
        print "Maximum of 2 numbers should be given."
        quit(1)

    # Initialization of some plot global parameters
    pl.figure(figsize=(figsize_x, figsize_y), dpi=100, facecolor='w', edgecolor='k')

    try:
        pl.ticklabel_format(useOffset=False)
    except AttributeError:
        pass

    pl.xlim(lbdmin0, lbdmax0)
    pl.xlabel('$\lambda$ [$\AA$]')

    if norm and not vshift:
        pl.ylabel('Normalized flux')
    else:
        pl.ylabel('Flux [Arbitrary Unit]')

    # Grid option
    if not nog:
        pl.grid()

    fluxmin, fluxmax = 0, 0
    p = False # available spectra to plot

    # Enlarge wavelength range for broadening spectra
    if broad:
        lbdmin = lbdmin0 - dft_dlbd_bdn
        lbdmax = lbdmax0 + dft_dlbd_bdn
    else:
        lbdmin = lbdmin0
        lbdmax = lbdmax0

    tracelist = []    # List of plot objects to modify properties if needed
    meanfluxlist = [] # To determine where plot the name of the observation on the plot

    gos = 1 # group of spectra

    # START OF THE LOOP ON EACH SPECTRUM

    for ind, itm in enumerate(obslist):

        # Read of spectrum data
        print "\nSpectrum", ind+1, "/", len(obslist), ":",
        print "read", itm, "...",
        a = get_data(itm)
        print "done"

        # Enlarge wavelength range for specific spectra
        if not broad:
            if broadlist[ind]:
                lbdmin = lbdmin0 - dft_dlbd_bdn
                lbdmax = lbdmax0 + dft_dlbd_bdn
            else:
                lbdmin = lbdmin0
                lbdmax = lbdmax0

        # Wavelength range selection
        print "Wavelength range selection"
        wave, flux, p = select_data(ind, a, lbdmin, lbdmax, obslist, lbdunitlist, hshiftlist)

        # Iterate for spectra not in the required wavelength range
        if not isinstance(wave, pl.ndarray) or not isinstance(flux, pl.ndarray):
            continue

        # Cosmics rejection (points with flux larger than n sigma from the mean
        #                    are rejected for the continuum determination)
        #if norm:
        print "Determination of cosmics and outliers points"
        msk, xrej, yrej = cosmic_rejection(wave, flux, dft_sig_rej)
        #else:
        #   msk = flux != pl.nan

        # Normalization option
        if norm or normlist[ind]:
            print "Selection of continuum points"

        if norm >= 2 or normlist[ind] >= 2:
            if norm > 2 or normlist[ind] > 2:
                print "Warning: norm = ", norm, "not implemented. Do it with norm = 2"
            a, b, xcont, ycont = linear_norm(wave, flux, msk)
            flux = flux/(a*wave + b)
            ycont = ycont/(a*xcont + b)
            # Plot continuum line if no grid
            if not nog:
                pl.plot(xcont, ycont+ind*vshift, '.', color='k')
        elif norm == 1 or normlist[ind] == 1:
            #print " Constant normalization, division by", max(flux[msk])
            flux = flux/max(flux[msk])
            fluxmax = 1.

        # Trace the continuum level for each star
        if norm or normlist[ind]:
            if nog:
                pl.plot(wave, pl.ones(len(wave))+ind*vshift, ':', color='grey')

        # Linear spacing option
        if nu:
            nu = int(nu)
            print "Constant interpolation of wavelength points"
            wave, flux = uniform_wave(wave, flux, lbdmin, lbdmax, n=nu)

        # Broadening option
        if broad:
            print "Change the resolution of the spectrum"
            if not nu:
                wave, flux = uniform_wave(wave, flux, lbdmin, lbdmax)
            flux = broadgauss(wave, flux, broad)
        elif broadlist[ind]:
        # in  [True, 'True', 'T', 'true']:
            print "Smooth the spectrum"
            flux = broadgauss(wave, flux, float(broadlist[ind]))

        # Shift option and group of spectra
        if vshift != 0.0:
            if grouplist[ind]:
                flux = flux + (int(grouplist[ind])-1)*vshift
                fluxmax = max(flux[msk])
                # Useful to know where plot legend entries
                meanfluxlist.append(pl.mean(flux))
            else:
                flux = flux + ind*vshift
                fluxmax = max(flux[msk])
                # Useful to know where plot legend entries
                meanfluxlist.append(pl.mean(flux))
        else:
            fluxmin = min([min(flux[msk]), fluxmin])
            fluxmax = max([max(flux[msk]), fluxmax])

        # Set of symbol for the current spectra
        if symlist[ind] != None:
            sym = str(symlist[ind]) + ''
        else:
            sym = ''

        # Black and white option
        if baw:
            l, = pl.plot(wave, flux, sym, label=namelist[ind], color='k')
        else:
            l, = pl.plot(wave, flux, sym, label=namelist[ind])

        # List of plot objects for modifiying properties
        tracelist.append(l)

        # Save option (post-processed spectrum in ASCII format)
        if savelist[ind] or wtf:
            outdata = pl.asarray([wave, flux]).T
            head = 'pspec.py output of '+itm
            ofname = namelist[ind] + '_' + str(int(lbdmin)) + '_' + str(int(lbdmax)) + '.dat'
            try:
                pl.savetxt(ofname, outdata, fmt='%12.4f %11.4e', header=head)
            except TypeError:
                pl.savetxt(ofname, outdata, fmt='%12.4f %11.4e')
            print "Output file required:", ofname


    # Available spectra to plot
    if not p:
        print ""
        print "No spectra available in the selected wavelength range."
        quit(1)

    #Plot section
    if ymin:
        fluxmin = ymin
    if ymax:
        fluxmax = ymax

    dflux = fluxmax-fluxmin

    if vshift:
        pl.ylim(fluxmin, fluxmax+0.1)
    else:
        pl.ylim(0, 1.2*fluxmax)

    if title:
        pl.title(title)

    # No center line option
    if not nocl:
        pl.axvline(x=lbd, linewidth=1.5, linestyle='-', color='k')

    # Manage of the legend
    if not nol:
        vlegend = 0.5
        #if vshift != 0.0: vlegend = vshift*1.5
        #For python 2.7

        if vshift:
            for ind, itm in enumerate(meanfluxlist):
                if lin:
                    x_shift = -0.1*lbdrange
                    #x_shift = -0.1*(max(wave)-min(wave))
                else:
                    #x_shift = +0.1*(max(wave)-min(wave))
                    x_shift = 0.01*lbdrange
                if alplist[ind]:
                    itm = alplist[ind]
                if symlist[ind]:
                    col = symlist[ind]
                    col = col[0]
                    pl.text(lbdmax+x_shift, itm+0.015, namelist[ind], fontsize=10, color=col)
                else:
                    pl.text(max(wave)+x_shift, itm+0.015, namelist[ind], fontsize=10)

        else:
            try:
                if lin:
                    pl.legend(tracelist[::-1], namelist[::-1], fontsize=8, frameon=False, \
                    loc='lower right',\
                    handletextpad=0.1, labelspacing=vlegend)
                else:
                    pl.legend(tracelist[::-1], namelist[::-1], fontsize=8, frameon=False, \
                    loc='center left', bbox_to_anchor=(1.0, 0.5),\
                    handletextpad=0.1, labelspacing=vlegend)
            #For python 2.6
            except:
                if lin:
                    pl.legend(tracelist[::-1], namelist[::-1], \
                    loc='lower right', bbox_to_anchor=(1.0, 0.5),\
                    handletextpad=0.1, labelspacing=vlegend)    
                else:
                    pl.legend(tracelist[::-1], namelist[::-1], \
                    loc='center left', bbox_to_anchor=(1.0, 0.5),\
                    handletextpad=0.1, labelspacing=vlegend)

    #Plot linelist
    ndiv_tot = 15

    if ll:
        print ""
        print "Overplot line identification:"
        print " Read", llname, "...",
        elt_ll, lbd_ll, dlbd_mean = select_ll(llname, lbdmin, lbdmax, lbdrange)
        print "done"

        lbd_old = 0
        fluxmax = max(flux[msk])

        if vshift:
            ylab = dflux/ndiv_tot
        else:
            ylab = (1+1./ndiv_tot)*fluxmax

        for i, elt in enumerate(elt_ll):
            if abs(lbd_ll[i]-lbd_old) < dlbd_mean:
                ylab = ylab + 1./ndiv_tot*fluxmax
            else:
                ylab = dflux/ndiv_tot

            pl.axvline(x=lbd_ll[i], linewidth=0.5, linestyle='-', color='0.75')
            pl.text(lbd_ll[i], ln, elt, rotation=90, fontsize=8, clip_on=True) #+str(lbd_ll_set[i])
            lbd_old = lbd_ll[i]

    print ""

    if nop:
        print "Computed spectra with:"
    else:
        print "Display spectra with:"

    if norm == 1:
        print " - constant normalization"
    if norm == 2:
        print " - linear normalization"
    if nu:
        print " - constant spacing for", nu, "wavelength points"
    if broad:
        print " - convolution by a gaussian function with sigma =", broad, "km/s"
    if vshift:
        print " - vertically shift of the spectra by ", vshift
    ngroup = len(set(grouplist))
    if ngroup > 1:
        print " -", ngroup, "groups of spectra"
    if ll:
        print " - line identification"
    if not norm and not broad and not vshift and not ll and not nu:
        print " - no options"
    print ""

    #Output section
    if output:
        ofname, ext = os.path.splitext(output)
        ext = ext[1:]
        if ext in ['png', 'eps', 'pdf']:
            ofname = ofname + '.' + ext
            if ext == 'eps':
                pl.savefig(ofname, dpi=100, format=ext, orientation='landscape', bbox_inches='tight')
            else:
                pl.savefig(ofname, dpi=100, format=ext, orientation='landscape')
    else:
        pl.savefig('plot.png', dpi=100, format='png', orientation='landscape', bbox_inches='tight')

    if nop:
        pl.draw()
    else:
        try:
            pl.show(not nop)
        except TypeError:
            print "nop option does not work for this version."
            pl.show()

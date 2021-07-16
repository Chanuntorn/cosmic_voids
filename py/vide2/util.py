

#various codes needed


import os,pylab, pdb
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc

#from .plotDefs
LIGHT_SPEED = 299792.458 

colorList = ['r', 'b', 'g', 'y', 'c', 'm', 
             'darkred', 'grey',
             'orange', 'darkblue',
             'indigo', 'lightseagreen', 'maroon', 'olive',
             'royalblue', 'palevioletred', 'seagreen', 'tomato',
             'aquamarine', 'darkslateblue',
             'khaki', 'lawngreen', 'mediumorchid',
             'orangered', 'thistle'
             'yellowgreen']

linewidth = 4
fontsize = 12


#from vide.plotUtil
def fill_between(x, y1, y2=0, ax=None, **kwargs):
    """Plot filled region between `y1` and `y2`.

    This function works exactly the same as matplotlib's fill_between, except
    that it also plots a proxy artist (specifically, a rectangle of 0 size)
    so that it can be added it appears on a legend.
    """
    ax = ax if ax is not None else pylab.gca()
    ax.fill_between(x, y1, y2, interpolate=True, **kwargs)
    p = pylab.Rectangle((0, 0), 0, 0, **kwargs)
    ax.add_patch(p)


#from vide.voidUtil
def getArray(objectList, attr):

  if hasattr(objectList[0], attr):
    ndim = np.shape( np.atleast_1d( getattr(objectList[0], attr) ) )[0]
    attrArr = np.zeros(( len(objectList), ndim ))

    for idim in range(ndim):
      attrArr[:,idim] = np.fromiter((np.atleast_1d(getattr(v, attr))[idim] \
                                    for v in objectList), float )

    if ndim == 1: attrArr = attrArr[:,0]

    return attrArr
  else:
    print(" Attribute", attr, "not found!")
    return -1


#from xcorlib.py

#CIC interpolation: cloud-in-cell
#first order (linear) weighting scheme
#moves cloud into a cell

def cic(x, Lbox, Lboxcut = 0, Nmesh = 128, weights = None):
    if weights == None: weights = 1
    wm = np.mean(weights)
    ws = np.mean(weights**2)

    #np.mod is division with remainder
    d = np.mod(x/(Lbox+2*Lboxcut)*Nmesh,1)

    box = ([Lboxcut,Lbox+Lboxcut],[Lboxcut,Lbox+Lboxcut],[Lboxcut,Lbox+Lboxcut])

    #histogramdd is histogram of smth multi-dimensional
    #np.roll shifts the value of an array
    rho = np.histogramdd(x, range = box, bins = Nmesh, weights = weights*(1-d[:,0])*(1-d[:,1])*(1-d[:,2]))[0] \
        + np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*d[:,0]*(1-d[:,1])*(1-d[:,2]))[0],1,0) \
        + np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*(1-d[:,0])*d[:,1]*(1-d[:,2]))[0],1,1) \
        + np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*(1-d[:,0])*(1-d[:,1])*d[:,2])[0],1,2) \
        + np.roll(np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*d[:,0]*d[:,1]*(1-d[:,2]))[0],1,0),1,1) \
        + np.roll(np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*d[:,0]*(1-d[:,1])*d[:,2])[0],1,0),1,2) \
        + np.roll(np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*(1-d[:,0])*d[:,1]*d[:,2])[0],1,1),1,2) \
        + np.roll(np.roll(np.roll(np.histogramdd(x, range = box, bins = Nmesh, weights = weights*d[:,0]*d[:,1]*d[:,2])[0],1,0),1,1),1,2)

    rho /= wm

    rho = rho/rho.mean() - 1.

    return (rho, wm, ws)

 
# Power spectra & correlation functions
def powcor(d1, d2, Lbox, Nbin = 10, scale = 'lin', cor = False, cic = True, dim = 1):

    Nmesh = len(d1)


    # CIC correction
    if cic:
        wid = np.indices(np.shape(d1))
        wid[np.where(wid >= Nmesh/2)] -= Nmesh
        wid = wid*np.pi/Nmesh + 1e-100
        wcic = np.prod(np.sin(wid)/wid,0)**2

    # Shell average power spectrum
    dk = 2*np.pi/Lbox
    #np.conj returns the complex conjugate
    Pk = np.conj(d1)*d2*(Lbox/Nmesh**2)**3
    if cic: Pk /= wcic**2

    #see function below
    (Nm, km, Pkm, SPkm) = shellavg(np.real(Pk), dk, Nmesh, Nbin = Nbin, xmin = 0., xmax = Nmesh*dk/2, scale = scale, dim = dim)

    # Inverse Fourier transform and shell average correlation function
    if cor:
        if cic: Pk *= wcic**2 # Undo cic-correction in correlation function
        dx = Lbox/Nmesh
        
        # inverse of the N-dimensional discrete Fourier Transform
        Xr = np.fft.irfftn(Pk)*(Nmesh/Lbox)**3

        (Nmx, rm, Xrm, SXrm) = shellavg(np.real(Xr), dx, Nmesh, Nbin = Nbin/2, xmin = dx, xmax = 140., scale = scale, dim = dim)

        return ((Nm, km, Pkm, SPkm),(Nmx, rm, Xrm, SXrm))


    else: return (Nm, km, Pkm, SPkm)




# Shell averaging
def shellavg(f, dx, Nmesh, Nbin = 10, xmin = 0., xmax = 1., scale = 'lin', dim = 1):
    x = np.indices(np.shape(f))
    x[np.where(x >= Nmesh/2)] -= Nmesh
    f = f.flatten()

    if scale == 'lin': bins = xmin+(xmax-xmin)* np.linspace(0,1,int(Nbin+1))
    if scale == 'log': bins = xmin*(xmax/xmin)**np.linspace(0,1,int(Nbin+1))

    #pdb.set_trace()
    
    if dim == 1: # 1D
        x = dx*np.sqrt(np.sum(x**2,0)).flatten()
        Nm = np.histogram(x, bins = bins)[0]
        xm = np.histogram(x, bins = bins, weights = x)[0]/Nm
        fm = np.histogram(x, bins = bins, weights = f)[0]/Nm
        fs = np.sqrt((np.histogram(x, bins = bins, weights = f**2)[0]/Nm - fm**2)/(Nm-1))
        return (Nm, xm, fm, fs)

    elif dim == 2: # 2D
        xper = dx*np.sqrt(x[0,:,:,:]**2 + x[1,:,:,:]**2 + 1e-100).flatten()
        xpar = dx*np.abs(x[2,:,:,:]).flatten()
        x = dx*np.sqrt(np.sum(x**2,0)).flatten()
        Nm = np.histogram2d(xper, xpar, bins = [bins,bins])[0]
        xmper = np.histogram2d(xper, xpar, bins = [bins,bins], weights = xper)[0]/Nm
        xmpar = np.histogram2d(xper, xpar, bins = [bins,bins], weights = xpar)[0]/Nm
        fm = np.histogram2d(xper, xpar, bins = [bins,bins], weights = f)[0]/Nm
        return (Nm, xmper, xmpar, fm)

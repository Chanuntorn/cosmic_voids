
#Editted Functions
import os, pdb
import numpy as np
import pylab as plt

from vide2.util import *


def plotnumfunction(catalogList,figDir="./",plotName="numberfunc",cumulative=True,binWidth=1):
    '''Plots a cumulative number function. This function takes in a list of void catalogs (or just one), an output directory for the figures made, a prefix for the plot titles, if cumulative is True it plots cumulative number function, and binwidth is the width of histograms in Mps/h. The function returns ellipDistList which is an array of len(catalogList), each element has array of size bins +/- 1 sigma.'''

    print("Plotting number function")

    #returns a list of arrays
    catalogList = np.atleast_1d(catalogList)


    plt.clf()
    plt.xlabel("$R_{eff}$ [$h^{-1}Mpc$]", fontsize=14)


    if cumulative:
        plt.ylabel(r"log ($n$ (> R) [$h^3$ Gpc$^{-3}$])", fontsize=14)

        ellipDistList = []
    else:
        plt.ylabel(r"log ($dn/dR$ [$h^3$ Gpc$^{-3}$])", fontsize=14)

        ellipDistList = []

    for (iSample,catalog) in enumerate(catalogList):
        sample = catalog.sampleInfo
        data = getArray(catalog.voids, 'radius')

        if sample.dataType == "observation":
            maskFile = sample.maskFile
            boxVol = vp.getSurveyProps(maskFile,sample.zBoundary[0],sample.zBoundary[1],sample.zRange[0],
                                       sample.zRange[1], "all",selectionFuncFilee=sample.selFunFile)[0]
            #print(boxVol)
            
        else:
            boxVol = sample.boxLen*sample.boxLen*(sample.zBoundaryMpc[1] - sample.zBoundaryMpc[0])
            #print(boxVol)


        boxVol *= 1.e-9 # Mpc->Gpc
        bins = int(100/binWidth)
        hist, binEdges = np.histogram(data, bins=bins, range=(0., 100.))
        binCenters = 0.5*(binEdges[1:] + binEdges[:-1])


        if cumulative:
            foundStart = False
            for iBin in range(len(hist)):
                if not foundStart and hist[iBin] == 0:
                    continue
                foundStart = True
                hist[iBin] = np.sum(hist[iBin:])

        nvoids = len(data)
        var = hist * (1. - hist/nvoids)
        sig = np.sqrt(var)
        lowerbound = hist - sig
        upperbound = hist + sig

        mean = np.log10(hist/boxVol)
        lowerbound = np.log10(lowerbound/boxVol)
        upperbound = np.log10(upperbound/boxVol)
        lineColor = colorList[iSample]
        lineTitle = sample.fullName
        

       # print(lowerbound)
       # print(upperbound)

        trim = (lowerbound > .01)
        mean = mean[trim]
        binCentersToUse = binCenters[trim]
        lower = lowerbound[trim]
        upper = upperbound[trim]

        alpha = 0.55
        fill_between(binCentersToUse, lower, upper,
                     color=lineColor, alpha=alpha,)
        lineStyle = '-'
        plt.plot(binCentersToUse, mean, lineStyle,
                 color=lineColor,label=lineTitle,
                linewidth=3)

        #print(lineTitle)
        
        #one, two = ax.get_legend_handles_labels()
             
        ellipDistList.append((binCentersToUse, mean, lower, upper))

    #pdb.set_trace()
    plt.legend(loc = "upper right", fancybox=True, prop={'size':14})
  
    #plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
    #plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
    plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

    return ellipDistList

#currently does not show plos in notebook
def computeXcor(catalog,figDir="./",Nmesh = 256,Nbin = 100):
    '''Computes and plots void-void and void-matter(galaxy) correlations. Catalog: catalog to analyze,
figDir: where to place plots, Nmesh: number of grid cells in cic mesh-interpolation, Nbin: number of bins in final plots'''
    # Parameters
    Lbox = catalog.boxLen[0] # Boxlength
    Lboxcut = 0
    Lbox -= 2*Lboxcut

    # Input particle arrays of shape (N,3)
    #particle position?
    xm = catalog.partPos # Halos / Galaxies / Dark matter
    xv = getArray(catalog.voids, 'macrocenter')

    # Interpolate to mesh
    dm, wm, ws = cic(xm, Lbox, Lboxcut = Lboxcut, Nmesh = Nmesh, weights = None)
    dv, wm, ws = cic(xv, Lbox, Lboxcut = Lboxcut, Nmesh = Nmesh, weights = None)

    # N-dimensional discrete Fourier transform
    dmk = np.fft.rfftn(dm)
    dvk = np.fft.rfftn(dv)

    # 1D Power spectra & correlation functions
    ((Nm, km, Pmm, SPmm),(Nmx, rm, Xmm, SXmm)) = powcor(dmk, dmk, Lbox, Nbin, 'lin', True, True, 1)
    ((Nm, km, Pvm, SPvm),(Nmx, rm, Xvm, SXvm)) = powcor(dvk, dmk, Lbox, Nbin, 'lin', True, True, 1)
    ((Nm, km, Pvv, SPvv),(Nmx, rm, Xvv, SXvv)) = powcor(dvk, dvk, Lbox, Nbin, 'lin', True, True, 1)

    #pdb.set_trace()

    # Number densities
    nm = np.empty(len(km)) #empty array
    nv = np.empty(len(km))
    nm[:] = len(xm)/Lbox**3
    nv[:] = len(xv)/Lbox**3


    # Plots
    mpl.rc('font', family='serif')
    ms = 2.5
    fs = 16
    mew = 0.1
    margin = 1.2
    kmin = km.min()/margin
    kmax = km.max()*margin
    rmin = rm.min()/margin
    rmax = rm.max()*margin

    # Density fields (projected)
    fig1 = plt.figure(1)
    ax1 = fig1.gca()    
    ax1.imshow(np.sum(dm[:,:,:]+1,2),extent=[0,Lbox,0,Lbox],aspect='equal',cmap='YlGnBu_r',interpolation='gaussian')
    ax1.set_xlabel(r'$x \;[h^{-1}\mathrm{Mpc}]$')
    ax1.set_ylabel(r'$y \;[h^{-1}\mathrm{Mpc}]$')
    ax1.set_title(r'Dark Matter')
    #fig1.savefig(figDir+'/dm.eps', bbox_inches="tight")
    #fig1.savefig(figDir+'/dm.pdf', bbox_inches="tight")
    fig1.savefig(figDir+'/dm.png', bbox_inches="tight")
    #plt.clf()

    #Voids
    fig2 = plt.figure(2)
    ax2 = fig2.gca()    
    ax2.imshow(np.sum(dv[:,:,:]+1,2)/Nmesh,extent=[0,Lbox,0,Lbox],aspect='equal',cmap='YlGnBu_r',interpolation='gaussian')
    ax2.set_xlabel(r'$x \;[h^{-1}\mathrm{Mpc}]$')
    ax2.set_ylabel(r'$y \;[h^{-1}\mathrm{Mpc}]$')
    ax2.set_title(r'Voids')
    #fig2.savefig(figDir+'/dv.eps', bbox_inches="tight") #, dpi=300
    #fig2.savefig(figDir+'/dv.pdf', bbox_inches="tight") #, dpi=300
    fig2.savefig(figDir+'/dv.png', bbox_inches="tight") #, dpi=300
    #plt.clf()

    # Power Spectra
    fig3 = plt.figure(3)
    ax3 = fig3.gca()
    pa ,= ax3.plot(km, Pmm, 'k-o', ms=0.8*ms, mew=mew, mec='k')
    #plt.plot(km, Pmm-1./nm, 'k--', ms=ms, mew=mew)
    ax3.fill_between(km, Pmm+SPmm, abs(Pmm-SPmm), color='k', alpha=0.2)
    pb ,= ax3.plot(km, Pvm, 'm-D', ms=ms, mew=mew, mec='k')
    ax3.plot(km, -Pvm, 'mD', ms=ms, mew=mew, mec='k')
    ax3.fill_between(km, abs(Pvm+SPvm), abs(Pvm-SPvm), color='m', alpha=0.2)
    pc ,= ax3.plot(km, Pvv, 'b-p', ms=1.3*ms, mew=mew, mec='k')

    #plt.plot(km, Pvv-1./nv, 'b--', ms=ms, mew=mew)
    ax3.fill_between(km, Pvv+SPvv, abs(Pvv-SPvv), color='b', alpha=0.2)
    ax3.set_xlabel(r'$k \;[h\mathrm{Mpc}^{-1}]$')
    ax3.set_ylabel(r'$P(k) \;[h^{-3}\mathrm{Mpc}^3]$')
    ax3.set_title(r'Power Spectra')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlim(kmin,kmax)
    #np.floor largest integer `i`, s.t. i <= x and np.ceil is the smallest int s.t. i >= x
    ax3.set_ylim(10**np.floor(np.log10(abs(Pvm[1:]).min()))/margin, max(10**np.ceil(np.log10(Pmm.max())),10**np.ceil(np.log10(Pvv.max())))*margin)
    ax3.legend([pa, pb, pc],['mm', 'vm', 'vv'],loc='best',prop={'size':12})
    #plt.savefig(figDir+'/power.eps', bbox_inches="tight")
    #plt.savefig(figDir+'/power.pdf', bbox_inches="tight")
    fig3.savefig(figDir+'/power.png', bbox_inches="tight")
    #plt.clf()

    #Correlation Function
    fig4 = plt.figure(4)
    ax4 = fig4.gca()
    pa ,= ax4.plot(rm, Xmm, 'k-o', ms=0.8*ms, mew=mew)
    ax4.fill_between(rm, abs(Xmm+SXmm), abs(Xmm-SXmm), color='k', alpha=0.2)
    pb ,= ax4.plot(rm, Xvm, 'm-D', ms=ms, mew=mew)
    ax4.plot(rm, -Xvm, 'mD', ms=ms, mew=mew)
    ax4.fill_between(rm, abs(Xvm+SXvm), abs(Xvm-SXvm), color='m', alpha=0.2)
    pc ,= ax4.plot(rm, Xvv, 'b-p', ms=1.3*ms, mew=mew)
    ax4.plot(rm, -Xvv, 'bp', ms=ms, mew=1.3*mew)
    ax4.fill_between(rm, abs(Xvv+SXvv), abs(Xvv-SXvv), color='b', alpha=0.2)
    ax4.set_xlabel(r'$r \;[h^{-1}\mathrm{Mpc}]$')
    ax4.set_ylabel(r'$\xi(r)$')
    ax4.set_title(r'Correlation Functions')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlim(rmin,rmax)
    ax4.set_ylim(min(10**np.floor(np.log10(abs(Xvm).min())),10**np.floor(np.log10(abs(Xmm).min())))/margin, max(10**np.ceil(np.log10(Xmm.max())),10**np.ceil(np.log10(Xvv.max())))*margin)
    ax4.legend([pa, pb, pc],['mm', 'vm', 'vv'],loc='best',prop={'size':12})
    #plt.savefig(figDir+'/correlation.eps', bbox_inches="tight") 
    #plt.savefig(figDir+'/correlation.pdf', bbox_inches="tight")
    fig4.savefig(figDir+'/correlation.png', bbox_inches="tight")
    #plt.clf()

    return 

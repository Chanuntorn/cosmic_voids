
#editted plotNumber Function
import os
import numpy as np
import pylab as plt

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
    ax = ax if ax is not None else plt.gca()
    ax.fill_between(x, y1, y2, interpolate=True, **kwargs)
    p = plt.Rectangle((0, 0), 0, 0, **kwargs)
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

        print(lowerbound)
        print(upperbound)

        trim = (lowerbound > .01)
        mean = mean[trim]
        binCentersToUse = binCenters[trim]
        lower = lowerbound[trim]
        upper = upperbound[trim]

        alpha = 0.55
        fill_between(binCentersToUse, lower, upper,
               label=lineTitle, color=lineColor,
               alpha=alpha,)
        lineStyle = '-'
        plt.plot(binCentersToUse, mean, lineStyle,
                color=lineColor,
                linewidth=3)
    
        ellipDistList.append((binCentersToUse, mean, lower, upper))

    plt.legend(loc = "upper right", fancybox=True, prop={'size':14})
  
    #plt.savefig(figDir+"/fig_"+plotName+".pdf", bbox_inches="tight")
    #plt.savefig(figDir+"/fig_"+plotName+".eps", bbox_inches="tight")
    plt.savefig(figDir+"/fig_"+plotName+".png", bbox_inches="tight")

    return ellipDistList


#editted plotNumber Function

def plotnumfunction(catalogList,figDir="./",plotName="numberfunc",cumulative=True,binWidth=1):
    '''Plots a cumulative number function. This function takes in a list of void catalogs (or just one), an output directory for the figures made, a prefix for the plot titles, if cumulative is True it plots cumulative number function, and binwidth is the width of histograms in Mps/h. The function returns ellipDistList which is an array of len(catalogList), each element has array of size bins +/- 1 sigma.'''

    print("Plotting number function")

    #returns a list of arrays
    catalogList = np.atleast_1d(catalogList)


    plt.clf()
    plt.xlabel("$R_{eff}$ [$h^{-1}Mpc$]", fontsize=14)


    if cumulative:
        plt.ylabel(r"log ($n$ (> R) [$h^3$ Gpc$^{-3}$])", fontsize=14)
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
        else:
            boxVol = sample.boxLen*sample.boxLen*(sample.zBoundaryMpc[1] - sample.zBoundaryMpc[0])


        boxVol *= 1.e-9 # Mpc->Gpc
        bins = 100/binWidth
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

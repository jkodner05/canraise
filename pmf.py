# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats as sp
from math import log, factorial, fsum
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import ternary

def tp(N):
    return N/log(N)

def combinatorial(n,k):
    return float(factorial(n)) / (factorial(k) * factorial(n-k))

def binom_point(p, n, k):
    nCk = combinatorial(n,k)
    return nCk * (p**k) * ((1-p)**(n-k))
    
def binom_pmf(p, n, minthresh, maxthresh):
    probs = []
    for k in range(0,n+1):
#        print binom_point(p,n,k), sp.binom.pmf(k,n,p)
        probs.append(sp.binom.pmf(k,n,p))
#        probs.append(binom_point(p,n,k))
#    exit()
#    print p, "\t", fsum(probs[0:kthresh+1]), "\t", fsum(probs)
    return fsum(probs[minthresh:maxthresh+1])

def binom_pmfs(n, minthresh, maxthresh):
    masses = []
    step = 0.01
    for p0 in np.arange(0,1+step,step):
        masses.append(binom_pmf(p0,n,minthresh,maxthresh))
    return masses
        

def pmf_test(N1, N2):
    R1 = binom_pmfs(N1,int(N1/log(N1)))
#    R2 = binom_pmfs(N2,int(N2/log(N2)))

    #print R1
    #print R2

#    for i, r1 in enumerate(R1):
#        r2 = R1[i]
#        if r1 > r2:
#            print r1, r2
    #    print r1 > r2


def can_g2(N1, N2):
    thresh1 = int(tp(N1))
    thresh2 = int(N2/tp(N2))
    print thresh1, N2-N1, thresh2
    print thresh1 + (N2-N1) > thresh2


def binom_cdfs(Ntr, Nfl, thetatr, thetafu):
    gtrans_bypnone = []
    gtransonly_bypnone = []
    gfull_bypnone = []
    gnone_bypnone = []

    step = 0.01
    for pnone in np.arange(0,1+step,step):
        gtrans = binom_pmf(pnone,Ntr,0,thetatr)
        gfull = binom_pmf(pnone,Ntr+Nfl,0,thetafu)
 

        gtrans_bypnone.append(gtrans)
        gfull_bypnone.append(gfull)

        gtransonly = 0.0
        for etr in range(thetafu+1 - Nfl, thetatr + 1):
            petr = sp.binom.pmf(etr,Ntr,pnone)
#            petr = binom_point(pnone,Ntr,etr)
            pefl = 0.0
            for efl in range(thetafu+1 - etr, Nfl + 1):
                pefl += sp.binom.pmf(efl,Nfl,pnone)
#                pefl += binom_point(pnone,Nfl,efl)
            gtransonly += (petr*pefl)
        gtransonly_bypnone.append(gtransonly)
        gnone = binom_pmf(pnone,Ntr+Nfl,thetafu+1,Ntr+Nfl) - gtransonly
        gnone_bypnone.append(gnone)
#        print gfull, gnone, gfull+gnone+gtransonly
#        print gfull, 1-(gnone+gtransonly)
#        print pnone, "\t", gtrans, "\t", gtransonly
    return gtrans_bypnone, gtransonly_bypnone, gfull_bypnone, gnone_bypnone



def binom_cdf(pnone, pfull, Ntr, Nfl, thetatr, thetafu):

    ptrans = 1-pnone-pfull
    gfull = 0.0
    gtransonly = 0.0

    for etr in range(0, thetafu +1): 
        petr = sp.binom.pmf(etr,Ntr,pnone)
        pefl = 0.0
        for efl in range(0, thetafu - etr +1):
            pefl += sp.binom.pmf(efl,Nfl,pnone+ptrans)
        gfull += (petr*pefl)

    for etr in range(thetafu+1 - Nfl, thetatr + 1):
        petr = sp.binom.pmf(etr,Ntr,pnone)
        pefl = 0.0
        for efl in range(thetafu+1 - etr, Nfl + 1):
            pefl += sp.binom.pmf(efl,Nfl,pnone+ptrans)
        gtransonly += (petr*pefl)

    gnone = 1- gfull - gtransonly
    return gnone, gfull, gtransonly


def dynamical(pnone_0, Ntr, Nfl, thetatr, thetafu):

    def mix(percentnew, old, new):
        return (percentnew*new) + ((1-percentnew)*old)

    pnone = float(pnone_0)/100
    ptrans = 0.0
    pfull = 1-pnone-ptrans

    updateratio = 0.1
    seq = [(pnone*100,pfull*100,ptrans*100)]

    for i in range(0,100):
        pnone_1, pfull_1, ptrans_1 = binom_cdf(pnone, pfull, Ntr, Nfl, thetatr, thetafu)
        pnone = mix(updateratio, pnone, pnone_1)
        pfull = mix(updateratio, pfull, pfull_1)
        ptrans = mix(updateratio, ptrans, ptrans_1)
        seq.append((pnone*100,pfull*100,ptrans*100))

    figure, tax = ternary.figure(scale=100)
    figure.set_dpi(300)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="black", multiple=5)
#    tax.gridlines(color="grey", multiple=1, linewidth=0.5)

    fontsize = 17
    tax.set_title("Outcomes from around Actuation Point", fontsize=fontsize+2, y=1.03)
    tax.clear_matplotlib_ticks()
    plt.axis('off')
    tax.left_axis_label("% Trans. Raising Input", fontsize=fontsize, offset = 0.14)
    tax.right_axis_label("% Canon. Raising Input", fontsize=fontsize, offset = 0.14)
    tax.bottom_axis_label("% Non-Productive Raising Input", fontsize=fontsize, offset = 0.10)
    tax.ticks(axis='lbr', linewidth=1, multiple=10, offset = 0.025, fontsize=14)
    
#    tax.plot(seq, linewidth=2.0, label="Curve")

    inits = []
    
    for offset in np.arange(-5,6,1):
        for offset2 in np.arange(8,-2,-2):
            pnone = float(pnone_0+offset)/100
            ptrans = 0.0+float(offset2)/100
            pfull = 1-pnone-ptrans
            print pnone, ptrans, pfull

            inits.append((pnone*100,pfull*100,ptrans*100))


            seq = [(pnone*100,pfull*100,ptrans*100)]

            for i in range(0,40):
                pnone_1, pfull_1, ptrans_1 = binom_cdf(pnone, pfull, Ntr, Nfl, thetatr, thetafu)
                pnone = mix(updateratio, pnone, pnone_1)
                pfull = mix(updateratio, pfull, pfull_1)
                ptrans = mix(updateratio, ptrans, ptrans_1)
                seq.append((pnone*100,pfull*100,ptrans*100))

                traj = tax.plot_colored_trajectory(seq, cmap="hsv", linewidth=1.5)

#    tax.scatter(inits, color="black", label = "initializations")

    tax.heatmap({}, scale=[0,40], cmap="gist_rainbow", vmax=40,vmin=0,cbarlabel="iteration number")
    tax._redraw_labels()
    ternary.plt.savefig("dynamical_%s.png" % int(updateratio*100))
    ternary.plt.show()


def binom_cdfs_3way(Ntr, Nfl, thetatr, thetafu):

    step = 0.01
    scale = 100
    
    nonedict = {}
    nonediffdict = {}
    fulldict = {}
    fulldiffdict = {}
    transdict = {}
    transonlydict = {}
    transdiffdict = {}

#    plot_3way(scale, "% Learning Transporent ai-Raising", nonedict, "Reds", "sample_%s.png" % scale, vmin=0)
#    exit()

    for pnone in np.arange(0,1+step,step):
        for pfull in np.arange(0,1+step,step):
            if pnone+pfull > 1.0:
                continue

            ptrans = 1-pnone-pfull
            gnone, gfull, gtransonly = binom_cdf(pnone, pfull, Ntr, Nfl, thetatr, thetafu)

            gtrans = binom_pmf(pnone,Ntr,0,thetatr)
            gfull_old = binom_pmf(pnone,Ntr+Nfl,0,thetafu)

            print pnone, pfull, ptrans, gtransonly
        
#            gnone = 0.0
#            for etr in range(thetatr - Nfl, Ntr + 1):
#                petr = sp.binom.pmf(etr,Ntr,pnone)
#                pefl = 0.0
#                for efl in range(thetatr+1 - etr, Nfl + 1):
#                    pefl += sp.binom.pmf(efl,Nfl,pnone+ptrans)
#                gnone += (petr*pefl)
#            print gnone, gnone_old
            

#            print gfull_new, gnone, gfull_new+gnone+gtransonly

            nonedict[(int(scale*pnone+0.1),int(scale*pfull+0.1))] = gnone*100
            nonediffdict[(int(scale*pnone+0.1),int(scale*pfull+0.1))] = (gnone-pnone)*100

            fulldict[(int(scale*pnone+0.1),int(scale*pfull+0.1))] = gfull*100
            fulldiffdict[(int(scale*pnone+0.1),int(scale*pfull+0.1))] = (gfull-pfull)*100

            transdict[(int(scale*pnone+0.1),int(scale*pfull+0.1))] = gtrans*100
            transonlydict[(int(scale*pnone+0.1),int(scale*pfull+0.1))] = gtransonly*100
            transdiffdict[(int(scale*pnone+0.1),int(scale*pfull+0.1))] = (gtransonly-ptrans)*100

    #        print gfull, 1-(gnone+gtransonly)
    #        print pnone, "\t", gtrans, "\t", gtransonly

    print "% learning gnone"
#    plot_3way(scale, "% Learning Non-Productive /ai/-Raising", nonedict, "Blues", "nonespace_%s.png" % scale, vmin=0)
    print "change learning gnone"
#    plot_3way(scale, "Absolute Change in Non-Productive /ai/-Raising", nonediffdict, "PRGn", "nonediffspace_%s.png" % scale, vmin=-100)

    print "% learning gfull"
#    plot_3way(scale, "% Learning Full ai-Raising", fulldict, "YlOrBr", "fullspace_%s.png" % scale, vmin=0)
    print "change learning gfull"
#    plot_3way(scale, "Absolute Change in Full /ai/-Raising", fulldiffdict, "PRGn", "fulldiffspace_%s.png" % scale, vmin=-100)

    print "% learning gtrans"
    plot_3way(scale, u"% Learning Transparent /aı/-Raising", transonlydict, "Reds", "transonlyspace_%s.png" % scale, vmin=0)
    print "change learning gtrans"
#    plot_3way(scale, u"Absolute Change in Transparent /aı/-Raising", transdiffdict, "PRGn", "transdiffspace_%s.png" % scale, vmin=-100)



def plot_3way(scale, title, data, color, fname, vmax=100, vmin=-100):

#    afig, ax = plt.subplots(figsize=(8,6))
#    figure, tax = ternary.figure(ax, scale=scale)
    figure, tax = ternary.figure(scale=scale)
    figure.set_dpi(300)
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="black", multiple=5)
#    tax.gridlines(color="grey", multiple=1, linewidth=0.5)

    fontsize = 17
    tax.set_title(title, fontsize=fontsize+2, y=1.03)
    tax.clear_matplotlib_ticks()
    plt.axis('off')
    tax.left_axis_label("% Trans. Raising Input", fontsize=fontsize, offset = 0.14)
    tax.right_axis_label("% Canon. Raising Input", fontsize=fontsize, offset = 0.14)
    tax.bottom_axis_label("% Non-Productive Raising Input", fontsize=fontsize, offset = 0.10)
    tax.ticks(axis='lbr', linewidth=1, multiple=10, offset = 0.025, fontsize=14)

    hax = tax.heatmap(data, scale=scale, cmap=color, vmax=vmax,vmin=vmin)
    
#    hax.cbar.ax.tick_params(labelsize=20)

    tax._redraw_labels()
    ternary.plt.savefig(fname)
    ternary.plt.show()




        
def plot_gprobs(gtransonly, gfull, gnone, title, corpussize, labels):

    def plotline(data, c, w):
        line = plt.plot([100*i for i in data])
        plt.setp(line, color=c, linewidth=w)

    def get_label(c, l):
        return mpatches.Patch(color=c, label=l)

    plt.cla()
    fig = plt.figure(dpi=300)
    axes = fig.add_subplot(1,1,1)
    axes.set_ylim([0,100])
    axes.set_xlim([0,len(gtransonly)])
#    plotline(gnone, 'deepskyblue', 3)
#    plotline(gfull, 'gold', 3)
#    plotline(gtransonly, 'red', 3)
    plotline(gnone, '#3d85c6', 3)
    plotline(gfull, '#f1c232', 3)
    plotline(gtransonly, '#e06666', 3)


    plt.ylabel('% Raising Type Learned', fontsize = 17)
    plt.xlabel('% Non-Raising Input', fontsize = 17)
    plt.legend(handles=[get_label("red",labels[0]),get_label("deepskyblue",labels[1]),get_label("gold",labels[2])], fontsize=15, loc="best")
    plt.legend(handles=[get_label("#e06666",labels[0]),get_label("#3d85c6",labels[1]),get_label("#f1c232",labels[2])], fontsize=15, loc="best")
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(title, fontsize = 19)
    plt.tight_layout()
    plt.savefig(title.replace("_","").replace(" ","_").replace(",","").replace("&","-") + corpussize + ".png")
    plt.show()




def calc_g2prob(Ntr, Nfu, title, corpussize, labels, do3way=False):
    Nfl = Nfu-Ntr
    thetatr = int(tp(Ntr))
    thetafu = int(tp(Nfu))
    print "Ntr", Ntr
    print "Nfl", Nfl
    print "thetatr", thetatr


    if do3way:
        print "Making ternary plots. Be patient"
        binom_cdfs_3way(Ntr, Nfl, thetatr, thetafu)
    print "Making 2D plots"
    gtrans_bypnone, gtransonly_bypnone, gfull_bypnone, gnone_bypnone = binom_cdfs(Ntr, Nfl, thetatr, thetafu)
    plot_gprobs(gtransonly_bypnone, gfull_bypnone, gnone_bypnone, title, corpussize, labels)
    print "Making dynamical plot"
    beststart = gtransonly_bypnone.index(max(gtransonly_bypnone))
    print "gtr best chance at", beststart, "% pnone"
    print max(gtransonly_bypnone)
    print gfull_bypnone[beststart]
    print gnone_bypnone[beststart]
    dynamical(beststart, Ntr, Nfl, thetatr, thetafu)

#####################
# ADD NEW EXPS HERE #
#####################
#Ntr equivalent
#Ntr+Nfl equivalent
#calc_g2prob(Ntr, Ntr+Nfu, "Title", "filename extension" % (Ntr, Ntr_Nfu), ["red line label","blue line label","gold line label"], do3way=False)


print "B>=1; ai vs tr"
Ntr = 103
Nai = 122
calc_g2prob(Ntr, Nai, "Probability of Learning Raising", "_Ntrans-%s_Nai-%s" % (Ntr, Nai), ["% learning transparent","% learning non-productive",u"% learning canonical /aı/"], do3way=False)

print "B>=5; ai vs tr"
Ntr = 45
Nai = 53
calc_g2prob(Ntr, Nai, "Probability of Learning Raising", "_Ntrans-%s_Nai-%s" % (Ntr, Nai), ["% learn transparent","% learn non-productive","% learn full ai"], do3way=True)

print "BB>=5; ai vs tr"
Ntr = 69
Nai = 82
calc_g2prob(Ntr, Nai, "Probability of Learning Raising", "_Ntrans-%s_Nai-%s" % (Ntr, Nai), ["% learn transparent","% learn non-productive","% learn full ai"], do3way=False)

print "BB>=1; ai vs tr"
Ntr = 155
Nai = 182
calc_g2prob(Ntr, Nai, "Probability of Learning Raising", "_Ntrans-%s_Nai-%s" % (Ntr, Nai), ["% learn transparent","% learn non-productive","% learn full ai"], do3way=False)







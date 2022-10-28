import requests
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from tess_stars2px import tess_stars2px_function_entry as tess2sec
from astroquery.mast import Tesscut
from astroquery.mast import Catalogs
from astroquery.gaia import Gaia
from astropy.wcs import WCS
import warnings
import sys, os
from astropy.io import fits
from requests.models import HTTPError
import time
import pycentess.maxsect as mxs
from pycentess.stackima.stackima import stackima
from pycentess.stackima.stackima import centroid
from pycentess.stackima.stackima import radec2dx
from pycentess.stackima.stackima import radec2dy
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib import container
import matplotlib


def writedatafits(pix, namefile, wcs):

    header = wcs.to_header()
    hdu = fits.PrimaryHDU(pix, header=header)
    hdu.writeto(namefile, overwrite=True)


def cutffi(ra,dec,sect):
    coord = SkyCoord(ra,dec,unit="deg")
    #if (mag>6.5):
    sizepx = (31,31)
    #else:
    #    sizepx = (99,51)

    hdulist = Tesscut.get_cutouts(coordinates=coord, sector=sect, size=sizepx) #[0]

    #cutout_hdu.info()
    cutout_table = hdulist[0][1].data
    cutout_table.columns
    imagette = np.array(cutout_table['FLUX'])
    quality  = np.array(cutout_table['QUALITY'])
    btjd     = np.array(cutout_table['TIME'])
    btjd_cor = np.array(cutout_table['TIMECORR'])


    wcs = WCS(header=hdulist[0][2].header)



    return wcs, imagette, quality, btjd, btjd_cor

def where_is_the_star(ra,dec,min_sector,max_sector):
    info   = tess2sec(0, ra, dec)
    cam    = info[4]
    ccd    = info[5]
    xraw   = info[6]
    yraw   = info[7]
    sector = info[3]
    if (sector[0]<0):
        print('The star in ({0},{1}) has not been observed in any sector ...'.format(ra,dec))
        return cam, ccd, xraw, yraw, sector

    flag = np.logical_and(sector>=min_sector,sector<=max_sector)
    cam = cam[flag]
    ccd = ccd[flag]
    xraw = xraw[flag]
    yraw = yraw[flag]
    sector = sector[flag]
    if (len(sector)<1):
        print('The star in ({0},{1}) has not been observed in the Sectors {2} -> {3}...'.format(ra,dec,min_sector,max_sector))
        return cam, ccd, xraw, yraw, sector


    return cam, ccd, xraw, yraw, sector

def blockPrint():
    sys.stdout = open(os.devnull, 'w')
# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def neighbors(rat,dect):
    warnings.filterwarnings("ignore")
    query = 'SELECT *, DISTANCE(POINT({0}, {1}), POINT(gaiadr3.gaia_source.ra, gaiadr3.gaia_source.dec)) AS ang_sep FROM gaiadr3.gaia_source WHERE 1 = CONTAINS(POINT({0}, {1}), CIRCLE(ra, dec, 0.0875)) AND gaiadr3.gaia_source.phot_g_mean_mag < 21. ORDER BY ang_sep ASC'.format(rat,dect)
    #print(query)
    blockPrint()
    job = Gaia.launch_job_async(query=query)
    results = job.get_results()
    enablePrint()
    ranei = results['ra'].data
    denei = results['dec'].data
    gnei  = results['phot_g_mean_mag'].data
    bpnei  = results['phot_bp_mean_mag'].data
    rpnei  = results['phot_rp_mean_mag'].data
    bpnei[np.isnan(bpnei)] = 0
    rpnei[np.isnan(rpnei)] = 0
    dsep  = results['ang_sep'].data*3600  
    tnei = gnei-0.00522555*(bpnei-rpnei)**3+0.0891337*(bpnei-rpnei)**2-0.633923*(bpnei-rpnei)+0.0324473
    ok = np.logical_or(bpnei==0, rpnei==0)
    tnei = np.where(ok,gnei-0.430,tnei)

    return ranei, denei, tnei



def wmean(x,w):
    num = np.sum(w*x)
    den = np.sum(w)
    meanw = num/den

    w1 = w/np.sum(w)
    sigma = np.std(x)
    sigmaw = sigma*np.sqrt(np.sum(w1**2))
    return meanw, sigmaw 


def stack(ratar,dectar,period,t0,dur,minSec, maxSec):
    matplotlib.rcParams['legend.handlelength'] = 0
    matplotlib.rcParams['legend.numpoints'] = 1

    cam, ccd, xraw, yraw, sect = where_is_the_star(ratar, dectar, minSec, maxSec)

    info_obs = np.column_stack((cam, ccd, xraw, yraw, sect))
    if (len(sect)<1): 
        return None

    if (sect[0]<0): 
        return None

    NSECTORs = len(sect)
    s = 0
    coord = SkyCoord(ratar, dectar, unit="deg")
    results = Catalogs.query_region(coordinates=coord, radius=0.005, catalog="TIC")
    tmag = results['Tmag'].data
    rat = results['ra'].data
    det = results['dec'].data
    idg = results['ID'].data
    print('Checking the good transits...')


    xcen = np.zeros(NSECTORs)
    ycen = np.zeros(NSECTORs)
    sxcen = np.zeros(NSECTORs)
    sycen = np.zeros(NSECTORs)


    for ls in range(NSECTORs):
        warnings.filterwarnings("ignore")
        wcs, imagette, quality, btjd, btjd_cor = cutffi(ratar,dectar,sect[ls])

        phase   = np.zeros(len(btjd))
        transit = np.zeros(len(btjd))
        ntr = 0
        if (sect[ls]<27):
            cadence = 0.020833043*24
        if (sect[ls]>26):
            cadence = 0.006944526*24
        if (sect[ls]>55):
            cadence = 0.002314814815*24

        for i in range(len(btjd)):
            if ((btjd[i]-t0)>=0):
                phase[i] = (btjd[i]-t0)/period-int((btjd[i]-t0)/period)
            else:
                phase[i] = (btjd[i]-t0)/period-int((btjd[i]-t0)/period) + 1
            if (phase[i]<=0.5):
                phase[i] = phase[i]+0.5
            else:
                phase[i] = phase[i]-0.5
            phase[i] = (phase[i]-0.5)*period*24
            if (abs(phase[i])<cadence/2): 
                ntr = ntr + 1
                npoints = int(dur/2/cadence)
                for k in range(i-2*npoints-1,i+2*npoints+1):
                    transit[k] = ntr
        print('In sector {0} there are {1} transits'.format(sect[ls],ntr))

        NUMTRANSIT = ntr
        ima_in = np.zeros((2*npoints+13,31,31))
        ima_out = np.zeros((2*npoints+13,31,31))
        q_in = np.zeros(2*npoints+13)
        q_in[:] = 99
        q_out = np.zeros(2*npoints+13)
        q_out[:] = 99
        #print(2*npoints+3)

        dx = np.zeros(NUMTRANSIT)
        dy = np.zeros(NUMTRANSIT)
        sdx = np.zeros(NUMTRANSIT)
        sdy = np.zeros(NUMTRANSIT)
        
        dx[:] = -999.9
        dy[:] = -999.9
        sdx[:] = -999.9
        sdy[:] = -999.9


        for nt in range(NUMTRANSIT):
            kout = -1
            kin = -1
            for i in range(len(transit)):
                if ((transit[i]==(nt+1) and abs(phase[i])<=dur/2 and quality[i]<1)):
#                if (np.logical_and(transit[i]==(nt+1),abs(phase[i])<dur/2)):
                    kin = kin + 1
                    ima_in[kin,:,:] = imagette[i,:,:]
                    q_in[kin] = quality[i]

                if ((transit[i]==(nt+1) and abs(phase[i])>=dur/2 and quality[i]<1)):
#                if (np.logical_and(transit[i]==(nt+1),abs(phase[i])>=dur/2)):
                    kout = kout + 1
                    ima_out[kout,:,:] = imagette[i,:,:]
                    q_out[kout] = quality[i]
            #print(kin, kout, sect[ls])
            if (kin>1 and kout>1):
                pixsin = stackima(pix=ima_in, flag=q_in)
                pixsout = stackima(pix=ima_out, flag=q_out)
                pixdiff = pixsout - pixsin
                #filename = 'DIFF{0}.{1}.fits'.format(nt,sect[ls])
                #writedatafits(pixdiff, filename, wcs)
                #filename = 'IN{0}.{1}.fits'.format(nt,sect[ls])
                #writedatafits(pixsin, filename, wcs)
                #filename = 'OUT{0}.{1}.fits'.format(nt,sect[ls])
                #writedatafits(pixsout, filename, wcs)
                xc, yc,sxc,syc = centroid(pix=pixdiff, mag=tmag[0])
                rac, dec = wcs.all_pix2world(xc,yc,1)
                #print('1',xc,yc,rac,dec,ratar,dectar,rat,det)
                dx[nt] = radec2dx(rac,dec,rat[0],det[0])
                dy[nt] = radec2dy(rac,dec,rat[0],det[0])
                sdx[nt] = sxc*21
                sdy[nt] = syc*21
                #print(dx[nt]*3600,dy[nt]*3600, sdx[nt], sdy[nt])
        ok = np.logical_and(dx>-900,dy>-900)
        dxf = dx[ok]*3600
        dyf = dy[ok]*3600
        sdxf = sdx[ok]
        sdyf = sdy[ok]
        
        #print(dxf,dyf,sdxf,sdyf)

        #wx = 1/sdxf**2
        #wy = 1/sdyf**2
        #xcen[ls], sxcen[ls] = wmean(dxf,wx)
        #ycen[ls], sycen[ls] = wmean(dyf,wy)
        
        ok = np.logical_and(abs(dxf)<900, abs(dyf)<900)
        for it in range(10):
            #wx = 1/sdxf**2
            #wy = 1/sdyf**2
            #xcen[ls], sxcen[ls] = wmean(dxf[ok],wx[ok])
            #ycen[ls], sycen[ls] = wmean(dyf[ok],wy[ok])
            xcen[ls] = np.mean(dxf[ok])
            ycen[ls] = np.mean(dyf[ok])
            sxcen[ls] = np.std(dxf[ok])
            sycen[ls] = np.std(dyf[ok])
            ok = np.logical_and(abs(dxf-xcen[ls])<3*sxcen[ls], abs(dyf-ycen[ls])<3*sycen[ls])

    
    ranei, decnei, tnei =  neighbors(rat[0],det[0])
    xnei = np.zeros(len(ranei))
    ynei = np.zeros(len(decnei))
    for i in range(len(ranei)):
        xnei[i] = radec2dx(ranei[i],decnei[i],rat[0],det[0])
        ynei[i] = radec2dy(ranei[i],decnei[i],rat[0],det[0])
    xnei = xnei*3600
    ynei = ynei*3600

    fig = plt.figure(figsize=(12, 12))
    plt.scatter(xnei,ynei,s=9*(22-tnei),color='black')
    plt.xlim([+120,-120])
    plt.ylim([-120,+120])
#    if (len(xcen)<2):
#        plt.errorbar(xcen,ycen,xerr=sxcen,yerr=sycen, fmt='x',color='red',elinewidth=1.0, capsize=3,
#                     label='SECTOR {0}'.format(sect[0]))
#        plt.xlabel("$\Delta \\alpha$* [arcsec]", size=24)
#        plt.ylabel("$\Delta \delta$ [arcsec]", size=24)
#        plt.legend(loc="best", fontsize=19)

    if (len(xcen)>0):
        usector = np.unique(sect)
        nsectors = len(usector)
        l = np.linspace(0.1, 0.9, endpoint=True, num=nsectors)
        colors = plt.get_cmap('rainbow',nsectors)
        newcolors = colors(l)
        for i_s, s in enumerate(usector):
            sel = sect == s
            xx = xcen[sel] 
            yy = ycen[sel] 
            ss = sect[sel][0]
            sxx = sxcen[sel] 
            syy = sycen[sel] 
            color = newcolors[i_s]
            plt.errorbar(xx,yy, xerr=sxx, yerr=syy, color=color, marker='D', ecolor=color, 
                         elinewidth=1.0, capsize=3, label='SECTOR {0}'.format(ss))

            ax = plt.gca()
            handles, labels = ax.get_legend_handles_labels()
            new_handles = []
            for h in handles:
            #only need to edit the errorbar legend entries
                if isinstance(h, container.ErrorbarContainer):
                    new_handles.append(h[0])
                else:
                    new_handles.append(h)

            ax.legend(new_handles, labels,numpoints=1,fontsize=19,loc='upper left')

            #plt.legend(loc="best", fontsize=19)

        plt.xlabel("$\Delta \\alpha$* [arcsec]", size=24)
        plt.ylabel("$\Delta \delta$ [arcsec]", size=24)
    ticid = int(idg[0])
    path = './vetting_tic-{0:011d}.pdf'.format(ticid)
    filepdf = os.path.abspath(path)
    plt.savefig(filepdf)
    plt.show()
#        colors = iter(cm.rainbow(np.linspace(0,1,len(xcen))))
#        c = next(colors)
#        print(xcen,'+/-',sxcen, ycen,'+/-',sycen)
#        plt.errorbar(xcen,ycen,xerr=sxcen,yerr=sycen, fmt='x',color=next(colors))
#        plt.show()


    






        
            


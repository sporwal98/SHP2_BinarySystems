import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit


#L0 = 3.0128#E-28
filts = ['Z','Y','J','H','K']
VegaToAB = [0.528,0.634,0.938,1.379,1.900]

def magToLum(mag, filt):
    idx = -1
    for i in range(len(filts)):
        if filts[i] == filt:
            idx = i
            
    return np.power(10,-0.4*(mag+VegaToAB[idx]-8.90))

def lumWithErrs(mag, magerr, filt):
    L = magToLum(mag,filt)

    idx = -1
    for i in range(len(filts)):
        if filts[i] == filt:
            idx = i

    dL = np.power(10,-0.4*(VegaToAB[idx]-8.90))*(-0.4)*np.power(10,-0.4*mag)*np.log(10)*magerr

    return L, dL

'''
def magToLumOLD(mag):
    return L0*np.power(10,-0.4*mag)
'''

def read(fname):
    df = pd.read_csv(fname)
    #print(df)
    for c in ['Z','Y','J','H','K']:
        df = df[df[c+'APERMAG3']>=0]
    return df


def plotLvT(df,wsa,wvsc,f):
    dft = df.copy(deep = True)
    dft = dft[dft['WSA'] == wsa]
    dft = dft[dft['WVSC'] == wvsc]
    L = magToLum(dft[f+'APERMAG3'])
    #err = -0.4*np.log(10)*L*dft[f+'APERMAG3ERR']
    
    plt.cla()
    plt.plot(dft['TIMEEQ'],L)
    #plt.plot(dft['TIMEEQ'],(L+err))
    #plt.plot(dft['TIMEEQ'],(L-err))

    plt.ylabel('L')
    plt.xlabel('t')
    plt.show()
    return

def plotVTime(df,var,wsa,wvsc):
    dft = df.copy(deep = True)
    dft = dft[dft['WSA'] == wsa]
    dft = dft[dft['WVSC'] == wvsc]

    
    plt.cla()
    plt.plot(dft['TIMEEQ'],dft[var])

    plt.show()
    
    return
    
def combPlot(df, wvsc,method=plt.scatter):
    dft = df.copy(deep  = True)
    dft = dft[dft['WVSC'] == wvsc]
    #dft.sort_values(by = ['SYNSEQNUM'])


    
    plt.cla()
    
    for c in ['Z','Y','J','H','K']:
        method(dft['XVAL'],dft[c+'APERMAG3'],label = c)
    plt.legend()
    plt.xlabel('Phase')
    plt.ylabel('Magnitude')
    #plt.show()
    plt.title('WVSC:'+str(wvsc))
    plt.savefig('plot_mag_wvsc_'+str(wvsc)+'.png')
    return

def combPlotL(df,wvsc,method=plt.scatter):
    dft = df.copy(deep  = True)
    dft = dft[dft['WVSC'] == wvsc]
    #dft.sort_values(by = ['SYNSEQNUM'])


    
    plt.cla()
    
    for c in ['Z','Y','J','H','K']:
        method(dft['XVAL'],magToLum(dft[c+'APERMAG3'],c),label = c)
    plt.legend()
    plt.xlabel('Phase')
    plt.ylabel('Luminosity')
    #plt.show()
    plt.title('WVSC:'+str(wvsc))
    plt.savefig('plot_lum_wvsc_'+str(wvsc)+'.png')
    return

def savegraphs(df):
    for wvsc in np.unique(df['WVSC']):
        print('WVSC: '+str(wvsc))
        combPlot(df,wvsc)
        combPlotL(df,wvsc)


def mlcheck(df, wvsc, spec, deg = 10):

    #DID NOT WORK
    from numpy.polynomial import polynomial
    dft = df.copy(deep = True)
    dft = dft[dft.WVSC == wvsc]



    c = polynomial.polyfit(dft['XVAL'],magToLum(dft[spec+'APERMAG3'],spec), deg)
    
    x = np.linspace(0.0, 1.0, 1000)
    mat = np.zeros((deg+1,len(x)))
    for i in range(deg+1):
        mat[i,:] = x**i



    y = c @ mat
        
    
    
    plt.cla()
    plt.scatter(dft['XVAL'],magToLum(dft[spec+'APERMAG3'],spec),label = spec)
    plt.plot(x,y,label = 'fit')
    plt.legend()
    plt.xlabel('Phase')
    plt.ylabel('Luminosity')
    #plt.show()
    plt.title('WVSC:'+str(wvsc))
    plt.savefig('plot_lum_wvsc_'+str(wvsc)+'.png')

def lumSplitTest(df,WVSC,spec):
    #31,95,98,107?,110,120,121,124?,137,146
    #INCLUDE ERRORS
    dft = df.copy(deep = True)
    dft = dft[dft.WVSC == WVSC]
    lums = magToLum(dft[spec+'APERMAG3'],spec)

    L1,L2 = np.max(lums),np.min(lums)
    L1 = L1-L2
    
    return L1,L2


def printLums(df,WVSC):
    specs = ['Z','Y','J','H','K']
    lambdas = [0.9917,1.0305,1.2483,1.6313,2.2010]#micrometers
    #width = [0.095,0.1,0.16,0.29,0.34]
    plt.cla()
    plt.figure()
    for i in range(5):
        L1,L2 = lumSplitTest(df,WVSC,specs[i])
        print(lumSplitTest(df,WVSC,specs[i]))

        plt.scatter(lambdas[i],L1, label = specs[i]+str(1))
        plt.scatter(lambdas[i],L2,label = specs[i]+str(2))

    plt.show()

def genLumPlots(df):
    #RUN FOR EB ONLY
    specs = ['Z','Y','J','H','K']
    lambdas = [0.9917e-6,1.0305e-6,1.2483e-6,1.6313e-6,2.2010e-6]#micrometers

    for wvsc in np.unique(df.WVSC):
        print(wvsc)
        lum1 = []
        lum2 = []
        for filt in specs:
            lums = lumSplitTest(df,wvsc,filt)
            lum1.append(lums[0])
            lum2.append(lums[1])
        plt.cla()
        plt.scatter(lambdas,lum1,label='Luminosity of Component 1')
        plt.scatter(lambdas,lum2,label='Luminosity of Component 2')
        plt.xlabel('Wavelength')
        plt.ylabel('Luminosity')
        plt.legend()
        plt.savefig('Source_'+str(wvsc)+'_LumVsWaveLen.png')
        
    return

def getLums(df,wvsc):
    #RUN FOR EB ONLY
    #Returns the luminosities of both components for each filter
    specs = filts
    lambdas = [0.9917e-6,1.0305e-6,1.2483e-6,1.6313e-6,2.2010e-6]#micrometers

    #print(wvsc)
    lum1 = []
    lum2 = []
    for filt in specs:
        lums = lumSplitTest(df,wvsc,filt)
        lum1.append(lums[0])
        lum2.append(lums[1])
        
    return lum1,lum2


def shiftlums(df,wvsc, filt):
    #dip fit
    
    dft = df[df.WVSC == wvsc]
    N = len(dft)
    print(N)

    
    lums = dft[filt+'LUM']
    lums = lums.to_numpy()

    phase = dft['XVAL']
    phase = phase.to_numpy()
    minarg= np.argmin(lums)
    phase = np.mod((phase - phase[minarg]+0.5),1)
    sorting = np.argsort(phase)
    lums = lums[sorting]
    phase = phase[sorting]
    minarg = np.argmin(lums)
    print(minarg)
    '''
    upper = minarg + 1
    for i in range(upper,N-1):
        if(lums[upper+1]>lums[upper]):
            upper = upper+1
        else:
            break

    lower = minarg - 1
    for i in range(lower,1,-1):
        if(lums[lower-1]>lums[lower]):
            lower = lower - 1
        else:
            break
    #'''
    #lower,upper = np.argmax(lums[0:minarg]), (minarg + np.argmax(lums[minarg + 1:N]))

    lower = minarg -10
    upper = minarg +12
    print('lower: '+str(lower)+',upper: '+str(upper))


    from numpy.polynomial import polynomial
    w = np.abs(phase - phase[minarg]+0.01)
    w = 1/w
    '''
    sig = np.sqrt(np.var(lums[lower:upper]))
    gauss = stats.norm(0.5,sig)
    '''
    
    coef = polynomial.polyfit(phase[lower:upper+1],lums[lower:upper+1],deg = 2,w = w[lower:upper+1])


    fitx = np.arange(phase[lower],phase[upper], 0.01)
    res = coef[0] + coef[1]*fitx + coef[2]*fitx**2

    
    plt.cla()

    plt.scatter(phase,lums)
    plt.plot(fitx,res,c = 'r')
    plt.show()




def createlumsfile(df,fname = 'LumSplitEB.csv'):
    systems = np.unique(df.WVSC)
    newDF = pd.DataFrame(columns = ['WVSC','Z1','Y1','J1','H1','K1','Z2','Y2','J2','H2','K2'])
    print(newDF)
    
    for wvsc in systems:
        lum1,lum2 = getLums(df,wvsc)
        new_entry = {'WVSC':int(wvsc),'Z1':lum1[0],'Y1':lum1[1],'J1':lum1[2],'H1':lum1[3],'K1':lum1[4],'Z2':lum2[0],'Y2':lum2[1],'J2':lum2[2],'H2':lum2[3],'K2':lum2[4]}
        newDF = newDF.append(new_entry,ignore_index = True)

    print(newDF)

    newDF.to_csv(fname)

#---------------------------------------#
#---------------------------------------#
#---------------------------------------#
#---------------------------------------#
#---------------------------------------#



def combPlotLwithErr(df,wvsc,method=plt.scatter):
    dft = df.copy(deep  = True)
    dft = dft[dft['WVSC'] == wvsc]
    #dft.sort_values(by = ['SYNSEQNUM'])


    
    plt.cla()
    
    for c in ['Z','Y','J','H','K']:
        mag = dft[c+'APERMAG3']
        magerr = dft[c+'APERMAG3ERR']
        L, dL = lumWithErrs(mag, magerr, c)
        method(dft['XVAL'],L,label = c)
        plt.errorbar(dft['XVAL'],L,yerr = dL,fmt = 'none')
    plt.legend()
    plt.xlabel('Phase')
    plt.ylabel('Luminosity / jansky')
    #plt.show()
    plt.title('WVSC:'+str(wvsc))
    plt.savefig('plot_lum_w_err_wvsc_'+str(wvsc)+'.png')
    return

def combPlotMwithErr(df,wvsc,method=plt.scatter):
    dft = df.copy(deep  = True)
    dft = dft[dft['WVSC'] == wvsc]
    #dft.sort_values(by = ['SYNSEQNUM'])


    
    plt.cla()
    
    for c in ['Z','Y','J','H','K']:
        mag = dft[c+'APERMAG3']
        magerr = dft[c+'APERMAG3ERR']
        method(dft['XVAL'],mag,label = c)
        plt.errorbar(dft['XVAL'],mag,yerr = magerr,fmt = 'none')
    plt.legend()
    plt.xlabel('Phase')
    plt.ylabel('Magnitude')
    #plt.show()
    plt.title('WVSC:'+str(wvsc))
    plt.savefig('plot_mag_w_err_wvsc_'+str(wvsc)+'.png')
    return

def savegraphswerr(df):
    for wvsc in np.unique(df['WVSC']):
        print('WVSC: '+str(wvsc))
        combPlotLwithErr(df,wvsc)
        combPlotMwithErr(df,wvsc)
    




#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#s




def lumSplitTestWErr(df,WVSC,spec):
    #31,95,98,107?,110,120,121,124?,137,146
    #INCLUDE ERRORS
    dft = df.copy(deep = True)
    dft = dft[dft.WVSC == WVSC]
    lums = dft[spec+'LUM']
    lums = lums.to_numpy()
    lumerr = dft[spec+'LUMERR']
    lumerr = lumerr.to_numpy()

    maxi, mini = np.argmax(lums),np.argmin(lums)
    L1,L2 = lums[maxi],lums[mini]
    L1 = L1-L2
    L1er = (lumerr[maxi]+lumerr[mini])#/2#CHANGE BACK
    L2er = lumerr[mini]#*2 CHANGE BACK
    
    return np.abs(L1),np.abs(L2),np.abs(L1er),np.abs(L2er)

def genLumPlotsWErr(df):
    #RUN FOR EB ONLY
    specs = ['Z','Y','J','H','K']
    lambdas = [0.9917e-6,1.0305e-6,1.2483e-6,1.6313e-6,2.2010e-6]

    for wvsc in np.unique(df.WVSC):
        print(wvsc)
        lum1 = []
        lum2 = []
        dL1 = []
        dL2 = []
        for filt in specs:
            lums = lumSplitTestWErr(df,wvsc,filt)
            lum1.append(lums[0])
            lum2.append(lums[1])
            dL1.append(lums[2])
            dL2.append(lums[3])
        plt.cla()
        plt.scatter(lambdas,lum1,label='Component 1')
        plt.errorbar(lambdas,lum1,yerr = dL1,fmt = 'none')
        plt.scatter(lambdas,lum2,label='Component 2')
        plt.errorbar(lambdas,lum2,yerr = dL2,fmt = 'none')
        plt.xlabel('Wavelength / m')
        plt.ylabel('Luminosity / jansky')
        plt.legend()
        plt.savefig('Source_'+str(wvsc)+'_LumVsWaveLen.png')

    return


def blackbodyfit(L,dL=None,Nini = 1.e-5,Tini = 15000,errs = False):
    '''input contains L and dL in each WL'''
    lambdas = np.array([0.9917e-6,1.0305e-6,1.2483e-6,1.6313e-6,2.2010e-6])
    params = [Nini,Tini]
    if(errs):
        fit,chiSq,covMat = LSF(modelBB, params, lambdas, L, dL = dL,errs = errs)
        return fit, chiSq,covMat
    else:
        fit, chiSq = LSF(modelBB, params, lambdas, L, dL = dL)
        return fit,chiSq
    #print('FIT: ',fit,chiSq)

    '''
    model = []
    for i in range(lambdas.size):
        expVal = modelBB(lambdas[i],fit[0],fit[1])
        #print('MODEL:',lambdas[i],expVal)
        model.append(expVal)

    
    fig = plt.figure()

    plt.errorbar(1.e6*lambdas, L, yerr=dL, fmt='o',label='both limits (default)')

    plt.plot(1.e6*lambdas,np.array(model),'-')

    plt.xlabel("Wavelength / micrometres")

    plt.show()
    '''
    

def modelBBold(wl, N, T):
    h=6.626e-34

    c=2.998e8

    k=1.380e-23

    expConst=h*c/(wl*k)

    return N*np.power(wl,-5)*1./( np.exp(expConst/T) -1.)

def modelBB(wl,N,T):
    h = 6.63e-34
    c = 3e8
    kB = 1.38e-23
    term1 = 2*h*c*c*np.power(wl,-5)
    term2 = 1/ (np.exp(h*c/(wl*kB*T))-1)
    
    return  N*term1*term2



def LSF(function, params, lamb,L,dL= None, errs = False):
    
    x = lamb
    y = L
    sigma = dL
    fit, _covMat,infoDict,_comment,_rank = curve_fit(function, x, y, p0 = params, sigma = sigma,full_output = True)
    chiSq = np.sum(infoDict['fvec'] **2)

    #print(_covMat)
    #print(fit)

    if(errs):
        return (fit, chiSq, _covMat)
    
    return (fit,chiSq)




import warnings
warnings.filterwarnings("ignore")
def getTemps(df,errs = False):
    lambs = np.array([0.9917e-6,1.0305e-6,1.2483e-6,1.6313e-6,2.2010e-6])
    IDs = np.unique(df.WVSC)
    T1 = []
    T2 = []
    a1 = []
    a2 = []
    chi1 = []
    chi2 = []
    T1er = []
    T2er = []
    for wvsc in IDs:
        #print('Processing wvsc : ',wvsc)
        L1 = []
        L2 = []
        L1er = []
        L2er = []
        for f in filts:
            l1,l2,l1err,l2err = lumSplitTestWErr(df,wvsc,f)
            L1.append(l1)
            L2.append(l2)
            L1er.append(l1err)
            L2er.append(l2err)


        L1 = np.array(L1)
        L2 = np.array(L2)
            
        if(errs):
            fit1,c1,cov1 = blackbodyfit(L1,dL = L1er,errs = errs)
            r1 = np.var(L1-modelBB(lambs,fit1[0],fit1[1]))
            fit2,c2,cov2 = blackbodyfit(L2,dL = L2er,errs = errs)
            r2 = np.var(L2-modelBB(lambs,fit2[0],fit2[1]))

            t1er = cov1*r1
            t2er = cov2*r2
            
            t1er = np.sqrt(np.diag(t1er))
            t2er = np.sqrt(np.diag(t2er))
            T1er.append(t1er[1])
            T2er.append(t2er[1])
        else:
            fit1,c1 = blackbodyfit(L1,dL = L1er)
            fit2,c2 = blackbodyfit(L2,dL = L2er)
        T1.append(fit1[1])
        a1.append(fit1[0])
        chi1.append(c1)
        T2.append(fit2[1])
        a2.append(fit2[0])
        chi2.append(c2)


    if(errs):
        return T1,T2,a1,a2,chi1,chi2,IDs,T1er,T2er
    return T1,T2,a1,a2,chi1,chi2,IDs

def fittedLumWLPlots(df):
    IDs = np.unique(df.WVSC)
    lambdas = np.array([0.9917e-6,1.0305e-6,1.2483e-6,1.6313e-6,2.2010e-6])
    for wvsc in IDs:
        print(wvsc)
        L1= []
        L1er = []
        L2 = []
        L2er = []
        for f in filts:
            l1,l2,l1er,l2er = lumSplitTestWErr(df,wvsc,f)
            L1.append(l1)
            L1er.append(l1er)
            L2.append(l2)
            L2er.append(l2er)

        fit1,c1 = blackbodyfit(L1,dL = L1er)
        fit2,c2 = blackbodyfit(L2,dL = L2er)
        
        model = []
        for i in range(lambdas.size):
            expVal = modelBB(lambdas[i],fit1[0],fit1[1])
            model.append(expVal)
        plt.cla()
        plt.errorbar(1.e6*lambdas, L1, yerr=L1er, fmt='o',label='both limits (default)')
        plt.plot(1.e6*lambdas,np.array(model),'-')
        plt.xlabel("Wavelength / micrometres")
        plt.ylabel('Luminosity / jansky')
        plt.savefig('fitted_lum_wl_plot_comp1_'+str(wvsc)+'.png')

        
        model = []
        for i in range(lambdas.size):
            expVal = modelBB(lambdas[i],fit2[0],fit2[1])
            model.append(expVal)
        plt.cla()
        plt.errorbar(1.e6*lambdas, L2, yerr=L2er, fmt='o',label='both limits (default)')
        plt.plot(1.e6*lambdas,np.array(model),'-')
        plt.xlabel("Wavelength / micrometres")
        plt.ylabel('Luminosity / jansky')
        plt.savefig('fitted_lum_wl_plot_comp2_'+str(wvsc)+'.png')


def fittedPlots(tdf,df):
    IDs = np.unique(tdf.WVSC)
    lambdas = np.array([0.9917e-6,1.0305e-6,1.2483e-6,1.6313e-6,2.2010e-6])
    plotx = np.arange(0.1e-6,5e-6,0.1e-6)
    N1s = tdf['N1']
    T1s = tdf['T1']
    N2s = tdf['N2']
    T2s = tdf['T2']
    
    for i in range(len(IDs)):
        wvsc = IDs[i]
        print(wvsc)
        L1= []
        L1er = []
        L2 = []
        L2er = []
        for f in filts:
            l1,l2,l1er,l2er = lumSplitTestWErr(df,wvsc,f)
            L1.append(l1)
            L1er.append(l1er)
            L2.append(l2)
            L2er.append(l2er)

        fit1 = (N1s[i],T1s[i])
        fit2 = (N2s[i],T2s[i])

        
        model = modelBB(plotx,fit1[0],fit1[1])
        '''
        model = []
        for i in range(lambdas.size):
            expVal = modelBB(lambdas[i],fit1[0],fit1[1])
            model.append(expVal)
        #'''
        plt.cla()
        plt.errorbar(1.e6*lambdas, L1, yerr=L1er, fmt='o',label='both limits (default)')
        plt.plot(1.e6*plotx,np.array(model),'-')
        plt.xlabel("Wavelength / micrometres")
        plt.ylabel('Luminosity / jansky')
        plt.savefig('fitted_lum_wl_plot_comp1_'+str(wvsc)+'.png')


        model = modelBB(plotx,fit2[0],fit2[1])
        '''
        model = []
        for i in range(lambdas.size):
            expVal = modelBB(lambdas[i],fit2[0],fit2[1])
            model.append(expVal)
        #'''
        plt.cla()
        plt.errorbar(1.e6*lambdas, L2, yerr=L2er, fmt='o',label='both limits (default)')
        plt.plot(1.e6*plotx,np.array(model),'-')
        plt.xlabel("Wavelength / micrometres")
        plt.ylabel('Luminosity / jansky')
        plt.savefig('fitted_lum_wl_plot_comp2_'+str(wvsc)+'.png')
    

def calcBoloLs(tdf,wlmin,wlmax, inc = 0.01e-6):
    BL1 = []
    BL2 = []
    N1 = tdf['N1'].to_numpy()
    N2 = tdf['N2'].to_numpy()
    T1 = tdf['T1'].to_numpy()
    T2 = tdf['T2'].to_numpy()

    
    for i in range(len(N1)):
        BL1.append(intLums(N1[i],T1[i],wlmin,wlmax,inc = inc))
        BL2.append(intLums(N2[i],T2[i],wlmin,wlmax,inc = inc))
    
    return BL1, BL2
        
        

def intLums(N, T, wlmin, wlmax,inc = 0.01e-5):
    wls = np.arange(wlmin,wlmax, inc)
    Ls = modelBB(wls,N,T)
        
    return np.sum(Ls)

def dist(p):
    return 1/p

def distErr(p,dp):
    return np.abs(dp/p**2)


def distances(gaia,tdf):
    dists = []
    distEr = []
    IDs = (tdf.WVSC).to_numpy()
    gaiaIDs = gaia.WVSC.to_numpy()
    p = gaia['PARALLAX'].to_numpy()
    dp = gaia['PARALLAX_ERROR'].to_numpy()
    d = dist(p)
    dd = distErr(p,dp)

    for i in range(len(IDs)):
        wvsc = IDs[i]
        if (wvsc in gaiaIDs):
            dists.append((d[gaiaIDs == wvsc])[0])
            distEr.append((dd[gaiaIDs == wvsc]))
        else:
            dists.append(-1.)
            distEr.append(-1.)

    Ds,DErs = [],[]

    for i in range(len(dists)):
        Ds.append(dists[i])
        DErs.append(distEr[i])

    
        
    return Ds, DErs, IDs

def absBolo(tdf):
    dist = tdf['DIST'].to_numpy()
    B1 = tdf['BOLO1'].to_numpy()
    B2 = tdf['BOLO2'].to_numpy()

    B1abs = B1*(dist/10)**2
    B2abs = B2*(dist/10)**2

    return B1abs, B2abs

def spec(T):
    if(T>=30000):
        return 'O'
    elif(10000<=T and T<30000):
        return 'B'
    elif(7500<=T and T<10000):
        return 'A'
    elif(6000<=T and T<7500):
        return 'F'
    elif(5200<=T and T<6000):
        return 'G'
    elif(3700<=T and T<5200):
        return 'K'
    elif(2400<=T and T<3700):
        return 'M'
    else:
        return 'COLD'

def spec_class(tdf):
    Cl1 = []
    Cl2 = []
    T1 = tdf.T1.to_numpy()
    T2 = tdf.T2.to_numpy()

    for i in range(len(T1)):
        Cl1.append(spec(T1[i]))
        Cl2.append(spec(T2[i]))

    return Cl1,Cl2

def msq(S, L):
    L0 = 4.75
    ratio = L/L0
    if(S == 'O'):
        return ratio>=30000
    elif(S == 'B'):
        return 25<=ratio and ratio<30000
    elif(S == 'A'):
        return 5<=ratio and ratio<25
    elif(S == 'F'):
        return 1.5<=ratio and ratio<5
    elif(S == 'G'):
        return 0.6<=ratio and ratio<1.5
    elif(S == 'K'):
        return 0.08<=ratio and ratio<0.6
    elif(S == 'M'):
        return ratio<=0.08
    else:
        return False

def mainsequence(tdf):
    S1 = tdf.SPEC1.to_numpy()
    S2 = tdf.SPEC2.to_numpy()
    L1 = tdf.BOLO1ABS.to_numpy()
    L2 = tdf.BOLO2ABS.to_numpy()

    MS1 = []
    MS2 = []

    for i in range(len(S1)):
        ms1 = msq(S1[i],L1[i])
        ms2 = msq(S2[i],L2[i])
        ms1 = 'MSQ' if ms1 else 'NOT_MSQ'
        ms2 = 'MSQ' if ms2 else 'NOT_MSQ'
        MS1.append(ms1)
        MS2.append(ms2)

    return MS1,MS2



#--------------------------------------------#
#--------------------------------------------#
#--------------------------------------------#
#--------------------------------------------#
#--------------------------------------------#





def getindividuallums(tdf,df):
    IDs = tdf.WVSC.to_numpy()
    Z1,Y1,J1,H1,K1 = [],[],[],[],[]
    Z1e,Y1e,J1e,H1e,K1e = [],[],[],[],[]
    Z2,Y2,J2,H2,K2 = [],[],[],[],[]
    Z2e,Y2e,J2e,H2e,K2e = [],[],[],[],[]
    

    for wvsc in IDs:
        L1s = []
        L2s = []
        L1ers = []
        L2ers = []
        for f in filts:
            L1,L2, L1er,L2er = lumSplitTestWErr(df,wvsc,f)
            L1s.append(L1)
            L2s.append(L2)
            L1ers.append(L1er)
            L2ers.append(L2er)
            
        Z1.append(L1s[0])
        Y1.append(L1s[1])
        J1.append(L1s[2])
        H1.append(L1s[3])
        K1.append(L1s[4])

        Z1e.append(L1ers[0])
        Y1e.append(L1ers[1])
        J1e.append(L1ers[2])
        H1e.append(L1ers[3])
        K1e.append(L1ers[4])

        Z2.append(L2s[0])
        Y2.append(L2s[1])
        J2.append(L2s[2])
        H2.append(L2s[3])
        K2.append(L2s[4])

        Z2e.append(L2ers[0])
        Y2e.append(L2ers[1])
        J2e.append(L2ers[2])
        H2e.append(L2ers[3])
        K2e.append(L2ers[4])

    
    tdf['Z1'] = Z1
    tdf['Y1'] = Y1
    tdf['J1'] = J1
    tdf['H1'] = H1
    tdf['K1'] = K1

    tdf['Z1ER'] = Z1e
    tdf['Y1ER'] = Y1e
    tdf['J1ER'] = J1e
    tdf['H1ER'] = H1e
    tdf['K1ER'] = K1e

    tdf['Z2'] = Z2
    tdf['Y2'] = Y2
    tdf['J2'] = J2
    tdf['H2'] = H2
    tdf['K2'] = K2

    tdf['Z2ER'] = Z2e
    tdf['Y2ER'] = Y2e
    tdf['J2ER'] = J2e
    tdf['H2ER'] = H2e
    tdf['K2ER'] = K2e

    tdf.to_csv('temps_new.csv')




def getTemps_scaled(df,gaia):

    gID = gaia['WVSC'].to_numpy()
    P = gaia['PARALLAX'].to_numpy()
    Ds = 1/P
    
    IDs = np.unique(df.WVSC)
    T1 = []
    T2 = []
    a1 = []
    a2 = []
    chi1 = []
    chi2 = []
    for wvsc in IDs:
        print('Processing wvsc : ',wvsc)
        D = Ds[gID == wvsc]
        if(D == np.inf or not(wvsc in gID)):
            D = 1.
        print(D)
        L1 = []
        L2 = []
        L1er = []
        L2er = []
        for f in filts:
            l1,l2,l1err,l2err = lumSplitTestWErr(df,wvsc,f)
            L1.append(l1)
            L2.append(l2)
            L1er.append(l1err)
            L2er.append(l2err)

        scale = (D/10)**2

        L1 = np.array(L1)*scale
        L2 = np.array(L2)*scale
        L1er = np.array(L1er)*scale
        L2er = np.array(L2er)*scale

        
        
        fit1,c1 = blackbodyfit(L1,dL = L1er)
        fit2,c2 = blackbodyfit(L2,dL = L2er)
        T1.append(fit1[1])
        a1.append(fit1[0])
        chi1.append(c1)
        T2.append(fit2[1])
        a2.append(fit2[0])
        chi2.append(c2)

    return T1,T2,a1,a2,chi1,chi2,IDs

def getTemps_script(scale=0):
    dft = read('data_with_errs.csv')
    dft = dft[dft.TYPE == 'EB']
    gaiat = pd.read_csv('parallax.csv')
    if(scale):
        T1,T2,a1,a2,chi1,chi2,ID = getTemps_scaled(dft,gaia)    
    else:
        T1,T2,a1,a2,chi1,chi2,ID,T1er,T2er = getTemps(dft,errs = True)

    tdft = pd.DataFrame.from_dict({'T1':T1,'T1ERR':T1er,'T2':T2,'T2ERR':T2er,'N1':a1,'N2':a2,'CHI1':chi1,'CHI2':chi2,'WVSC':ID})


    tdft['BOLO1'],tdft['BOLO2'] = calcBoloLs(tdft,0.1e-6,100e-6)
    D,DE, _= distances(gaia,tdft)
    tdft['DIST'],tdft['DISTERR'] = D,DE

    B1A,B2A = absBolo(tdft)
    tdft['BOLO1ABS'],tdft['BOLO2ABS'] = B1A,B2A
    tdft['SPEC1'], tdft['SPEC2'] = spec_class(tdft)
    tdft['MSQ1'],tdft['MSQ2'] = mainsequence(tdft)
    
    return tdft            

def readall(dft = 'data_with_errs.csv',gaiat = 'parallax.csv',tdft = 'temps.csv'):
    dft = read(dft)
    dft = dft[dft.TYPE == 'EB']
    gaiat = pd.read_csv(gaiat)
    tdft = pd.read_csv(tdft)

    return dft,gaiat,tdft


def HR(tdf):
    ID = np.append(tdf.WVSC.to_numpy,tdf.WVSC.to_numpy)
    T1 = tdf.T1.to_numpy()
    T2 = tdf.T2.to_numpy()
    T = np.append(T1,T2,axis = 0)
    S1 = tdf.SPEC1.to_numpy()
    S2 = tdf.SPEC2.to_numpy()
    S = np.append(S1,S2)
    M1 = tdf.MSQ1.to_numpy()
    M2 = tdf.MSQ2.to_numpy()
    M = np.append(M1,M2)

    B1ABS = tdf.BOLO1ABS.to_numpy()
    B2ABS = tdf.BOLO2ABS.to_numpy()
    BABS = np.append(B1ABS,B2ABS)

    L0 = 4.5
    BABS = BABS/L0



    co = np.unique(S)
    cols = []
    for s in S:
        cols.append(np.where(co == s))
    plt.cla()
    plt.scatter(T,BABS,c = cols)
    plt.xlabel('Temperature / K')
    plt.ylabel('Absolute Bolometric Luminosity / L0(Solar Luminosity)')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('SPEC_COL_LUM_VS_T.png')




    cols = []
    for i in range(len(M)):
        cols.append(int(M[i]=='MSQ'))
    
    plt.cla()
    plt.scatter(T,BABS,c = cols)
    plt.xlabel('Temperature / K')
    plt.ylabel('Absolute Bolometric Luminosity / L0(Solar Luminosity)')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('MSQ_COL_LUM_VS_T.png')


    plt.cla()
    plt.hist(T,bins = int(np.sqrt(len(T))))
    plt.xlabel('T / K')
    plt.ylabel('Frequency')
    plt.savefig('T_DIST.png')

    plt.cla()
    plt.hist(S,bins = int(np.sqrt(len(S))))
    plt.xlabel('Spectral Class')
    plt.ylabel('Frequency')
    plt.savefig('S_DIST.png')

    plt.cla()
    plt.hist(M,bins = int(np.sqrt(len(M))))
    plt.xlabel('MSQ or NOT')
    plt.ylabel('Frequency')
    plt.savefig('M_DIST.png')

    

    
    
    
    

    
    
    
#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#

def surface(tdf):
    return






    
    



import numpy as np
from scipy.integrate import cumtrapz as inttrapz

def bin_data(fluxdata, lamdata, errordata, binsize):

    binsize=int(binsize)

    steps=range(0,len(fluxdata),int(binsize))
    #inds=range(0,len(steps)+1,1)
    
    Flux=np.array([])
    Lambda=np.array([])
    Error=np.array([])

    for i in enumerate(steps):
        i=i[0]
    
        
        fpoint=sum(fluxdata[(binsize)*i:(binsize)*(i+1)])/binsize
        Flux=np.append(Flux,fpoint)

        lpoint=sum(lamdata[(binsize)*i:(binsize)*(i+1)])/binsize
        Lambda=np.append(Lambda,lpoint)

        errorpoints=errordata[(binsize)*i:(binsize)*(i+1)]
        epoint=(sum(errorpoints**2)**.5)/binsize
        Error=np.append(Error,epoint)

    difference=(len(fluxdata)/binsize-len(fluxdata)/float(binsize))

    if difference!=0:
        Flux=Flux[0:len(Flux)-2]
        Lambda=Lambda[0:len(Lambda)-2]
        Error=Error[0:len(Error)-2]


    return Flux, Lambda, Error


def spectrum_cut(fluxdata, lamdata, errordata, minl, maxl):

    Flux=[]
    Lambda=[]
    Error=[]
 
    if len(fluxdata)==len(lamdata):

        for i in range(0,len(fluxdata)):

            if minl <= lamdata[i] <= maxl:
                Flux.append(fluxdata[i])
                Lambda.append(lamdata[i])
                Error.append(errordata[i])

        Flux=np.array(Flux)
        Lambda=np.array(Lambda)
        Error=np.array(Error)

    else:
         raise Exception('Flux and Lambda arrays need to be the same length')


    return Flux, Lambda, Error


def column_density(I0, IM, lamdata, lam0, f, star_v):
#column density due to a single ion
#I0 intensity of the star (w/o gas)
#IM intensity measured
#lamdata range over which the column density should be calculated
#lam0 - line of interest
#f - oscillator strength
#star_v - star velocity

    pi = np.pi          #unitless
    m_e = 9.109e-28     #grams
    c = 2.998e10        #cm/s
    echarge = 4.803e-10 #g^.5 cm^(1.5) s^-1
    AtoCM=1e-8

    tau0 = (pi * echarge**2) / (m_e * c) 

    velocity = ((lamdata - lam0) / lam0) * c
    

    #Column Density cm^-2
    N = ((tau0**-1) / (AtoCM * lam0 * f) * inttrapz(np.log(I0/IM), velocity))
    N = N[-1:] #last element of N 

    return velocity, N


def column_density_diff(Flux1, Error1, Flux2, Error2, lamdata, lam0, f, star_v):
#column density due to a single ion
#Flux1 - Intensity day 1
#Flux2 - Intensity day 2
#lamdata range over which the column density should be calculated
#lam0 - line of interest
#f - oscillator strength
#star_v - speed of star relative to Earth. away is positive km/s

#outputs
#velocity cm/s
#column density and error (N) atoms/cm^2


    km2cm = 1e5
    pi = np.pi          #unitless
    m_e = 9.109e-28     #grams
    c = 2.998e10        #cm/s
    echarge = 4.803e-10 #g^.5 cm^(1.5) s^-1
    AtoCM=1e-8

    tau0 = (pi * echarge**2) / (m_e * c) 

    velocity = (((lamdata - lam0) / lam0) * c) - (star_v * km2cm) 


    #Column Density Difference  cm^-2
    N = ((tau0**-1) / (AtoCM * lam0 * f) * inttrapz(np.log(Flux2/Flux1), velocity))
    N = N[-1:] #last element of N 

    #ERROR
    #Error is one sigma
    log_f_error = log_flux_error(Flux1, Error1, Flux2, Error2)
    integral_error = trapz_error(log_f_error , velocity)
    
    N_error = integral_error * ((tau0**-1) / (AtoCM * lam0 * f))
    N_error = np.array(N_error) #There is no reason for this I just like seeing exponents

    return N, N_error, velocity





def column_density_no_approx(I0, IM, lamdata, lam0, f, star_v):
#column density due to a single ion
#I0 intensity of the star (w/o gas)
#IM intensity measured
#lamdata range over which the column density should be calculated
#lam0 - line of interest
#f - oscillator strength
#star_v - star velocity

    pi = np.pi          #unitless
    m_e = 9.109e-28     #grams
    c = 2.998e10        #cm/s
    echarge = 4.803e-10 #g^.5 cm^(1.5) s^-1
    AtoCM=1e-8

    tau0 = (pi * echarge**2) / (m_e * c) 

    velocity = ((lamdata - lam0) / lam0) * c

    #Column Density cm^-2

    #NEED TO FIND MIDPOINTS OF TRAPEZOIDAL RULE
    N = ((tau0**-1) / (AtoCM * lamdata * f) * inttrapz(np.log(I0/IM), velocity)) 

    return velocity, N

def log_flux_error(Flux1, Error1, Flux2, Error2):
#CALCULATES error of that wierd ln(I0/IM) function

    log_f_err = np.array([])
    for j in np.arange(0, len(Flux1), 1):
    
        log_f_errpt = (( Error1[j] / Flux1[j] )**2 + ( Error2[j] / Flux2[j] )**2 )**.5

        log_f_err = np.append(log_f_err, log_f_errpt)

    return log_f_err


def trapz_error(yerror , x):
# Calculating the error through the trapezoidal RULEEEEEZ
# no x error for STIS spectrograph data.
   
    error = 0

    for j in np.arange(0 , len(yerror)-1, 1):

        sigmapoint = ((x[j+1] - x[j]) / 2.0) * ((yerror[j]**2 + yerror[j+1]**2)**.5)
    
        error+=sigmapoint

    return error

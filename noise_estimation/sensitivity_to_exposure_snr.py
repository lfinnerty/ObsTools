from astropy.modeling.physical_models import BlackBody
import astropy.units as u
import numpy as np
import sys

if __name__ == '__main__':
    band = sys.argv[1]

    ### Star info
    sttemp = 5050. * u.K
    strad = 0.82 * u.R_sun
    dist = 19.764*u.pc

    # sttemp = 10000. * u.K
    # strad = 2.7 * u.R_sun
    # dist = 7.68*u.pc

    ### Band info
    if band == 'L':
        wv = 3.8*u.micron
        refflux = 0.15*u.mJy
        reftime = 1.*u.hour
        refsnr = 10.
    elif band == 'M':
        wv = 4.8*u.micron
        refflux = 1.2*u.mJy
        reftime = 1.*u.hour
        refsnr = 10.
    elif band == 'K':
        wv = 2.2*u.micron
        refflux = 0.5*u.mJy
        reftime = 1.*u.hour
        refsnr = 10.
    ### Desired integration info
    obstime = 120*u.second
    
    
    bb = BlackBody(temperature=sttemp)
    flux = bb(wv)

    obsflux = flux * u.sr*np.pi*(strad.to(u.cm))**2/((dist.to(u.cm))**2)
    obsmagAB = -2.5*np.log10(obsflux.to(u.Jy)/(3631*u.Jy))
    if band == 'K':
        obsmagK = -2.5*np.log10(obsflux.to(u.Jy)/(670*u.Jy))
        print('K mag', obsmagK)

    print('Source flux = ',obsflux.to(u.Jy))

    ### Calculate the SNR given reftime and refflux
    snr = refsnr * np.sqrt(obstime.to(u.second)/reftime.to(u.second)) * np.sqrt((obsflux.to(u.Jy)/refflux.to(u.Jy)))
    print('SNR per integration = ', snr)
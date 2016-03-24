# My_IDL_Routines
Some of my routines in IDL 

## The IGM.pro:
Monte Carlo Simulation of Intergalactic Medium absorption along several lines of sight. I followed Inoue and Iwata (2008), with the updated parameters from Inoue et al. (2011). 
I used the provided probability density functions and generated IGM absorbers along each sightline with variety of neutral hydrogen column density NHi, doppler broadening parameter b and redshift z. 
The results are saved in separate files. The routine "tra23.pro" reads these files and creates absorption spectra for each sightline. 
Another routine "av3.35.pro" takes into account all the simulated sightlines and calculates the mean and the median IGM absorption spectrum of a source at redshift of z=zs.


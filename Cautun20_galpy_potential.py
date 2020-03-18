# Cautun (2020) potential 
# Thanks to Thomas Callingham (Durham University, UK) which implemented the potential within galpy

import numpy as np
from galpy.potential import NFWPotential
from galpy.potential import DiskSCFPotential
from galpy.potential import SCFPotential
from galpy.potential import scf_compute_coeffs_axi
from galpy.potential import scf_compute_coeffs_spherical
from galpy.potential import mwpot_helpers
from galpy.util import bovy_conversion
# Suppress the np floating-point warnings that this code generates...
old_error_settings= np.seterr(all='ignore')




ro= 8.122
vo= 229
sigo = bovy_conversion.surfdens_in_msolpc2(vo=vo,ro=ro)
rhoo = bovy_conversion.dens_in_msolpc3(vo=vo,ro=ro)

#Cautun DM halo
fb   = 4.825 / 30.7 # Planck 1 baryon fraction
m200 = 0.97e12  # the DM halo mass
conc = 9.4

#Cautun Bulge
r0_bulge  = 0.075/ro
rcut_bulge= 2.1/ro
rho0_bulge= 103/rhoo

#Cautun Stellar Discs
zd_thin    = 0.3/ro
Rd_thin    =2.63/ro
Sigma0_thin= 731./sigo
zd_thick    = 0.9/ro
Rd_thick    = 3.80/ro
Sigma0_thick= 101./sigo

#Cautun Gas Discs
Rd_HI= 7./ro
Rm_HI= 4./ro
zd_HI= 0.085/ro
Sigma0_HI= 53/sigo
Rd_H2= 1.5/ro
Rm_H2= 12./ro
zd_H2= 0.045/ro
Sigma0_H2= 2200/sigo

# Cautun CGM
A = 0.19
Beta = -1.46
critz0 = ((127.5/(1e9))/rhoo)
R200   = 219/ro #R200 for cgm


def gas_dens(R,z):
    return mwpot_helpers.expsech2_dens_with_hole(R,z,Rd_HI,Rm_HI,zd_HI,Sigma0_HI) + mwpot_helpers.expsech2_dens_with_hole(R,z,Rd_H2,Rm_H2,zd_H2,Sigma0_H2)

def stellar_dens(R,z):
    return mwpot_helpers.expexp_dens(R,z,Rd_thin,zd_thin,Sigma0_thin) + mwpot_helpers.expexp_dens(R,z,Rd_thick,zd_thick,Sigma0_thick)

def bulge_dens(R,z):
    return mwpot_helpers.core_pow_dens_with_cut(R,z,1.8,r0_bulge,rcut_bulge,
                                                rho0_bulge,0.5)

def cgm_dens(R,z):
    r = np.sqrt(R**2+(z**2))
    dens_cgm = 200 * critz0 * A * fb * (r/R200)**Beta 
    if r>R200:
        dens_cgm*=np.exp(1-r/R200)
    return dens_cgm



#dicts used in DiskSCFPotential 
sigmadict = [{'type':'exp','h':Rd_HI,'amp':Sigma0_HI, 'Rhole':Rm_HI},
             {'type':'exp','h':Rd_H2,'amp':Sigma0_H2, 'Rhole':Rm_H2},
             {'type':'exp','h':Rd_thin,'amp':Sigma0_thin, 'Rhole':0.},
             {'type':'exp','h':Rd_thick,'amp':Sigma0_thick, 'Rhole':0.}]

hzdict = [{'type':'sech2', 'h':zd_HI},
          {'type':'sech2', 'h':zd_H2},
          {'type':'exp', 'h':zd_thin},
          {'type':'exp', 'h':zd_thick}]


#generate separate disk and halo potential - and combined potential
Cautun_bulge= SCFPotential(\
    Acos=scf_compute_coeffs_axi(bulge_dens,20,10,a=0.1)[0], a=0.1, ro=ro, vo=vo )

Cautun_cgm= SCFPotential(\
    Acos=scf_compute_coeffs_spherical(cgm_dens,20,a=20)[0], a=20, ro=ro, vo=vo )

Cautun_disk= DiskSCFPotential( dens=lambda R,z: gas_dens(R,z) + stellar_dens(R,z), Sigma=sigmadict, \
                              hz=hzdict, a=2.5, N=30, L=30, ro=ro, vo=vo )


Cautun_unContracted = NFWPotential( conc=conc, mvir=m200/1.e12, vo=vo, ro=ro, H=67.77, Om=0.307, overdens=200.0 * (1.-fb), wrtcrit=True )



# functions for calculating the contraction of the DM halo given the baryonic mass distribution
def Uncontracted_DMDens(R,z): 
    return Cautun_unContracted.dens( R, z, use_physical=False )

def Baryon_Dens(R,z): #Total baryon profile in array friendly format
    TotalDens = gas_dens(R,z) + stellar_dens(R,z) + bulge_dens(R,z) + cgm_dens(R,z)
    return TotalDens


from Cautun20_contraction import potential_contract_mass_profile
Contracted_rho_dm, rspace, MCum_bar, MCum_DM, MCum_DM_contracted, rspace_MCum = \
        potential_contract_mass_profile( Uncontracted_DMDens, Baryon_Dens, f_bar=fb )


def Contracted_DM_dens(R,z):
    r = np.sqrt( (R**2.) + (z**2.) )
    #Ensure intepolation doesn't return errors
    if r<rspace[0]:
        return Contracted_rho_dm[0]
    if r>rspace[-1]: 
        return Contracted_rho_dm[-1] * np.exp( 1 - ((r/rspace[-1])**2) )
    contracted_dm_dens = np.interp( r, rspace, Contracted_rho_dm )
    return contracted_dm_dens

Cautun_halo= SCFPotential(\
    Acos=scf_compute_coeffs_spherical( Contracted_DM_dens,60,a=50 )[0], a=50, ro=ro, vo=vo )
# Go back to old floating-point warnings settings
np.seterr(**old_error_settings)


Cautun20 = Cautun_halo + Cautun_disk + Cautun_bulge + Cautun_cgm

    
    
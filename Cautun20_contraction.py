import numpy as np
from scipy import integrate

def contract_enclosed_mass( mass_DM, mass_bar, f_bar=0.157 ):
    """ Returns the contracted DM profile given the 'uncontracted' profile and that of the baryonic distribution.
   
   Args:
      mass_DM       : enclosed mass in the DM component in the absence of baryons. 
                          It corresponds to '(1-baryon_fraction) * enclosed mass in
                          DMO (dark matter only) simulations'.
      mass_bar      : enclosed baryonic mass for which to calculate the DM profile.
      f_bar         : optional cosmic baryonic fraction.
   Returns:
      Array of 'contracted' enclosed masses.
   """
    eta_bar = mass_bar / mass_DM * (1.-f_bar) / f_bar  # the last two terms account for transforming the DM mass into the corresponding baryonic mass in DMO simulations
    increase_factor = 0.45 + 0.38 * (eta_bar + 1.)**0.53
    return mass_DM * increase_factor


def NFW_enclosed_mass( r, M200, conc ):
    """ Returns the enclosed mass for an NFW profile.
   
   Args:
      r             : radial distance in 'kpc'.
      M200          : halo mass in 'M_solar'.
      conc          : halo concentration w.r.t. halo radius R_200 = average enclosed 
                          density is 200 times the critical one.
   Returns:
      Array of enclosed masses for each 'r' value in units of 'M_solar'.
   """
    G_newton = 4.302e-6        # Newton constant in kpc Msun^-1 (km/s)^2
    H0 = 1.e-1 * h             # Hubble constant in (km/s) kpc^-1
    rho_critical = 3.*H0*H0 / (8.*np.pi * G_newton)
    R200 = (M200 / (4.*np.pi/3. * 200. * rho_critical ))**(1./3.)
    xsi = r / R200 * conc
    return M200 * (np.log(1.+xsi) - xsi/(1.+xsi)) / (np.log(1.+conc) - conc/(1.+conc))


def density_from_enclosed_mass( r_bins, enclosed_mass, out_r_bins ):
    """ Converts an array of enclosed masses to 3D densities.
   
   Args:
      r_bins             : the radial bins at which the enclosed mass is defined.
      enclosed_mass      : array of enclosed masses.
      out_r_bins         : array of radial distances at which the density will be interpolated.
   Returns:
      Array of densities at the location of the output radial distances.
   """
    bin_mass = enclosed_mass[1:] - enclosed_mass[:-1]
    shell_vol= 4.*np.pi/3. * (r_bins[1:]**3 - r_bins[:-1]**3)
    bin_dens = bin_mass / shell_vol
    r_vals = np.sqrt( r_bins[1:] * r_bins[:-1] )
    return np.interp( out_r_bins, r_vals, bin_dens )



def enclosedMass( Func, rbins, N=10000 ):  #Gives the enclosed mass in r for (R,z) func
    r1 = np.logspace( -3, np.log10(rbins[-1]*1.1), N+1 )  # the bins used for the enclosed mass calculation
    r  = np.sqrt( r1[1:] * r1[:-1] )
    dr = r1[1:] - r1[:-1]
    
    I = lambda x, xr: Func( R=xr * np.cos(x), z=xr * np.sin(x) ) * 4*np.pi * xr**2 * np.cos(x) 
    
    shellMass = np.zeros( r.shape[0] )
    for i in range( r.shape[0] ):
        shellMass[i] = integrate.quad( I, 0., np.pi/2, args=( r[i], ) )[0] * dr[i]
    return np.interp( rbins, r, shellMass.cumsum() )

def enclosedMass_spherical( Func, rbins, N=10000 ):  #Gives the enclosed mass in r for spherically symmetric (R,z) func
    r1 = np.logspace( -3, np.log10(rbins[-1]*1.1), N+1 )  # the bins used for the enclosed mass calculation
    r  = np.sqrt( r1[1:] * r1[:-1] )
    dr = r1[1:] - r1[:-1]
    dV = 4*np.pi * r**2 * dr
    
    shellMass = Func( r, 0 ) * dV
    return np.interp( rbins, r, shellMass.cumsum() )



def potential_contract_mass_profile( rho_DMO_func, rho_Baryon_func, f_bar=0.157 ):
    """ Returns the contracted DM density in galpy units for the given baryonic profile.
   
   Args:
      rho_DM_func        : function that gives the 'uncontracted' DM density at coordinates (R,z).
      rho_Baryon_func    : function that gives the baryonic density at coordinates (R,z).
      f_bar              : optional cosmic baryonic fraction.
   Returns:
      Contracted_rho_dm  : contracted DM density at a set of radial distances.
      rvals              : array of radial distances for which the density was calculated.
      MCum_bar           : the enclosed baryonic mass (at a different set of distances from the density).
      MCum_DM            : the enclosed input DM mass. 
      MCum_DM_contracted : the enclosed contracted DM mass.
      rspace             : the radial values for which the enclosed mass was calculated.
   """
    # Create logarithmic r grid
    rspace = np.logspace( -2, 2, 201 )
    rvals  = np.sqrt( rspace[1:] * rspace[:-1] )
    
    # Load Densities on r grid
    MCum_DM  = enclosedMass_spherical( rho_DMO_func, rspace )      # the enclosed mass in DM
    MCum_bar = enclosedMass( rho_Baryon_func, rspace )             # the enclosed mass in baryons
    MCum_DM_contracted = contract_enclosed_mass( MCum_DM, MCum_bar, f_bar=f_bar )
        
    Contracted_rho_dm = density_from_enclosed_mass( rspace, MCum_DM_contracted, rvals )
    return Contracted_rho_dm, rvals, MCum_bar, MCum_DM, MCum_DM_contracted, rspace



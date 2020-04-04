import numpy as np
from scipy import integrate

def contract_enclosed_mass( mass_DM, mass_bar, f_bar=0.157 ):
    """ Returns the contracted DM enclosed mass given the 'uncontracted' profile and that of the baryonic distribution.
   
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
    increase_factor = 0.45 + 0.38 * (eta_bar + 1.16)**0.53
    return mass_DM * increase_factor

def enclosed_mass( r, density ):
    """ Calculates the enclosed mass.
   
   Args:
      r             : the radial distances for which we have density values
      density       : array of density values
   Returns:
      Array of enclosed masses.
   """
    # calculate the bin edges
    r_edges = np.zeros( len(density)+1 )
    r_edges[0]    = 0.
    r_edges[1:-1] = 0.5 * (r[:-1] + r[1:])
    r_edges[-1]   = r[-1]
    
    # calculate the volume of each bin
    delta_V = 4./3. * np.pi * (r_edges[1:]**3 - r_edges[:-1]**3)
    
    return np.interp( r, r_edges[1:], (delta_V * density).cumsum() )

def contract_density( density_DM, density_bar, mass_DM, mass_bar, f_bar=0.157 ):
    """ Returns the contracted DM density profile given the 'uncontracted' density and that of the baryonic distribution.
    It uses the differential (d/dr) form of Eq. (11) from Cautun et al (2020).
   
   Args:
      density_DM    : array of DM densities. 
                          It corresponds to '(1-baryon_fraction) * density in
                          DMO (dark matter only) simulations'.
      density_bar   : array of baryonic densities.
      mass_DM       : enclosed mass in the DM component in the absence of baryons. 
                          It corresponds to '(1-baryon_fraction) * enclosed mass in
                          DMO (dark matter only) simulations'.
      mass_bar      : enclosed baryonic mass for which to calculate the DM profile.
      f_bar         : optional cosmic baryonic fraction.
   Returns:
      Array of 'contracted' DM densities.
   """
        
    eta_bar = mass_bar / mass_DM * (1.-f_bar) / f_bar  # the last two terms account for transforming the DM mass into the corresponding baryonic mass in DMO simulations
    first_factor = 0.45 + 0.38 * (eta_bar + 1.16)**0.53
    temp         = density_bar - eta_bar * density_DM * f_bar / (1.-f_bar)
    const_term   = 0.38 * 0.53 * (eta_bar + 1.16)**(0.53-1.) * (1.-f_bar) / f_bar * temp
    
    return density_DM * first_factor + const_term


def NFW_enclosed_mass( r, M200, conc, hFactor=0.6777 ):
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
    H0 = 1.e-1 * hFactor       # Hubble constant in (km/s) kpc^-1
    rho_critical = 3.*H0*H0 / (8.*np.pi * G_newton)
    
    R200 = (M200 / (4.*np.pi/3. * 200. * rho_critical ))**(1./3.)
    xsi = r / R200 * conc
    return M200 * (np.log(1.+xsi) - xsi/(1.+xsi)) / (np.log(1.+conc) - conc/(1.+conc))

def NFW_density( r, M200, conc, hFactor=0.6777 ):
    """ Returns the density for an NFW profile.
   
   Args:
      r             : radial distance in 'kpc'.
      M200          : halo mass in 'M_solar'.
      conc          : halo concentration w.r.t. halo radius R_200 = average enclosed 
                          density is 200 times the critical one.
   Returns:
      Array of enclosed masses for each 'r' value in units of 'M_solar'.
   """
    G_newton = 4.302e-6        # Newton constant in kpc Msun^-1 (km/s)^2
    H0 = 1.e-1 * hFactor       # Hubble constant in (km/s) kpc^-1
    rho_critical = 3.*H0*H0 / (8.*np.pi * G_newton)
    
    R200 = (M200 / (4.*np.pi/3. * 200. * rho_critical ))**(1./3.)
    Rs  = R200 / conc
    xsi = r / Rs
    rho_0  = M200 / (4.*np.pi * Rs**3 * (np.log(1.+conc) - conc/(1.+conc)) )
    return rho_0 / xsi / (1+xsi)**2


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



def enclosedMass( Func, rbins, N=2000 ):  #Gives the enclosed mass in r for (R,z) func
    r1 = np.logspace( np.log10(rbins[0]*1.e-3), np.log10(rbins[-1]*1.1), N+1 )  # the bins used for the enclosed mass calculation
    r  = np.sqrt( r1[1:] * r1[:-1] )
    dr = r1[1:] - r1[:-1]
    
    I = lambda x, xr: Func( R=xr * np.cos(x), z=xr * np.sin(x) ) * 4*np.pi * xr**2 * np.cos(x) 
    
    shellMass = np.zeros( r.shape[0] )
    for i in range( r.shape[0] ):
        shellMass[i] = integrate.quad( I, 0., np.pi/2, args=( r[i], ) )[0] * dr[i]
    return np.interp( rbins, r1[1:], shellMass.cumsum() )

def enclosedMass_spherical( Func, rbins ):  #Gives the enclosed mass in r for spherically symmetric (R,z) func
    I = lambda x: Func(x,0) * 4.* np.pi * x*x
    r_range = np.column_stack( ( np.hstack( (1.e-3*rbins[0],rbins[:-1]) ), rbins ) )
    out = np.zeros( len(rbins) )
    for i in range( len(rbins) ):
        out[i] = integrate.quad( I, r_range[i,0], r_range[i,1] )[0]
    return out.cumsum()

def sphericalAverage( Func, rbins ):  #Gives the spherical average in r for (R,z) func
    I = lambda x, xr: Func( R=xr * np.cos(x), z=xr * np.sin(x) ) * np.cos(x) 
    meanValue = np.zeros( rbins.shape[0] )
    for i in range( rbins.shape[0] ):
        meanValue[i] = integrate.quad( I, 0., np.pi/2, args=( rbins[i], ) )[0]
    return meanValue



def potential_contract_DM_halo( rho_DMO_func, rho_Baryon_func, f_bar=0.157 ):
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
    
    # calculate masses and densities on r grid
    MCum_DM  = enclosedMass_spherical( rho_DMO_func, rspace )      # the enclosed mass in DM
    MCum_bar = enclosedMass( rho_Baryon_func, rspace )             # the enclosed mass in baryons
    rho_DM   = rho_DMO_func( rspace, 0 )                           # DM density at each bin position
    rho_bar  = sphericalAverage( rho_Baryon_func, rspace )         # baryonic density at each bin position
    
    # contract the DM density profile
    rho_dm_contracted = contract_density( rho_DM, rho_bar, MCum_DM, MCum_bar, f_bar=f_bar )
    MCum_DM_contracted = contract_enclosed_mass( MCum_DM, MCum_bar, f_bar=f_bar ) 
    
    return rho_dm_contracted, MCum_DM_contracted, MCum_bar, MCum_DM, rspace



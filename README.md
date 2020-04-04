# The Milky Way potential according to Cautun et al (2020)
**Last reviewed:** v1.0

Contains a table (as an ASCII file) that gives the Milky Way radial mass profile and a python code to: 

i) calculate the contraction of the dark matter (DM) halo induced by the condensation of baryons at the centre of Milky Way mass haloes. These results are based on studying state-of-the-art hydrodynamical simulations such as [EAGLE](https://ui.adsabs.harvard.edu/abs/2015MNRAS.446..521S/abstract), [Apostle](https://ui.adsabs.harvard.edu/abs/2016MNRAS.457.1931S/abstract) and [Auriga](https://ui.adsabs.harvard.edu/abs/2017MNRAS.467..179G/abstract).

ii) implement the Milky Way mass profile into [galpy](https://docs.galpy.org/en/v1.5.0/), where it can be used to perform various Galactic dynamic calculations. The mass profile corresponds to the best fitting model given the [Eilers et al (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...871..120E/abstract) Gaia DR2 rotation curve and the [Callingham et al (2019)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.5453C/abstract) total mass measurement of our galaxy.

The code included here is based on the results and the best fitting Milky Way model that have been presented in [Cautun et al (2020)](https://arxiv.org/abs/1911.04557).  



## Getting Started

This is a stand alone python code that only makes use of a few common python modules: numpy and scipy. For the potential, the code uses [galpy](https://docs.galpy.org/en/v1.5.0/), which can be downloaded [here](https://github.com/jobovy/galpy).
 
You can just download the two python codes and either place them in the working directory, or point python towards the directory where the files are located. The latter can be done by adding the directory path to the environment variable $PYTHONPATH or by adding the path directly within python, for example using:
```
import sys
dir_path = "absolute_path_to_file_directory"
if dir_path not in sys.path:
    sys.path.append( dir_path )
```
 
## The Milky Way radial mass profile

If you are interested only in the spherically averaged mass profile of our galaxy, just download the [MW_enclosed_mass_profile.txt](./MW_enclosed_mass_profile.txt) ASCII file. It gives the enclosed mass and its associated 68 percentile confidence interval for a set of radial distances for: stars, all baryons, dark matter, and total matter. 

Examples of how to read in the data files and run various tasks are provided in the jupyter notebook file [example_notebook.ipynb](./example_notebook.ipynb) (this is saved as a pdf file at found [here](./example_notebook.pdf)). Some of the examples from the notebook are discussed below. 


## Calculating the contracted DM profile for a distribution of baryons

Start by importing the corresponding function:

```
from Cautun20_contraction import contract_enclosed_mass
```
This takes two arguments: 1) the enclosed DM mass and 2) the enclosed baryonic mass. Optionally, you can also specify the cosmic baryonic fraction. A quick example on how to calculate the DM profile for a galactic-mass halo that does not contain any baryons:
```
from Cautun20_contraction import NFW_enclosed_mass
import numpy as np

r = np.logspace( -1, 2, 61 )
mass_total = NFW_enclosed_mass( r, M200=1.e12, conc=9 )      # total mass
f_bar      = 0.157   # cosmic baryon fraction
mass_DM    = mass_total * (1.-f_bar)  # DM mass
mass_bar   = 0.                       # no baryons
mass_DM_contracted = contract_enclosed_mass( mass_DM, mass_bar, f_bar=f_bar )
```
To see the difference, you can compare the original and contracted DM masses:
```
# plot the enclosed mass (actually plots mass / r^2 to decrease the dynamical range)

import matplotlib.pyplot as plt

plt.loglog( r, mass_DM / r**2, 'k-', label="original" )
plt.loglog( r, mass_DM_contracted / r**2, 'r--', label="contracted" )
plt.legend()
plt.show()

plt.semilogx( r, mass_DM_contracted / mass_DM, label="contracted / original" )
plt.legend()
```
This is a trivial example and the difference is a constant multiplication factor. When using realistic baryonic distributions, the difference between the original and the contracted DM profiles is more complex (see other examples in the [jupyter notebook](./example_notebook.ipynb)).

The code also comes with a simple function that calculates the contracted DM density. For best results, you need to give to the function the spherically averaged DM and baryon densities, as well as the enclosed DM and baryon masses at the same radial distances. The following code calculate the contracted DM density for the same case as above, that is for the case of zero baryonic mass.
```
from Cautun20_contraction import contract_density, NFW_density

# First calculate the NFW radial density profile for the DM component for the same halo used for the enclosed mass calculation.
density_NFW = NFW_density( r, M200=1.e12, conc=9 ) * (1.-f_bar)

density_DM_contracted = contract_density( density_NFW, 0., mass_DM=mass_DM, mass_bar=mass_bar, f_bar=f_bar )
```


## Calculating orbits with the Cautun et al (2020) Galactic mass profile

The package also comes with a set of python functions that implement the best fitting mass profile in the galpy package for galactic dynamics. The mass profile contains 7 components: a thin and a thick stellar disc, an HI and a molecular disc, a stellar bulge, a circumgalactic medium (CGM) component, and a contracted DM halo. To use this mass profile, just load the module included here (this will take about one minute since it performs some calculations when loading the module).

```
from Cautun20_galpy_potential import Cautun20
```
Now Cautun20 is a potential including all the components described above. It can be used to calculate various galactic quantities such as circular rotation curve, stellar and satellite orbits. For examples, see the [galpy](https://docs.galpy.org/en/v1.5.0/#tutorials) documentation, such as the one associated to the various [Milky Way potentials](https://docs.galpy.org/en/v1.5.0/reference/potential.html#new-in-v1-5-milky-way-like-potentials) implemented within galpy.


You can also access the various components of the potential, such as the halo or the bulge:
```
Cautun_halo, Cautun_Discs, Cautun_Bulge, Cautun_cgm = Cautun20
```
Alternatively, if you need the spherically averaged enclosed DM or baryonic mass profiles, these have been calculated when loading the module and can be accessed as:
```
from Cautun20_galpy_potential import rspace, rho_DM_contracted, MassCum_DM_contracted, MassCum_bar, MassCum_DM_uncontracted
```
where rspace gives the radial values from the Galactic Centre for which the enclosed masses and density were calculated. Note that these quantities are expressed in internal galpy units.


## Contributors
* **Marius Cautun (Leiden University)** 
* **Thomas Callingham (Durham University)** - *Implemented the Galactic mass profile in galpy.*


## Reference
This code and accompanying input data are freely available. If using this code,
a derivative work or results thereof, please cite:
[Cautun et al (2020)](https://arxiv.org/abs/1911.04557)

If you have any questions or would like help in using the code, please email:
> marius 'dot' cautun 'at' gmail 'dot' com


## License

This project is licensed under GNU GENERAL PUBLIC LICENSE Version 3 - see the [LICENSE.md](LICENSE.md) file for details.

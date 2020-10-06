# Perseus Cluster Gamma Ray Calibration

This repository contain the notebooks and scripts used for the Perseus cluster model calibration in the context of CTA observations.


## Relevant biliographie

-  Hanisch and Erickson (1980) [https://ui.adsabs.harvard.edu/abs/1980AJ.....85..183H/abstract] \
Non detection of the Perseus halo

-  Noordam and  de Bruyn (1982) [https://ui.adsabs.harvard.edu/abs/1982Natur.299..597N/abstract] \
WSRT processing applyied to Perseus

- Pedlar et al. (1990) [https://ui.adsabs.harvard.edu/abs/1990MNRAS.246..477P/abstract, https://www.researchgate.net/publication/234294211_The_Radio_Structure_of_NGC1275] \
Radio structure around NGC1275

- Burns et al. (1992) [https://ui.adsabs.harvard.edu/abs/1992ApJ...388L..49B/abstract] \
Confirmation of the Perseus mini-halo

- Sijbring (1993) [https://ui.adsabs.harvard.edu/abs/1993rchl.book.....S/abstract] \
WSRT observations of Perseus

- Ettori et al. (1999) [https://ui.adsabs.harvard.edu/abs/1998MNRAS.300..837E/abstract] \
ROSAT PSPC observations of Perseus \
H0 = 50 km/s/Mpc

- Jones & Forman (1999) [https://ui.adsabs.harvard.edu/abs/1999ApJ...511...65J/abstract] \
Einstein telescope obesrvation of nearby clusters and groups \
See [https://cdsarc.unistra.fr/viz-bin/cat/J/ApJ/511/65] for the full table (Table 4) \
H0 = 50 km/s/Mpc

- Churazov et al. (2003) [https://ui.adsabs.harvard.edu/abs/2003ApJ...590..225C/abstract] \
XMM observations of Perseus: Sx and T \
H0 = 50 km/s/Mpc ?

- Churazov et al. (2004) [https://ui.adsabs.harvard.edu/abs/2004MNRAS.347...29C/abstract] \
XMM observations of Perseus: gas motion \
H0 = 70 km/s/Mpc

- Sanders et al. (2004) [https://ui.adsabs.harvard.edu/abs/2004MNRAS.349..952S/abstract] \
Small scale thermodynamics and abundance in Perseus \
H0 = 70 km/s/Mpc

- Gitti et al. (2002) [https://ui.adsabs.harvard.edu/abs/2002A%26A...386..456G/abstract] \
Modeling of the radio mini-halo in a re0-acceleration framework. Profile, 3 point spectrum and spectral index profile are shown \
H0 = 50 km/s/Mpc

- Sanders et al. (2004) [https://ui.adsabs.harvard.edu/abs/2004MNRAS.349..952S/abstract] \
Chandra observations in the core, showing temperature and density profile

- Taylor et al. (2006) [https://ui.adsabs.harvard.edu/abs/2006MNRAS.368.1500T/abstract] \
Magnetic field in the center of Perseus

- Simionescu et al. (2011) [https://ui.adsabs.harvard.edu/abs/2011Sci...331.1576S/abstract] \
Baryon physics to the edge of Perseus in 2 specific directions \
H0 = 70 km/s/Mpc

- Aleksic et al. (2010) [https://ui.adsabs.harvard.edu/abs/2010ApJ...710..634A/abstract] \
MAGIC observations of Perseus

- Aleksic et al. (2012) [https://ui.adsabs.harvard.edu/abs/2012A%26A...541A..99A/abstract] \
MAGIC observations of Perseus

- Werner et al. (2013) [https://ui.adsabs.harvard.edu/abs/2013Natur.502..656W/abstract] \
Metalicity profile of Perseus

- Zandanel et al. (2014) [http://adsabs.harvard.edu/abs/2014MNRAS.438..124Z] \
Hadronic model for Perseus compared to Pedlar et al. (1990) data

- Zhuravleva et al. (2014) [https://ui.adsabs.harvard.edu/abs/2015MNRAS.450.4184Z/abstract] \
Density fluctuation in Perseus

- Urban et al. (2014) [https://ui.adsabs.harvard.edu/abs/2014MNRAS.437.3939U/abstract] \
Azimuthaly resolved X-ray spectroscopy of Perseus

- Ahnen et al. (2016) [http://adsabs.harvard.edu/abs/2016A%26A...589A..33A] \
MAGIC observations of Perseus

- Gendron-Marsolais et al. (2017) [http://adsabs.harvard.edu/abs/2017MNRAS.469.3872G]\
Deep VLA observations of Perseus

- Hitomi Collaboration (2018) [https://ui.adsabs.harvard.edu/abs/2018PASJ...70....9H/abstract] \
Gas motion in Perseus

- Donnert et al. (2018) [https://ui.adsabs.harvard.edu/abs/2018SSRv..214..122D/abstract] \
Magnetic field review

- Van Weeren et al. (2018) [https://ui.adsabs.harvard.edu/abs/2019SSRv..215...16V/abstract] \
Diffuse cluster scale radio emission review

- Bykov et al. (2018) [https://ui.adsabs.harvard.edu/abs/2019SSRv..215...14B/abstract] \
Shocks and nonthermal particles in clusters, a review

- Johnson et al. (2020) [https://ui.adsabs.harvard.edu/abs/2020ApJ...888..101J/abstract] \
Determination of magnetic field from RM

## Gas thermal density

- X-ray properties (Jones & Forman 1999): \
L_x = 14.948 ergs/s (0.5-4.5 keV) \
R_c = 0.28 Mpc \
β = 0.58 \
n_0 = 4.050 cm^-3 \
M_gas = 1.10 10^14 Msun (1 Mpc) \
M_gas = 1.93 10^14 Msun (5 Rc) \
M_tot = 3.59 10^14 Msun (1 Mpc) \
M_tot = 5.31 10^14 Msun (5 Rc)\
H0 = 50 km/s/ Mpc

- double beta-model (Churazov 2004, rescaled from Churazov 2003, from Jones & Forman 1999): \
n_e(r) = 4.6e-2 / (1+(r/57)**2)**1.8 + 4.8e-3 / (1+(r/200)**2)**0.87 cm^-3 \
3 β1 / 2 = 1.8 => β1 = 1.2, 3 β2 / 2 = 0.87 => β2 = 0.58 \
H0 = 70 km/s/ Mpc

- Ettori et al. (1999): \
R_200 = 2.7 Mpc, M_tot (2.3 Mpc) = 1.19 ± 0.12 10^15 M_sun (80% C.L.) \
Best-fit β-model β = 0.63 ± 0.01, r_c = 9.3 ± 0.2 arcmin  \
H0 = 50 km/s/Mpc

- Simionescu et al. (2011):\
Slope of the density profile α = 1.68 ± 0.04 above 0.7 Mpc, for the East and Northwest arms. This is consistent with previous results from ROSAT (Ettori et al. 1999) data extending out to ∼ 1.4 Mpc.\
This implies β2 = 1.68/3 = 0.56 ± 0.02 \
However, looking at Ettori et al. (1999), Table 5, it is not consistent with Simionescu et al.: E: β = 0.88 or powerlaw index = 2.3 ; N: beta = 0.67 or powerlaw index = 1.87 ; W: beta = 0.61 or powerlaw index = 1.87 \
Correcting for clumping leads to a slope of 2.5, i.e. beta = 0.83, in good agreement, this time, with Ettori et al. (1999), but who did not correct for clumping. \
R200 = 1.79  ± 0.04 Mpc, M = 6.65+0.43-0.46 × 10^14 Msun \
The best-fitting NFW function has a concentration parameter c = 5.0 ± 0.5
H0 = 70 km/s/Mpc

- Zhuravleva et al. (2014) : \
Single beta-model fit with r_c = 1.26 arcmin and β=0.53, ok to 300 kpc \
H0 = 72 km/s/Mpc

- Urban et al. (2014): \
Density beta model at r > 10 arcmin: \
r_c = 13.18±2.31 and beta = 0.71±0.05, n0 is not provided
Pressure profile: \
Planck V (2013) with r500=59.7arcmin±0.4 arcmin gives a good fit
H0 = 70 km/s/Mpc


## Gas temperature

- From Ettori et al. (1999):\
The temperature profile is flat in the outskirt (0.3-2 Mpc), reaching 6.3 keV \
! Pressure in the outermost bin fixed so that the deprojected temperature match 6.3 keV ! \
H0 = 50 km/s/Mpc

- From Churazov et al. (2003), valid in the core: \
T_e (r) = 7 * (1+(r/71)**3) / (2.3+(r/71)**3) keV
H0 = 50 km/s/Mpc

- The projected temperature profile is provided in Simionescu et al. (2011)
H0 = 70 km/s/Mpc

- Urban et al. (2014): \
Deprojected temperature profile beyond 100 kpc: T(r) = T0 [ (r / r_cool)^ac +T_min/T0] / [1+(r / r_cool)^ac] * (r / r_t)^(−a) / [1 + (r / r_t) ^ b]^(c/b) \
[T0, T_min/T0, r_cool (kpc), ac, r_t (Mpc), a, b, c]=[4.06, 0.72, 294, 6.72, 1.6, 0.33, 16.24, 2.36] \


## Metalicity

- Constant metalicity for Perseus: Z  = 0.306 ± 0.012 (Werner et al. 2013) \
This is true down to 400 kpc

- Core metalicity gets to Z = 0.5 at r < 100 kpc (Sanders et al. 2004)


## Magnetic field

- B = 25 uG at 10 kpc, very crude approximation (Taylor et al. 2006)

- B(r) \propto n_e^[0.4-0.7] with central magnetic field in [3.9-5.4] uG

- B_0 cannot be determined better than in a range of 3 even under ideal conditions (Johnson et al. 2020)

- If the magnetic flux is conserved, B = B_star (rho/<rho>)^2/3 (e.g. Donnert et al. 2018 for a review)

## Cosmic rays

In MAGIC paper: "The reported radial spectral steepening of the radio mini-halo emission (Gitti et al. 2002) could be an observational artifact owing to a poor signal-to-noise ratio in the outer core of the cluster and the ambiguity in determining the large scale Fourier components owing to nonuniform coverage of the Fourier plane and missing short-baseline information: the so-called “missing zero spacing”-problem of interferometric radio observations. By comparing the spectral index distribution of the three radio maps (at 92 cm, 49 cm, and 21 cm, Sijbring 1993), radial spectral flattening depending on the chosen radial direction seems possible."


## Redshift

z = 0.017284 (Hitomi colaboration 2018)

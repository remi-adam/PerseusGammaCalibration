# PerseusGammaCalibration
This repository contain the notebooks and scripts used for the Perseus cluster model calibration in the context of CTA observations.


## Relevant biliographie

- Ettori et al. (1999) [https://ui.adsabs.harvard.edu/abs/1998MNRAS.300..837E/abstract]
ROSAT PSPC observations of Perseus
H0 = 50 km/s/Mpc

- Jones & Forman (1999) []
Einstein telescope obesrvation of nearby clusters and groups
See [https://cdsarc.unistra.fr/viz-bin/cat/J/ApJ/511/65] for the full table (Table 4)
H0 = 50 km/s/Mpc

- Churazov et al. (2003) [https://ui.adsabs.harvard.edu/abs/2003ApJ...590..225C/abstract]
XMM observations of Perseus: Sx and T
H0 = 50 km/s/Mpc ?

- Churazov et al. (2004) [https://ui.adsabs.harvard.edu/abs/2004MNRAS.347...29C/abstract]
XMM observations of Perseus: gas motion
H0 = 70 km/s/Mpc

- Gitti et al. (2002) [https://ui.adsabs.harvard.edu/abs/2002A%26A...386..456G/abstract]
Modeling of the radio mini-halo in a re0-acceleration framework. Profile, 3 point spectrum and spectral index profile are shown
H0 = 50 km/s/Mpc

- Simionescu et al. (2011) [https://ui.adsabs.harvard.edu/abs/2011Sci...331.1576S/abstract]
Baryon physics to the edge of Perseus in 2 specific directions
H0 = 70 km/s/Mpc

- Aleksic et al. (2012) []
MAGIC observations of Perseus.


## Gas thermal density

- beta model



- double beta-model (Churazov 2004, rescaled from Churazov 2003, from Jones & Forman 1999):
n_e(r) = 4.6e-2 / (1+(r/57)**2)**1.8 + 4.8e-3 / (1+(r/200)**2)**0.87 cm**-3
Probably disputable for the large scales.
H0 = 70 km/s/ Mpc


- Ettori et al. (1999)




## Gas temperature

- 

## Magnetic field

- B = 25 uG at 10 kpc

- B(r) \propto n_e**0.5 ?
for Coma

## Cosmic rays

In MAGIC paper: "The reported radial spectral steepening of the radio mini-halo emission (Gitti et al. 2002) could be an observational artifact owing to a poor signal-to-noise ratio in the outer core of the cluster and the ambiguity in determining the large scale Fourier components owing to nonuniform coverage of the Fourier plane and missing short-baseline information: the so-called “missing zero spacing”-problem of interferometric radio observations. By comparing the spectral index distribution of the three radio maps (at 92 cm, 49 cm, and 21 cm, Sijbring 1993), radial spectral flattening depending on the chosen radial direction seems possible."

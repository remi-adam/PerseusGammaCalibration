#!/bin/bash

source activate ctaphys
source '/pbs/home/r/radam/Project/CTA/Phys/Software/Config/ctaphys_current_bash'

python /pbs/home/r/radam/Project/CTA/Phys/Software/PerseusGammaCalibration/cr_constrain_distribution_from_radio.py

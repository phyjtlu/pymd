# pymd
This is a set of scripts to run semi-classical Langevin molecular dynamics on junctions connecting to several electron or phonon baths. The details of the method are presented in Prog. Surf. Sci. [https://doi.org/10.1016/j.progsurf.2018.07.002]. 

To do the molecular dynamics, we need a force driver. Currently, we implemented Siesta and Brenner potential. 

The important feature of this script is that, the statistics of the bath degrees of freeom is quantum. For example, the zero point fluctuations are included; they fulfill the quantum-mechanical Bose-Einstein distribution. Moreover, the electron bath is allowed to be in a nonequilibrium steady state (non-thermal).

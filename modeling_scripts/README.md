# Modeling scripts for integrative structure determination of Aβ40-sTREM2 complex

Directory for the production scripts for the structure of the Aβ40fibril-sTREM2 complex.

## List of files and directories:

`modeling_scripts`   	 contains all scripts used for modeling
- run_modeling.sh is the relevant submission scripts for a SGE cluster (Wynton at UCSF / QB3 Cluster at UCSF).
- modeling.py is the modeling script used by PMI to generate the samples

## Input for all the scripts

All inputs are located in the [data](../data) directory.

## Scripts used to generate the solution structures of the Abet40-sTREM2 complex: 

- `run_modeling.sh` and `modeling.py` are used in combination to generate the solution structure of Aβ40-sTREM2 complex using a dataset of 15 crosslinks, excluded volume, connectivity restraints, and roto-translational symmetry contraints.

## Information

_Author(s)_: Dibyendu Mondal

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
Soluble TREM2 inhibits secondary nucleation of Aβ fibrillization and enhances cellular uptake of fibrillar Aβ.
Ketaki D. Belsare#, Haifan Wu#, Dibyendu Mondal#, Annalise Bond, Jia Jin, Hyunil Jo, Addison E. Roush, Kala Bharath Pilla, Andrej Sali*, Carlo Condello*, and William DeGrado*
Proceedings of the National Academy of Sciences DOI:10.1073/pnas.2114486119.  

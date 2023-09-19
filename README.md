# v2f
An OpenFOAM package for an anisotropic k-epsilon-v2-f turbulent FENE-P viscoelastic model.

Full publication found in Physics of Fluids. https://doi.org/10.1063/5.0159668

The code was developed within foam-extend 4.0 and requires an update to OpenFOAM v2212+.

turbulenceModels/...

-- Contains the source code for the Viscoelastic model embedded within the RAS class structure, with additional polymeric fields solved.

Channel/Viscoelastic/...

-- Contains the case file for quick simulations. 

-- RUNCASE-v2.sh is a bash script which can choose the case with specific rheological parameters and run, plotting the data afterwards.

-- The coefficients can also be varied within RASProperties 

For more details or queries, please email me.

Thanks :)

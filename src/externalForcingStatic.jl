function [Fexternal, Fdof] = externalForcingStatic()

#owens externalForcingStatic function for the OWENS toolkit
# **********************************************************************
# *                   Part of the SNL OWENS toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [Fexternal, Fdof] = externalForcingStatic()
#
#   This function specifies external forcing for a static analysis.
#   Fexternal is a vector of loads and Fdof is a corresponding vector of
#   degrees of freedom the concentrated loads in Fexternal correspond to.
#   The global degree of freedom number corresponding with the local degree
#   of freedom of a node may be calculated by:
#   globalDOFNumber = (nodeNumber-1)*6 + localDOFnumber
#   The localDOFnumber may range from 1 to 6 such that 1 corresponds to a
#   force in "x direction" of the co-rotating hub frame. 2 and 3
#   corresponds to a force in the "y" and "z directions" respectively. 4,
#   5, and 6 correspond to a moment about the "x", "y", and "z" directions
#   respectively.

#
#      input: NONE

#      output:
#      Fexternal     = vector of external loads (forces/moments)
#      Fdof          = vector of corresponding DOF numbers to apply loads to

Fexternal = [];
Fdof = [];

end
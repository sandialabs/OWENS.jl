function setInitialConditions(initCond,u,numDOFPerNode)
    #setInitialConditions sets initial conditions
    # **********************************************************************
    # *                   Part of the SNL OWENS Toolkit                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
    #   [u] =  setInitialConditions(initCond,u,numDOFPerNode)
    #
    #   This function reads initial conditions from file
    #
    #   input:
    #   initCond      = array containing initial conditions
    #                     initCond(i,1) = node number for init cond i
    #                     initCond(i,2) = local DOF number for init cond i
    #                     initCond(i,3) = value for init cond i
    #   u             = displacement vector
    #   numDOFPerNode = number of degrees of freedom per node
    #
    #   output:
    #    u             = displacement vector modified for initial conditions

    len=size(initCond) #get number of specified initial conditions
    #unspecified initial conditions are assumed to
    #be zero

    for i=1:len[1] #loop over initial conditions
        if (initCond[i,2]>numDOFPerNode) #error check
            error("setInitalConditios:: DOF greater than numDOFPerNode")
        end
        index = (initCond[i,1]-1)*numDOFPerNode + initCond[i,2] #calculate global DOF number for initial condition
        u[index] = initCond[i,3] #specify initial condition using global DOF number
    end

    return u

end

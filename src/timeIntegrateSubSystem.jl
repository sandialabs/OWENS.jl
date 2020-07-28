function timeIntegrateSubSystem(M,K,C,F,delta_t,u,udot,uddot)
    #timeIntegrateSubSystem integrates a system using Newmark-Beta method
    # **********************************************************************
    # *                   Part of the SNL OWENS Toolkit                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
    #
    #   [unp1,udotnp1,uddotnp1] = timeIntegrateSubSystem(M,K,C,F,delta_t,u,udot,uddot)
    #
    #   #This function perform integration of a system using the Newmark-Beta
    #   method(constant-average acceleration sceheme).
    #
    #   input:
    #   M        = system mass matrix
    #   K        = system sttiffness matrix
    #   C        = system damping matrix
    #   F        = system force vector
    #   delta_t  = time step
    #   u        = displacement at beginning of time step
    #   udot     = velocity at beginning of time step
    #   uddot    = acceleration at beginning of time step
    #
    #
    #   output:
    #   unp1        = displacement at end of time step
    #   udotnp1     = velocity at end of time step
    #   uddotnp1    = acceleration at end of time step
    #
    alpha = 0.5 #constant avg accel scheme
    gamma = 0.5
    beta = 0.5*gamma

    a1 = alpha*delta_t
    a2 = (1.0-alpha)*delta_t
    a3 = 1.0/(beta*delta_t*delta_t)
    a4 = a3*delta_t
    a5 = 1.0/gamma-1.0
    a6 = alpha/(beta*delta_t)
    a7 = alpha/beta - 1.0
    a8 = delta_t*(alpha/gamma-1.0)

    A = a3*u + a4*udot + a5*uddot
    B = a6*u + a7*udot + a8*uddot

    Khat = K + a3.*M + a6.*C
    Fhat = F + M*(A') + C*(B')

    unp1 = Khat\Fhat

    uddotnp1 = a3*(unp1-u) - a4*udot - a5*uddot
    udotnp1 =  udot + a2*uddot + a1*uddotnp1

    return unp1,udotnp1,uddotnp1

end


struct Mesh
    nodeNum
    numEl
    numNodes
    x
    y
    z
    elNum
    conn
    type
end

mutable struct Ort
    Psi_d
    Theta_d
    Twist_d
    Length
    elNum
    Offset
end

mutable struct BC_struct
    numpBC
    pBC
    numsBC
    nummBC
    isConstrained
    map
    redVectorMap
end

mutable struct BladeData
    numBlades
    bladeNum
    h
    nodeNum
    elementNum
    remaining
end

mutable struct SectionPropsArray
    ac
    twist
    rhoA
    EIyy
    EIzz
    GJ
    EA
    rhoIyy
    rhoIzz
    rhoJ
    zcm
    ycm
    a
    EIyz
    alpha1
    alpha2
    alpha3
    alpha4
    alpha5
    alpha6
    rhoIyz
    b
    a0
    aeroCenterOffset
end

mutable struct El
    props
    elLen
    psi
    theta
    roll
    rotationalEffects
end

mutable struct NodalTerms
    concLoad
    concStiff
    concMass
    concStiffGen
    concMassGen
    concDampGen
end

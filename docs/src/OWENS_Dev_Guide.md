---
author:
- Kevin R. Moore
bibliography:
- ac_sources.bib
title: Offshore Wind Energy Simulator (OWENS) Glue Code Basic Theory,
  Frames of Reference, and Inter-Module Coupling Methods
---

# Introduction

The Offshore Wind Energy Simulator (OWENS) Code is a floating
turbine-platform code capable of steady, unsteady, and modal analysis.
It is comprised of, or coupled to many modules for aerodynamics,
hydrodynamics, mooring dynamics, and finite element formulations.

*How* many of these items are used is not straightforward or a-priori
and can easily lead to confusion or mistake (like how the aerodynamics
uses inflow data; is it just simple based on the mean hub velocity, or
in the rotating blade frame of reference? If the latter, what is that
frame of reference so that we make sure we don't incorrectly map the
vectors?) This document is intended as a developers guide for the intent
of clearing up the relative ambiguity between the many moving and
modular pieces. It is not a full theory document or users guide. For
theory and use, please see the theory and user guides respectively.
Within this document, we give enough detail about the operations of the
\"glue\" code, the frames of reference, and coupling methods to aid
further development.

# Contributing

Please make all feature changes and bug fixes as branches and then
create pull requests against the dev branch. The dev branch will be
periodically pulled into master.

# Glue Code Basic Theory

Due to the nonlinearity and non-closed-form nature of the floating
turbine problem in both the design and temporal space, OWENS uses many
types of iterative solution methods, specific to the physics involved.
Loose or Gauss-Seidel iteration is used for the coupling between the
modules. gives an overview of the general code structure for transient
analysis, with the main equations being solved and the relative location
of the iteration loops.

::: algorithm
::: spacing
1.18

::: algorithmic
Initialize simulation (all module parameters, meshing, composite
properties, etc) Initialize displacements and forces as 0

Get driveshaft displacement from generator torque Get rotor speed,
rotational displacement from generator torque, aero torque, and turbine
inertia Get aero forces from rotor speed, position Get structural
displacements and reaction forces from aero forces Update
displacements$_\text{guess}$ with the calculated displacements
:::
:::

[]{#alg:Transient label="alg:Transient"}
:::

# WindIO Usage
The majority of the WindIO standard can be used for the structural and aerodynamic design.  Airfoil polars are forthcoming.  However, there are the following additional and limitations:
  - Gravity single value or array of 3, is in the global frame of reference, so positive gravity would be pulling the structure up.
  - Start and end nd_arc values are the only input supported.  Otherwise, there are as many as 12 different ways to specify this input.
  - Grids must be defined from end to end, set number of plies to 1e-6 at unused spanwise sections.
  - The blade object is reused for the tower, struts, cables, etc. since these elements can be modeled at the same fidelity as the blade and can be simplified via simpler inputs, but the other WindIO definitions don't allow for the more complex inputs.

# Global Frame of Reference

Thanks to some legacy VAWT frames of reference where wind was originally
coming from the top down, and then from the left right, the convention
was to simply rotate the axes instead of properly redefining the global
frame of reference to match standard math convention again. Therefore,
an azimuth angle of 0, or yaw angle for standard VAWT frame of reference
is top dead center. Therefore, blade 1 of a vawt should also be top dead
center, i.e. aligned with the y+ axis as shown below.

![Global Frame of Reference. Wind comes in from the left in the
direction of the positive X-axis, the positive Y-axis is 90 degrees
counter-clockwise to the X-axis. Z-axis is vertical. Positive rotations
follow right hand rule. Note that the turbine is modeled in the rotating
frame of reference. Therefore, the turbine stays stationary in its
initial meshed position, excepting deformations, and inertial effects
are calculated. The developer must use the calculated rotational
position to map to the other models as
necessary.](./figs/global_FOR3.png){#fig:ac_velocities width="60%"}

# InflowWind and TurbSim Velocity Inflow

This is how it is defined in the manual, despite inverting the positive
direction of rotation compared to standard notation.

![Wind frame of reference is the same as global. However, the angle
follows TurbSim and [InflowWind](https://github.com/old-NWTC/InflowWind)
notation Wind propagation angle is zero when aligned with the positive
X-axis and clockwise positive, in the direction of the negative negative
Y-axis.](./figs/inflow_wind.pdf){#fig:ac_velocities width="50%"}

# Airfoils

Airfoil geometry is typically defined starting at the trailing edge (1,0), looping around the bottom/pressure side to the leading edge (0,0), and then back along the top/suction side to the trailing edge (1,0) again.  There are several automatic corrections to handle files that slightly ill formatted, with output warnings, but this is limited and should not be relied on.

Airfoil polars are currently in different formats for OWENSAero and AeroDyn, see examples.  The airfoil name in the structural definition is assumed to be the same airfoil name as the geometry and the polar data.  The path to these files is automatically resolved as relative to the run location, so it expects an airfoils folder there.

Future work is planned to include enabling OWENSAero to read in the same airfoil files as AeroDyn, and to automatically generate polars from xfoil.jl.

# OWENSAero Aerodynamics (AC and DMS)

![VAWT 2D section looking downwards with induced velocity $w$ vector
broken into $u$ and $v$ components depicted by arrows. Airfoils show
example blade locations with dots aligning to the circumferential
discretization. Aero frame of reference is the same as global, however a
blade is at 0 degrees azimuth when it is aligned with the Y-axis. If an
aero module is 2D, it is made quasi-3D by stacking slices from lower to
higher in the Z-axis.](./figs/vawt_slice.pdf){#fig:ac_velocities
width="50%"}

Note: CACTUS, does not fully follow this scheme and differs from the
global frame of reference by switching the Y and Z axes for a VAWT.
Also, be careful with the geometry inputs; if the blade 1 starts out at
the \"south\" position, as opposed to the north, then it will behave as
though it were a north starting blade rotating clockwise, and a
symmetric gust velocity will match (if the simple iec uniform gust is
used). All else for CACTUS follows the description above.

Aerodynamics forces are radial inward positive and tangential in direction of rotation positive.
So, the blade local accelerations in the radial will be outward positive and in the tangential 
will be in the opposite direction of rotation for a CCW rotor.  CW rotors are still in work and
would need to have the local frame of reference rotation verified.  Therefor, for the added mass
 equations, the forces should be negated. 

## OWENSFEA Structures

Yaw is equivalent to the global azimuth angle as above, and is denoted
in the code as Psi. Zero is top dead center, right hand rule positive
rotation. The subscript D indicates degrees. Theta is the element slope,
or delta, and

The rotation sequence is a Roll (Rx) - Pitch (Ry) - Yaw (Rz) sequence in
OWENSFEA (R = Rz \* Ry \*Rx in multiplication order). This corresponds
to Twist - Theta - Psi angle ordering. In the matrix multiplication this
is R = $R_{psi}$ \* $R_{theta}$ \* $R_{twist}$, which is a \[3,2,1\]
matrix ordering.

::: center
   Name    Rotation Axis   OWENS Angle Names
  ------- --------------- -------------------
   Roll       about X            Twist
   Pitch      about Y            Theta
    Yaw       about Z             Psi

  : Module Level of Detail in This Document
:::

We apply these in order of Roll, then Pitch, then Yaw, which is R = Rz
\* Ry \* Rx (application is from right to left in rotations) for a
rotation of vector x by R. i.e R(x) = Rz \* Ry \* Rx \* x.

In the code I now have ang = \[ort.Psid\[idx\],
ort.Thetad\[idx\],ort.Twistd\[idx\]+180,-90\]; DCM = CH2G \*
createGeneralTransformationMatrix(ang,\[3,2,1,2\]); Where I have an
additional rotation by -90 about the Y axis to align Z with X before
applying the other rotation sequences. This is due to AD15 using the IEC
conventions with +Z along the blade, but the beam solver using the more
typical beam formulation of +X along the beam.

This is almost consistent between blades and struts (have a +90 rotation
on struts for some reason that I haven't sorted out yet).

Rotating Frame of Reference, 6 DOF where 1 = turbine vertical force, 2 =
turbine 2D slice tangential force, 3 = turbine 2D slice normal force, 4
= blade M25 twisting moment, 5 = blade curvature twisting moment, 6 =
blade sweep moment.

The Mesh matches the global frame of reference of x, y, and z. Element
length is the length along the element.

The mesh itself is comprised only of components. For example, a tower,
two blades, and four struts. The element number is sequential. There are
overlapping points where each component connects. The mesh has a
connectivity vector, which has rows corresponding to each element,
column 1 corresponds to the \"master\" node, and column 2 corresponds to
the \"slave\" node. The element connection in the mesh is only
intra-component. I.e. there is no connectivity between components - that
is defined in the joint matrix, which has columns for: Joint Number,
Joint Master Node, Joint Slave Node, Joint Type, Joint Mass, Not Used,
PsiD, ThetaD. The D indicates angle in degrees. Joint types are: (0 =
weld(fixed), 1=pinned, 2 = hinge joint with axis about slave node
element's e2 axis, 3 = hinge joint axis about slave node element's e1
axis, 4 = hinge joint axis about slave node element's e3 axis). The Psi
and Theta are of the slave node (or its closest neighbor of the same
component due to the gaps in element mesh connectivity). The not used
column is just filled with zeros. The \"flapwise\" normal vector of an
element is forced to be away from the machine for consistency. During
the meshing process, the component type need to be known in order to get
this right: Mesh Type: 0-blade 1-tower 2-strut.

## Marine Considerations
When we turn on AddedMass_Coeff_Ca>0.0, this causes a few changes.  The added mass is included in the GXBeam 6x6 inertia terms (2,2 and 3,3) since a force coupling is unstable.  This changes the more smooth gravitational and centrifugal loads, so they are offloaded to the aero model, and are handled in the buoyancy and centrifugal force respectively.  Note, that for the GX implementation, this is required, but
for the Timoshenko implementation, added mass is not coupled to the internal centrifugal force, so one
could use TNB or ROM with added mass, turn gravity on, turn it off for the aero, and turn the centrifugal
forcing off for aero as well.  That's a lot to change, so it is recommended to keep with the current
implementation until development occurs otherwise. The example cases and/or tests contain verification cases of the added mass and buoyancy implementations.

## Composites

## Hydrodynamics
Hydrodynamics are applied into OWENS using the HydroDyn module native to [OpenFAST](https://github.com/OpenFAST/openfast), but applied here through a standalone library linked to OWENS. The physics contained in the HydroDyn standalone library is virtually unchanged from its native application in OpenFAST. HydroDyn user documentation can be found [here](https://openfast.readthedocs.io/en/main/source/user/hydrodyn/index.html).

### Platform Structure
To enable capturing the hydrodynamics in OWENS for floating turbines, a second mesh is created to represent the platform.
OWENS assumes the platform is a rigid body, so the mesh is composed of a single element with effectively infinite stiffness.
The two nodes forming this element are:

* The platform reference point (i.e., the platform at the still water level along the vertical centerline of its center of gravity) is the bottom node
* The tower base point is the top node

The structures composing the platform mesh itself are the same as the OWENSFEA structures documented above, though most of the definitions to define the mesh can be set to zero or are very simple.
There are two main exceptions:

* The platform center of gravity must be defined, since neither node is at this location.
  This is defined in the `zcm` parameter in the `SectionPropsArray` for the platform mesh.
* The 6x6 rigid body mass matrix is defined as a concentrated input at the platform reference point node

A full example of these definitions can be found in the floating platform example notebook (TODO).

### Solve Procedure
When applying hydrodynamics, OWENS adds additional Gauss-Seidel loops to solve for the structural motions of both the platform mesh and the "topside" mesh (i.e. everything above the platform) in a way that generates global convergence. This process is:

1. Solve via a Gauss-Seidel loop for the topside mesh using platform motions from the *previous* time step
  * Setting platform motions to zero is how OWENS runs for land-based turbines
2. Solve via an outer Gauss-Seidel loop and an inner Newton-Raphson loop for the platform mesh, using updated external hydrodynamic loads from HydroDyn and tower base loads from step 1
  * The Newton-Raphson loop is contained in `Unsteady.OWENS_HD_Coupled_Solve`, and is necessary to properly account for added mass impacting both HydroDyn and OWENS in a tightly coupled manner.
  This solve process is complex, but is documented in detail in Section 2.4.1 of Devin et al. (2022) [here](https://www.mdpi.com/1996-1073/16/5/2462).
  This is a very similar procedure used to coupled HydroDyn to the ElastoDyn module in OpenFAST.
3. Solve via a Gauss-Seidel loop for the topside mesh using platform motions from the *current* time step, i.e. the motions calculated in step 2.

This flow chart visualizes this solve procedure, and which step corresponds to which modules connecting to the core OWENS glue code:
![Hydrodynamics Solve Procedure Flowchart](figs/hydro_solve_flowchart.pdf)

### Reference Frame Conversions
Since HydroDyn operates in the global reference frame, the hydrodynamics capabilities in OWENS deviate from the general convention in OWENS of operating in the hub reference frame.
As a result, operations manipulating the platform mesh (e.g. running `OWENSFEA.structuralDynamicsTransient`) are also applied in the global reference frame.
Functionally, this means that the `rbData` vector input is set to zero when using the platform mesh in `OWENSFEA.structuralDynamicsTransient`.
**Care must be taken for I/O procedures between the platform and topside meshes!**
The `Unsteady_utilities.frame_convert` function is intended to make this process simpler.

## Mooring
Mooring dynamics for floating platforms in OWENS are calculated using a standalone library linking to [MoorDyn](https://moordyn.readthedocs.io/en/latest/), another module native to OpenFAST.
The mooring loads calculated in MoorDyn are applied as an external load to the same platform mesh as the hydrodynamics, as described above.
Since the mooring dynamics do not need to account for added mass, the coupling between MoorDyn and OWENS is much simpler, and the mooring dynamics are applied in the Gauss-Seidel solve procedure without being modified by the additional Newton-Raphson solve in `Unsteady.OWENS_HD_Coupled_Solve`.
Note that MoorDyn operates in the global reference frame, so be careful with reference frame conversions for MoorDyn I/O.

# Coupling Methods

## Inflow - Aero

Direct: Aero module supplies an x-y-z and time coordinate, inflow
returns x-y-z velocity. This is repeated for all blade discrete points
as per the aero formulation.

## Aero - Turbine Structure

Loose Iteration: Structure provides blade local radius, twist, sweep,
and 6 DOF velocities, aero returns forces, moments. This is iterated on
until convergence. It is preferred to change this to a N-dimensional
root solver and pass gradients to the root solver to increase
performance.

Specifically, from the meshing process, the starting and ending node
numbers for the blades are known and the aerodynamic loads mapped to the
elements between those nodes.

## Turbine Structure - Platform Structure

Same as Aero-Structure

## Hydro - Platform Structure - Mooring

Same as Aero-Structure

## Structures - Composites

Initialization is Direct: Structures provides macro geometry, Composites
provide sectional properties. Composite Failure is Direct: Structure
provides strains, composites provides failure. Buckling is also
calculated.

## Controllers - Control Elements

Direct: Controllers provide reactionary inputs to control inputs in real
time based on dynamics.

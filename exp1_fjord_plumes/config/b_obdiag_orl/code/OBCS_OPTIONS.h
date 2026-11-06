C CPP options file for OBCS package
C Use this file for selecting options within the OBCS package
C Exp1 b_twoplume3d: only the eastern OB, Stevens (1990) (pkg/obcs/obcs_calc_stevens.F) with
C prescribed depth-mean normal velocity OBEu (needs ALLOW_OBCS_PRESCRIBE) + U/V sponge.
C Orlanski stays compiled (as in the thesis config verification/BayIW_Rectan) but is NOT used
C at run time (useOrlanskiEast=.FALSE.): in 3-D it left the depth-mean u free and drew inflow.
C No OBCS balance (it ignores addMass, pkg/obcs/obcs_balance_flow.F:300-315).

#ifndef OBCS_OPTIONS_H
#define OBCS_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_OBCS
C Package-specific Options & Macros go here

C Enable individual open boundaries
#undef ALLOW_OBCS_NORTH
#undef ALLOW_OBCS_SOUTH
#define ALLOW_OBCS_EAST
#undef ALLOW_OBCS_WEST

C Orlanski radiation code: compiled but not selected at run time in (b)
#define ALLOW_ORLANSKI

C Enable OB values to be prescribed via external fields that are read
C from a file
#define ALLOW_OBCS_PRESCRIBE

C Enable OB conditions following Stevens (1990)
#define ALLOW_OBCS_STEVENS

C Allow sponge layer treatment of open boundary conditions
#define ALLOW_OBCS_SPONGE

C Include hooks to sponge layer treatment of pkg/seaice variables
#undef ALLOW_OBCS_SEAICE_SPONGE

C balance barotropic velocity
C OB diagnosis (2026-11-06, notes/pd1_b_twoplume3d.md): balance for V6
C (pkg/obcs/obcs_balance_flow.F) and Stevens baroclinic u copied from the
C interior for V3 (pkg/obcs/obcs_calc_stevens.F:39-42, 193-207).
#define ALLOW_OBCS_BALANCE
#define OBCS_STEVENS_USE_INTERIOR_VELOCITY

C Allow to add barotropic tidal contributions to OB velocity
#undef ALLOW_OBCS_TIDES

C Use older implementation of obcs in seaice-dynamics
C note: most of the "experimental" options listed below have not yet
C       been implementated in new version.
#undef OBCS_UVICE_OLD

#ifdef OBCS_UVICE_OLD
C     The following five CPP options are experimental and aim to deal
C     with artifacts due to the low-frequency specification of sea-ice
C     boundary conditions compared to the model forcing frequency.
C     Ice convergence at edges can cause model to blow up.  The
C     following CPP option fixes this problem at the expense of less
C     accurate boundary conditions.
#undef OBCS_SEAICE_AVOID_CONVERGENCE

C     Smooth the component of sea-ice velocity perpendicular to the edge.
#undef OBCS_SEAICE_SMOOTH_UVICE_PERP

C     Smooth the component of sea ice velocity parallel to the edge.
#undef OBCS_SEAICE_SMOOTH_UVICE_PAR

C     Compute rather than specify seaice velocities at the edges.
#undef OBCS_SEAICE_COMPUTE_UVICE
#endif /* OBCS_UVICE_OLD */

C     Smooth the tracer sea-ice variables near the edges.
#undef OBCS_SEAICE_SMOOTH_EDGE

#endif /* ALLOW_OBCS */
#endif /* OBCS_OPTIONS_H */

C CPP options file for OBCS package
C Use this file for selecting options within the OBCS package
C Exp1 c_fjord3d: only the eastern OB (fjord mouth), Stevens (1990)
C (pkg/obcs/obcs_calc_stevens.F) with prescribed depth-mean normal velocity
C OBEu (needs ALLOW_OBCS_PRESCRIBE) and T/S relaxed toward the PRESCRIBED
C time series OBEt/OBEs (shelf pycnocline oscillation); sponge for U/V and
C T/S. Orlanski stays compiled (as in a_line2d, b_twoplume3d) but is NOT
C used at run time (useOrlanskiEast=.FALSE.): it drained the domain (eta
C -11 m in 24 h). No OBCS balance (it ignores addMass,
C pkg/obcs/obcs_balance_flow.F:300-315).

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

C Orlanski radiation code: compiled but not selected at run time in (c)
#define ALLOW_ORLANSKI

C Enable OB values to be prescribed via external fields that are read
C from a file
C  Exp1 c_fjord3d: prescribed OB u (Stevens depth-mean) and T/S time series (shelf forcing)
#define ALLOW_OBCS_PRESCRIBE

C Enable OB conditions following Stevens (1990)
#define ALLOW_OBCS_STEVENS
C R37 (2026-11-06, notes/pd1_b_twoplume3d.md V3): the Stevens baroclinic OB velocity
C is copied from the interior column (pkg/obcs/obcs_calc_stevens.F:39-42, 193-207);
C the default (the OB's own previous-step u, OBEuStevens) grew a spurious
C overturning cell at the OB in (b) and diverged in the shelf-forced (c) runs.
#define OBCS_STEVENS_USE_INTERIOR_VELOCITY

C Allow sponge layer treatment of open boundary conditions
#define ALLOW_OBCS_SPONGE

C Include hooks to sponge layer treatment of pkg/seaice variables
#undef ALLOW_OBCS_SEAICE_SPONGE

C balance barotropic velocity
#undef ALLOW_OBCS_BALANCE

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

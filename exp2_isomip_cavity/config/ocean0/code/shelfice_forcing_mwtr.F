#include "SHELFICE_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: SHELFICE_FORCING_MWTR
C     !INTERFACE:
      SUBROUTINE SHELFICE_FORCING_MWTR(
     U                    gTr_arr,
     I                    iMin,iMax,jMin,jMax, kLev, bi, bj,
     I                    myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R SHELFICE_FORCING_MWTR                (Exp. 2, ISOMIP+)
C     | o Meltwater tracer source under the ice shelf.
C     *==========================================================*
C     | Adds to the tendency gTr_arr of a passive tracer C the
C     | virtual source of meltwater volume fraction
C     |   F_C = -shelfIceFreshWaterFlux / rhoConstFresh   [m/s]
C     | (shelfIceFreshWaterFlux < 0 = melting, kg/m^2/s;
C     |  pkg/shelfice/shelfice_thermodynamics.F:645-660), so that
C     | the tracer inventory sum(C dV) is the meltwater volume
C     | (freshwater density rhoConstFresh = 1000 kg/m^3, AD16
C     | Table 4) that has entered the ocean.
C     | With useRealFreshWaterFlux=.FALSE. there is no volume flux
C     | and pkg/ptracers' PTRACERS_EvPrRn route
C     | (ptracers_forcing_surf.F:113-132) is inactive, so the
C     | source is added here, distributed in the vertical exactly
C     | like the salt forcing of S/R SHELFICE_FORCING_S
C     | (pkg/shelfice/shelfice_forcing.F:122-224): over the top
C     | wet cell kTopC and, with SHELFICEboundaryLayer, the part
C     | of cell kTopC+1 that completes a full cell thickness
C     | (Losch 2008, sect. 2.2).
C     | Called from PTRACERS_APPLY_FORCING (Exp. 2 copy).
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SHELFICE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gTr_arr   :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     kLev      :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gTr_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER kLev, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

#ifdef ALLOW_SHELFICE
C     !LOCAL VARIABLES:
C     i,j     :: Loop counters
C     kp1,km1 :: index of next/previous level
C     drLoc   :: fractional cell width of boundary layer
C     fC      :: meltwater volume flux into the ocean [m/s]
C     gCloc   :: local tendency in boundary layer
      INTEGER i, j
      INTEGER Kp1, Km1
      _RS     drLoc
      _RL     fC, gCloc
CEOP

      IF ( SHELFICEboundaryLayer ) THEN
C--   same weights as SHELFICE_FORCING_S (shelfice_forcing.F:180-206)
       DO j=1,sNy
        DO i=1,sNx
         IF ( kLev .LT. Nr .AND. kLev .EQ. kTopC(I,J,bi,bj) ) THEN
          kp1 = MIN(kLev+1,Nr)
          drLoc = drF(kLev)*( 1. _d 0 - _hFacC(I,J,kLev,bi,bj) )
          drLoc = MIN( drLoc, drF(Kp1) * _hFacC(I,J,Kp1,bi,bj) )
          drLoc = MAX( drLoc, 0. _d 0)
          fC    = -shelfIceFreshWaterFlux(i,j,bi,bj)/rhoConstFresh
          gCloc = fC/( drF(kLev)*_hFacC(I,J,kLev,bi,bj)+drLoc )
          gTr_arr(i,j) = gTr_arr(i,j) + gCloc
         ELSEIF ( kLev .GT. 1 .AND. kLev-1 .EQ. kTopC(I,J,bi,bj) ) THEN
          km1 = MAX(kLev-1,1)
          drLoc = drF(km1)*( 1. _d 0 - _hFacC(I,J,km1,bi,bj) )
          drLoc = MIN( drLoc, drF(kLev) * _hFacC(I,J,kLev,bi,bj) )
          drLoc = MAX( drLoc, 0. _d 0)
          fC    = -shelfIceFreshWaterFlux(i,j,bi,bj)/rhoConstFresh
          gCloc = fC/( drF(km1)*_hFacC(I,J,km1,bi,bj)+drLoc )
          gTr_arr(i,j) = gTr_arr(i,j) + gCloc
     &         * drLoc*recip_drF(kLev)* _recip_hFacC(i,j,kLev,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ELSE
C--   no boundary layer: all in the top wet cell, as the salt forcing
C     through surfaceForcingS (shelfice_forcing_surf.F:88-98, applied at
C     kSurfC with recip_drF*recip_hFacC)
       DO j=1,sNy
        DO i=1,sNx
         IF ( kTopC(I,J,bi,bj).NE.0 .AND.
     &        kLev .EQ. kTopC(I,J,bi,bj) ) THEN
          fC = -shelfIceFreshWaterFlux(i,j,bi,bj)/rhoConstFresh
          gTr_arr(i,j) = gTr_arr(i,j)
     &         + fC*recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

#endif /* ALLOW_SHELFICE */
      RETURN
      END

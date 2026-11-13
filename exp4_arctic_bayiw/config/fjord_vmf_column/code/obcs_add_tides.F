C--- Exp. 4 header (fjord_vmf, Task 4 PD1, Esteban 2026-11-12)
C     Modified copy of pkg/obcs/obcs_add_tides.F, MITgcm checkpoint69j, git 8e1c4076 (2026-01-12)
C     sha256 of the original: bb3ca4686a292f46de2a2a16f2c5728435efd449d62ab5e89456a05614245ea8
C     Changes, all between 'C--- Exp. 4 BEGIN' / 'C--- Exp. 4 END' lines:
C     1. with EXP4_TIDE_RAMP_S defined (OBCS_OPTIONS.h), the tidal velocity is
C        multiplied by r(t) = (1 - cos(pi t / EXP4_TIDE_RAMP_S)) / 2 for
C        t < EXP4_TIDE_RAMP_S and 1 after (one M2 cycle; V3, PD1 decision 4).
C     Everything else is the original file, line for line.
C--- Exp. 4 header end
#include "OBCS_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: OBCS_ADD_TIDES

C     !INTERFACE:
      SUBROUTINE OBCS_ADD_TIDES( myTime, myIter, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE OBCS_ADD_TIDES
C     | o Add barotropic tidal velocity to OB value
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "OBCS_PARAMS.h"
#include "OBCS_GRID.h"
#include "OBCS_FIELDS.h"

C     !INPUT/OUTPUT PARAMETERS:
      _RL myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_OBCS
#ifdef ALLOW_OBCS_TIDES
C     !FUNCTIONS:

C     !LOCAL VARIABLES:
C     bi, bj       :: tile indices
C     i,j,k        :: loop indices
C     iB, jB       :: local index of open boundary
C     msgBuf       :: Informational/error message buffer
      INTEGER bi, bj
      INTEGER i, j, k, iB, jB
      INTEGER td
      _RL  timePhase, timeCos, timeSin
C--- Exp. 4 BEGIN (ramp factor)
#ifdef EXP4_TIDE_RAMP_S
      _RL  rampFac
#endif
C--- Exp. 4 END
      _RL  uTide, vTide

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('OBCS_ADD_TIDES',myThid)
#endif

C--   Add tidal currents:
      DO td=1,OBCS_nTidalComp
       IF ( OBCS_tidalPeriod(td) .NE. zeroRL ) THEN

        timePhase = myTime * twoRL * PI / OBCS_tidalPeriod(td)
        timeCos = COS( timePhase )
        timeSin = SIN( timePhase )
C--- Exp. 4 BEGIN (V3: ramp of the tide from 0 to full amplitude)
#ifdef EXP4_TIDE_RAMP_S
C     half-cosine ramp of the tidal velocity over EXP4_TIDE_RAMP_S seconds
        rampFac = oneRL
        IF ( myTime .LT. EXP4_TIDE_RAMP_S ) rampFac =
     &       halfRL*( oneRL - COS( PI*myTime/EXP4_TIDE_RAMP_S ) )
        timeCos = timeCos*rampFac
        timeSin = timeSin*rampFac
#endif
C--- Exp. 4 END
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)

#ifdef ALLOW_OBCS_EAST
          IF ( tileHasOBE(bi,bj) ) THEN
           DO j=1-OLy,sNy+OLy
            iB = OB_Ie(j,bi,bj)
            IF ( iB.NE.OB_indexNone ) THEN
             uTide = OBE_uTideCs(j,td,bi,bj) * timeCos
     &             + OBE_uTideSn(j,td,bi,bj) * timeSin
             vTide = OBE_vTideCs(j,td,bi,bj) * timeCos
     &             + OBE_vTideSn(j,td,bi,bj) * timeSin
             DO k=1,Nr
              OBEu(j,k,bi,bj) = OBEu(j,k,bi,bj)
     &                        + uTide * maskW(iB,j,k,bi,bj)
              OBEv(j,k,bi,bj) = OBEv(j,k,bi,bj)
     &                        + vTide * maskS(iB,j,k,bi,bj)
             ENDDO
            ENDIF
           ENDDO
          ENDIF
#endif /* ALLOW_OBCS_EAST */

#ifdef ALLOW_OBCS_WEST
          IF ( tileHasOBW(bi,bj) ) THEN
           DO j=1-OLy,sNy+OLy
            iB = OB_Iw(j,bi,bj)
            IF ( iB.NE.OB_indexNone ) THEN
             uTide = OBW_uTideCs(j,td,bi,bj) * timeCos
     &             + OBW_uTideSn(j,td,bi,bj) * timeSin
             vTide = OBW_vTideCs(j,td,bi,bj) * timeCos
     &             + OBW_vTideSn(j,td,bi,bj) * timeSin
             DO k=1,Nr
              OBWu(j,k,bi,bj) = OBWu(j,k,bi,bj)
     &                        + uTide * maskW(1+iB,j,k,bi,bj)
              OBWv(j,k,bi,bj) = OBWv(j,k,bi,bj)
     &                        + vTide * maskS(iB,j,k,bi,bj)
             ENDDO
            ENDIF
           ENDDO
          ENDIF
#endif /* ALLOW_OBCS_WEST */

#ifdef ALLOW_OBCS_NORTH
          IF ( tileHasOBN(bi,bj) ) THEN
           DO i=1-OLx,sNx+OLx
            jB = OB_Jn(i,bi,bj)
            IF ( jB.NE.OB_indexNone ) THEN
             uTide = OBN_uTideCs(i,td,bi,bj) * timeCos
     &             + OBN_uTideSn(i,td,bi,bj) * timeSin
             vTide = OBN_vTideCs(i,td,bi,bj) * timeCos
     &             + OBN_vTideSn(i,td,bi,bj) * timeSin
             DO k=1,Nr
              OBNu(i,k,bi,bj) = OBNu(i,k,bi,bj)
     &                        + uTide * maskW(i,jB,k,bi,bj)
              OBNv(i,k,bi,bj) = OBNv(i,k,bi,bj)
     &                        + vTide * maskS(i,jB,k,bi,bj)
             ENDDO
            ENDIF
           ENDDO
          ENDIF
#endif /* ALLOW_OBCS_NORTH */

#ifdef ALLOW_OBCS_SOUTH
          IF ( tileHasOBS(bi,bj) ) THEN
           DO i=1-OLx,sNx+OLx
            jB = OB_Js(i,bi,bj)
            IF ( jB.NE.OB_indexNone ) THEN
             uTide = OBS_uTideCs(i,td,bi,bj) * timeCos
     &             + OBS_uTideSn(i,td,bi,bj) * timeSin
             vTide = OBS_vTideCs(i,td,bi,bj) * timeCos
     &             + OBS_vTideSn(i,td,bi,bj) * timeSin
             DO k=1,Nr
              OBSu(i,k,bi,bj) = OBSu(i,k,bi,bj)
     &                        + uTide * maskW(i,jB,k,bi,bj)
              OBSv(i,k,bi,bj) = OBSv(i,k,bi,bj)
     &                        + vTide * maskS(i,1+jB,k,bi,bj)
             ENDDO
            ENDIF
           ENDDO
          ENDIF
#endif /* ALLOW_OBCS_SOUTH */

C    end bi,bj loops
         ENDDO
        ENDDO

C    end if (tidalPeriod <> 0) and end loop on td (= tidal comp. index)
       ENDIF
      ENDDO

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('OBCS_ADD_TIDES',myThid)
#endif

#endif /* ALLOW_OBCS_TIDES */
#endif /* ALLOW_OBCS */

      RETURN
      END

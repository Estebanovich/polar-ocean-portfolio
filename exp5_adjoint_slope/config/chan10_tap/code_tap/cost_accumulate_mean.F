#include "COST_OPTIONS.h"
#ifdef ALLOW_GMREDI
# include "GMREDI_OPTIONS.h"
#endif

C--- Exp. 5: replaces pkg/cost/cost_accumulate_mean.F. pkg/cost/cost_tile.F calls it once per
C    time step while myTime > endTime - lastinterval (ALLOW_COST_TEST), after the step, with the
C    state that model/src/forward_step.F then hands to DO_THE_MODEL_IO. Instead of the time means
C    of the original, it accumulates J, the onshore heat transport across the v face of global
C    row jqb_glob (spec Sec. 5; option A, notes/v4_revision_posthoc.md), into objf_test:
C      J = -(rho0 cp / N) sum_steps sum_i sum_k [ (vVel hFacS drF + PsiY(k+1) - PsiY(k))
C            * (thetabar - thetaF) dxG + Fgm ] * maskS * jScale         [TW]
C    thetabar = (theta(j) + theta(j-1))/2; the bolus transport per unit width is
C    PsiY(k+1) - PsiY(k) (pkg/gmredi/gmredi_residual_flow.F:76-81, flipSign4LHCoord = +1), used
C    only with GM_AdvForm; Fgm is the GM/Redi tensor (Redi, in the advective form) meridional
C    theta flux through the face, as pkg/gmredi/gmredi_ytransport.F:128-186 computes it:
C      Fgm = -yA [ Kvy (theta(j) - theta(j-1)) / dyC + Kvz dTdz ],  yA = dxG drF hFacS,
C    the Kvz term only with GM_ExtraDiag. J > 0 is southward (onshore).
C    Python twin: analysis/heat_transport.py.
      SUBROUTINE COST_ACCUMULATE_MEAN( myThid )
      IMPLICIT NONE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"
#ifdef ALLOW_GMREDI
# include "GMREDI.h"
#endif
#include "cost.h"
#include "EXP5_COST.h"

      INTEGER myThid

#ifdef ALLOW_COST_TEST
      INTEGER i, j, k, km1, kp1, jg
      INTEGER bi, bj
      _RL deltaTfrac, flux, trans, tbar, yA, fgm, dTdz, maskp1
      _RL psiBelow

      deltaTfrac = deltaTClock/lastinterval

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j = 1, sNy
         jg = myYGlobalLo - 1 + (bj-1)*sNy + j
         IF ( jg .EQ. jqb_glob ) THEN
          DO k = 1, Nr
           km1 = MAX(k-1, 1)
           kp1 = MIN(k+1, Nr)
           maskp1 = 1. _d 0
           IF ( k .GE. Nr ) maskp1 = 0. _d 0
           DO i = 1, sNx
            tbar = halfRL*( theta(i,j,k,bi,bj) + theta(i,j-1,k,bi,bj) )
            trans = vVel(i,j,k,bi,bj)*hFacS(i,j,k,bi,bj)*drF(k)
#if (defined ALLOW_GMREDI && defined GM_BOLUS_ADVEC)
            IF ( useGMRedi .AND. GM_AdvForm ) THEN
             psiBelow = 0. _d 0
             IF ( k .LT. Nr ) psiBelow = GM_PsiY(i,j,kp1,bi,bj)
             trans = trans + psiBelow - GM_PsiY(i,j,k,bi,bj)
            ENDIF
#endif
            flux = trans*( tbar - thetaF )*dxG(i,j,bi,bj)
#ifdef ALLOW_GMREDI
            IF ( useGMRedi ) THEN
             yA = dxG(i,j,bi,bj)*drF(k)*hFacS(i,j,k,bi,bj)
             fgm = -yA*Kvy(i,j,k,bi,bj)*recip_dyC(i,j,bi,bj)
     &           *( theta(i,j,k,bi,bj) - theta(i,j-1,k,bi,bj) )
# ifdef GM_EXTRA_DIAGONAL
             IF ( GM_ExtraDiag ) THEN
              dTdz = halfRL*(
     &          halfRL*recip_drC(k)*(
     &            maskC(i,j-1,km1,bi,bj)*maskC(i,j-1,k,bi,bj)
     &             *( theta(i,j-1,km1,bi,bj) - theta(i,j-1,k,bi,bj) )
     &           +maskC(i,j,km1,bi,bj)*maskC(i,j,k,bi,bj)
     &             *( theta(i,j,km1,bi,bj) - theta(i,j,k,bi,bj) ) )
     &         +halfRL*recip_drC(kp1)*(
     &            maskC(i,j-1,k,bi,bj)*maskC(i,j-1,kp1,bi,bj)*maskp1
     &             *( theta(i,j-1,k,bi,bj) - theta(i,j-1,kp1,bi,bj) )
     &           +maskC(i,j,k,bi,bj)*maskC(i,j,kp1,bi,bj)*maskp1
     &             *( theta(i,j,k,bi,bj) - theta(i,j,kp1,bi,bj) ) ) )
              fgm = fgm - yA*Kvz(i,j,k,bi,bj)*dTdz
             ENDIF
# endif
             flux = flux + fgm
            ENDIF
#endif
            objf_test(bi,bj) = objf_test(bi,bj)
     &         - rhoConst*HeatCapacity_Cp*flux*maskS(i,j,k,bi,bj)
     &           *jScale*deltaTfrac
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDDO
#endif /* ALLOW_COST_TEST */

      RETURN
      END

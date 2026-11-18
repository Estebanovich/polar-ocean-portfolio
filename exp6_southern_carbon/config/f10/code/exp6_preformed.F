#include "GCHEM_OPTIONS.h"

CBOP
C !ROUTINE: EXP6_PREFORMED
C !INTERFACE: ==========================================================
      SUBROUTINE EXP6_PREFORMED( bi, bj, myThid )
C !DESCRIPTION:
C     Exp. 6 (portfolio): preformed phosphate. Ptracer 6 is set equal to
C     PO4 (ptracer 3) in the top model level every time step and has no
C     biological source below, so PO4 - PO4pre is the regenerated
C     phosphate (Lauderdale et al. 2013; Munday et al. 2014).
C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
C !INPUT PARAMETERS: ===================================================
      INTEGER bi, bj, myThid
CEOP
      INTEGER i, j
      DO j=1,sNy
       DO i=1,sNx
        pTracer(i,j,1,bi,bj,6) = pTracer(i,j,1,bi,bj,3)
       ENDDO
      ENDDO
      RETURN
      END

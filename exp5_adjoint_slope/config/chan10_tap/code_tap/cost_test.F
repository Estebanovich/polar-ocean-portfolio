#include "COST_OPTIONS.h"

C--- Exp. 5: replaces pkg/cost/cost_test.F. pkg/cost/cost_final.F calls COST_TEST before summing
C    mult_test*objf_test into fc. The original overwrites objf_test with a point temperature;
C    here objf_test already holds J (accumulated in cost_accumulate_mean.F), so this routine
C    leaves it unchanged.
      SUBROUTINE COST_TEST( myThid )
      IMPLICIT NONE
      INTEGER myThid
      RETURN
      END

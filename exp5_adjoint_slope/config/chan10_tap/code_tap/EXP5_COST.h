C     Exp. 5 cost function J (spec Sec. 5): the v-face row of the J line (global, 1-based;
C     = j_qb + 1 of cases.yaml / case.json), the freezing reference temperature (degC) and
C     the unit factor (W -> TW).
      INTEGER jqb_glob
      _RL thetaF, jScale
      PARAMETER ( jqb_glob = 18 )
      PARAMETER ( thetaF = -1.9 _d 0, jScale = 1. _d -12 )

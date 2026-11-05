CBOP
C    !ROUTINE: SIZE.h
C    !INTERFACE:
C    include SIZE.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | SIZE.h Declare size of underlying computational grid.
C     | Exp1 PD2 convergence case conv_dx2p5 (ruling R20): 2-D x-z, serial.
C     | Nx = 216 = len(stretched_x(dx0=2.5, dxmax=200, r=1.03,
C     |   L=20000)); Nr = H/dz = 500/2.5 = 200. Other headers are
C     | symlinks to config/a_line2d/code. run_campaign.py checks
C     | case.json (nx, ny, nz) against this file before running.
C     *==========================================================*
C     \ev
CEOP
C     Voodoo numbers controlling data layout.
C     sNx :: No. X points in sub-grid.
C     sNy :: No. Y points in sub-grid.
C     OLx :: Overlap extent in X.
C     OLy :: Overlat extent in Y.
C     nSx :: No. sub-grids in X.
C     nSy :: No. sub-grids in Y.
C     nPx :: No. of processes to use in X.
C     nPy :: No. of processes to use in Y.
C     Nx  :: No. points in X for the total domain.
C     Ny  :: No. points in Y for the total domain.
C     Nr  :: No. points in Z for full process domain.
      INTEGER sNx
      INTEGER sNy
      INTEGER OLx
      INTEGER OLy
      INTEGER nSx
      INTEGER nSy
      INTEGER nPx
      INTEGER nPy
      INTEGER Nx
      INTEGER Ny
      INTEGER Nr
      PARAMETER (
     &           sNx = 216,
     &           sNy =   1,
     &           OLx =   3,
     &           OLy =   3,
     &           nSx =   1,
     &           nSy =   1,
     &           nPx =   1,
     &           nPy =   1,
     &           Nx  = sNx*nSx*nPx,
     &           Ny  = sNy*nSy*nPy,
     &           Nr  = 200)

C     MAX_OLX :: Set to the maximum overlap region size of any array
C     MAX_OLY    that will be exchanged. Due to the way packages
C                employ the exchange routines, these need to be set
C                appropriately.
      INTEGER MAX_OLX
      INTEGER MAX_OLY
      PARAMETER ( MAX_OLX = OLx,
     &            MAX_OLY = OLy )


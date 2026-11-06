CBOP
C    !ROUTINE: SIZE.h
C    !INTERFACE:
C    include SIZE.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | SIZE.h Declare size of underlying computational grid.
C     | Exp1 b_twoplume3d: 3-D box, 5 MPI ranks split in y.
C     | Nx = 100 (1.0 km / dx 10 m), Ny = 120 (1.2 km / dy 10 m),
C     | Nr = H/dz = 300/10 = 30, from exp1_fjord_plumes/cases.yaml
C     | (b_twoplume3d). Ny includes the two land rows j = 0 and
C     | Ny-1 (side walls) and Nx the land column Nx-1 (cases.yaml
C     | land_j_rows / land_i_cols). D2 (2026-11-06, dx = dy = dz = 5 m): sNx = 200,
C     | sNy = 240/5 = 48, Nr = 60. (10 m production: sNx = 100, sNy = 120/5)
C     | = 24; rank 0 holds land row j = 0, rank 4 holds j = Ny-1.
C     | run_campaign.py checks case.json (nx, ny, nz) against
C     | this file before running; np must be nPx*nPy = 5.
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
     &           sNx = 200,
     &           sNy =  48,
     &           OLx =   3,
     &           OLy =   3,
     &           nSx =   1,
     &           nSy =   1,
     &           nPx =   1,
     &           nPy =   5,
     &           Nx  = sNx*nSx*nPx,
     &           Ny  = sNy*nSy*nPy,
     &           Nr  =  60)

C     MAX_OLX :: Set to the maximum overlap region size of any array
C     MAX_OLY    that will be exchanged. Due to the way packages
C                employ the exchange routines, these need to be set
C                appropriately.
      INTEGER MAX_OLX
      INTEGER MAX_OLY
      PARAMETER ( MAX_OLX = OLx,
     &            MAX_OLY = OLy )


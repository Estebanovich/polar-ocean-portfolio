CBOP
C    !ROUTINE: SIZE.h
C    !INTERFACE:
C    include SIZE.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | SIZE.h Declare size of underlying computational grid.
C     | Exp1 c_fjord3d: 3-D hydrostatic fjord, 5 MPI ranks in x.
C     | Nx = 200 (40 km / dx 200 m), Ny = 30 (6 km / dy 200 m),
C     | Nr = H/dz = 500/10 = 50, from exp1_fjord_plumes/cases.yaml
C     | (c_fjord3d). Ny includes the two land rows j = 0 and Ny-1
C     | (side walls) and Nx the land column Nx-1 (cases.yaml
C     | land_j_rows / land_i_cols). Tiles: sNx = 200/5 = 40,
C     | sNy = 30; rank 0 holds the ice face, rank 4 the mouth (OB).
C     | np must be nPx*nPy = 5 (test_c_grid_matches_size_h_and_5_ranks).
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
     &           sNx = 40,
     &           sNy =  30,
     &           OLx =   3,
     &           OLy =   3,
     &           nSx =   1,
     &           nSy =   1,
     &           nPx =   5,
     &           nPy =   1,
     &           Nx  = sNx*nSx*nPx,
     &           Ny  = sNy*nSy*nPy,
     &           Nr  =  50)

C     MAX_OLX :: Set to the maximum overlap region size of any array
C     MAX_OLY    that will be exchanged. Due to the way packages
C                employ the exchange routines, these need to be set
C                appropriately.
      INTEGER MAX_OLX
      INTEGER MAX_OLY
      PARAMETER ( MAX_OLX = OLx,
     &            MAX_OLY = OLy )


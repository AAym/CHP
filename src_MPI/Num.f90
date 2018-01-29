MODULE Num

  ! Description :
  !   Ce module est créé afin d'implémenter la transition entre trois indicages.
  !      * L'indicage DIRECT : (x,y)
  !      * L'indicage LOCAL : (i,j)
  !      * L'indicage GLOBAL : (k)
  !   Schéma : cf NoteProjet.txt

  !Use MPI

  Use Data

  IMPLICIT NONE


CONTAINS

  ! Local -> Global (Fonctionnel)
  Function Global(i,j) Result(k)

    Integer, Intent(in) :: i,j
    Integer :: k

    k = Nx*(j-1)+i

  End Function


  ! Global -> Local (Fonctionnel)
  Function Local(k) Result(IJ)

    Integer, Intent(in) :: k
    Integer, Dimension(2) :: IJ

    If (mod(k,Nx)==0) Then
      IJ(1) = Nx
      IJ(2) = k/Nx
    Else
      IJ(1) = mod(k,Nx)
      IJ(2) = k/Nx+1
    End If

  End Function


  ! Global -> Direct (A verifier)
  Function Direct(k) Result(XY)

    Integer, Intent(in) :: k
    Real, Dimension(2) :: IJ,XY

    IJ = Local(k)
    XY(1) = Lx/(Nx+1)*IJ(1)
    XY(2) = Ly/(Ny+1)*IJ(2)

  End Function

END MODULE Num

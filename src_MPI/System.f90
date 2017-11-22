MODULE System

  ! Description :
  !   Ce module permet d'implémenter le système matriciel qu'on cherchera ensuite à résoudre
  !   Il implémente également la solution initiale


  !Use MPI

  Use Data
  Use Num
  Use Func

  IMPLICIT NONE



CONTAINS

  ! Definit la solution du initiale du système : U0
  Subroutine Init(U)

    Real, Dimension(Nx*Ny), Intent(out):: U
    Integer :: k
    Real, Dimension(2) :: XY
    Integer, Dimension(2) :: IJ

    U = 1

  End Subroutine


  ! Definit le second membre de notre système matriciel : F + U + CdB (Fonctionnel)
  Subroutine BuiltSecondMembre(SecondMembre,U,t)

    Real, Dimension(Nx*Ny), Intent(out):: SecondMembre
    Real, Dimension(Nx*Ny), Intent(in):: U
    Real :: t
    Integer :: k
    Real, Dimension(2) :: XY
    Integer, Dimension(2) :: IJ

    Do k = 1,(Nx*Ny)
      XY = Direct(k)
      IJ = Local(k)
      SecondMembre(k) = dt*F(XY(1),XY(2),t) + U(k)
      If (IJ(1)==1) Then
        SecondMembre(k) = SecondMembre(k) + D*dt/(dx*dx)*G(0.,XY(2))
      End If
      If (IJ(2)==1) Then
        SecondMembre(k) = SecondMembre(k) + D*dt/(dy*dy)*H(XY(1),0.)
      End If
      If (IJ(1)==Nx) Then
        SecondMembre(k) = SecondMembre(k) + D*dt/(dx*dx)*G(Lx,XY(2))
      End If
      If (IJ(2)==Ny) Then
        SecondMembre(k) = SecondMembre(k) + D*dt/(dy*dy)*H(XY(1),Ly)
      End If
!      Print *, "k = ", k, "SecondMembre(k) = ", SecondMembre(k)
    End Do

    !Print *, "SecondMembre = ", SecondMembre

  End Subroutine


  ! Définit le produit matriciel A*U (Fonctionnel)
  Function ProdMat(U) Result(AU)

    Real, Dimension(Nx*Ny), Intent(in):: U
    Real, Dimension(Nx*Ny) :: AU
    Integer :: k
    Real, Dimension(2) :: XY
    Integer, Dimension(2) :: IJ

    Do k = 1,(Nx*Ny)
      IJ = Local(k)
      AU(k) = ( 1 + 2*D*dt/(dx*dx) + 2*D*dt/(dy*dy) )* U(k)
      If (IJ(1)/=1) Then
        AU(k) = AU(k) - D*dt/(dx*dx)*U(Global(IJ(1)-1,IJ(2)))
      End If
      If (IJ(2)/=1) Then
        AU(k) = AU(k) - D*dt/(dy*dy)*U(Global(IJ(1),IJ(2)-1))
      End If
      If (IJ(1)/=Nx) Then
        AU(k) = AU(k) - D*dt/(dx*dx)*U(Global(IJ(1)+1,IJ(2)))
      End If
      If (IJ(2)/=Ny) Then
        AU(k) = AU(k) - D*dt/(dy*dy)*U(Global(IJ(1),IJ(2)+1))
      End If
!      Print *, "k = ", k, "AU(k) = ", Au(k)
    End Do

    !Print *, "AU = ", AU

  End Function


END MODULE System

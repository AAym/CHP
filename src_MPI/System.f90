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

  ! Definit le second membre de notre système matriciel : F + U + CdB (Fonctionnel)
  Subroutine BuiltSecondMembre(SecondMembre,U,t,i1,iN,Bord_inf,Bord_sup)
	

    Integer,intent(in)::i1, iN
    Real, Dimension(i1*Nx+1:iN*Nx), Intent(out):: SecondMembre
    Real, Dimension(i1*Nx+1:iN*Nx), Intent(in):: U
    Real, Dimension(Nx),Intent(in)::Bord_inf,Bord_sup
    Real, intent(in) :: t
    Integer :: k
    Real, Dimension(2) :: XY
    Integer, Dimension(2) :: IJ

	
    Do k = i1*Nx+1,iN*Nx
      !print*, "k= ", k
      !XY = Direct(k)
      IJ = Local(k)
      XY(1)=IJ(1)/dx
      XY(2)=IJ(2)/dy
      SecondMembre(k) = dt*F(XY(1),XY(2),t) + U(k)
      If (IJ(1)==1) Then
        SecondMembre(k) = SecondMembre(k) + D*dt/(dx*dx)*H(0.,XY(2))
      End If
      If (IJ(2)==i1+1) Then
        SecondMembre(k) = SecondMembre(k) + D*dt/(dy*dy)*Bord_inf(IJ(1))
      End If
      If (IJ(1)==Nx) Then
        SecondMembre(k) = SecondMembre(k) + D*dt/(dx*dx)*H(Lx,XY(2))
      End If
      If (IJ(2)==iN+1) Then
        SecondMembre(k) = SecondMembre(k) + D*dt/(dy*dy)*Bord_sup(IJ(1))
      End If
!      Print *, "k = ", k, "SecondMembre(k) = ", SecondMembre(k)
    End Do

    !Print *, "SecondMembre = ", SecondMembre

  End Subroutine


  ! Définit le produit matriciel A*U (Fonctionnel)
  Function ProdMat(U, i1, iN) Result(AU)

    Real, Dimension(i1*Nx+1:iN*Nx+Nx), Intent(in):: U
    Integer, intent(in)::i1, iN
    Real, Dimension(i1*Nx+1:iN*Nx+Nx) :: AU
    Integer :: k
    Real, Dimension(2) :: XY
    Integer, Dimension(2) :: IJ

    Do k = i1*Nx+1, iN*Nx+Nx
      IJ = Local(k)
      AU(k) = ( 1 + 2*D*dt/(dx*dx) + 2*D*dt/(dy*dy) )* U(k)
      If (IJ(1)/=1) Then
        AU(k) = AU(k) - D*dt/(dx*dx)*U(Global(IJ(1)-1,IJ(2)))
      End If
      If (IJ(2)/=i1+1) Then
        AU(k) = AU(k) - D*dt/(dy*dy)*U(Global(IJ(1),IJ(2)-1))
      End If
      If (IJ(1)/=Nx) Then
        AU(k) = AU(k) - D*dt/(dx*dx)*U(Global(IJ(1)+1,IJ(2)))
      End If
      If (IJ(2)/=iN+1) Then
        AU(k) = AU(k) - D*dt/(dy*dy)*U(Global(IJ(1),IJ(2)+1))
      End If
!      Print *, "k = ", k, "AU(k) = ", Au(k)
    End Do

    !Print *, "AU = ", AU

  End Function


END MODULE System

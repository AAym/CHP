MODULE Func

  ! Description :
  !   Ce module est créé afin d'implémenter les fonctions second membre et de bord
  !   ainsi que la solution exacte si elle existe

  ! Note :
  !   Penser à vérifier les fonctions


  !Use MPI

  Use Data

  IMPLICIT NONE


CONTAINS

  !Fonction f
  Real Function F(x,y,t)

    Real, Intent(in) :: x, y, t

    If (SystType == "Stationnaire") Then
      f = 2*(y-y*y + x-x*x)
    Else If (SystType == "Sinusoidal") Then
      f = sin(x) + cos(y)
    Else If (SystType == "Exponentiel") Then
      f = exp(-(x-Lx/2)**2) * exp(-(y-Ly/2)**2) * cos(3.14/2*t)
    End If

  End Function


  !Fonction g
  Real Function G(x,y)

    Real, Intent(in) :: x,y

    If (SystType == "Stationnaire") Then
      g = 0
    Else If (SystType == "Sinusoidal") Then
      g = sin(x) + cos(y)
    Else If (SystType == "Exponentiel") Then
      g = 0
    End If

  End Function


  !Fonction h
  Real Function H(x,y)

    Real, Intent(in) :: x,y

    If (SystType == "Stationnaire") Then
      h = 0
    Else If (SystType == "Sinusoidal") Then
      h = sin(x) + cos(y)
    Else If (SystType == "Exponentiel") Then
      h = 1
    End If

  End Function


  !Solution exacte
  Function Exact(x,y) Result(Uexact)

    Real, Intent(in) :: x,y
    Real :: Uexact

    If (SystType == "Stationnaire") Then
      Uexact = x*(1-x)*(1-y)*y
    Else If (SystType == "Sinusoidal") Then
      Uexact = sin(x) + cos(y)
    Else If (SystType == "Exponentiel") Then
      Print*, "Pas de solution exacte"
    End If

  End Function

  subroutine charge(me,N,Np,i1,in)
    
    implicit none

    integer,intent(in)::me,N,Np
    integer,intent(inout)::i1,in
    integer::diveuc,reste

    diveuc=N/Np
    reste=N-Np*diveuc

    if(me<reste) then
       i1=me*(diveuc+1)+1
       in=(1+me)*(diveuc+1)
    else
       i1=1+reste+me*diveuc
       in=i1+diveuc
    end if
    i1=i1-1
    iN=iN-1
  end subroutine charge


END MODULE Func

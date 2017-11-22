PROGRAM Main

  !USE MPI

  Use Data
  Use Num
  Use Func
  Use System
  Use GradientConjugue

  IMPLICIT NONE

  ! Variables test
    Integer :: k1, k2, k3, k4
    Real, Dimension(:), Allocatable  :: Test
    Real, Dimension(:), Allocatable  :: Test2
  ! Variables résolution
    Real, Dimension(:), Allocatable :: U , SecondMembre, Uexact, Uprev
    Integer, Dimension(2) :: IJ
    Real, Dimension(2) :: XY
    Integer :: k
    Real :: t


  ! Lecture variables
    Call ReadData()
    Print *, "----------------Variables-----------------"
    Print *, "Lx = ", Lx
    Print *, "Ly = ", Ly
    Print *, "Nx = ", Nx
    Print *, "Ny = ", Ny
    Print *, "SystType = ", SystType
    Print *, "epsilon = ", epsilon
    Print *, "dx = ", dx
    Print *, "dy = ", dy
    Print *, "dt= ", dt
    Print *, "------------------------------------------"

  !Allocation
    Allocate(U(Nx*Ny), SecondMembre(Nx*Ny), Uexact(Nx*Ny), Uprev(Nx*Ny))
    Allocate(Test(Nx*Ny), Test2(Nx*Ny))


  ! Tests
    Print *, "------------------Tests-------------------"
  ! Test Numérotation
    k1 = 13
    K2 = Nx*Ny
    IJ = Local(k1)
    k3 = Global(IJ(1),IJ(2))
    IJ = Local(k2)
    k4 = Global(IJ(1),IJ(2))
    If (k1==k3 .and. k2==k4 ) Then
      Print*, "Numérotation ok"
    Else
      Print*, "Problème de numérotation"
    End If
    Print *, "------------------------------------------"

  !Résolution
    t = 0
    Call Init(U)
    UPrev = 0

    Do while (dot_product(U-UPrev,U-UPrev) > 1E-10)
      Print *, "-------------------------------------"
      Print *, "t = ", t
      Call BuiltSecondMembre(SecondMembre,U,t)
      UPrev = U
      Call GC(U,SecondMembre)
      Print *, "-------------------------------------"
      t = t + dt
    End Do

    If (SystType == "Stationnaire" .or. SystType == "Sinusoidal") Then
      Do k=1,(Nx*Ny)
        XY = Direct(k)
        Uexact(k) = Exact(XY(1),XY(2))
      End Do
      If (dot_product(U-Uexact,U-Uexact) < 1E-5) Then
          Print *, "L'algorithme a convergé vers la solution exacte"
      End if
    End If

    Print *, U
  !Deallocation
    Deallocate(U, SecondMembre, Uexact, Uprev)



END PROGRAM Main


! Tout droits réservés LG

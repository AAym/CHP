PROGRAM Main

  !USE MPI

  Use Data
  Use Num
  Use Func
  Use System
  Use GradientConjugue
  include "mpif.h"

  IMPLICIT NONE

    integer,dimension(MPI_STATUS_SIZE)::status
    integer,parameter::tag=100
    integer::statinfo

    call MPI_INIT(statinfo)
    call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)

  ! Variables test
    Integer :: k1, k2, k3, k4
    Real, Dimension(:), Allocatable  :: Test
    Real, Dimension(:), Allocatable  :: Test2
  ! Variables résolution
    Real, Dimension(:), Allocatable :: U , SecondMembre, Uprev
    Integer, Dimension(2) :: IJ
    Real, Dimension(2) :: XY
    Integer :: k
    Real :: t, max



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
  
    call charge(me, Np, Ny, i1, in)

  !Allocation
    if (me==0) then
       
       Allocate(U(Nx*(in-i1)), SecondMembre(Nx*(in-i1)), Uprev(Nx*(in-i1)))
       !Allocate(Test(Nx*Ny), Test2(Nx*Ny)) utilité ?

    else

       Allocate(U(Nx*(in-i1+1)), SecondMembre(Nx*(in-i1+1)), Uprev(Nx*(in-i1+1)))

    end if


  !Résolution
    t = 0
    Call Init(U)
    UPrev = 0

    Do while (max > 1E-1)
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

    call MPI_FINALIZE(statinfo)

END PROGRAM Main


! Tout droits réservés LG

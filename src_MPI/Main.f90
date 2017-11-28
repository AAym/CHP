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
    Real, Dimension(:), Allocatable  :: Test, Test2, Bord_inf, Bord_sup
  ! Variables résolution
    Real, Dimension(:), Allocatable :: U , SecondMembre, Uprev
    Integer, Dimension(2) :: IJ, nb_lignes
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
    Print *, "Recouvrement : ", R
    Print *, "------------------------------------------"
  
    call charge(me, Np, Ny, i1, in)

    !Allocation
    nb_lignes = in-i1+1+(R-1)
    Allocate(U(Nx*nb_lignes), SecondMembre(Nx*nb_lignes), Uprev(Nx*nb_lignes))
    Allocate(Bord_inf(Nx), Bord_sup(Nx))
  !Résolution
    t = 0
    U = 0
    UPrev = 0
    Bord_inf = 0
    Bord_sup = 0

    if(me==0) then
	do i=1,Nx
		Bord_inf(i)=h(i*dx,0)
	end do
    end if

    if (me==Np-1)
	do i=1,Nx
		Bord_sup(i)=h(i*dx,Ly)
	end do
    end if

    
    Do while (max > 1E-1)
       ! Ne pas oublier la construction/actualisation de bord_inf et bord_sup
       Call BuiltSecondMembre(SecondMembre,U,t,nb_lignes,i1,Bord_inf,Bord_sup)
       UPrev = U
       Call GC(U,SecondMembre)
       
       !! Communication
       do i=1,Nx
          Bord_inf = U( (nb_lignes-2*(R-1))*Nx + 1 + i)
          Bord_sup = U( (1+2*(R-1)-1)*Nx + 1 + i)
       end do

       Allocate(Tmp(Nx))
       if (me==0) then    
          Tmp = Bord_sup
          Call MPI_Recv(Bord_sup,Nx,MPI_REAL,me+1,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)   
          Call MPI_Send(Tmp, Nx,MPI_REAL,me+1,MPI_ANY_TAG,MPI_COMM_WORLD,statinfo)
       else if (me==Np-1) then
          Tmp = Bord_infw
          Call MPI_Recv(Bord_inf,Nx,MPI_REAL,me-1,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)
          Call MPI_Send(Tmp, Nx,MPI_REAL,me-1,MPI_ANY_TAG,MPI_COMM_WORLD,statinfo)

       else 

          Call MPI_Send(Bord_sup, Nx,MPI_REAL,me+1,MPI_ANY_TAG,MPI_COMM_WORLD,statinfo)
          Call MPI_Send(Bord_inf, Nx,MPI_REAL,me-1,MPI_ANY_TAG,MPI_COMM_WORLD,statinfo)  
          Call MPI_Recv(Bord_sup,Nx,MPI_REAL,me+1,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)   
          Call MPI_Recv(Bord_inf,Nx,MPI_REAL,me-1,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)   
       end if
     
       
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
    Deallocate(U, SecondMembre, Uprev, Bord_inf, Bord_sup, Test, Test2)

    call MPI_FINALIZE(statinfo)

END PROGRAM Main


! Tout droits réservés LG

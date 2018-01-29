PROGRAM Main

  Use MPI

  Use Data
  Use Num
  Use Func
  Use System
  Use GradientConjugue
  
  IMPLICIT NONE

    integer,dimension(MPI_STATUS_SIZE)::status
    integer,parameter::tag=100
    integer::statinfo

    

  ! Variables test
    Integer :: Np, me, i1, iN
    Real, Dimension(:), Allocatable  ::  Bord_inf, Bord_sup
  ! Variables résolution
    Real, Dimension(:), Allocatable :: U , SecondMembre, Uprev, Tmp,Uexact
    Integer, Dimension(2) :: IJ
    Real, Dimension(2) :: XY
    Integer :: k, i, j
    Real :: t, max, maxtemp


    call MPI_INIT(statinfo)
    call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)


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
  
    call charge(me, Ny, Np, i1, in)
    print*, "me : ", me
    print*, "i1 : ", i1
    print*, "iN : ", iN
    !Allocation
    !nb_lignes = in-i1+1
    !Allocate(U(Nx*nb_lignes), SecondMembre(Nx*nb_lignes), Uprev(Nx*nb_lignes),Uexact(Nx*nb_lignes))
    Allocate(U(i1*Nx+1:iN*Nx+Nx), SecondMembre(i1*Nx+1:iN*Nx+Nx), Uprev(i1*Nx+1:iN*Nx+Nx),Uexact(i1*Nx+1:iN*Nx+Nx))
    Allocate(Bord_inf(Nx), Bord_sup(Nx),Tmp(Nx))
 
    !Résolution
    t = 0
    U = 0
    UPrev = 0
    ! print*,"J'initialise les bords"
    Bord_inf = 0
    Bord_sup = 0
    max=3

    if(me==0) then
	do i=1,Nx
		Bord_inf(i)=G(i*dx,0.0)
	end do
    end if

    if (me==Np-1) then
	do i=1,Nx
		Bord_sup(i)=G(i*dx,Ly)
	end do
    end if

    
    !Do while (max > 1E-1)
    Do i=1,1
       print*, max
       ! Ne pas oublier la construction/actualisation de bord_inf et bord_sup
       !print*, "Ca tourne !!!"
       Call BuiltSecondMembre(SecondMembre,U,t,i1,iN,Bord_inf,Bord_sup)
       print*,SecondMembre
       UPrev = U

       ! Code vérifié jusqu'ici
       Call GC(U,SecondMembre, i1, iN)

       !! Communication

       !Bord_inf = U(i1*Nx+1:i1*Nx+Nx)   !if i1==0 ?
       !Bord_sup = U(iN*Nx+1:iN*Nx+Nx)   !if iN==Ny-1  ?

       if (Np /= 1) then

          if (me==0) then    
             Tmp = Bord_sup
             Call MPI_Recv(Bord_sup,Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)   
             Call MPI_Send(Tmp, Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,statinfo)

          else if (me==Np-1) then
             Tmp = Bord_inf
             Call MPI_Send(Tmp, Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,statinfo)
             Call MPI_Recv(Bord_inf,Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)
          
          else 
             Call MPI_Send(Bord_sup, Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,statinfo)
             Call MPI_Send(Bord_inf, Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,statinfo)  
             Call MPI_Recv(Bord_sup,Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)   
             Call MPI_Recv(Bord_inf,Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)   
          end if

       end if
             
       !Print *, "-------------------------------------"
       t = t + dt
       max=dot_product(U-UPrev,U-UPrev)
       if (me/=0) then
          Call MPI_Send(max, 1,MPI_REAL,0,tag,MPI_COMM_WORLD,statinfo)
       end if
       if (me==0) then
          do j=1,Np-1
             Call MPI_Recv(maxtemp,1,MPI_REAL,j,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)
             max=max+maxtemp
          end do
          do j=1,Np-1
             Call MPI_Send(max, 1,MPI_REAL,j,tag,MPI_COMM_WORLD,statinfo)
          end do
       end if
       if (me/=0) then
          Call MPI_Recv(max,1,MPI_REAL,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)
       end if
    End Do

    if (me==0) then
       If (SystType == "Stationnaire" .or. SystType == "Sinusoidal") Then
          Do k=1,iN*Nx+Nx
             XY = Local(k)
             print*,k,XY(1),XY(2)
             Uexact(k) = Exact(XY(1)/Nx,XY(2)/Ny)
          End Do
          If (dot_product(U-Uexact,U-Uexact) < 1E-5) Then
             Print *, "L'algorithme a convergé vers la solution exacte"
          End if
          print*, "Solution exacte : "
          print*,Uexact
       End If
       print*, "Solution calculée : "
       Print *, U
    end if

  !Deallocation
    Deallocate(U, SecondMembre, Uprev, Bord_inf, Bord_sup, Tmp, Uexact)

    call MPI_FINALIZE(statinfo)

END PROGRAM Main

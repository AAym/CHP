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
    Integer :: k, i
    Real :: t, max


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
    Allocate(U(i1*Nx+1:iN*Nx), SecondMembre(i1*Nx+1:iN*Nx), Uprev(i1*Nx+1:iN*Nx),Uexact(i1*Nx+1:iN*Nx))
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

    
    Do while (max > 1E-1)
       ! Ne pas oublier la construction/actualisation de bord_inf et bord_sup
       print*, "Ca tourne !!!"
       Call BuiltSecondMembre(SecondMembre,U,t,i1,iN,Bord_inf,Bord_sup)
       
       UPrev = U
       ! Code vérifié jusqu'ici
       Call GC(U,SecondMembre, i1, iN)

       !! Communication

       Bord_inf = U(i1*Nx+1:i1*Nx+Nx)
       Bord_sup = U((iN-1)*Nx+1:iN*Nx)

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
       !max=dot_product(U-UPrev,U-UPrev) ! Ce n'est pas ça qui fait planter, mais la norme n'est calculée que par proc et pas pour la matrice entière
    End Do

    !If (SystType == "Stationnaire" .or. SystType == "Sinusoidal") Then
    !  Do k=1,(Nx*Ny)
    !    XY = Direct(k)
    !    Uexact(k) = Exact(XY(1),XY(2))
    !  End Do
    !  If (dot_product(U-Uexact,U-Uexact) < 1E-5) Then
    !      Print *, "L'algorithme a convergé vers la solution exacte"
    !  End if
    !End If

    !Print *, U
  !Deallocation
    Deallocate(U, SecondMembre, Uprev, Bord_inf, Bord_sup, Tmp, Uexact)

    call MPI_FINALIZE(statinfo)

END PROGRAM Main

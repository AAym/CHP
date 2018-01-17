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
    Integer :: k1, k2, k3, k4, Np, me, i1, in,nb_lignes
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
    
    !Allocation
    print*, "ligne", in, i1, Ny
    nb_lignes = in-i1+1
    Allocate(U(Nx*nb_lignes), SecondMembre(Nx*nb_lignes), Uprev(Nx*nb_lignes),Uexact(Nx*nb_lignes))
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
		Bord_inf(i)=h(i*dx,0.0)
	end do
    end if

    if (me==Np-1) then
	do i=1,Nx
		Bord_sup(i)=h(i*dx,Ly)
	end do
    end if
    ! print*,"je rentre dans la boucle"
    
    Do while (max > 1E-1)
       ! Ne pas oublier la construction/actualisation de bord_inf et bord_sup

       !print*,"J'ai contruit 2nd membre 1"
       Call BuiltSecondMembre(SecondMembre,U,t,nb_lignes,i1,Bord_inf,Bord_sup)
       !print*,"J'ai contruit 2nd membre"
       
       UPrev = U
       !print*, "taille de U", size(U),Nx,nb_lignes, size(Uprev),size(Bord_inf), size(Bord_sup)
       Call GC(U,SecondMembre)
       ! print*,"j'ai fait le gradient conjugué"
       !! Communication
       do i=1,Nx
          !print*, " taille de Nx", Nx, size(U), (nb_lignes-2*(R-1))*Nx + 1 + i, (1+2*(R-1)-1)*Nx + 1 + i,R 
          !Bord_inf = U( (nb_lignes-2*(R-1))*Nx + 1 + i)
          !Bord_sup = U( (1+2*(R-1)-1)*Nx + 1 + i)

          Bord_inf = U(1:Nx)
          Bord_sup = U((nb_lignes-1)*Nx+1:nb_lignes*Nx)
       end do
       if (Np /= 1) then

          if (me==2) then
             Print*, "--------------------------------------------------"
             Print*, "--------------------------------------------------"
             Print*, Bord_sup
             Print*, "--------------------------------------------------"
             Print*, Bord_inf
          end if
          if (me==0) then    
             Tmp = Bord_sup
        !     print*,"Je suis le proc", me,me+1
             Call MPI_Recv(Bord_sup,Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)   
             Call MPI_Send(Tmp, Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,statinfo)
             !print*, "j'ai reçu"
          else if (me==Np-1) then
             Tmp = Bord_inf
             
       !      print*,"Je suis le proc", me,me-1
             Call MPI_Send(Tmp, Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,statinfo)
             !print*, "j'ai envoyé"
             Call MPI_Recv(Bord_inf,Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)
          
          else 
      !       print*, "Je suis perdu"
             Call MPI_Send(Bord_sup, Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,statinfo)
             Call MPI_Send(Bord_inf, Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,statinfo)  
             Call MPI_Recv(Bord_sup,Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)   
             Call MPI_Recv(Bord_inf,Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,statinfo)   
          end if

          if (me==2) then
             Print*, "--------------------------------------------------"
             Print*, "--------------------------------------------------"
             Print*, Bord_sup
             Print*, "--------------------------------------------------"
             Print*, Bord_inf
          end if
       end if
             
       Print *, "-------------------------------------"
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
    Deallocate(U, SecondMembre, Uprev, Bord_inf, Bord_sup,Tmp)

    call MPI_FINALIZE(statinfo)

END PROGRAM Main

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
    Real :: t1, t2, t_exec, t_exec_loc
    Real :: t1_com, t2_com, t_com_loc, t_com
    t_com = 0
    call cpu_time(t1)
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

    
    Do while (max > 1E-8)
       !print*, "max : ", max

       !print*, "Ca tourne !!!"
       Call BuiltSecondMembre(SecondMembre,U,t,i1,iN,Bord_inf,Bord_sup)
       !print*,SecondMembre
       UPrev = U

       ! Code vérifié jusqu'ici
       Call GC(U,SecondMembre, i1, iN)

       !! Communication
       if (me/=0)then
          Bord_inf = U(i1*Nx+1:i1*Nx+Nx)
       end if

       if (me/=Np-1)then
          Bord_sup = U(iN*Nx+1:iN*Nx+Nx)
       end if

       call cpu_time(t1_com)
   
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
       
       call cpu_time(t2_com)
       t_com_loc = t_com_loc + t2_com-t1_com
       !Print *, "-------------------------------------"
       t = t + dt
       max=dot_product(U-UPrev,U-UPrev)
       call cpu_time(t1_com)
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
    call cpu_time(t2_com)
    t_com_loc =  t_com_loc + t2_com-t1_com

    
    if (me==0) then
       If (SystType == "Stationnaire" .or. SystType == "Sinusoidal") Then
          Do k=1,iN*Nx+Nx
             XY = Local(k)
             !print*,k,XY(1),XY(2)
             Uexact(k) = Exact(XY(1)*dx,XY(2)*dy)
          End Do

          If (dot_product(U-Uexact,U-Uexact) < 1E-5) Then
             Print *, "L'algorithme a convergé vers la solution exacte"
          End if
          !print*, "Solution exacte : "
          !print*,Uexact
       End If
       !print*, "Solution calculée : "
       !Print *, U
    end if

    ! Ecriture dans un fichier
    
    open(13,file="Solution_calculee.txt",form="formatted")
    open(14,file="Solution_exact.txt",form="formatted")
    do k=i1*Nx+1,iN*Nx+Nx
       XY = Local(k)
       write(13,'(F10.5,F10.5,F10.5)')XY(1)*dx,XY(2)*dy,U(k)
       write(14,'(F10.5,F10.5,F10.5)')XY(1)*dx,XY(2)*dy,Uexact(k)
    end do
    close(13)
    close(14)
       
    !Deallocation

    Deallocate(U, SecondMembre, Uprev, Bord_inf, Bord_sup, Tmp, Uexact)

    
    call cpu_time(t2)
    t_exec_loc = t2 -t1
    print*, "Temps local proc", me,": ", t_exec_loc
    print*, "Temps local comm proc", me,": ", t_com_loc
     
    if (Np > 1) then
       call MPI_REDUCE(t_exec_loc,t_exec,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD)
       call MPI_REDUCE(t_com_loc,t_com,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD)
    
       if ( me==0) then
          print*, "Temps global : "
          print*, t_exec
          print*, "Temps comm global : "
          print*, t_com
       endif
    endif
     call MPI_FINALIZE(statinfo)

END PROGRAM Main

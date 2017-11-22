module Arret


Use Data
Use Func


implicit none

CONTAINS

  FUNCTION STOP(U,line,tps)

    implicit none

    Real,dimension(:),intent(in)::U
    integer,intent(in)::line
    real,intent(in)::tps
    integer::i
    real::test,STOP,zero

    zero=0.0d0

    STOP=0

    do i=Nx*(line-1)+2,(Nx)*(line)-1
       test = ((2*U(i)-U(i+1)-U(i-1))/dx**2+(2*U(i)-U(i+Nx)-U(i-Nx))/dy**2)-f(dx*(i-(Nx*(line-1)+1)),dy*line,tps)
       if(abs(test)>STOP) then
          STOP = abs(test)
       end if
    end do


    !! Si bord gauche
    test = (2*U((line-1)*Nx+1)-U((line-1)*Nx+2)-g(zero,dy*line))/dx**2+(2*U((line-1)*Nx+1)-U(line*Nx+1)-U((line-2)*Nx+1))/dy**2&
         &-f(dx,dy*line,tps)

    if(abs(test)>STOP) then
      STOP=abs(test)
    end if



    !! Si bord droit
    test=((2*U(line*Nx)-U(line*Nx-1)-g(Lx,dy*line))/dx**2+(2*U(line*Nx)-U(line*Nx+Nx)-U((line-1)*Nx))/dy**2)-f(Lx-dx,dy*line,tps)

    if(abs(test)>STOP) then
      STOP=abs(test)
    end if




  end FUNCTION STOP





end module Arret

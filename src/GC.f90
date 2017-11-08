MODULE GradientConjugue

  Use Data
  Use Num
  Use Func
  Use System

  IMPLICIT NONE



CONTAINS

  Subroutine GC(U,SecondMembre)

    !In
      Real, Dimension(Nx*Ny), Intent(inout):: U
      Real, Dimension(Nx*Ny), Intent(in):: SecondMembre
    !Num
      Integer :: k
      Real, Dimension(2) :: XY
      Integer, Dimension(2) :: IJ
    !VarGC
      Real, Dimension(Nx*Ny) :: R, RPrev, P, Z
      Real :: alpha, beta

    !Initilialisation
      R = SecondMembre - ProdMat(U)
      Rprev = R
      P = R

      Do while ( dot_product(R,R) > epsilon)
        Print *, "ProdScal = ", dot_product(R,R)
        Z = ProdMat(P)
        Print *, "Z = ", Z
        alpha = dot_product(R,R)/dot_product(P,Z)
        U = U + alpha*P
        Rprev = R
        R = R - alpha*Z
        beta = dot_product(R,R)/dot_product(RPrev,RPrev)
        P = R + beta*P
        Print *, "||R|| = ", dot_product(R,R)
      !  Call SLEEP(1)
      End Do

  End Subroutine



END MODULE GradientConjugue

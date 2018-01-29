MODULE GradientConjugue

  Use Data
  Use Num
  Use Func
  Use System

  IMPLICIT NONE



CONTAINS

  Subroutine GC(U,SecondMembre, i1, iN)

    !In
      Real, Dimension(i1*Nx+1:iN*Nx), Intent(inout):: U
      Real, Dimension(i1*Nx+1:iN*Nx), Intent(in):: SecondMembre
      Integer, intent(in)::i1, iN
    !Num
      Integer :: k
      Real, Dimension(2) :: XY
      Integer, Dimension(2) :: IJ
    !VarGC
      Real, Dimension(i1*Nx+1:iN*Nx) :: R, RPrev, P, Z
      Real :: alpha, beta

    !Initilialisation
      R = SecondMembre - ProdMat(U, i1, iN)
      Rprev = R
      P = R

      Do while ( dot_product(R,R) > epsilon)
        !Print *, "ProdScal = ", dot_product(R,R)
        Z = ProdMat(P, i1, iN)
        !Print *, "Z = ", Z
        alpha = dot_product(R,R)/dot_product(P,Z)
        U = U + alpha*P
        Rprev = R
        R = R - alpha*Z
        beta = dot_product(R,R)/dot_product(RPrev,RPrev)
        P = R + beta*P
        !Print *, "||R|| = ", dot_product(R,R)
      !  Call SLEEP(1)
      End Do

  End Subroutine



END MODULE GradientConjugue

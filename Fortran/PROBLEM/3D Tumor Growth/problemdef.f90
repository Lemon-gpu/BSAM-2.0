MODULE ProblemDef
USE NodeInfoDef, ONLY: r8
IMPLICIT NONE
!
REAL(KIND=r8), PUBLIC:: alpha, eps, eps2, gamma, mu, shearvel
!
REAL(KIND=r8), PUBLIC::   pp, AA, Zeta, epep,epep2,DtDh,MD
REAL(KIND=r8) :: x0, x1, x2, x3,xx1,xx2,xx3,yy1,yy2,yy3, y0, y1, y2, y3, sigmax, sigmay
REAL(KIND=r8) :: xx11,xx12,yy11,yy12
REAL(KIND=r8) :: xx21,xx22,yy21,yy22
REAL(KIND=r8) :: xx31,xx32,yy31,yy32
REAL(KIND=r8) :: xx41,xx42,yy41,yy42
REAL(KIND=r8), PUBLIC:: DD, FFF, Tau, Rou00, Rou01, Gama, kk1, kk2, dede1, dede2
REAL(KIND=r8), PUBLIC:: kkm1, kkm2, kkp1, kkp2

END MODULE ProblemDef

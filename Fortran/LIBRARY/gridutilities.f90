!==========================================================================
! BSAM 2.0: Block-Structured Adaptive Multigrid Solver
!==========================================================================
!
! WiseSoft: Innovators, Brothers.
!
! (c) Copyright Steven M. Wise, 2006
! Department of Mathematics
! University of California at Irvine
! swise@math.uci.edu
!
! (c) Copyright Steven M. Wise, 2007
! Department of Mathematics
! University of Tennessee
! swise@math.utk.edu
!
! (c) Copyright Steven M. Wise, 2015
! Department of Mathematics
! University of Tennessee
! swise@math.utk.edu
!
! (c) Copyright Zhenlin Guo, 2016
! Department of Mathematics
! University of California, Irvine
! zhenling@math.uci.edu
! -----------------------------------------------------------------------
! This software is made available for research and instructional use only.
! You may copy and use this software without charge for these
! non-commercial purposes, provided that the copyright notice and
! associated text is reproduced on all copies. For all other uses,
! including distribution of modified versions, please contact the authors.
!
! Commercial use is strictly forbidden without permission.
!
! This software is made available "as is" without any assurance that it
! will work for your purposes. The software may in fact have defects,
! so use the software at your own risk.
!
! -----------------------------------------------------------------------
! File:             gridutilities.f90
! Purpose:          BSAM grid utilities module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 2.0 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------
MODULE GridUtilities
!
! Contains the commonly-used 2 and 3D grid utilities.  These only operate on
! uniform, rectangular grid patches.
!
CONTAINS
!
   FUNCTION ULap2D(a) RESULT(ulapresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! 2D UNDIVIDED laplacian of a(0:mx(1)+1,0:mx(2)+1). The result is stored in
! ulapresult(1:mx(1),1:mx(2)).
!
      REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2, &
         1:SIZE(a,2)-2):: ulapresult
!
      INTEGER, DIMENSION(1:2):: mx
!
      mx(1) = SIZE(a,1)-2
      mx(2) = SIZE(a,2)-2
!
      ulapresult(1:mx(1),1:mx(2)) =        a(2:mx(1)+1,1:mx(2)  ) &
         +        a(0:mx(1)-1,1:mx(2)  ) &
         +        a(1:mx(1)  ,2:mx(2)+1) &
         +        a(1:mx(1)  ,0:mx(2)-1) &
         - 4.0_r8*a(1:mx(1)  ,1:mx(2)  )
!
   END FUNCTION ULap2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION UDiv2D(f1,f2) RESULT(udivresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! UNDIVIDED divergence of the 2D flux function
!
!  f = (f1(0:mx(1),1:mx(2)),f2(1:mx(1),0:mx(2))).
!
! The result is stored in udivresult(1:mx(1),1:mx(2)).
!
      REAL(KIND=r8), DIMENSION(0:,1:), INTENT(IN):: f1
      REAL(KIND=r8), DIMENSION(1:,0:), INTENT(IN):: f2
      REAL(KIND=r8), DIMENSION(1:SIZE(f2,1), &
         1:SIZE(f1,2)):: udivresult
!
      INTEGER, DIMENSION(1:2):: mx
!
      mx(1) = SIZE(f2,1)
      mx(2) = SIZE(f1,2)
!
      udivresult(1:mx(1),1:mx(2)) = f1(1:mx(1)  ,1:mx(2)  ) &
         - f1(0:mx(1)-1,1:mx(2)  ) &
         + f2(1:mx(1)  ,1:mx(2)  ) &
         - f2(1:mx(1)  ,0:mx(2)-1)
!
   END FUNCTION UDiv2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Restriction2D(a) RESULT(restrictionresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(:,:,:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(1:SIZE(a,1)/2, &
         1:SIZE(a,2)/2, &
         1:SIZE(a,3)):: restrictionresult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:2):: cmx, mx
!
      mx(1)  = SIZE(a,1)
      mx(2)  = SIZE(a,2)
      cmx(1:2) = mx(1:2)/2
      ncv = SIZE(a,3)
!
      restrictionresult(1:cmx(1)    ,1:cmx(2)    ,1:ncv) &
         =    0.25_r8*(a(2: mx(1)  :2,2: mx(2)  :2,1:ncv) &
         +             a(1: mx(1)-1:2,2: mx(2)  :2,1:ncv) &
         +             a(2: mx(1)  :2,1: mx(2)-1:2,1:ncv) &
         +             a(1: mx(1)-1:2,1: mx(2)-1:2,1:ncv))
!
   END FUNCTION Restriction2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION UDivPhi3D(p,v1,v2,v3) RESULT(udivresult)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(0:,0:,0:), INTENT(IN):: p
      REAL(KIND=r8), DIMENSION(0:,1:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(1:,0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION(1:,1:,0:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION(1:SIZE(p,1)-2,1:SIZE(p,2)-2,1:SIZE(p,3)-2):: udivresult
!
      INTEGER, DIMENSION(1:3):: mx
      REAL(KIND=r8), DIMENSION(0:SIZE(p,1)-2,1:SIZE(p,2)-2,1:SIZE(p,3)-2):: f1
      REAL(KIND=r8), DIMENSION(1:SIZE(p,1)-2,0:SIZE(p,2)-2,1:SIZE(p,3)-2):: f2
      REAL(KIND=r8), DIMENSION(1:SIZE(p,1)-2,1:SIZE(p,2)-2,0:SIZE(p,3)-2):: f3
!
      mx(1) = SIZE(p,1)-2
      mx(2) = SIZE(p,2)-2
      mx(3) = SIZE(p,3)-2
!
! Calculate the UNDIVIDED 3D flux function:
      f1(0:mx(1)  ,1:mx(2)  ,1:mx(3)  )   &
         = 0.5_r8*(p(1:mx(1)+1,1:mx(2)  ,1:mx(3)  )   &
         + p(0:mx(1)  ,1:mx(2)  ,1:mx(3)  ))  &
         *v1(0:mx(1)  ,1:mx(2)  ,1:mx(3)  )
!
      f2(1:mx(1)  ,0:mx(2)  ,1:mx(3)  )   &
         = 0.5_r8*(p(1:mx(1)  ,1:mx(2)+1,1:mx(3)  )   &
         + p(1:mx(1)  ,0:mx(2)  ,1:mx(3)  ))  &
         *v2(1:mx(1)  ,0:mx(2)  ,1:mx(3))
!
      f3(1:mx(1)  ,1:mx(2)  ,0:mx(3)  )   &
         = 0.5_r8*(p(1:mx(1)  ,1:mx(2)  ,1:mx(3)+1)   &
         + p(1:mx(1)  ,1:mx(2)  ,0:mx(3)  ))  &
         *v3(1:mx(1)  ,1:mx(2)  ,0:mx(3))
!
! Calculate the UNDIVIDED divergence of the flux:
      udivresult(1:mx(1)  ,1:mx(2)  ,1:mx(3))     &
         = UDiv3D(f1(0:mx(1)  ,1:mx(2)  ,1:mx(3)),    &
         f2(1:mx(1)  ,0:mx(2)  ,1:mx(3)),    &
         f3(1:mx(1)  ,1:mx(2)  ,0:mx(3)))
!
   END FUNCTION UDivPhi3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Prolongation2D(a) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! This prolongation algorithm assumes no ghost layers are present.
!
      REAL(KIND=r8), DIMENSION(:,:,:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(1:SIZE(a,1)*2, &
         1:SIZE(a,2)*2, &
         1:SIZE(a,3)):: presult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:2):: cmx, mx
!
      cmx(1) = SIZE(a,1)
      cmx(2) = SIZE(a,2)
      mx(1:2) = cmx(1:2)*2
      ncv = SIZE(a,3)
!
      presult(2:mx(1)  :2,2:mx(2)  :2,1:ncv) = a(1:cmx(1),1:cmx(2),1:ncv)
      presult(1:mx(1)-1:2,2:mx(2)  :2,1:ncv) = a(1:cmx(1),1:cmx(2),1:ncv)
      presult(2:mx(1)  :2,1:mx(2)-1:2,1:ncv) = a(1:cmx(1),1:cmx(2),1:ncv)
      presult(1:mx(1)-1:2,1:mx(2)-1:2,1:ncv) = a(1:cmx(1),1:cmx(2),1:ncv)
!
   END FUNCTION Prolongation2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION BiLinProlongationP1(a) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! This bilinear prolongation algorithm assumes exactly one ghost layer, i.e.,
! mbc = 1.
!
      REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
         0:(SIZE(a,2)-2)*2+1, &
         1: SIZE(a,3)       ):: presult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:2):: cmx, mx
      REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
         0: SIZE(a,2)-2   +1, &
         1: SIZE(a,3)       ):: b
!
      cmx(1) = SIZE(a,1)-2
      cmx(2) = SIZE(a,2)-2
      mx(1:2) = cmx(1:2)*2
      ncv = SIZE(a,3)
!
! Linear interpolation in the x-direction first:
      b(0: mx(1)  :2,0:cmx(2)+1  ,1:ncv) &
         =  3.0_r8*a(0:cmx(1)    ,0:cmx(2)+1  ,1:ncv) &
         +         a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:ncv)
      b(1: mx(1)+1:2,0:cmx(2)+1  ,1:ncv) &
         =         a(0:cmx(1)    ,0:cmx(2)+1  ,1:ncv) &
         +  3.0_r8*a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:ncv)
!
! Linear interpolation in the y-direction second:
      presult(0: mx(1)+1  ,0: mx(2)  :2,1:ncv) &
         = (3.0_r8*b(0: mx(1)+1  ,0:cmx(2)    ,1:ncv) &
         +         b(0: mx(1)+1  ,1:cmx(2)+1  ,1:ncv))/16.0_r8
      presult(0: mx(1)+1  ,1: mx(2)+1:2,1:ncv) &
         = (       b(0: mx(1)+1  ,0:cmx(2)    ,1:ncv) &
         +  3.0_r8*b(0: mx(1)+1  ,1:cmx(2)+1  ,1:ncv))/16.0_r8
!
   END FUNCTION BiLinProlongationP1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION BiLinProlongationP2(a) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! This bilinear prolongation algorithm assumes exactly two ghost layers, i.e.,
! mbc = 2.
!
      REAL(KIND=r8), DIMENSION(-1:,-1:,1:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
         -1:(SIZE(a,2)-4)*2+2, &
         1: SIZE(a,3)       ):: presult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:2):: cmx, mx
      REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
         -1: SIZE(a,2)-4   +2, &
         1: SIZE(a,3)       ):: b
!
      cmx(1) = SIZE(a,1)-4
      cmx(2) = SIZE(a,2)-4
      mx(1:2) = cmx(1:2)*2
      ncv = SIZE(a,3)
!
! Linear interpolation in the x-direction first:
      b(-1: mx(1)+1:2,-1:cmx(2)+2  ,1:ncv) &
         =  3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:ncv) &
         +         a(-1:cmx(1)    ,-1:cmx(2)+2  ,1:ncv)
!
      b( 0: mx(1)+2:2,-1:cmx(2)+2  ,1:ncv) &
         =  3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:ncv) &
         +         a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,1:ncv)
!
! Linear interpolation in the y-direction second:
      presult(-1: mx(1)+2  ,-1: mx(2)+1:2,1:ncv) &
         = (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:ncv) &
         +         b(-1: mx(1)+2  ,-1:cmx(2)    ,1:ncv))/16.0_r8
!
      presult(-1: mx(1)+2  , 0: mx(2)+2:2,1:ncv) &
         = (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:ncv) &
         +         b(-1: mx(1)+2  , 1:cmx(2)+2  ,1:ncv))/16.0_r8
!
   END FUNCTION BiLinProlongationP2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION BiLinProlongationP1MC(a) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! This mass-corrected bilinear prolongation algorithm assumes exactly one ghost
! layer, i.e., mbc = 1.
!
      REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
         0:(SIZE(a,2)-2)*2+1, &
         1: SIZE(a,3)       ):: presult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:2):: cmx, mx
      REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
         0: SIZE(a,2)-2   +1, &
         1: SIZE(a,3)       ):: b
      REAL(KIND=r8), DIMENSION(1: SIZE(a,1)-2     , &
         1: SIZE(a,2)-2     , &
         1: SIZE(a,3)       ):: cor
!
      cmx(1) = SIZE(a,1)-2
      cmx(2) = SIZE(a,2)-2
      mx(1:2) = cmx(1:2)*2
      ncv = SIZE(a,3)
!
! Linear interpolation in the x-direction first:
      b(0: mx(1)  :2,0:cmx(2)+1  ,1:ncv) &
         =  3.0_r8*a(0:cmx(1)    ,0:cmx(2)+1  ,1:ncv) &
         +         a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:ncv)
      b(1: mx(1)+1:2,0:cmx(2)+1  ,1:ncv) &
         =         a(0:cmx(1)    ,0:cmx(2)+1  ,1:ncv) &
         +  3.0_r8*a(1:cmx(1)+1  ,0:cmx(2)+1  ,1:ncv)
!
! Linear interpolation in the y-direction second:
      presult(0: mx(1)+1  ,0: mx(2)  :2,1:ncv) &
         = (3.0_r8*b(0: mx(1)+1  ,0:cmx(2)    ,1:ncv) &
         +         b(0: mx(1)+1  ,1:cmx(2)+1  ,1:ncv))/16.0_r8
      presult(0: mx(1)+1  ,1: mx(2)+1:2,1:ncv) &
         = (       b(0: mx(1)+1  ,0:cmx(2)    ,1:ncv) &
         +  3.0_r8*b(0: mx(1)+1  ,1:cmx(2)+1  ,1:ncv))/16.0_r8
!
! The mass correction is computed only for regular cells; not for ghost cells.
! Compute mass correction:
      cor(1:cmx(1),1:cmx(2),1:ncv) &
         = a(1:cmx(1),1:cmx(2),1:ncv) &
         - 0.25_r8*(presult(2:mx(1)  :2,2:mx(2)  :2,1:ncv) &
         +          presult(1:mx(1)-1:2,2:mx(2)  :2,1:ncv) &
         +          presult(2:mx(1)  :2,1:mx(2)-1:2,1:ncv) &
         +          presult(1:mx(1)-1:2,1:mx(2)-1:2,1:ncv))
!
! Add the mass correction:
      presult(2:mx(1)  :2,2:mx(2)  :2,1:ncv) &
         = presult(2:mx(1)  :2,2:mx(2)  :2,1:ncv)+cor(1:cmx(1),1:cmx(2),1:ncv)
      presult(1:mx(1)-1:2,2:mx(2)  :2,1:ncv) &
         = presult(1:mx(1)-1:2,2:mx(2)  :2,1:ncv)+cor(1:cmx(1),1:cmx(2),1:ncv)
      presult(2:mx(1)  :2,1:mx(2)-1:2,1:ncv) &
         = presult(2:mx(1)  :2,1:mx(2)-1:2,1:ncv)+cor(1:cmx(1),1:cmx(2),1:ncv)
      presult(1:mx(1)-1:2,1:mx(2)-1:2,1:ncv) &
         = presult(1:mx(1)-1:2,1:mx(2)-1:2,1:ncv)+cor(1:cmx(1),1:cmx(2),1:ncv)
!
   END FUNCTION BiLinProlongationP1MC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION BiLinProlongationP2MC(a) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! This mass-corrected bilinear prolongation algorithm assumes exactly two ghost
! layers, i.e., mbc = 2.
!
      REAL(KIND=r8), DIMENSION(-1:,-1:,1:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
         -1:(SIZE(a,2)-4)*2+2, &
         1: SIZE(a,3)       ):: presult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:2):: cmx, mx
      REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
         -1: SIZE(a,2)-4   +2, &
         1: SIZE(a,3)       ):: b
      REAL(KIND=r8), DIMENSION( 1: SIZE(a,1)-4     , &
         1: SIZE(a,2)-4     , &
         1: SIZE(a,3)       ):: cor
!
      cmx(1) = SIZE(a,1)-4
      cmx(2) = SIZE(a,2)-4
      mx(1:2) = cmx(1:2)*2
      ncv = SIZE(a,3)
!
! Linear interpolation in the x-direction first:
      b(-1: mx(1)+1:2,-1:cmx(2)+2  ,1:ncv) &
         =  3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:ncv) &
         +         a(-1:cmx(1)    ,-1:cmx(2)+2  ,1:ncv)
!
      b( 0: mx(1)+2:2,-1:cmx(2)+2  ,1:ncv) &
         =  3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,1:ncv) &
         +         a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,1:ncv)
!
! Linear interpolation in the y-direction second:
      presult(-1: mx(1)+2  ,-1: mx(2)+1:2,1:ncv) &
         = (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:ncv) &
         +         b(-1: mx(1)+2  ,-1:cmx(2)    ,1:ncv))/16.0_r8
!
      presult(-1: mx(1)+2  , 0: mx(2)+2:2,1:ncv) &
         = (3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,1:ncv) &
         +         b(-1: mx(1)+2  , 1:cmx(2)+2  ,1:ncv))/16.0_r8
!
! The mass correction is computed only for regular cells; not for ghost cells.
! Compute mass correction:
      cor(1:cmx(1),1:cmx(2),1:ncv) &
         = a(1:cmx(1),1:cmx(2),1:ncv) &
         - 0.25_r8*(presult(2:mx(1)  :2,2:mx(2)  :2,1:ncv) &
         +          presult(1:mx(1)-1:2,2:mx(2)  :2,1:ncv) &
         +          presult(2:mx(1)  :2,1:mx(2)-1:2,1:ncv) &
         +          presult(1:mx(1)-1:2,1:mx(2)-1:2,1:ncv))
!
! Add the mass correction:
      presult(2:mx(1)  :2,2:mx(2)  :2,1:ncv) &
         = presult(2:mx(1)  :2,2:mx(2)  :2,1:ncv)+cor(1:cmx(1),1:cmx(2),1:ncv)
      presult(1:mx(1)-1:2,2:mx(2)  :2,1:ncv) &
         = presult(1:mx(1)-1:2,2:mx(2)  :2,1:ncv)+cor(1:cmx(1),1:cmx(2),1:ncv)
      presult(2:mx(1)  :2,1:mx(2)-1:2,1:ncv) &
         = presult(2:mx(1)  :2,1:mx(2)-1:2,1:ncv)+cor(1:cmx(1),1:cmx(2),1:ncv)
      presult(1:mx(1)-1:2,1:mx(2)-1:2,1:ncv) &
         = presult(1:mx(1)-1:2,1:mx(2)-1:2,1:ncv)+cor(1:cmx(1),1:cmx(2),1:ncv)
!
   END FUNCTION BiLinProlongationP2MC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 2D Face/Edge Variable Grid Utilities:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Restriction2DV1(w1) RESULT(restrictionresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Restriction of v1 (east-west) edge functions by a simple average. This
! method is flux preserving:
!
! Input:                 w1(0: mx(1),1: mx(2),1:)
!
! Output: restrictionresult(0:cmx(1),1:cmx(2),1:)
!
      REAL(KIND=r8), DIMENSION(0:,1:,1:), INTENT(IN):: w1
      REAL(KIND=r8), DIMENSION(0:(SIZE(w1,1)-1)/2, &
         1: SIZE(w1,2)   /2, &
         1: SIZE(w1,3)      ):: restrictionresult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:2):: cmx, mx
!
      mx(1)  = SIZE(w1,1)-1
      mx(2)  = SIZE(w1,2)
      cmx(1:2) = mx(1:2)/2
      nfv = SIZE(w1,3)
!
      restrictionresult(0:cmx(1)  ,1:cmx(2)    ,1:nfv) &
         = 0.5_r8*(w1(0: mx(1):2,1: mx(2)-1:2,1:nfv) &
         +         w1(0: mx(1):2,2: mx(2)  :2,1:nfv))
!
   END FUNCTION Restriction2DV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Restriction2DV2(w2) RESULT(restrictionresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Restriction of v2 (north-south) edge functions by a simple average. This
! method is flux preserving:
!
! Input:                 w2(1: mx(1),0: mx(2),1:)
!
! Output: restrictionresult(1:cmx(1),0:cmx(2),1:)
!
      REAL(KIND=r8), DIMENSION(1:,0:,1:), INTENT(IN):: w2
      REAL(KIND=r8), DIMENSION(1: SIZE(w2,1)   /2, &
         0:(SIZE(w2,2)-1)/2, &
         1: SIZE(w2,3)      ):: restrictionresult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:2):: cmx, mx
!
      mx(1)  = SIZE(w2,1)
      mx(2)  = SIZE(w2,2)-1
      cmx(1:2) = mx(1:2)/2
      nfv = SIZE(w2,3)
!
      restrictionresult(1:cmx(1)    ,0:cmx(2)  ,1:nfv) &
         = 0.5_r8*(w2(1: mx(1)-1:2,0: mx(2):2,1:nfv) &
         +         w2(2: mx(1)  :2,0: mx(2):2,1:nfv))
!
   END FUNCTION Restriction2DV2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Prolongation2DV1Ex(v1) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Prolongation of v1 (east-west) edge functions by bi-linear interpolation.
!
! Input:       v1(-1:cmx(1)+1,0:cmx(2)+1,1:) (with x1 ghost cells)
!
! Output: presult(-1: mx(1)+1,1: mx(2)  ,1:)
!
      REAL(KIND=r8), DIMENSION(-1:,0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:(SIZE(v1,1)-3)*2+1, &
         0:(SIZE(v1,2)-2)*2+1, &
         1: SIZE(v1,3)        ):: presult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:2):: cmx, mx
      REAL(KIND=r8), DIMENSION(-2:(SIZE(v1,1)-3)*2+2, &
         0:(SIZE(v1,2)-2)*2+1, &
         1: SIZE(v1,3)        ):: tmp
!
      cmx(1) = SIZE(v1,1)-3
      cmx(2) = SIZE(v1,2)-2
      mx(1:2) = cmx(1:2)*2
      nfv = SIZE(v1,3)
!
! Not optimized. Can save some memory here by eliminating tmp.
!
! A points:
      tmp(-2: mx(1)+2:2,0: mx(2)  :2,1:nfv) &
         = 0.25_r8*(3.0_r8*v1(-1:cmx(1)+1  ,0:cmx(2)    ,1:nfv) &
         +                 v1(-1:cmx(1)+1  ,1:cmx(2)+1  ,1:nfv))
!
! B points:
      tmp(-2: mx(1)+2:2,1: mx(2)+1:2,1:nfv) &
         = 0.25_r8*(       v1(-1:cmx(1)+1  ,0:cmx(2)    ,1:nfv) &
         +          3.0_r8*v1(-1:cmx(1)+1  ,1:cmx(2)+1  ,1:nfv))
!
! C points:
      tmp(-1: mx(1)+1:2,0: mx(2)  :2,1:nfv) &
         = 0.50_r8*(      tmp(-2: mx(1)  :2,0: mx(2)  :2,1:nfv) &
         +                tmp( 0: mx(1)+2:2,0: mx(2)  :2,1:nfv))
!
! D points:
      tmp(-1: mx(1)+1:2,1: mx(2)+1:2,1:nfv) &
         = 0.50_r8*(      tmp(-2: mx(1)  :2,1: mx(2)+1:2,1:nfv) &
         +                tmp( 0: mx(1)+2:2,1: mx(2)+1:2,1:nfv))
!
      presult(-1:mx(1)+1,0:mx(2)+1,1:nfv) = tmp(-1:mx(1)+1,0:mx(2)+1,1:nfv)
!
   END FUNCTION Prolongation2DV1Ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Prolongation2DV1(v1) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Prolongation of v1 (east-west) edge functions by bi-linear interpolation.
!
! Input:       v1(0:cmx(1),0:cmx(2)+1,1:) (no x1 ghost cells)
!
! Output: presult(0: mx(1),1: mx(2)  ,1:)
!
      REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(0:(SIZE(v1,1)-1)*2, &
         1:(SIZE(v1,2)-2)*2, &
         1: SIZE(v1,3)      ):: presult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:2):: cmx, mx
!
      cmx(1) = SIZE(v1,1)-1
      cmx(2) = SIZE(v1,2)-2
      mx(1:2) = cmx(1:2)*2
      nfv = SIZE(v1,3)
!
! A points:
      presult(0: mx(1)  :2,2: mx(2)  :2,1:nfv) &
         = 0.25_r8*(3.0_r8*v1(0:cmx(1)    ,1:cmx(2)    ,1:nfv) &
         +                 v1(0:cmx(1)    ,2:cmx(2)+1  ,1:nfv))
!
! B points:
      presult(0: mx(1)  :2,1: mx(2)-1:2,1:nfv) &
         = 0.25_r8*(       v1(0:cmx(1)    ,0:cmx(2)-1  ,1:nfv) &
         +          3.0_r8*v1(0:cmx(1)    ,1:cmx(2)    ,1:nfv))
!
! C points:
      presult(1:mx(1)-1:2,2: mx(2)  :2,1:nfv) &
         = 0.50_r8*(   presult(0:mx(1)-2:2,2: mx(2)  :2,1:nfv) &
         +             presult(2:mx(1)  :2,2: mx(2)  :2,1:nfv))
!
! D points:
      presult(1:mx(1)-1:2,1: mx(2)-1:2,1:nfv) &
         = 0.50_r8*(   presult(0:mx(1)-2:2,1: mx(2)-1:2,1:nfv) &
         +             presult(2:mx(1)  :2,1: mx(2)-1:2,1:nfv))
!
   END FUNCTION Prolongation2DV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Prolongation2DV2Ex(v2) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Prolongation of v2 (north-south) edge functions by bi-linear interpolation.
!
! Input:       v2(0:cmx(1)+1,-1:cmx(2)+1,1:) (with x2 ghost cells)
!
! Output: presult(0: mx(1)+1,-1: mx(2)+1,1:)
!
      REAL(KIND=r8), DIMENSION(0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:(SIZE(v2,1)-2)*2+1, &
         -1:(SIZE(v2,2)-3)*2+1, &
         1: SIZE(v2,3)        ):: presult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:2):: cmx, mx
      REAL(KIND=r8), DIMENSION( 0:(SIZE(v2,1)-2)*2+1, &
         -2:(SIZE(v2,2)-3)*2+2, &
         1: SIZE(v2,3)        ):: tmp
!
      cmx(1) = SIZE(v2,1)-2
      cmx(2) = SIZE(v2,2)-3
      mx(1:2) = cmx(1:2)*2
      nfv = SIZE(v2,3)
!
! Not optimized. Can save some memory here by eliminating tmp.
!
! A points:
      tmp(0: mx(1)  :2,-2: mx(2)+2:2,1:nfv) &
         = 0.25_r8*(3.0_r8*v2(0:cmx(1)    ,-1:cmx(2)+1  ,1:nfv) &
         +                 v2(1:cmx(1)+1  ,-1:cmx(2)+1  ,1:nfv))
!
! B points:
      tmp(1: mx(1)+1:2,-2: mx(2)+2:2,1:nfv) &
         = 0.25_r8*(       v2(0:cmx(1)    ,-1:cmx(2)+1  ,1:nfv) &
         +          3.0_r8*v2(1:cmx(1)+1  ,-1:cmx(2)+1  ,1:nfv))
!
! C points:
      tmp(0: mx(1)  :2,-1: mx(2)+1:2,1:nfv) &
         = 0.50_r8*(      tmp(0: mx(1)  :2,-2: mx(2)  :2,1:nfv) &
         +                tmp(0: mx(1)  :2, 0: mx(2)+2:2,1:nfv))
!
! D points:
      tmp(1: mx(1)+1:2,-1: mx(2)+1:2,1:nfv) &
         = 0.50_r8*(      tmp(1: mx(1)+1:2,-2: mx(2)  :2,1:nfv) &
         +                tmp(1: mx(1)+1:2, 0: mx(2)+2:2,1:nfv))
!
      presult(0:mx(1)+1,-1:mx(2)+1,1:nfv) = tmp(0:mx(1)+1,-1:mx(2)+1,1:nfv)
!
   END FUNCTION Prolongation2DV2Ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Prolongation2DV2(v2) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Prolongation of v2 (north-south) edge functions by bi-linear interpolation.
!
! Input:       v2(0:cmx(1)+1,0:cmx(2),1:) (no x2 ghost cells)
!
! Output: presult(1: mx(1)  ,0: mx(2),1:)
!
      REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION(1:(SIZE(v2,1)-2)*2, &
         0:(SIZE(v2,2)-1)*2, &
         1: SIZE(v2,3)      ):: presult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:2):: cmx, mx
!
      cmx(1) = SIZE(v2,1)-2
      cmx(2) = SIZE(v2,2)-1
      mx(1:2) = cmx(1:2)*2
      nfv = SIZE(v2,3)
!
! A points:
      presult(2: mx(1)  :2,0: mx(2)  :2,1:nfv) &
         = 0.25_r8*(3.0_r8*v2(1:cmx(1)    ,0:cmx(2)    ,1:nfv) &
         +                 v2(2:cmx(1)+1  ,0:cmx(2)    ,1:nfv))
!
! B points:
      presult(1: mx(1)-1:2,0: mx(2)  :2,1:nfv) &
         = 0.25_r8*(       v2(0:cmx(1)-1  ,0:cmx(2)    ,1:nfv) &
         +          3.0_r8*v2(1:cmx(1)    ,0:cmx(2)    ,1:nfv))
!
! C points:
      presult(2: mx(1)  :2,1: mx(2)-1:2,1:nfv) &
         = 0.50_r8*(  presult(2: mx(1)  :2,0: mx(2)-2:2,1:nfv) &
         +            presult(2: mx(1)  :2,2: mx(2)  :2,1:nfv))
!
! D points:
      presult(1: mx(1)-1:2,1: mx(2)-1:2,1:nfv) &
         = 0.50_r8*(  presult(1: mx(1)-1:2,0: mx(2)-2:2,1:nfv) &
         +            presult(1: mx(1)-1:2,2: mx(2)  :2,1:nfv))
!
   END FUNCTION Prolongation2DV2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION ULap2DV1(v1) RESULT(ulapresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! 2D undivided laplacian of v1 (east-west) edge functions.
!
! Input:          v1(-1:mx(1)+1,0:mx(2)+1)
!
! Output: ulapresult( 0:mx(1)  ,1:mx(2)  )
!
! Intermediate
! Fluxes:         f1( 0:mx(1)+1,1:mx(2)  ) (  cell centered)
!                 f2( 0:mx(1)  ,0:mx(2)  ) (vertex centered)
!
      REAL(KIND=r8), DIMENSION(-1:,0:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(0:SIZE(v1,1)-3+1,1:SIZE(v1,2)-2):: f1
      REAL(KIND=r8), DIMENSION(0:SIZE(v1,1)-3  ,0:SIZE(v1,2)-2):: f2
      REAL(KIND=r8), DIMENSION(0:SIZE(v1,1)-3  ,1:SIZE(v1,2)-2):: ulapresult
!
      INTEGER, DIMENSION(1:2):: mx
!
      mx(1) = SIZE(v1,1)-3
      mx(2) = SIZE(v1,2)-2
!
      f1(0:mx(1)+1,1:mx(2)) = v1( 0:mx(1)+1,1:mx(2)  ) &
         - v1(-1:mx(1)  ,1:mx(2)  )
      f2(0:mx(1)  ,0:mx(2)) = v1( 0:mx(1)  ,1:mx(2)+1) &
         - v1( 0:mx(1)  ,0:mx(2)  )
!
      ulapresult(0:mx(1),1:mx(2)) = f1(1:mx(1)+1,1:mx(2)  ) &
         - f1(0:mx(1)  ,1:mx(2)  ) &
         + f2(0:mx(1)  ,1:mx(2)  ) &
         - f2(0:mx(1)  ,0:mx(2)-1)
!
   END FUNCTION ULap2DV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION ULap2DV2(v2) RESULT(ulapresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! 2D undivided laplacian of v2 (north-south) edge functions.
!
! Input:          v1(0:mx(1)+1,-1:mx(2)+1)
!
! Output: ulapresult(1:mx(1)  , 0:mx(2)  )
!
! Intermediate
! Fluxes:        f1( 0:mx(1)  , 0:mx(2)  ) (vertex centered)
!                f2( 1:mx(1)  , 0:mx(2)+1) (cell   centered)
!
      REAL(KIND=r8), DIMENSION(0:,-1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION(0:SIZE(v2,1)-2,0:SIZE(v2,2)-3  ):: f1
      REAL(KIND=r8), DIMENSION(1:SIZE(v2,1)-2,0:SIZE(v2,2)-3+1):: f2
      REAL(KIND=r8), DIMENSION(1:SIZE(v2,1)-2,0:SIZE(v2,2)-3  ):: ulapresult
!
      INTEGER, DIMENSION(1:2):: mx
!
      mx(1) = SIZE(v2,1)-2
      mx(2) = SIZE(v2,2)-3
!
      f1(0:mx(1),0:mx(2)  ) = v2(1:mx(1)+1, 0:mx(2)  ) &
         - v2(0:mx(1)  , 0:mx(2)  )
      f2(1:mx(1),0:mx(2)+1) = v2(1:mx(1)  , 0:mx(2)+1) &
         - v2(1:mx(1)  ,-1:mx(2)  )
!
      ulapresult(1:mx(1),0:mx(2)) = f1(1:mx(1)  ,0:mx(2)  ) &
         - f1(0:mx(1)-1,0:mx(2)  ) &
         + f2(1:mx(1)  ,1:mx(2)+1) &
         - f2(1:mx(1)  ,0:mx(2)  )
!
   END FUNCTION ULap2DV2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 3D CC Grid Utilities:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION ULap3D(a) RESULT(ulapresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! 3D UNDIVIDED laplacian of a(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1). The result is
! stored in ulapresult(1:mx(1),1:mx(2)).
!
      REAL(KIND=r8), DIMENSION(0:,0:,0:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2, &
         1:SIZE(a,2)-2, &
         1:SIZE(a,3)-2):: ulapresult
!
      INTEGER, DIMENSION(1:3):: mx
!
      mx(1) = SIZE(a,1)-2
      mx(2) = SIZE(a,2)-2
      mx(3) = SIZE(a,3)-2
!
      ulapresult(1:mx(1),1:mx(2),1:mx(3)) &
         =        a(2:mx(1)+1,1:mx(2)  ,1:mx(3)  ) &
         +        a(0:mx(1)-1,1:mx(2)  ,1:mx(3)  ) &
         +        a(1:mx(1)  ,2:mx(2)+1,1:mx(3)  ) &
         +        a(1:mx(1)  ,0:mx(2)-1,1:mx(3)  ) &
         +        a(1:mx(1)  ,1:mx(2)  ,2:mx(3)+1) &
         +        a(1:mx(1)  ,1:mx(2)  ,0:mx(3)-1) &
         - 6.0_r8*a(1:mx(1)  ,1:mx(2)  ,1:mx(3)  )
!
   END FUNCTION ULap3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION ULap3DV1(v1) RESULT(ulapresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! 3D undivided laplacian of v1 (east-west) edge functions.
!
! Input:          v1(-1:mx(1)+1,0:mx(2)+1,0:mx(3)+1)
!
! Output: ulapresult( 0:mx(1)  ,1:mx(2)  ,1:mx(3)  )
!
! Intermediate
!
      REAL(KIND=r8), DIMENSION(-1:,0:,0:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(0:SIZE(v1,1)-3  ,1:SIZE(v1,2)-2  ,1:SIZE(v1,3)-2):: ulapresult
!
      INTEGER, DIMENSION(1:3):: mx
!
      mx(1) = SIZE(v1,1)-3
      mx(2) = SIZE(v1,2)-2
      mx(3) = SIZE(v1,3)-2
!
      ulapresult(0:mx(1),1:mx(2),1:mx(3)) =        v1( 1:mx(1)+1,1:mx(2)  ,1:mx(3)  ) &
         +        v1(-1:mx(1)-1,1:mx(2)  ,1:mx(3)  ) &
         +        v1( 0:mx(1)  ,2:mx(2)+1,1:mx(3)  ) &
         +        v1( 0:mx(1)  ,0:mx(2)-1,1:mx(3)  ) &
         +        v1( 0:mx(1)  ,1:mx(2)  ,2:mx(3)+1) &
         +        v1( 0:mx(1)  ,1:mx(2)  ,0:mx(3)-1) &
         - 6.0_r8*v1( 0:mx(1)  ,1:mx(2)  ,1:mx(3)  )
!
   END FUNCTION ULap3DV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION ULap3DV2(v2) RESULT(ulapresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! 3D undivided laplacian of v2 (east-west) edge functions.
!
! Input:          v2( 0:mx(1)+1,-1:mx(2)+1,0:mx(3)+1)
!
! Output: ulapresult( 1:mx(1)  , 0:mx(2)  ,1:mx(3)  )
!
! Intermediate
!
      REAL(KIND=r8), DIMENSION(0:,-1:,0:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION(1:SIZE(v2,1)-2  ,0:SIZE(v2,2)-3  ,1:SIZE(v2,3)-2):: ulapresult
!
      INTEGER, DIMENSION(1:3):: mx
!
      mx(1) = SIZE(v2,1)-2
      mx(2) = SIZE(v2,2)-3
      mx(3) = SIZE(v2,3)-2
!
      ulapresult(1:mx(1),0:mx(2),1:mx(3)) =        v2( 2:mx(1)+1, 0:mx(2)  ,1:mx(3)  ) &
         +        v2( 0:mx(1)-1, 0:mx(2)  ,1:mx(3)  ) &
         +        v2( 1:mx(1)  , 1:mx(2)+1,1:mx(3)  ) &
         +        v2( 1:mx(1)  ,-1:mx(2)-1,1:mx(3)  ) &
         +        v2( 1:mx(1)  , 0:mx(2)  ,2:mx(3)+1) &
         +        v2( 1:mx(1)  , 0:mx(2)  ,0:mx(3)-1) &
         - 6.0_r8*v2( 1:mx(1)  , 0:mx(2)  ,1:mx(3)  )
!
   END FUNCTION ULap3DV2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION ULap3DV3(v3) RESULT(ulapresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! 3D undivided laplacian of v3 (east-west) edge functions.
!
! Input:          v3( 0:mx(1)+1, 0:mx(2)+1,-1:mx(3)+1)
!
! Output: ulapresult( 1:mx(1)  , 1:mx(2)  , 0:mx(3)  )
!
! Intermediate
!
      REAL(KIND=r8), DIMENSION(0:,0:,-1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION(1:SIZE(v3,1)-2  ,1:SIZE(v3,2)-2  ,0:SIZE(v3,3)-3):: ulapresult
!
      INTEGER, DIMENSION(1:3):: mx
!
      mx(1) = SIZE(v3,1)-2
      mx(2) = SIZE(v3,2)-2
      mx(3) = SIZE(v3,3)-3
!
      ulapresult(1:mx(1),1:mx(2),0:mx(3)) =        v3( 2:mx(1)+1, 1:mx(2)  , 0:mx(3)  ) &
         +        v3( 0:mx(1)-1, 1:mx(2)  , 0:mx(3)  ) &
         +        v3( 1:mx(1)  , 2:mx(2)+1, 0:mx(3)  ) &
         +        v3( 1:mx(1)  , 0:mx(2)-1, 0:mx(3)  ) &
         +        v3( 1:mx(1)  , 1:mx(2)  , 1:mx(3)+1) &
         +        v3( 1:mx(1)  , 1:mx(2)  ,-1:mx(3)-1) &
         - 6.0_r8*v3( 1:mx(1)  , 1:mx(2)  , 0:mx(3)  )
!
   END FUNCTION ULap3DV3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION UDiv3D(f1,f2,f3) RESULT(udivresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! UNDIVIDED divergence of the 3D flux function
!
!  (f1(0:mx(1),1:mx(2),1:mx(3)),
!   f2(1:mx(1),0:mx(2),1:mx(3)),
!   f3(1:mx(1),1:mx(2),0:mx(3))).
!
! The result is stored in udivresult(1:mx(1),1:mx(2),1:mx(3)).
!
      REAL(KIND=r8), DIMENSION(0:,1:,1:), INTENT(IN):: f1
      REAL(KIND=r8), DIMENSION(1:,0:,1:), INTENT(IN):: f2
      REAL(KIND=r8), DIMENSION(1:,1:,0:), INTENT(IN):: f3
      REAL(KIND=r8), DIMENSION(1:SIZE(f2,1), &
         1:SIZE(f3,2), &
         1:SIZE(f1,3)):: udivresult
!
      INTEGER, DIMENSION(1:3):: mx
!
      mx(1) = SIZE(f2,1)
      mx(2) = SIZE(f3,2)
      mx(3) = SIZE(f1,3)
!
      udivresult(1:mx(1),1:mx(2),:mx(3)) = f1(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ) &
         - f1(0:mx(1)-1,1:mx(2)  ,1:mx(3)  ) &
         + f2(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ) &
         - f2(1:mx(1)  ,0:mx(2)-1,1:mx(3)  ) &
         + f3(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ) &
         - f3(1:mx(1)  ,1:mx(2)  ,0:mx(3)-1)
!
   END FUNCTION UDiv3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Restriction3D(a) RESULT(restrictionresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(1:SIZE(a,1)/2, &
         1:SIZE(a,2)/2, &
         1:SIZE(a,3)/2, &
         1:SIZE(a,4)   ):: restrictionresult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:3):: cmx, mx
!
      mx(1)  = SIZE(a,1)
      mx(2)  = SIZE(a,2)
      mx(3)  = SIZE(a,3)
      cmx(1:3) = mx(1:3)/2
      ncv = SIZE(a,4)
!
      restrictionresult(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv) &
         = 0.125_r8*(a(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         +           a(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         +           a(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         +           a(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         +           a(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         +           a(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         +           a(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         +           a(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv))
!
   END FUNCTION Restriction3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Restriction3DV1(w1) RESULT(restrictionresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Restriction of v1 (east-west) edge functions by a simple average. This
! method is flux preserving:
!
! Input:                 w1(0: mx(1),1: mx(2),1:)
!
! Output: restrictionresult(0:cmx(1),1:cmx(2),1:)
!
      REAL(KIND=r8), DIMENSION(0:,1:,1:,1:), INTENT(IN):: w1
      REAL(KIND=r8), DIMENSION(0:(SIZE(w1,1)-1)/2, &
         1: SIZE(w1,2)   /2, &
         1: SIZE(w1,3)   /2, &
         1: SIZE(w1,4)          ):: restrictionresult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:3):: cmx, mx
!
      mx(1)  = SIZE(w1,1)-1
      mx(2)  = SIZE(w1,2)
      mx(3)  = SIZE(w1,3)
      cmx(1:3) = mx(1:3)/2
      nfv = SIZE(w1,4)
!
      restrictionresult(0:cmx(1)  ,1:cmx(2)    ,1:cmx(3)    ,1:nfv) &
         = 0.25_r8*(w1(0: mx(1):2,1: mx(2)-1:2,1: mx(3)-1:2,1:nfv) &
         +          w1(0: mx(1):2,2: mx(2)  :2,1: mx(3)-1:2,1:nfv) &
         +          w1(0: mx(1):2,1: mx(2)-1:2,2: mx(3)  :2,1:nfv) &
         +          w1(0: mx(1):2,2: mx(2)  :2,2: mx(3)  :2,1:nfv))
!
   END FUNCTION Restriction3DV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Restriction3DV2(w2) RESULT(restrictionresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Restriction of v2 (north-south) edge functions by a simple average. This
! method is flux preserving:
!
! Input:                 w2(1: mx(1),0: mx(2),1:)
!
! Output: restrictionresult(1:cmx(1),0:cmx(2),1:)
!
      REAL(KIND=r8), DIMENSION(1:,0:,1:,1:), INTENT(IN):: w2
      REAL(KIND=r8), DIMENSION(1: SIZE(w2,1)   /2, &
         0:(SIZE(w2,2)-1)/2, &
         1: SIZE(w2,3)   /2, &
         1: SIZE(w2,4)          ):: restrictionresult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:3):: cmx, mx
!
      mx(1)  = SIZE(w2,1)
      mx(2)  = SIZE(w2,2)-1
      mx(3)  = SIZE(w2,3)
      cmx(1:3) = mx(1:3)/2
      nfv = SIZE(w2,4)
!
      restrictionresult(1:cmx(1)    ,0:cmx(2)  ,1:cmx(3)    ,1:nfv) &
         = 0.25_r8*(w2(1: mx(1)-1:2,0: mx(2):2,1: mx(3)-1:2,1:nfv) &
         +          w2(2: mx(1)  :2,0: mx(2):2,1: mx(3)-1:2,1:nfv) &
         +          w2(1: mx(1)-1:2,0: mx(2):2,2: mx(3)  :2,1:nfv) &
         +          w2(2: mx(1)  :2,0: mx(2):2,2: mx(3)  :2,1:nfv))
!
   END FUNCTION Restriction3DV2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Restriction3DV3(w3) RESULT(restrictionresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Restriction of v2 (north-south) edge functions by a simple average. This
! method is flux preserving:
!
! Input:                 w3(1: mx(1),0: mx(2),1:)
!
! Output: restrictionresult(1:cmx(1),0:cmx(2),1:)
!
      REAL(KIND=r8), DIMENSION(1:,1:,0:,1:), INTENT(IN):: w3
      REAL(KIND=r8), DIMENSION(1: SIZE(w3,1)   /2, &
         1: SIZE(w3,2)   /2, &
         0:(SIZE(w3,3)-1)/2, &
         1: SIZE(w3,4)          ):: restrictionresult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:3):: cmx, mx
!
      mx(1)  = SIZE(w3,1)
      mx(2)  = SIZE(w3,2)
      mx(3)  = SIZE(w3,3)-1
      cmx(1:3) = mx(1:3)/2
      nfv = SIZE(w3,4)
!
      restrictionresult(1:cmx(1)    ,1:cmx(2)    ,0:cmx(3)  ,1:nfv) &
         = 0.25_r8*(w3(1: mx(1)-1:2,1: mx(2)-1:2,0: mx(3):2,1:nfv) &
         +          w3(2: mx(1)  :2,1: mx(2)-1:2,0: mx(3):2,1:nfv) &
         +          w3(1: mx(1)-1:2,2: mx(2)  :2,0: mx(3):2,1:nfv) &
         +          w3(2: mx(1)  :2,2: mx(2)  :2,0: mx(3):2,1:nfv))
!
   END FUNCTION Restriction3DV3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Prolongation3D(a) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! This prolongation algorithm assumes no ghost layers are present.
!
      REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(1:SIZE(a,1)*2, &
         1:SIZE(a,2)*2, &
         1:SIZE(a,3)*2, &
         1:SIZE(a,4)   ):: presult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:3):: cmx, mx
!
      cmx(1) = SIZE(a,1)
      cmx(2) = SIZE(a,2)
      cmx(3) = SIZE(a,3)
      mx(1:3) = cmx(1:3)*2
      ncv = SIZE(a,4)
!
      presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         =   a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
!
   END FUNCTION Prolongation3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Prolongation3DV1Ex(v1) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Prolongation of v1 (east-west) edge functions by bi-linear interpolation.
!
! Input:       v1(-1:cmx(1)+1,0:cmx(2)+1,1:) (with x1 ghost cells)
!
! Output: presult(-1: mx(1)+1,1: mx(2)  ,1:)
!
      REAL(KIND=r8), DIMENSION(-1:,0:,0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:(SIZE(v1,1)-3)*2+1, &
         0:(SIZE(v1,2)-2)*2+1, &
         0:(SIZE(v1,3)-2)*2+1, &
         1: SIZE(v1,4)        ):: presult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:3):: cmx, mx
      REAL(KIND=r8), DIMENSION(-2:(SIZE(v1,1)-3)*2+2, &
         0:(SIZE(v1,2)-2)*2+1, &
         0:(SIZE(v1,3)-2)*2+1, &
         1: SIZE(v1,4)        ):: tmp
!
      cmx(1) = SIZE(v1,1)-3
      cmx(2) = SIZE(v1,2)-2
      cmx(3) = SIZE(v1,3)-2
      mx(1:3) = cmx(1:3)*2
      nfv = SIZE(v1,4)
!
! Not optimized. Can save some memory here by eliminating tmp.
!
! A points:
      tmp(-2: mx(1)+2:2, 0: mx(2)  :2, 0: mx(3)  :2,1:nfv) &
         = 0.25_r8*(3.0_r8*v1(-1:cmx(1)+1  , 0:cmx(2)    , 0:cmx(3)    ,1:nfv) &
         +                 v1(-1:cmx(1)+1  , 1:cmx(2)+1  , 0:cmx(3)    ,1:nfv))
!
! B points:
      tmp(-2: mx(1)+2:2, 0: mx(2)  :2, 1: mx(3)+1:2,1:nfv) &
         = 0.25_r8*(3.0_r8*v1(-1:cmx(1)+1  , 0:cmx(2)    , 1:cmx(3)+1  ,1:nfv) &
         +                 v1(-1:cmx(1)+1  , 1:cmx(2)+1  , 1:cmx(3)+1  ,1:nfv))
!
! C points:
      tmp(-2: mx(1)+2:2, 1: mx(2)+1:2, 0: mx(3)  :2,1:nfv) &
         = 0.25_r8*(       v1(-1:cmx(1)+1  , 0:cmx(2)    , 0:cmx(3)    ,1:nfv) &
         +          3.0_r8*v1(-1:cmx(1)+1  , 1:cmx(2)+1  , 0:cmx(3)    ,1:nfv))
!
! D points:
      tmp(-2: mx(1)+2:2, 1: mx(2)+1:2, 1: mx(3)+1:2,1:nfv) &
         = 0.25_r8*(       v1(-1:cmx(1)+1  , 0:cmx(2)    , 1:cmx(3)+1  ,1:nfv) &
         +          3.0_r8*v1(-1:cmx(1)+1  , 1:cmx(2)+1  , 1:cmx(3)+1  ,1:nfv))
!
! A1 points:
      tmp(-2: mx(1)+2:2, 0: mx(2)  :2, 0: mx(3)  :2,1:nfv) &
         = 0.25_r8*(       tmp(-2: mx(1)+2:2, 0: mx(2)  :2, 1: mx(3)+1:2,1:nfv) &
         +          3.0_r8*tmp(-2: mx(1)+2:2, 0: mx(2)  :2, 0: mx(3)  :2,1:nfv))
!
! B1 points:
      tmp(-2: mx(1)+2:2, 0: mx(2)  :2, 1: mx(3)+1:2,1:nfv) &
         = 0.25_r8*(3.0_r8*tmp(-2: mx(1)+2:2, 0: mx(2)  :2, 1: mx(3)+1:2,1:nfv) &
         +                 tmp(-2: mx(1)+2:2, 0: mx(2)  :2, 0: mx(3)  :2,1:nfv))
!

! C1 points:
      tmp(-2: mx(1)+2:2, 1: mx(2)+1:2, 0: mx(3)  :2,1:nfv) &
         = 0.25_r8*(       tmp(-2: mx(1)+2:2, 1: mx(2)+1:2, 1: mx(3)+1:2,1:nfv) &
         +          3.0_r8*tmp(-2: mx(1)+2:2, 1: mx(2)+1:2, 0: mx(3)  :2,1:nfv))
!
! D1 points:
      tmp(-2: mx(1)+2:2, 1: mx(2)+1:2, 1: mx(3)+1:2,1:nfv) &
         = 0.25_r8*(3.0_r8*tmp(-2: mx(1)+2:2, 1: mx(2)+1:2, 1: mx(3)+1:2,1:nfv) &
         +                 tmp(-2: mx(1)+2:2, 1: mx(2)+1:2, 0: mx(3)  :2,1:nfv))
!
! E1 points:
      tmp(-1: mx(1)+1:2, 0: mx(2)  :2, 1: mx(3)+1:2,1:nfv) &
         = 0.50_r8*(       tmp(-2: mx(1)  :2, 0: mx(2)  :2, 1: mx(3)+1:2,1:nfv) &
         +                 tmp( 0: mx(1)+2:2, 0: mx(2)  :2, 1: mx(3)+1:2,1:nfv))
!
! F1 points:
      tmp(-1: mx(1)+1:2, 1: mx(2)+1:2, 1: mx(3)+1:2,1:nfv) &
         = 0.50_r8*(       tmp(-2: mx(1)  :2, 1: mx(2)+1:2, 1: mx(3)+1:2,1:nfv) &
         +                 tmp( 0: mx(1)+2:2, 1: mx(2)+1:2, 1: mx(3)+1:2,1:nfv))
!
! E2 points:
      tmp(-1: mx(1)+1:2, 0: mx(2)  :2, 0: mx(3)  :2,1:nfv) &
         = 0.50_r8*(       tmp(-2: mx(1)  :2, 0: mx(2)  :2, 0: mx(3)  :2,1:nfv) &
         +                 tmp( 0: mx(1)+2:2, 0: mx(2)  :2, 0: mx(3)  :2,1:nfv))
!
! F2 points:
      tmp(-1: mx(1)+1:2, 1: mx(2)+1:2, 0: mx(3)  :2,1:nfv) &
         = 0.50_r8*(       tmp(-2: mx(1)  :2, 1: mx(2)+1:2, 0: mx(3)  :2,1:nfv) &
         +                 tmp( 0: mx(1)+2:2, 1: mx(2)+1:2, 0: mx(3)  :2,1:nfv))
!
      presult(-1: mx(1)+1  , 0: mx(2)+1  , 0: mx(3)+1  ,1:nfv) &
         =                 tmp(-1: mx(1)+1  , 0: mx(2)+1  , 0: mx(3)+1  ,1:nfv)
!
   END FUNCTION Prolongation3DV1Ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Prolongation3DV2Ex(v2) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Prolongation of v1 (east-west) edge functions by bi-linear interpolation.
!
! Input:       v1(-1:cmx(1)+1,0:cmx(2)+1,1:) (with x1 ghost cells)
!
! Output: presult(-1: mx(1)+1,1: mx(2)  ,1:)
!
      REAL(KIND=r8), DIMENSION( 0:,-1:,0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:(SIZE(v2,1)-2)*2+1, &
         -1:(SIZE(v2,2)-3)*2+1, &
         0:(SIZE(v2,3)-2)*2+1, &
         1: SIZE(v2,4)        ):: presult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:3):: cmx, mx
      REAL(KIND=r8), DIMENSION( 0:(SIZE(v2,1)-2)*2+1, &
         -2:(SIZE(v2,2)-3)*2+2, &
         0:(SIZE(v2,3)-2)*2+1, &
         1: SIZE(v2,4)        ):: tmp
!
      cmx(1) = SIZE(v2,1)-2
      cmx(2) = SIZE(v2,2)-3
      cmx(3) = SIZE(v2,3)-2
      mx(1:3) = cmx(1:3)*2
      nfv = SIZE(v2,4)
!
! Not optimized. Can save some memory here by eliminating tmp.
!
! A points:
      tmp( 0: mx(1)  :2,-2: mx(2)+2:2, 0: mx(3)  :2,1:nfv) &
         = 0.25_r8*( 3.0_r8*v2( 0:cmx(1)    ,-1:cmx(2)+1  , 0:cmx(3)    ,1:nfv) &
         +                  v2( 1:cmx(1)+1  ,-1:cmx(2)+1  , 0:cmx(3)    ,1:nfv))
!
! B points:
      tmp( 0: mx(1)  :2,-2: mx(2)+2:2, 1: mx(3)+1:2,1:nfv) &
         = 0.25_r8*( 3.0_r8*v2( 0:cmx(1)    ,-1:cmx(2)+1  , 1:cmx(3)+1  ,1:nfv) &
         +                  v2( 1:cmx(1)+1  ,-1:cmx(2)+1  , 1:cmx(3)+1  ,1:nfv))
!
! C points:
      tmp( 1: mx(1)+1:2,-2: mx(2)+2:2, 0: mx(3)  :2,1:nfv) &
         = 0.25_r8*(        v2( 0:cmx(1)    ,-1:cmx(2)+1  , 0:cmx(3)    ,1:nfv) &
         +           3.0_r8*v2( 1:cmx(1)+1  ,-1:cmx(2)+1  , 0:cmx(3)    ,1:nfv))
!
! D points:
      tmp( 1: mx(1)+1:2,-2: mx(2)+2:2, 1: mx(3)+1:2,1:nfv) &
         = 0.25_r8*(        v2( 0:cmx(1)    ,-1:cmx(2)+1  , 1:cmx(3)+1  ,1:nfv) &
         +           3.0_r8*v2( 1:cmx(1)+1  ,-1:cmx(2)+1  , 1:cmx(3)+1  ,1:nfv))
!
! A1 points:
      tmp( 0: mx(1)  :2,-2: mx(2)+2:2, 0: mx(3)  :2,1:nfv) &
         = 0.25_r8*(       tmp( 0: mx(1)  :2,-2: mx(2)+2:2, 1: mx(3)+1:2,1:nfv) &
         +          3.0_r8*tmp( 0: mx(1)  :2,-2: mx(2)+2:2, 0: mx(3)  :2,1:nfv))
!
! B1 points:
      tmp( 0: mx(1)  :2,-2: mx(2)+2:2, 1: mx(3)+1:2,1:nfv) &
         = 0.25_r8*(3.0_r8*tmp( 0: mx(1)  :2,-2: mx(2)+2:2, 1: mx(3)+1:2,1:nfv) &
         +                 tmp( 0: mx(1)  :2,-2: mx(2)+2:2, 0: mx(3)  :2,1:nfv))
!

! C1 points:
      tmp( 1: mx(1)+1:2,-2: mx(2)+2:2, 0: mx(3)  :2,1:nfv) &
         = 0.25_r8*(       tmp( 1: mx(1)+1:2,-2: mx(2)+2:2, 1: mx(3)+1:2,1:nfv) &
         +          3.0_r8*tmp( 1: mx(1)+1:2,-2: mx(2)+2:2, 0: mx(3)  :2,1:nfv))
!
! D1 points:
      tmp( 1: mx(1)+1:2,-2: mx(2)+2:2, 1: mx(3)+1:2,1:nfv) &
         = 0.25_r8*(3.0_r8*tmp( 1: mx(1)+1:2,-2: mx(2)+2:2, 1: mx(3)+1:2,1:nfv) &
         +                 tmp( 1: mx(1)+1:2,-2: mx(2)+2:2, 0: mx(3)  :2,1:nfv))
!
! E1 points:
      tmp( 0: mx(1)  :2,-1: mx(2)+1:2, 1: mx(3)+1:2,1:nfv) &
         = 0.50_r8*(       tmp( 0: mx(1)  :2,-2: mx(2)  :2, 1: mx(3)+1:2,1:nfv) &
         +                 tmp( 0: mx(1)  :2, 0: mx(2)+2:2, 1: mx(3)+1:2,1:nfv))
!
! F1 points:
      tmp( 1: mx(1)+1:2,-1: mx(2)+1:2, 1: mx(3)+1:2,1:nfv) &
         = 0.50_r8*(       tmp( 1: mx(1)+1:2,-2: mx(2)  :2, 1: mx(3)+1:2,1:nfv) &
         +                 tmp( 1: mx(1)+1:2, 0: mx(2)+2:2, 1: mx(3)+1:2,1:nfv))
!
! E2 points:
      tmp( 0: mx(1)  :2,-1: mx(2)+1:2, 0: mx(3)  :2,1:nfv) &
         = 0.50_r8*(       tmp( 0: mx(1)  :2,-2: mx(2)  :2, 0: mx(3)  :2,1:nfv) &
         +                 tmp( 0: mx(1)  :2, 0: mx(2)+2:2, 0: mx(3)  :2,1:nfv))
!
! F2 points:
      tmp( 1: mx(1)+1:2,-1: mx(2)+1:2, 0: mx(3)  :2,1:nfv) &
         = 0.50_r8*(       tmp( 1: mx(1)+1:2,-2: mx(2)  :2, 0: mx(3)  :2,1:nfv) &
         +                 tmp( 1: mx(1)+1:2, 0: mx(2)+2:2, 0: mx(3)  :2,1:nfv))
!
      presult( 0: mx(1)+1  ,-1: mx(2)+1  , 0: mx(3)+1  ,1:nfv) &
         =                 tmp( 0: mx(1)+1  ,-1: mx(2)+1  , 0: mx(3)+1  ,1:nfv)
!
   END FUNCTION Prolongation3DV2Ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Prolongation3DV3Ex(v3) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Prolongation of v1 (east-west) edge functions by bi-linear interpolation.
!
! Input:       v1(-1:cmx(1)+1,0:cmx(2)+1,1:) (with x1 ghost cells)
!
! Output: presult(-1: mx(1)+1,1: mx(2)  ,1:)
!
      REAL(KIND=r8), DIMENSION( 0:,0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:(SIZE(v3,1)-2)*2+1, &
         0:(SIZE(v3,2)-2)*2+1, &
         -1:(SIZE(v3,3)-3)*2+1, &
         1: SIZE(v3,4)        ):: presult
!
      INTEGER:: nfv
      INTEGER, DIMENSION(1:3):: cmx, mx
      REAL(KIND=r8), DIMENSION( 0:(SIZE(v3,1)-2)*2+1, &
         0:(SIZE(v3,2)-2)*2+1, &
         -2:(SIZE(v3,3)-3)*2+2, &
         1: SIZE(v3,4)        ):: tmp
!
      cmx(1) = SIZE(v3,1)-2
      cmx(2) = SIZE(v3,2)-2
      cmx(3) = SIZE(v3,3)-3
      mx(1:3) = cmx(1:3)*2
      nfv = SIZE(v3,4)
!
! Not optimized. Can save some memory here by eliminating tmp.
!
! A points:
      tmp( 0: mx(1)  :2, 0: mx(2)  :2,-2: mx(3)+2:2,1:nfv) &
         = 0.25_r8*( 3.0_r8*v3( 0:cmx(1)    , 0:cmx(2)    ,-1:cmx(3)+1  ,1:nfv) &
         +                  v3( 1:cmx(1)+1  , 0:cmx(2)    ,-1:cmx(3)+1  ,1:nfv))
!
! B points:
      tmp( 0: mx(1)  :2, 1: mx(2)+1:2,-2: mx(3)+2:2,1:nfv) &
         = 0.25_r8*( 3.0_r8*v3( 0:cmx(1)    , 1:cmx(2)+1  ,-1:cmx(3)+1  ,1:nfv) &
         +                  v3( 1:cmx(1)+1  , 1:cmx(2)+1  ,-1:cmx(3)+1  ,1:nfv))
!
! C points:
      tmp( 1: mx(1)+1:2, 0: mx(2)  :2,-2: mx(3)+2:2,1:nfv) &
         = 0.25_r8*(        v3( 0:cmx(1)    , 0:cmx(2)    ,-1:cmx(3)+1  ,1:nfv) &
         +           3.0_r8*v3( 1:cmx(1)+1  , 0:cmx(2)    ,-1:cmx(3)+1  ,1:nfv))
!
! D points:
      tmp( 1: mx(1)+1:2, 1: mx(2)+1:2,-2: mx(3)+2:2,1:nfv) &
         = 0.25_r8*(        v3( 0:cmx(1)    , 1:cmx(2)+1  ,-1:cmx(3)+1  ,1:nfv) &
         +           3.0_r8*v3( 1:cmx(1)+1  , 1:cmx(2)+1  ,-1:cmx(3)+1  ,1:nfv))
!
! A1 points:
      tmp( 0: mx(1)  :2, 0: mx(2)  :2,-2: mx(3)+2:2,1:nfv) &
         = 0.25_r8*(       tmp( 0: mx(1)  :2, 1: mx(2)+1:2,-2: mx(3)+2:2,1:nfv) &
         +          3.0_r8*tmp( 0: mx(1)  :2, 0: mx(2)  :2,-2: mx(3)+2:2,1:nfv))
!
! B1 points:
      tmp( 0: mx(1)  :2, 1: mx(2)+1:2,-2: mx(3)+2:2,1:nfv) &
         = 0.25_r8*(3.0_r8*tmp( 0: mx(1)  :2, 1: mx(2)+1:2,-2: mx(3)+2:2,1:nfv) &
         +                 tmp( 0: mx(1)  :2, 0: mx(2)  :2,-2: mx(3)+2:2,1:nfv))
!

! C1 points:
      tmp( 1: mx(1)+1:2, 0: mx(2)  :2,-2: mx(3)+2:2,1:nfv) &
         = 0.25_r8*(       tmp( 1: mx(1)+1:2, 1: mx(2)+1:2,-2: mx(3)+2:2,1:nfv) &
         +          3.0_r8*tmp( 1: mx(1)+1:2, 0: mx(2)  :2,-2: mx(3)+2:2,1:nfv))
!
! D1 points:
      tmp( 1: mx(1)+1:2, 1: mx(2)+1:2,-2: mx(3)+2:2,1:nfv) &
         = 0.25_r8*(3.0_r8*tmp( 1: mx(1)+1:2, 1: mx(2)+1:2,-2: mx(3)+2:2,1:nfv) &
         +                 tmp( 1: mx(1)+1:2, 0: mx(2)  :2,-2: mx(3)+2:2,1:nfv))
!
! E1 points:
      tmp( 0: mx(1)  :2, 1: mx(2)+1:2,-1: mx(3)+1:2,1:nfv) &
         = 0.50_r8*(       tmp( 0: mx(1)  :2, 1: mx(2)+1:2,-2: mx(3)  :2,1:nfv) &
         +                 tmp( 0: mx(1)  :2, 1: mx(2)+1:2, 0: mx(3)+2:2,1:nfv))
!
! F1 points:
      tmp( 1: mx(1)+1:2, 1: mx(2)+1:2,-1: mx(3)+1:2,1:nfv) &
         = 0.50_r8*(       tmp( 1: mx(1)+1:2, 1: mx(2)+1:2,-2: mx(3)  :2,1:nfv) &
         +                 tmp( 1: mx(1)+1:2, 1: mx(2)+1:2, 0: mx(3)+2:2,1:nfv))
!
! E2 points:
      tmp( 0: mx(1)  :2, 0: mx(2)  :2,-1: mx(3)+1:2,1:nfv) &
         = 0.50_r8*(       tmp( 0: mx(1)  :2, 0: mx(2)  :2,-2: mx(3)  :2,1:nfv) &
         +                 tmp( 0: mx(1)  :2, 0: mx(2)  :2, 0: mx(3)+2:2,1:nfv))
!
! F2 points:
      tmp( 1: mx(1)+1:2, 0: mx(2)  :2,-1: mx(3)+1:2,1:nfv) &
         = 0.50_r8*(       tmp( 1: mx(1)+1:2, 0: mx(2)  :2,-2: mx(3)  :2,1:nfv) &
         +                 tmp( 1: mx(1)+1:2, 0: mx(2)  :2, 0: mx(3)+2:2,1:nfv))
!
      presult( 0: mx(1)+1  , 0: mx(2)+1  ,-1: mx(3)+1  ,1:nfv) &
         =                 tmp( 0: mx(1)+1  , 0: mx(2)+1  ,-1: mx(3)+1  ,1:nfv)
!
   END FUNCTION Prolongation3DV3Ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION TriLinProlongationP1(a) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! This trilinear prolongation algorithm assumes exactly one ghost layer, i.e.,
! mbc = 1.
!
      REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
         0:(SIZE(a,2)-2)*2+1, &
         0:(SIZE(a,3)-2)*2+1, &
         1: SIZE(a,4)):: presult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:3):: cmx, mx
      REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
         0: SIZE(a,2)-2   +1, &
         0: SIZE(a,3)-2   +1, &
         1: SIZE(a,4)       ):: b
      REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
         0:(SIZE(a,2)-2)*2+1, &
         0: SIZE(a,3)-2   +1, &
         1: SIZE(a,4)       ):: c
!
      cmx(1) = SIZE(a,1)-2
      cmx(2) = SIZE(a,2)-2
      cmx(3) = SIZE(a,3)-2
      mx(1:3) = cmx(1:3)*2
      ncv = SIZE(a,4)
!
! Later, use presult in place of b to save storage space:
!
! Linear interpolation in the x-direction first:
      b(0: mx(1)  :2,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv) &
         = 3.0_r8*a(0:cmx(1)    ,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv) &
         +        a(1:cmx(1)+1  ,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv)
      b(1: mx(1)+1:2,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv) &
         =        a(0:cmx(1)    ,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv) &
         + 3.0_r8*a(1:cmx(1)+1  ,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv)
!
! Linear interpolation in the y-direction second:
      c(0: mx(1)+1  ,0: mx(2)  :2,0:cmx(3)+1  ,1:ncv) &
         = 3.0_r8*b(0: mx(1)+1  ,0:cmx(2)    ,0:cmx(3)+1  ,1:ncv) &
         +        b(0: mx(1)+1  ,1:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv)
      c(0: mx(1)+1  ,1: mx(2)+1:2,0:cmx(3)+1  ,1:ncv) &
         =        b(0: mx(1)+1  ,0:cmx(2)    ,0:cmx(3)+1  ,1:ncv) &
         + 3.0_r8*b(0: mx(1)+1  ,1:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv)
!
! Linear interpolation in the z-direction third:
      presult(0: mx(1)+1  ,0: mx(2)+1  ,0: mx(3)  :2,1:ncv) &
         = (3.0_r8*c(0: mx(1)+1  ,0: mx(2)+1  ,0:cmx(3)    ,1:ncv) &
         +         c(0: mx(1)+1  ,0: mx(2)+1  ,1:cmx(3)+1  ,1:ncv))/64.0_r8
      presult(0: mx(1)+1  ,0: mx(2)+1  ,1: mx(3)+1:2,1:ncv) &
         = (       c(0: mx(1)+1  ,0: mx(2)+1  ,0:cmx(3)    ,1:ncv) &
         +  3.0_r8*c(0: mx(1)+1  ,0: mx(2)+1  ,1:cmx(3)+1  ,1:ncv))/64.0_r8
!
   END FUNCTION TriLinProlongationP1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION TriLinProlongationP2(a) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! This trilinear prolongation algorithm assumes exactly two ghost layers, i.e.,
! mbc = 2.
!
      REAL(KIND=r8), DIMENSION(-1:,-1:,-1:,1:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
         -1:(SIZE(a,2)-4)*2+2, &
         -1:(SIZE(a,3)-4)*2+2, &
         1: SIZE(a,4)       ):: presult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:3):: cmx, mx
      REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
         -1: SIZE(a,2)-4   +2, &
         -1: SIZE(a,3)-4   +2, &
         1: SIZE(a,4)       ):: b
      REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
         -1:(SIZE(a,2)-4)*2+2, &
         -1: SIZE(a,3)-4   +2, &
         1: SIZE(a,4)       ):: c
!
      cmx(1) = SIZE(a,1)-4
      cmx(2) = SIZE(a,2)-4
      cmx(3) = SIZE(a,3)-4
      mx(1:3) = cmx(1:3)*2
      ncv = SIZE(a,4)
!
! Linear interpolation in the x-direction first:
      b(-1: mx(1)+1:2,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv) &
         = 3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv) &
         +        a(-1:cmx(1)    ,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv)
      b( 0: mx(1)+2:2,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv) &
         = 3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv) &
         +        a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv)
!
! Linear interpolation in the y-direction second:
      c(-1: mx(1)+2  ,-1: mx(2)+1:2,-1:cmx(3)+2  ,1:ncv) &
         = 3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2  ,1:ncv) &
         +        b(-1: mx(1)+2  ,-1:cmx(2)    ,-1:cmx(3)+2  ,1:ncv)
      c(-1: mx(1)+2  , 0: mx(2)+2:2,-1:cmx(3)+2  ,1:ncv) &
         = 3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2  ,1:ncv) &
         +        b(-1: mx(1)+2  , 1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv)
!
! Linear interpolation in the z-direction second:
      presult(-1: mx(1)+2  ,-1: mx(2)+2  ,-1: mx(3)+1:2,1:ncv) &
         = (3.0_r8*c(-1: mx(1)+2  ,-1: mx(2)+2  , 0:cmx(3)+1  ,1:ncv) &
         +         c(-1: mx(1)+2  ,-1: mx(2)+2  ,-1:cmx(3)    ,1:ncv))/64.0_r8
      presult(-1: mx(1)+2  ,-1: mx(2)+2  , 0: mx(3)+2:2,1:ncv) &
         = (3.0_r8*c(-1: mx(1)+2  ,-1: mx(2)+2  , 0:cmx(3)+1  ,1:ncv) &
         +         c(-1: mx(1)+2  ,-1: mx(2)+2  , 1:cmx(3)+2  ,1:ncv))/64.0_r8
!
   END FUNCTION TriLinProlongationP2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION TriLinProlongationP1MC(a) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! This trilinear prolongation algorithm assumes exactly one ghost layer, i.e.,
! mbc = 1.
!
      REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
         0:(SIZE(a,2)-2)*2+1, &
         0:(SIZE(a,3)-2)*2+1, &
         1: SIZE(a,4)):: presult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:3):: cmx, mx
      REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
         0: SIZE(a,2)-2   +1, &
         0: SIZE(a,3)-2   +1, &
         1: SIZE(a,4)       ):: b
      REAL(KIND=r8), DIMENSION(0:(SIZE(a,1)-2)*2+1, &
         0:(SIZE(a,2)-2)*2+1, &
         0: SIZE(a,3)-2   +1, &
         1: SIZE(a,4)       ):: c
      REAL(KIND=r8), DIMENSION(1: SIZE(a,1)-2     , &
         1: SIZE(a,2)-2     , &
         1: SIZE(a,3)-2     , &
         1: SIZE(a,4)       ):: cor
!
      cmx(1) = SIZE(a,1)-2
      cmx(2) = SIZE(a,2)-2
      cmx(3) = SIZE(a,3)-2
      mx(1:3) = cmx(1:3)*2
      ncv = SIZE(a,4)
!
! Later, use presult in place of b to save storage space:
!
! Linear interpolation in the x-direction first:
      b(0: mx(1)  :2,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv) &
         = 3.0_r8*a(0:cmx(1)    ,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv) &
         +        a(1:cmx(1)+1  ,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv)
      b(1: mx(1)+1:2,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv) &
         =        a(0:cmx(1)    ,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv) &
         + 3.0_r8*a(1:cmx(1)+1  ,0:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv)
!
! Linear interpolation in the y-direction second:
      c(0: mx(1)+1  ,0: mx(2)  :2,0:cmx(3)+1  ,1:ncv) &
         = 3.0_r8*b(0: mx(1)+1  ,0:cmx(2)    ,0:cmx(3)+1  ,1:ncv) &
         +        b(0: mx(1)+1  ,1:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv)
      c(0: mx(1)+1  ,1: mx(2)+1:2,0:cmx(3)+1  ,1:ncv) &
         =        b(0: mx(1)+1  ,0:cmx(2)    ,0:cmx(3)+1  ,1:ncv) &
         + 3.0_r8*b(0: mx(1)+1  ,1:cmx(2)+1  ,0:cmx(3)+1  ,1:ncv)
!
! Linear interpolation in the z-direction third:
      presult(0: mx(1)+1  ,0: mx(2)+1  ,0: mx(3)  :2,1:ncv) &
         = (3.0_r8*c(0: mx(1)+1  ,0: mx(2)+1  ,0:cmx(3)    ,1:ncv) &
         +         c(0: mx(1)+1  ,0: mx(2)+1  ,1:cmx(3)+1  ,1:ncv))/64.0_r8
      presult(0: mx(1)+1  ,0: mx(2)+1  ,1: mx(3)+1:2,1:ncv) &
         = (       c(0: mx(1)+1  ,0: mx(2)+1  ,0:cmx(3)    ,1:ncv) &
         +  3.0_r8*c(0: mx(1)+1  ,0: mx(2)+1  ,1:cmx(3)+1  ,1:ncv))/64.0_r8
!
! The mass correction is computed only for regular cells; not for ghost cells.
! Compute mass correction:
      cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv) &
         =                 a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv) &
         - 0.125_r8*(presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         +           presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         +           presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         +           presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         +           presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         +           presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         +           presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         +           presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv))
!
! Add the mass correction:
      presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         = presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         = presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         = presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         = presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         = presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         = presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         = presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         = presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
!
   END FUNCTION TriLinProlongationP1MC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION TriLinProlongationP2MC(a) RESULT(presult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! This trilinear prolongation algorithm assumes exactly two ghost layers, i.e.,
! mbc = 2.
!
      REAL(KIND=r8), DIMENSION(-1:,-1:,-1:,1:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
         -1:(SIZE(a,2)-4)*2+2, &
         -1:(SIZE(a,3)-4)*2+2, &
         1: SIZE(a,4)       ):: presult
!
      INTEGER:: ncv
      INTEGER, DIMENSION(1:3):: cmx, mx
      REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
         -1: SIZE(a,2)-4   +2, &
         -1: SIZE(a,3)-4   +2, &
         1: SIZE(a,4)       ):: b
      REAL(KIND=r8), DIMENSION(-1:(SIZE(a,1)-4)*2+2, &
         -1:(SIZE(a,2)-4)*2+2, &
         -1: SIZE(a,3)-4   +2, &
         1: SIZE(a,4)       ):: c
      REAL(KIND=r8), DIMENSION( 1: SIZE(a,1)-4     , &
         1: SIZE(a,2)-4     , &
         1: SIZE(a,3)-4     , &
         1: SIZE(a,4)       ):: cor
!
      cmx(1) = SIZE(a,1)-4
      cmx(2) = SIZE(a,2)-4
      cmx(3) = SIZE(a,3)-4
      mx(1:3) = cmx(1:3)*2
      ncv = SIZE(a,4)
!
! Linear interpolation in the x-direction first:
      b(-1: mx(1)+1:2,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv) &
         = 3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv) &
         +        a(-1:cmx(1)    ,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv)
      b( 0: mx(1)+2:2,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv) &
         = 3.0_r8*a( 0:cmx(1)+1  ,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv) &
         +        a( 1:cmx(1)+2  ,-1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv)
!
! Linear interpolation in the y-direction second:
      c(-1: mx(1)+2  ,-1: mx(2)+1:2,-1:cmx(3)+2  ,1:ncv) &
         = 3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2  ,1:ncv) &
         +        b(-1: mx(1)+2  ,-1:cmx(2)    ,-1:cmx(3)+2  ,1:ncv)
      c(-1: mx(1)+2  , 0: mx(2)+2:2,-1:cmx(3)+2  ,1:ncv) &
         = 3.0_r8*b(-1: mx(1)+2  , 0:cmx(2)+1  ,-1:cmx(3)+2  ,1:ncv) &
         +        b(-1: mx(1)+2  , 1:cmx(2)+2  ,-1:cmx(3)+2  ,1:ncv)
!
! Linear interpolation in the z-direction third:
      presult(-1: mx(1)+2  ,-1: mx(2)+2  ,-1: mx(3)+1:2,1:ncv) &
         = (3.0_r8*c(-1: mx(1)+2  ,-1: mx(2)+2  , 0:cmx(3)+1  ,1:ncv) &
         +         c(-1: mx(1)+2  ,-1: mx(2)+2  ,-1:cmx(3)    ,1:ncv))/64.0_r8
      presult(-1: mx(1)+2  ,-1: mx(2)+2  , 0: mx(3)+2:2,1:ncv) &
         = (3.0_r8*c(-1: mx(1)+2  ,-1: mx(2)+2  , 0:cmx(3)+1  ,1:ncv) &
         +         c(-1: mx(1)+2  ,-1: mx(2)+2  , 1:cmx(3)+2  ,1:ncv))/64.0_r8
!
! The mass correction is computed only for regular cells; not for ghost cells.
! Compute mass correction:
      cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv) &
         = a(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv) &
         - 0.125_r8*(presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         +           presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         +           presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         +           presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         +           presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         +           presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         +           presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         +           presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv))
!
! Add the mass correction:
      presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         = presult(2: mx(1)  :2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         = presult(1: mx(1)-1:2,2: mx(2)  :2,2: mx(3)  :2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         = presult(2: mx(1)  :2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         = presult(1: mx(1)-1:2,1: mx(2)-1:2,2: mx(3)  :2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         = presult(2: mx(1)  :2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         = presult(1: mx(1)-1:2,2: mx(2)  :2,1: mx(3)-1:2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         = presult(2: mx(1)  :2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
      presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         = presult(1: mx(1)-1:2,1: mx(2)-1:2,1: mx(3)-1:2,1:ncv) &
         +     cor(1:cmx(1)    ,1:cmx(2)    ,1:cmx(3)    ,1:ncv)
!
   END FUNCTION TriLinProlongationP2MC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Advecx2D(v1,v2,h) RESULT(udivresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! UNDIVIDED divergence of the 2D flux function
!
!  f = (f1(0:mx(1),1:mx(2)),f2(1:mx(1),0:mx(2))).
!
! The result is stored in udivresult(1:mx(1),1:mx(2)).
!
      REAL(KIND=r8), DIMENSION(-1:, 0:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION( 0:,-1:), INTENT(IN):: v2
      real(KIND=r8), intent(in) :: h
      REAL(KIND=r8), DIMENSION(0:SIZE(v1,1)-3, &
         1:SIZE(v1,2)-2):: udivresult
!
      INTEGER, DIMENSION(1:2):: mx
!
      mx(1) = SIZE(v1,1)-3
      mx(2) = SIZE(v1,2)-2
!
      udivresult(0:mx(1),1:mx(2)) =          v1( 0:mx(1)  ,1:mx(2)  )             &
         *         (v1( 1:mx(1)+1,1:mx(2)  )             &
         -          v1(-1:mx(1)-1,1:mx(2)  ))/(2.0_r8*h) &
         + 0.25_r8*(v2( 1:mx(1)+1,1:mx(2)  )             &
         +          v2( 0:mx(1)  ,1:mx(2)  ))            &
         *         (v1( 0:mx(1)  ,2:mx(2)+1)             &
         -          v1( 0:mx(1)  ,1:mx(2)  ))/h          &
         + 0.25_r8*(v2( 1:mx(1)+1,0:mx(2)-1)             &
         +          v2( 0:mx(1)  ,0:mx(2)-1))            &
         *         (v1( 0:mx(1)  ,1:mx(2)  )             &
         -          v1( 0:mx(1)  ,0:mx(2)-1))/h
!
   END FUNCTION Advecx2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Advecy2D(v1,v2,h) RESULT(udivresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! UNDIVIDED divergence of the 2D flux function
!
!  f = (f1(0:mx(1),1:mx(2)),f2(1:mx(1),0:mx(2))).
!
! The result is stored in udivresult(1:mx(1),1:mx(2)).
!
      REAL(KIND=r8), DIMENSION(-1:, 0:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION( 0:,-1:), INTENT(IN):: v2
      real(KIND=r8), intent(in) :: h
      REAL(KIND=r8), DIMENSION(1:SIZE(v1,1)-3, &
         0:SIZE(v1,2)-2):: udivresult
!
      INTEGER, DIMENSION(1:2):: mx
!
      mx(1) = SIZE(v1,1)-3
      mx(2) = SIZE(v1,2)-2
!
      udivresult(1:mx(1),0:mx(2)) =          v2(1:mx(1)  , 0:mx(2)  )             &
         *         (v2(1:mx(1)  , 1:mx(2)+1)             &
         -          v2(1:mx(1)  ,-1:mx(2)-1))/(2.0_r8*h) &
         + 0.25_r8*(v1(1:mx(1)  , 1:mx(2)+1)             &
         +          v1(1:mx(1)  , 0:mx(2)  ))            &
         *         (v2(2:mx(1)+1, 0:mx(2)  )             &
         -          v2(1:mx(1)  , 0:mx(2)  ))/h          &
         + 0.25_r8*(v1(0:mx(1)-1, 1:mx(2)+1)             &
         +          v1(0:mx(1)-1, 0:mx(2)  ))            &
         *         (v2(1:mx(1)  , 0:mx(2)  )             &
         -          v2(0:mx(1)-1, 0:mx(2)  ))/h
!
   END FUNCTION Advecy2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function ULapVar3D(c,a) result(lapvarresult)
!use GridUtilities, only: UDiv2D
      use ProblemDef

      implicit none
!
! Level independent, 2D UNDIVIDED laplacian operator for non-constant
! diffusivity, e.g div*(c*grad(a)).
!
      real(kind=r8), dimension(0:,0:,0:), intent(in):: c
      real(kind=r8), dimension(0:,0:,0:), intent(in):: a
      real(kind=r8), dimension(1:size(a,1)-2,1:size(a,2)-2,1:size(a,3)-2):: lapvarresult
!
      INTEGER, DIMENSION(1:3):: mx
      real(kind=r8), dimension(0:size(a,1)-2  ,1:size(a,2)-2  ,1:size(a,3)-2  ):: f1
      real(kind=r8), dimension(1:size(a,1)-2  ,0:size(a,2)-2  ,1:size(a,3)-2  ):: f2
      real(kind=r8), dimension(1:size(a,1)-2  ,1:size(a,2)-2  ,0:size(a,3)-2  ):: f3
      real(kind=r8), dimension(0:size(a,1)-2+1,0:size(a,2)-2+1,0:size(a,3)-2+1):: m
!
      mx(1) = size(a,1)-2; mx(2) = size(a,2)-2; mx(3) = size(a,3)-2
!
      m(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1) = 0.5_r8*(c(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1))
!
! Calculate the UNDIVIDED 2D flux function:
      f1(0:mx(1),1:mx(2),1:mx(3)) &
         = (m(1:mx(1)+1,1:mx(2),1:mx(3))+m(0:mx(1),1:mx(2),1:mx(3))) &
         * (a(1:mx(1)+1,1:mx(2),1:mx(3))-a(0:mx(1),1:mx(2),1:mx(3)))
!
      f2(1:mx(1),0:mx(2),1:mx(3)) &
         = (m(1:mx(1),1:mx(2)+1,1:mx(3))+m(1:mx(1),0:mx(2),1:mx(3))) &
         * (a(1:mx(1),1:mx(2)+1,1:mx(3))-a(1:mx(1),0:mx(2),1:mx(3)))
!
      f3(1:mx(1),1:mx(2),0:mx(3)) &
         = (m(1:mx(1),1:mx(2),1:mx(3)+1)+m(1:mx(1),1:mx(2),0:mx(3))) &
         * (a(1:mx(1),1:mx(2),1:mx(3)+1)-a(1:mx(1),1:mx(2),0:mx(3)))
!
! Calculate the UNDIVIDED divergence of the flux:
      lapvarresult(1:mx(1),1:mx(2),1:mx(3)) = UDiv3D(f1(0:mx(1),1:mx(2),1:mx(3)), &
         f2(1:mx(1),0:mx(2),1:mx(3)), &
         f3(1:mx(1),1:mx(2),0:mx(3)))

   end function ULapVar3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   elemental function Df(c) result(dfresult)
! Derivative of f-function in equation for chemical potential
! f=0.25*c^2*(1-c)^2
! f' = c*(c-1)*(c-1/2)
      use ProblemDef
      implicit none
!
      real(KIND=r8), intent(in):: c
      real(KIND=r8):: dfresult
!
      dfresult = 72.0_r8*c*(c-1.0_r8)*(c-0.5_r8)
!
   end function Df
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   elemental function D2f(c) result(dfresult)
! Derivative of f-function in equation for chemical potential
! f=0.25*c^2*(1-c)^2
! f' = c*(c-1)*(c-1/2)
      use ProblemDef
      implicit none
!
      real(KIND=r8), intent(in):: c
      real(KIND=r8):: dfresult
!
      dfresult = 36.0_r8*(6.0_r8*c*c - 6.0_r8*c + 1.0_r8)
!
   end function D2f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elemental function Gf(c) result(dfresult)
! Derivative of f-function in equation for chemical potential
! f=0.25*c^2*(1-c)^2
! f' = c*(c-1)*(c-1/2)
      use ProblemDef
      implicit none
!
      real(KIND=r8), intent(in):: c
      real(KIND=r8):: dfresult
!
      dfresult = 36.0_r8*c**2*(c-1.0_r8)**2
!
   end function Gf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function ULapVar2D(c,a) result(lapvarresult)
!use GridUtilities, only: UDiv2D
      use ProblemDef

      implicit none
!
! Level independent, 2D UNDIVIDED laplacian operator for non-constant
! diffusivity, e.g div*(c*grad(a)).
!
      real(kind=r8), dimension(0:,0:), intent(in):: c
      real(kind=r8), dimension(0:,0:), intent(in):: a
      real(kind=r8), dimension(1:size(a,1)-2,1:size(a,2)-2):: lapvarresult
!
      INTEGER, DIMENSION(1:2):: mx
      real(kind=r8), dimension(0:size(a,1)-2  ,1:size(a,2)-2  ):: f1
      real(kind=r8), dimension(1:size(a,1)-2  ,0:size(a,2)-2  ):: f2
      real(kind=r8), dimension(0:size(a,1)-2+1,0:size(a,2)-2+1):: m
!
      mx(1) = size(a,1)-2; mx(2) = size(a,2)-2
!
      m(0:mx(1)+1,0:mx(2)+1) = 0.5_r8*(c(0:mx(1)+1,0:mx(2)+1))
!
! Calculate the UNDIVIDED 2D flux function:
      f1(0:mx(1),1:mx(2)) &
         = (m(1:mx(1)+1,1:mx(2))+m(0:mx(1),1:mx(2))) &
         * (a(1:mx(1)+1,1:mx(2))-a(0:mx(1),1:mx(2)))
!
      f2(1:mx(1),0:mx(2)) &
         = (m(1:mx(1),1:mx(2)+1)+m(1:mx(1),0:mx(2))) &
         * (a(1:mx(1),1:mx(2)+1)-a(1:mx(1),0:mx(2)))
!
! Calculate the UNDIVIDED divergence of the flux:
      lapvarresult(1:mx(1),1:mx(2)) = UDiv2D(f1(0:mx(1),1:mx(2)), &
         f2(1:mx(1),0:mx(2)))

   end function ULapVar2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION AbsGrad2D(f1) RESULT(agresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! UNDIVIDED divergence of the 2D flux function
!
! The result is stored in udivresult(1:mx(1),1:mx(2)).
!
      REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: f1
      REAL(KIND=r8), DIMENSION(1:SIZE(f1,1)-2, &
         1:SIZE(f1,2)-2):: agresult
!
      INTEGER, DIMENSION(1:2):: mx
!
      mx(1) = SIZE(f1,1)-2
      mx(2) = SIZE(f1,2)-2
!
      agresult(1:mx(1),1:mx(2)) = sqrt((f1(2:mx(1)+1,1:mx(2)  ) &
         -  f1(0:mx(1)-1,1:mx(2)  ))**2 &
         + (f1(1:mx(1)  ,2:mx(2)+1) &
         -  f1(1:mx(1)  ,0:mx(2)-1))**2)
!
   END FUNCTION AbsGrad2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION AbsGrad(f1,h) RESULT(udivresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! UNDIVIDED divergence of the 2D flux function
!
!  f = (f1(0:mx(1),1:mx(2)),f2(1:mx(1),0:mx(2))).
!
! The result is stored in udivresult(1:mx(1),1:mx(2)).
!
      REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: f1
      REAL(KIND=r8),                   INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:SIZE(f1,1)-2, &
         1:SIZE(f1,2)-2):: udivresult
!
      INTEGER, DIMENSION(1:2):: mx
!
      mx(1) = SIZE(f1,1)-2
      mx(2) = SIZE(f1,2)-2
!
      udivresult(1:mx(1),1:mx(2)) = sqrt(  (f1(2:mx(1)+1,1:mx(2)  ) &
         -  f1(0:mx(1)-1,1:mx(2)  ))**2 &
         + (f1(1:mx(1)  ,2:mx(2)+1) &
         -  f1(1:mx(1)  ,0:mx(2)-1))**2)/(2.0_r8*h)
!
   END FUNCTION AbsGrad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   FUNCTION LapVar2D(c,a,h) RESULT(lapvarresult)
!use NodeInfoDef
      USE ProblemDef
      IMPLICIT none
!
! Level independent, 2D laplacian operator for non-constant
! diffusivity, e.g div*(c*grad(a)).
!
      REAL(KIND=r8), dimension(0:,0:), intent(in):: c
      REAL(KIND=r8), dimension(0:,0:), intent(in):: a
      REAL(KIND=r8), intent(in) :: h
      REAL(KIND=r8), dimension(1:size(a,1)-2,1:size(a,2)-2):: lapvarresult
!
      INTEGER, DIMENSION(1:2):: mx
      real(kind=r8), dimension(1:size(a,1)-2,1:size(a,2)-2):: xs
      real(kind=r8), dimension(1:size(a,1)-2) :: x
      integer :: i,j
!
      mx(1) = size(a,1)-2; mx(2) = size(a,2)-2

      lapvarresult = ULapVar2D(c,a)/(h*h)

   end function LapVar2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE GridUtilities

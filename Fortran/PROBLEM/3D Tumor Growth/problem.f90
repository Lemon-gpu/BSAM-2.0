MODULE Problem
!
! This MODULE CONTAINS routines for initializing the problem and for driving
! the Multigrid routines in the MODULE AFASRoutines.  In particular, the
! routines in AFASRoutines USE RelaxGrid#D, Operator#D, Source#D, and
! SourceCorrection#D where # = 2 and 3.
!
! Templates for the Cahn-Hilliard Equation.  We use Eyre's gradient stable
! scheme for time stepping.
!
CONTAINS
!
   SUBROUTINE SetupProblem
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      INTEGER:: ierror
      NAMELIST/problemdata/alpha, eps, gamma, shearvel,pp, AA, Zeta, epep,DtDh,MD, &
         xx11, xx12, yy11, yy12, &
         xx21, xx22, yy21, yy22, &
         xx31, xx32, yy31, yy32, &
         xx41, xx42, yy41, yy42, &
         x0,y0,x1,y1,x2,y2,x3,y3, sigmax, sigmay, &
         DD, FFF, Tau, Rou00, Rou01, Gama, kk1, kk2, dede1, dede2,&
         kkp1, kkp2, kkm1, kkm2
!
      OPEN(UNIT=75,FILE='problemdata.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
      IF(ierror/=0) THEN
         PRINT *,'Error opening input file problemdata.dat. Program stop.'
         STOP
      END IF
      READ(75,NML=problemdata)
      CLOSE(75)
      OPEN(UNIT=76,FILE='output.dat',STATUS='OLD',ACTION='WRITE',FORM='FORMATTED', &
         POSITION='APPEND')
      WRITE(76,NML=problemdata)
      CLOSE(76)
!
      epep2 = epep*epep
!
   END SUBROUTINE SetupProblem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE AfterRun
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
   END SUBROUTINE AfterRun
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE Initialize2D(q,v1,v2,mx,h,xlower)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap2D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(OUT):: q
      REAL(KIND=r8), DIMENSION(0:,1:,1:), INTENT(OUT):: v1
      REAL(KIND=r8), DIMENSION(1:,0:,1:), INTENT(OUT):: v2
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
!
      INTEGER:: i, j
      REAL(KIND=r8):: pi, r1, r2, r3, r4, tmp, x, y, z, RR, ttt, Rinf, R0
!
!delta = epep*SQRT(2.0_r8)
!PI = acos(-1.0_r8)
!
      ttt= 0.001
      Rinf = 1.8_r8
      R0 = 0.5_r8
      RR = SQRT(R0**2+FFF*Rinf**2*ttt)
!
      DO i = 1, mx(1)
         x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
         DO j = 1, mx(2)
            y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
            r1 = SQRT(((x-12.8_r8)/2.1_r8)**2+((y-12.8_r8)/1.9)**2)

!
            r1 = SQRT(((x-12.0_r8)/2.0_r8)**2+((y-12.8_r8)/2.0_r8)**2)
            r2 = SQRT(((x-12.8_r8)/2.0_r8)**2+((y-12.0_r8)/2.0_r8)**2)
            r3 = SQRT(((x-13.6_r8)/2.0_r8)**2+((y-12.8_r8)/2.0_r8)**2)

            !
            tmp=	 -0.5_r8*(1.0_r8+TANH((r1-0.5_r8)/(2.0_r8*SQRT(2.0_r8)*epep))) &
               *0.5_r8*(1.0_r8+TANH((r2-0.5_r8)/(2.0_r8*SQRT(2.0_r8)*epep)))&
               *0.5_r8*(1.0_r8+TANH((r3-0.6_r8)/(2.0_r8*SQRT(2.0_r8)*epep)))&
               +1.0_r8
            q(i,j,1) = tmp
!
            q(i,j,2) = 0.0_r8
!
            q(i,j,3) = 1.0_r8                                                 ! mu
!
            q(i,j,4) = max(1.0_r8, 2.0_r8-(x-12.8)**2/12.8_r8**2)                                                 ! mu
!
         END DO
      END DO
!
!
   END SUBROUTINE Initialize2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE Initialize3D(q,v1,v2,v3,mx,h,xlower)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(OUT):: q
      REAL(KIND=r8), DIMENSION(0:,1:,1:,1:), INTENT(OUT):: v1
      REAL(KIND=r8), DIMENSION(1:,0:,1:,1:), INTENT(OUT):: v2
      REAL(KIND=r8), DIMENSION(1:,1:,0:,1:), INTENT(OUT):: v3
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
!
      INTEGER:: i, j, k
      REAL(KIND=r8):: r1, r2, r3, r4, tmp, x, y, z, randdd
!
      v1 = 0.1_r8
      v2 = 0.1_r8
      v3 = 0.1_r8

!
      DO i = 1, mx(1)
         x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
         DO j = 1, mx(2)
            y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
            DO k = 1, mx(3)
               z = (REAL(k,KIND=r8)-0.5_r8)*h+xlower(3)
!
               r1 = SQRT(((x-12.8_r8)/2.1_r8)**2+((y-12.8_r8)/1.9)**2+((z-12.8_r8)/1.9)**2)
!
               r1 = SQRT(((x-12.8_r8)/1.5_r8)**2+((y-12.8_r8)/0.5)**2+((z-12.8_r8)/0.5)**2)
               r2 = SQRT(((x-12.8_r8)/0.5_r8)**2+((y-12.8_r8)/1.5)**2+((z-12.8_r8)/0.5)**2)
               r3 = SQRT(((x-12.8_r8)/0.5_r8)**2+((y-12.8_r8)/0.5)**2+((z-12.8_r8)/1.5)**2)
!
!    r1 = SQRT((x-6.4_r8)**2+(y-3.2_r8)**2)
!    q(i,j,1) =  -TANH((r1-1.6)/(2.0*SQRT(2.0)*epep))

!kissing bubble

               tmp = 1.0_r8 - (1.0_r8 - 0.5_r8*(1.0_r8-TANH((r1-2.0_r8)/(2.0_r8*SQRT(2.0_r8)*epep)))) &
                  * (1.0_r8 - 0.5_r8*(1.0_r8-TANH((r2-2.0_r8)/(2.0_r8*SQRT(2.0_r8)*epep)))) &
                  * (1.0_r8 - 0.5_r8*(1.0_r8-TANH((r3-2.0_r8)/(2.0_r8*SQRT(2.0_r8)*epep))))
!
               q(i,j,k,1) = tmp!
               q(i,j,k,2) = 0.0!
               Call random_number(randdd)
               q(i,j,k,3) = randdd
            END DO
         END DO
      END DO
!
   END SUBROUTINE Initialize3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE AfterStep2D(q,qc,mx,h,xlower,level)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: q
      REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qc
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      INTEGER, INTENT(IN):: level
!
      INTEGER:: tmp
      INTEGER, DIMENSION(1:2):: cmx
      REAL(KIND=r8):: ch, ch2, h2
!
      cmx(1:2) = mx(1:2)/2
      h2 = h*h; ch = h*2.0_r8; ch2 = ch*ch
!
      IF(level==0) THEN
         integralresult(1) &
            = integralresult(1)+h2*SUM(q(1:mx(1),1:mx(2),1))
      ELSE
         integralresult(1) &
            = integralresult(1)+h2*SUM(q( 1: mx(1),1: mx(2),1)) &
            -                  ch2*SUM(qc(1:cmx(1),1:cmx(2),1))
      END IF
!
   END SUBROUTINE AfterStep2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE AfterStep3D(q,qc,mx,h,xlower,level)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: q
      REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qc
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      INTEGER, INTENT(IN):: level
!
      INTEGER:: tmp
      INTEGER, DIMENSION(1:3):: cmx
      REAL(KIND=r8):: ch, ch3, h3
!
      cmx(1:3) = mx(1:3)/2
      h3 = h*h*h; ch = h*2.0_r8; ch3 = ch*ch*ch
!
      IF(level==0) THEN
         integralresult(1) &
            = integralresult(1)+h3*SUM(q(1:mx(1),1:mx(2),1:mx(3),1))
      ELSE
         integralresult(1) &
            = integralresult(1)+h3*SUM(q( 1: mx(1),1: mx(2),1: mx(3),1)) &
            -                  ch3*SUM(qc(1:cmx(1),1:cmx(2),1:cmx(3),1))
      END IF
!
   END SUBROUTINE AfterStep3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE SetAux(info)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
!
      INTEGER:: i, j, k
      REAL(KIND=r8):: time, x, y, z
!
!time = info%time
!
   END SUBROUTINE SetAux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE SetSrc(info)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
!
      INTEGER:: i, j, k
      REAL(KIND=r8):: time, x, y, z
!
!time = info%time
!
   END SUBROUTINE SetSrc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE SetErrFlagsUser2D(qrte,qnew,errorflags,mx,cmx,h,xlower,level)
      USE NodeInfoDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: qrte
      REAL(KIND=r8), DIMENSION(0:,0:,1:), INTENT(IN):: qnew
      INTEGER, DIMENSION(1:,1:), INTENT(OUT):: errorflags
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      INTEGER, DIMENSION(1:2), INTENT(IN):: cmx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      INTEGER, INTENT(IN):: level
!
      INTEGER:: i, j
      REAL(KIND=r8):: x, y
!
      errorflags(1:mx(1),1:mx(2)) = 0
!
! Error tagging based on the second differences of the composition:
      WHERE(  ABS(       qnew(2:mx(1)+1,1:mx(2)  ,1)  &
         -     2.0_r8*qnew(1:mx(1)  ,1:mx(2)  ,1)  &
         +            qnew(0:mx(1)-1,1:mx(2)  ,1)) &
         + ABS(       qnew(1:mx(1)  ,2:mx(2)+1,1)  &
         -     2.0_r8*qnew(1:mx(1)  ,1:mx(2)  ,1)  &
         +            qnew(1:mx(1)  ,0:mx(2)-1,1)) > 2.0_r8*qtolerance(level))
         errorflags(1:mx(1),1:mx(2)) = 1
      END WHERE
!
   END SUBROUTINE SetErrFlagsUser2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE SetErrFlagsUser3D(qrte,qnew,errorflags,mx,cmx,h,xlower,level)
      USE NodeInfoDef
      USE ProblemDef

      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(IN):: qrte
      REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: qnew
      INTEGER, DIMENSION(1:,1:,1:), INTENT(OUT):: errorflags
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      INTEGER, DIMENSION(1:3), INTENT(IN):: cmx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      INTEGER, INTENT(IN):: level
!
      integer:: i, j, k, mmm, nnn
      REAL(KIND=r8):: pi, r1, r2, r3, r4, tmp, x, y, z, A, h2, inte
!
      errorflags(1:mx(1),1:mx(2),1:mx(3)) = 0

!Error tagging based on the second differences of the phase variable:
      WHERE(  ABS(       qnew(2:mx(1)+1,1:mx(2)  ,1:mx(3)  ,1)  &
         -     2.0_r8*qnew(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)  &
         +            qnew(0:mx(1)-1,1:mx(2)  ,1:mx(3)  ,1)) &
         + ABS(       qnew(1:mx(1)  ,2:mx(2)+1,1:mx(3)  ,1)  &
         -     2.0_r8*qnew(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)  &
         +            qnew(1:mx(1)  ,0:mx(2)-1,1:mx(3)  ,1)) &
         + ABS(       qnew(1:mx(1)  ,1:mx(2)  ,2:mx(3)+1,1)  &
         -     2.0_r8*qnew(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)  &
         +            qnew(1:mx(1)  ,1:mx(2)  ,0:mx(3)-1,1)) &
         > 3.0_r8*qtolerance(level))
         errorflags(1:mx(1),1:mx(2),1:mx(3)) = 1
      END WHERE



!h2 = h*h
!!
!errorflags(1:mx(1),1:mx(2),1:mx(3)) = 0
!!
!!mmm = mx(1)
!!nnn = mx(2)
!!
!DO i = 1, mx(1)
!  x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
!  DO j = 1, mx(2)
!    y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
!	DO k = 1, mx(3)
!		z = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(3)
!!!
!   if      ((level==0).AND.(x>xx11).AND.(x<xx12).AND.(y>yy11).AND.(y<yy12)) then
!!!
!   errorflags(i,j,k) = 1
!!!
!   end if
!!!
!		END DO
!	END DO
!END DO


   END SUBROUTINE SetErrFlagsUser3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User-Supplied AFAS Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 2D Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE RelaxGrid2D(qnew,qold,v1,v1old,v2,v2old,aux,f,mx,h,xlower,ipass)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN OUT):: qnew
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN):: aux
      REAL(KIND=r8), DIMENSION( 1:, 1:,1:), INTENT(IN):: f
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      INTEGER, INTENT(IN):: ipass
!
      INTEGER:: i, j
      REAL(KIND=r8):: det, h2, tmp1, tmp2, tmp3
      REAL(KIND=r8), DIMENSION(1:2,1:2):: a
      REAL(KIND=r8), DIMENSION(1:2):: b
      INTEGER:: ptb_flag
      REAL(KIND=r8):: omega1, omega2, tmp, tmpp, dn1, dn2, dn3, dn4,x,y, m1, m2, m3, m4, mobij
!
      h2 = h*h
!
      tmp1 = 0.5_r8*dt/h2
      tmp2 = epep*epep/h2
      tmp3 = 4.0_r8*epep*epep/h2
!
      DO j = 1, mx(2)
         y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
         DO i = 1+MODULO(j+ipass,2), mx(1), 2
            x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
!
            mobij =   Mob(qold(i,j,1))
            m1 = 0.5_r8*(Mob(qnew(i-1,j,1))+Mob(qnew(i,j,1)))
            m2 = 0.5_r8*(Mob(qnew(i+1,j,1))+Mob(qnew(i,j,1)))
            m3 = 0.5_r8*(Mob(qnew(i,j-1,1))+Mob(qnew(i,j,1)))
            m4 = 0.5_r8*(Mob(qnew(i,j+1,1))+Mob(qnew(i,j,1)))
!
            a(1,1) = 1.0_r8
            a(1,2) = tmp1*(m1+m2+m3+m4)
            a(2,2) = 1.0_r8
            a(2,1) = - fof2(qnew(i,j,1))-tmp3
!
            b(1) = f(i,j,1)+tmp1*(m1*qnew(i-1,j,2)+m2*qnew(i+1,j,2) &
               +                m3*qnew(i,j-1,2)+m4*qnew(i,j+1,2))
!
            b(2) = f(i,j,2)-tmp2*(qnew(i+1,j  ,1)        +qnew(i-1,j  ,1)     &
               +                qnew(i  ,j+1,1)        +qnew(i  ,j-1,1))    &
               +           fof1(qnew(i  ,j  ,1)) - fof2(qnew(i  ,j  ,1))    &
               *                qnew(i  ,j  ,1)
!
! Solve using Cramer's rule:
            det = 1.0_r8-a(1,2)*a(2,1)
!
            qnew(i,j,1) = (b(1)-a(1,2)*b(2))/det
            qnew(i,j,2) = (b(2)-a(2,1)*b(1))/det
!
            mobij =    MobD(qnew(i  ,j  ,1))
            m1 = 0.5_r8*(MobD(qnew(i-1,j  ,1))+mobij)
            m2 = 0.5_r8*(MobD(qnew(i+1,j  ,1))+mobij)
            m3 = 0.5_r8*(MobD(qnew(i  ,j-1,1))+mobij)
            m4 = 0.5_r8*(MobD(qnew(i  ,j+1,1))+mobij)
!
            qnew(i,j,3) = (f(i,j,3)+(m1*qnew(i-1,j,3)+m2*qnew(i+1,j,3)+m3*qnew(i,j-1,3)+m4*qnew(i,j+1,3))/h2) &
               /((qold(i,j,1)+epep)+ (m1+m2+m3+m4)/h2)
!
            qnew(i,j,4) = f(i,j,4)
!
            qnew(i,j,5) = f(i,j,5)
!
         END DO
      END DO
   END SUBROUTINE RelaxGrid2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE UpdateAuxVcycle2D(qnew,qold,v1,v1old,v2,v2old,aux,ll,mx,h,xlower)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN OUT):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
!
!aux(0:mx(1)+1,0:mx(2)+1,1) =
!
   END SUBROUTINE UpdateAuxVcycle2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Operator2D(qnew,qold,v1,v1old,v2,v2old,aux,mx,h,xlower) RESULT(opresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap2D, UDiv2D, Df, ULapVar2D, Gf, AbsGrad2D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN):: aux
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:nccv):: opresult
!
      REAL(KIND=r8):: h2, tmp1, tmp2
!
      h2 = h*h
      tmp1 = dt/h2/2.0_r8
      tmp2 = epep*epep/h2
!
      opresult(1:mx(1),1:mx(2),1) =  qnew(1:mx(1)  ,1:mx(2)  ,1)            &
         - tmp1*ULapVar2D(Mob(qnew(0:mx(1)+1,0:mx(2)+1,1)),          &
         qnew(0:mx(1)+1,0:mx(2)+1,2))
!
      opresult(1:mx(1),1:mx(2),2) =  qnew(1:mx(1)  ,1:mx(2)  ,2)  &
         - Fof(qnew(1:mx(1)  ,1:mx(2)  ,1)) &
         + tmp2*ULap2D(qnew(0:mx(1)+1,0:mx(2)+1,1))

      opresult(1:mx(1)  ,1:mx(2)  ,3) = (qold(1:mx(1)  ,1:mx(2)  ,1) +epep)      &
         * qnew(1:mx(1)  ,1:mx(2)  ,3)       &
         - ULapVar2D(MobD(qnew(0:mx(1)+1,0:mx(2)+1,1)),          &
         qnew(0:mx(1)+1,0:mx(2)+1,3))/h2
!
      opresult(1:mx(1),1:mx(2),4) =  qnew(1:mx(1)  ,1:mx(2)  ,4)
!
      opresult(1:mx(1),1:mx(2),5) =  qnew(1:mx(1)  ,1:mx(2)  ,5)
!
   END FUNCTION Operator2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Source2D(qnew,qold,v1,v1old,v2,v2old,aux,ll,mx,h,xlower) &
      RESULT(sourceresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap2D, AbsGrad, LapVar2D, AbsGrad2D, ULapVar2D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:nccv):: sourceresult
      REAL(KIND=r8):: h2, tmp, tmp1

      h2 = h*h
      tmp1 = dt/2.0_r8
      sourceresult(1:mx(1),1:mx(2),1)  =          qold(1:mx(1)  ,1:mx(2)  ,1) &
         +	tmp1*(ULapVar2D(Mob(qold(0:mx(1)+1,0:mx(2)+1,1)),          &
         qold(0:mx(1)+1,0:mx(2)+1,2))/h2 &
         + pp*qold(1:mx(1)  ,1:mx(2)  ,1) &
         *       qold(1:mx(1)  ,1:mx(2)  ,3) &
         - AA*qold(1:mx(1)  ,1:mx(2)  ,1))

      sourceresult(1:mx(1),1:mx(2),2) =  -epep* Zeta*qold(1:mx(1)  ,1:mx(2)  ,3)/MD

      sourceresult(1:mx(1),1:mx(2),3) = 0.0_r8
!
      sourceresult(1:mx(1),1:mx(2),4)  =         qold(1:mx(1)  ,1:mx(2)  ,4)
!
      sourceresult(1:mx(1),1:mx(2),5) = 0.5_r8* (ULapVar2D(Mob(qold(0:mx(1)+1,0:mx(2)+1,1)),          &
         qold(0:mx(1)+1,0:mx(2)+1,2))/h2         &
         + pp*qold(1:mx(1)  ,1:mx(2)  ,1) &
         *       qold(1:mx(1)  ,1:mx(2)  ,3) &
         - AA*qold(1:mx(1)  ,1:mx(2)  ,1))
!
   END FUNCTION Source2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION SourceUpdate2D(qnew,qold,v1,v1old,v2,v2old,aux,ll,mx,h,xlower) &
      RESULT(updateresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: UDiv2D, ULap2D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:nccv):: updateresult
!
      REAL(KIND=r8):: h2, tmp, tmp1
!
      h2 = h*h
      tmp = dt/h2
      tmp1 = dt /2.0_r8
!!   x < R
!
      updateresult(1:mx(1),1:mx(2),1) = tmp1*(pp*qnew(1:mx(1)  ,1:mx(2)  ,1) &
         *      qnew(1:mx(1)  ,1:mx(2)  ,3) &
         -   AA*qnew(1:mx(1)  ,1:mx(2)  ,1))
!
      updateresult(1:mx(1),1:mx(2),2) = 0.0_r8
!
      updateresult(1:mx(1),1:mx(2),3) = 0.0_r8
!
      updateresult(1:mx(1),1:mx(2),4) = 0.0
!
      updateresult(1:mx(1),1:mx(2),5) =      0.5_r8*(pp*qnew(1:mx(1)  ,1:mx(2)  ,1) &
         *      qnew(1:mx(1)  ,1:mx(2)  ,3) &
         -   AA*qnew(1:mx(1)  ,1:mx(2)  ,1))
!
   END FUNCTION SourceUpdate2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 2D Edge Problem Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE RelaxGrid2DEdge(qnew,qold,v1,v1old,v2,v2old,aux,f,f1,f2,level,mx,h,xlower,ipass)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN OUT):: qnew
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN)    :: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN OUT):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN)    :: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN OUT):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN)    :: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN)    :: aux
      REAL(KIND=r8), DIMENSION( 1:, 1:,1:), INTENT(IN)    :: f
      REAL(KIND=r8), DIMENSION( 0:, 1:,1:), INTENT(IN)    :: f1
      REAL(KIND=r8), DIMENSION( 1:, 0:,1:), INTENT(IN)    :: f2
      INTEGER, INTENT(IN):: level
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      INTEGER, INTENT(IN):: ipass
!
      INTEGER:: i, j
      REAL(KIND=r8):: h2, a11, a12, a13, a14, a15
      REAL(KIND=r8), DIMENSION(1:5):: b, c
!
      h2 = h*h

!
   END SUBROUTINE RelaxGrid2DEdge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Operator2DV1(qnew,qold,v1,v1old,v2,v2old,aux,mx,h,xlower) RESULT(opresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap2DV1
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN):: aux
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(0:mx(1),1:mx(2),1:nfcv):: opresult
!
      REAL(KIND=r8):: h2
!
      h2 = h*h
!
   END FUNCTION Operator2DV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Source2DV1(qnew,qold,v1,v1old,v2,v2old,aux,ll,mx,h,xlower) RESULT(sourceresult)
      USE NodeInfoDef
      USE GridUtilities
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(0:mx(1),1:mx(2),1:nfcv):: sourceresult
!
      REAL(KIND=r8):: tmp
!
   END FUNCTION Source2DV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION SourceUpdate2DV1(qnew,qold,v1,v1old,v2,v2old,aux,ll,mx,h,xlower) &
      RESULT(updateresult)
      USE NodeInfoDef
      USE GridUtilities
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(0:mx(1),1:mx(2),1:nfcv):: updateresult
!
      REAL(kind=r8):: tmp
!

   END FUNCTION SourceUpdate2DV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Operator2DV2(qnew,qold,v1,v1old,v2,v2old,aux,mx,h,xlower) RESULT(opresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap2DV2
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,1:), INTENT(IN):: aux
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),0:mx(2),1:nfcv):: opresult
!
      REAL(KIND=r8):: h2
!
      h2 = h*h
!
   END FUNCTION Operator2DV2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Source2DV2(qnew,qold,v1,v1old,v2,v2old,aux,ll,mx,h,xlower) RESULT(sourceresult)
      USE NodeInfoDef
      USE GridUtilities
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),0:mx(2),1:nfcv):: sourceresult
!
      REAL(KIND=r8):: tmp
!
   END FUNCTION Source2DV2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION SourceUpdate2DV2(qnew,qold,v1,v1old,v2,v2old,aux,ll,mx,h,xlower) &
      RESULT(updateresult)
      USE NodeInfoDef
      USE GridUtilities
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),0:mx(2),1:nfcv):: updateresult
!
      REAL(KIND=r8):: tmp
!
      tmp = dt*gamma/(2.0_r8*h)
!

   END FUNCTION SourceUpdate2DV2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION UDivPhi2D(p,v1,v2) RESULT(udivresult)
      USE NodeInfoDef
      USE ProblemDef
      USE GridUtilities, ONLY: UDiv2D
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: p
      REAL(KIND=r8), DIMENSION(0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(1:,0:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION(1:SIZE(p,1)-2,1:SIZE(p,2)-2):: udivresult
!
      INTEGER, DIMENSION(1:2):: mx
      REAL(KIND=r8), DIMENSION(0:SIZE(p,1)-2,1:SIZE(p,2)-2):: f1
      REAL(KIND=r8), DIMENSION(1:SIZE(p,1)-2,0:SIZE(p,2)-2):: f2
!
      mx(1) = SIZE(p,1)-2; mx(2) = SIZE(p,2)-2
!
! Calculate the UNDIVIDED 2D flux function:
      f1(0:mx(1),1:mx(2)) &
         = 0.5_r8*(p(1:mx(1)+1,1:mx(2))+p(0:mx(1),1:mx(2)))*v1(0:mx(1),1:mx(2))
!
      f2(1:mx(1),0:mx(2)) &
         = 0.5_r8*(p(1:mx(1),1:mx(2)+1)+p(1:mx(1),0:mx(2)))*v2(1:mx(1),0:mx(2))
!
! Calculate the UNDIVIDED divergence of the flux:
      udivresult(1:mx(1),1:mx(2)) = UDiv2D(f1(0:mx(1),1:mx(2)), &
         f2(1:mx(1),0:mx(2)))
!
   END FUNCTION UDivPhi2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 3D Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE RelaxGrid3D(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,f,mx,h,xlower,ipass)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN OUT):: qnew
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION(0:,0:,0:,1:), INTENT(IN):: aux
      REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(IN):: f
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      INTEGER, INTENT(IN):: ipass
!
      INTEGER:: i, j, k
      REAL(KIND=r8):: d1, d2, d3, d4, d5, d6, det, h2, tmp1, tmp2, tmp3, tmp4
      REAL(KIND=r8):: m1, m2, m3, m4, m5, m6, mobij
      REAL(KIND=r8), DIMENSION(1:2,1:2):: a
      REAL(KIND=r8), DIMENSION(1:2):: b
!
      h2 = h*h
!
      tmp1 = 0.5_r8*dt/h2
      tmp2 = epep2/h2
      tmp3 = 6.0_r8*epep2/h2
!
      a(1,1) = 1.0_r8
      a(2,2) = 1.0_r8
!
! These loops must always be structured this way, i.e., so that i,j,k = 1 is a
! black cell. (The first cell relaxed is always i = 2, j,k = 1, since ipass = 1
! is called first.)
!
      DO k = 1, mx(3)
         DO j = 1, mx(2)
            DO i = 1+MODULO(k+j+ipass,2), mx(1), 2
!
!
               mobij =   Mob(qnew(i,j,k,1))
               m1 = 0.5_r8*(Mob(qnew(i-1,j,k,1)) + mobij)
               m2 = 0.5_r8*(Mob(qnew(i+1,j,k,1)) + mobij)
               m3 = 0.5_r8*(Mob(qnew(i,j-1,k,1)) + mobij)
               m4 = 0.5_r8*(Mob(qnew(i,j+1,k,1)) + mobij)
               m5 = 0.5_r8*(Mob(qnew(i,j,k-1,1)) + mobij)
               m6 = 0.5_r8*(Mob(qnew(i,j,k+1,1)) + mobij)
!
               a(1,2) = tmp1*(m1+m2+m3+m4+m5+m6)
!
               a(2,1) = -fof2(qnew(i,j,k,1))-tmp3
!
               b(1) = f(i,j,k,1)+tmp1*(m1*qnew(i-1,j,k,2)+m2*qnew(i+1,j,k,2) &
                  +                  m3*qnew(i,j-1,k,2)+m4*qnew(i,j+1,k,2) &
                  +                  m5*qnew(i,j,k-1,2)+m6*qnew(i,j,k+1,2))
!
               b(2) = f(i,j,k,2)- tmp2*(qnew(i+1,j,k,1)+qnew(i-1,j,k,1) &
                  +                   qnew(i,j+1,k,1)+qnew(i,j-1,k,1) &
                  +                   qnew(i,j,k+1,1)+qnew(i,j,k-1,1))&
                  +           		fof1(qnew(i  ,j ,k ,1)) - fof2(qnew(i  ,j,k  ,1))&
                  *                qnew(i  ,j ,k ,1)
!
! Solve for Cahn-Hilliard using Cramer's rule:
               det = 1.0_r8-a(1,2)*a(2,1)
!
               qnew(i,j,k,1) = (b(1)-a(1,2)*b(2))/det
               qnew(i,j,k,2) = (b(2)-a(2,1)*b(1))/det
!
               mobij =      MobD(qnew(i  ,j  ,k  ,1))
               m1 = 0.5_r8*(MobD(qnew(i-1,j  ,k  ,1)) + mobij)
               m2 = 0.5_r8*(MobD(qnew(i+1,j  ,k  ,1)) + mobij)
               m3 = 0.5_r8*(MobD(qnew(i  ,j-1,k  ,1)) + mobij)
               m4 = 0.5_r8*(MobD(qnew(i  ,j+1,k  ,1)) + mobij)
               m5 = 0.5_r8*(MobD(qnew(i  ,j  ,k-1,1)) + mobij)
               m6 = 0.5_r8*(MobD(qnew(i  ,j  ,k+1,1)) + mobij)
!
               qnew(i,j,k,3) = (f(i,j,k,3)          &
                  +(m1*qnew(i-1,j,k,3)   &
                  + m2*qnew(i+1,j,k,3)   &
                  + m3*qnew(i,j-1,k,3)   &
                  + m4*qnew(i,j+1,k,3)   &
                  + m5*qnew(i,j,k-1,3)   &
                  + m6*qnew(i,j,k+1,3))/h2)/(qold(i,j,k,1)+ (m1+m2+m3+m4+m5+m6)/h2)
!
            END DO
         END DO
      END DO
!
   END SUBROUTINE RelaxGrid3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE UpdateAuxVcycle3D(qnew,qold,aux,ll,mx,h,xlower)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN OUT):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
!
!aux(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1) =
!
   END SUBROUTINE UpdateAuxVcycle3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE RelaxGrid3DEdge(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,f,f1,f2,f3,level,mx,h,xlower,ipass)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN OUT):: qnew
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN    ):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN OUT):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN    ):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN OUT):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN    ):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN OUT):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN    ):: v3old
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN    ):: aux
      REAL(KIND=r8), DIMENSION( 1:, 1:, 1:,1:), INTENT(IN    ):: f
      REAL(KIND=r8), DIMENSION( 0:, 1:, 1:,1:), INTENT(IN    ):: f1
      REAL(KIND=r8), DIMENSION( 1:, 0:, 1:,1:), INTENT(IN    ):: f2
      REAL(KIND=r8), DIMENSION( 1:, 1:, 0:,1:), INTENT(IN    ):: f3
      INTEGER, INTENT(IN):: level
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      INTEGER, INTENT(IN):: ipass
!
      INTEGER:: i, j, k
      REAL(KIND=r8):: d1, d2, d3, d4, d5, d6, det, h2, tmp1, tmp2, tmp3, tmp4
      REAL(KIND=r8), DIMENSION(1:2,1:2):: a
      REAL(KIND=r8):: a11, a12, a13, a14, a15
      REAL(KIND=r8), DIMENSION(1:7):: b, c
!
      h2 = h*h
!
      a11  = 1.0_r8 + 6.0_r8 * mu * dt/h2
      a12  = dt/h
      a13  = dt*mu/h2
!
      tmp1 = dt/h2
      tmp2 = epep2/h2
      tmp3 = 6.0_r8*epep2/h2
!
      a(1,1) = 1.0_r8
      a(1,2) = 6.0_r8*tmp1
      a(2,2) = 1.0_r8
!
! These loops must always be structured this way, i.e., so that i,j,k = 1 is a
! black cell. (The first cell relaxed is always i = 2, j,k = 1, since ipass = 1
! is called first.)
!
      DO k = 1, mx(3)
         DO j = 1, mx(2)
            DO i = 1+MODULO(j+ipass,2), mx(1), 2
!
               b(1) =      (f1(i  ,j  ,k  ,1)      &
                  + a13* (v1(i-1,j  ,k  ,1)      &
                  +       v1(i+1,j  ,k  ,1)      &
                  +       v1(i  ,j-1,k  ,1)      &
                  +       v1(i  ,j+1,k  ,1)      &
                  +       v1(i  ,j  ,k-1,1)      &
                  +       v1(i  ,j  ,k+1,1))     &
               !- a12*qnew(i+1,j  ,k  ,3)
                  )/a11
!
               b(2) =      (f1(i-1,j  ,k  ,1)      &
                  + a13* (v1(i-2,j  ,k  ,1)      &
                  +       v1(i  ,j  ,k  ,1)      &
                  +       v1(i-1,j-1,k  ,1)      &
                  +       v1(i-1,j+1,k  ,1)      &
                  +       v1(i-1,j  ,k-1,1)      &
                  +       v1(i-1,j  ,k+1,1))     &
               !+ a12*qnew(i-1,j  ,k  ,3)
                  )/a11
!
               b(3) =      (f2(i  ,j  ,k  ,1)      &
                  + a13* (v2(i-1,j  ,k  ,1)      &
                  +       v2(i+1,j  ,k  ,1)      &
                  +       v2(i  ,j-1,k  ,1)      &
                  +       v2(i  ,j+1,k  ,1)      &
                  +       v2(i  ,j  ,k-1,1)      &
                  +       v2(i  ,j  ,k+1,1))     &
!			     !- a12*qnew(i  ,j+1,k  ,3)
                  )/a11
!
               b(4) =      (f2(i  ,j-1,k  ,1)      &
                  + a13* (v2(i-1,j-1,k  ,1)      &
                  +       v2(i+1,j-1,k  ,1)      &
                  +       v2(i  ,j-2,k  ,1)      &
                  +       v2(i  ,j  ,k  ,1)      &
                  +       v2(i  ,j-1,k-1,1)      &
                  +       v2(i  ,j-1,k+1,1))     &
!			     !+ a12*qnew(i  ,j-1,k  ,3)
                  )/a11
!
               b(5) =      (f3(i  ,j  ,k  ,1)      &
                  + a13* (v3(i-1,j  ,k  ,1)      &
                  +       v3(i+1,j  ,k  ,1)      &
                  +       v3(i  ,j-1,k  ,1)      &
                  +       v3(i  ,j+1,k  ,1)      &
                  +       v3(i  ,j  ,k-1,1)      &
                  +       v3(i  ,j  ,k+1,1))     &
!			     !- a12*qnew(i  ,j  ,k+1,3)
                  )/a11
!
               b(6) =      (f3(i  ,j  ,k-1,1)      &
                  + a13* (v3(i-1,j  ,k-1,1)      &
                  +       v3(i+1,j  ,k-1,1)      &
                  +       v3(i  ,j-1,k-1,1)      &
                  +       v3(i  ,j+1,k-1,1)      &
                  +       v3(i  ,j  ,k-2,1)      &
                  +       v3(i  ,j  ,k  ,1))     &
!			     !+ a12*qnew(i  ,j  ,k-2,3)
                  )/a11
!
!			b(7) =        f(i  ,j  ,k  ,3)
!
!			c(1) = 0.0!(-b(1)+b(2)-b(3)+b(4)-b(5)+b(6)-b(7))*(a11/(6.0_r8*a12))
!
               c(2) = (b(1)+(a12/a11)*c(1))
               c(3) = (b(2)-(a12/a11)*c(1))
               c(4) = (b(3)+(a12/a11)*c(1))
               c(5) = (b(4)-(a12/a11)*c(1))
               c(6) = (b(5)+(a12/a11)*c(1))
               c(7) = (b(6)-(a12/a11)*c(1))
!
!			qnew(i,j,k,3) = qnew(i,j,k,3)+omega*(c(1)-qnew(i,j,k,3))
!
               v1(i  ,j  ,k  ,1) = f1(i  ,j  ,k  ,1)!b(1)
               v1(i-1,j  ,k  ,1) = f1(i-1,j  ,k  ,1)!b(2)
               v2(i  ,j  ,k  ,1) = f2(i  ,j  ,k  ,1)!b(3)
               v2(i  ,j-1,k  ,1) = f2(i  ,j-1,k  ,1)!b(4)
               v3(i  ,j  ,k  ,1) = f3(i  ,j  ,k  ,1)!b(5)
               v3(i  ,j  ,k-1,1) = f3(i  ,j  ,k-1,1)!b(6)
!
               IF(level<=rootlevel) THEN
                  IF(i==1) THEN
                     v1(mx(1)  ,j      ,k,1) = v1(0        ,j     ,k,1)
                     v1(mx(1)+1,j      ,k,1) = v1(1        ,j     ,k,1)
                  ELSE IF(i==mx(1)) THEN
                     v1( 0     ,j      ,k,1) =  v1(mx(1)   ,j     ,k,1)
                     v1(-1     ,j      ,k,1) =  v1(mx(1)-1 ,j     ,k,1)
                  END IF
!
                  IF(j==1) THEN
                     v2( i     ,mx(2)  ,k,1) =  0.0_r8
                     v2( i     ,mx(2)+1,k,1) = -v2(i      ,mx(2)-1,k,1)
                  ELSE IF(j==mx(2)) THEN
                     v2( i     , 0     ,k,1) =  0.0_r8
                     v2( i     ,-1     ,k,1) = -v2(i      ,1      ,k,1)
                  END IF
!
                  IF(k==1) THEN
                     v3( i     ,j      ,mx(3)  ,1) =  0.0_r8
                     v3( i     ,j      ,mx(3)+1,1) = -v3(i      ,j,mx(3)-1,1)
                  ELSE IF(k==mx(3)) THEN
                     v3( i     ,j      , 0     ,1) =  0.0_r8
                     v3( i     ,j      ,-1     ,1) = -v3(i      ,j,1      ,1)
                  END IF
               END IF
!
            END DO
         END DO
      END DO
!
   END SUBROUTINE RelaxGrid3DEdge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Operator3D(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,mx,h,xlower) RESULT(opresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D, ULapVar3D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: aux
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:mx(3),1:nccv):: opresult
!
      REAL(KIND=r8):: h2, tmp1, tmp2
!
      h2 = h*h
      tmp1 = 0.5_r8*dt/h2
      tmp2 = epep2/h2
!
      opresult(1:mx(1),1:mx(2),1:mx(3),1) &
         = qnew(1:mx(1),1:mx(2),1:mx(3),1) &
         - tmp1*ULapVar3D(Mob(qnew(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1)),   &
         qnew(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,2))
!
      opresult(1:mx(1),1:mx(2),1:mx(3),2) &
         =             qnew(1:mx(1),1:mx(2),1:mx(3),2) &
         -        fof1(qnew(1:mx(1),1:mx(2),1:mx(3),1)) &
         + tmp2*ULap3D(qnew(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1))
!
      opresult(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,3)      &
         =                qold(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)*  	&
         qnew(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,3)      &
         - ULapVar3D(MobD(qnew(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1)),    &
         qnew(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,3)) /h2
!
   END FUNCTION Operator3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Source3D(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,ll,mx,h,xlower) RESULT(sourceresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D, ULapVar3D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:mx(3),1:nccv):: sourceresult
!
      REAL(KIND=r8):: h2, tmp, tmp1
!
      h2 = h*h
      tmp = dt/h2
      tmp1 = 0.5_r8*dt
!
      sourceresult(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1) =   &
         qold(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)     &
         +	tmp1*(ULapVar3D(Mob(qold(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1)),   &
         qold(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,2))/h2 &
         + 					 pp*qold(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)     &
         *    				    qold(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,3)     &
         - 					 AA*qold(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1))
!
      sourceresult(1:mx(1),1:mx(2),1:mx(3),2) = -epep* Zeta*qold(1:mx(1)  ,1:mx(2) ,1:mx(3) ,3)
!
      sourceresult(1:mx(1),1:mx(2),1:mx(3),3) = 0.0_r8
!
   END FUNCTION Source3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION SourceUpdate3D(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,ll,mx,h,xlower) RESULT(updateresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: UDiv3D, ULap3D, UDivPhi3D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),1:mx(3),1:nccv):: updateresult
!
      REAL(KIND=r8):: tmp
!
      tmp = 0.5_r8*dt
!
      updateresult(1:mx(1),1:mx(2),1:mx(3),1) =  tmp*( pp*qnew(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)     &
         *    				    qnew(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,3)     &
         - 					 AA*qnew(1:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1))
!
!
      updateresult(1:mx(1),1:mx(2),1:mx(3),2) = 0.0_r8
!
      updateresult(1:mx(1),1:mx(2),1:mx(3),3) = 0.0_r8
!
   END FUNCTION SourceUpdate3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Operator3DV1(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,mx,h,xlower) RESULT(opresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D, ULap3DV1
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: aux
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(0:mx(1),1:mx(2),1:mx(3),1:nfcv):: opresult
!
      REAL(KIND=r8):: h2, tmp1, tmp2
!
      h2 = h*h
      tmp1 = dt/h2
      tmp2 = epep2/h2
!
      opresult( 0:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)        &
         = v1( 0:mx(1)  ,1:mx(2)  ,1:mx(3)  ,1)        !&
      !-mu*ULap3DV1(v1(-1:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1))*dt/h2 !&
      !+ (qnew( 1:mx(1)+1,1:mx(2)  ,1:mx(3)  ,3)        &
      !-  qnew( 0:mx(1)  ,1:mx(2)  ,1:mx(3)  ,3))*dt/h
!
   END FUNCTION Operator3DV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Source3DV1(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,ll,mx,h,xlower) RESULT(sourceresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(0:mx(1),1:mx(2),1:mx(3),1:nfcv):: sourceresult
!
      REAL(KIND=r8):: h2, tmp
!
      h2 = h*h
      tmp = dt/h2
!
      sourceresult(0:mx(1),1:mx(2),1:mx(3),1) &
         = v1old(0:mx(1),1:mx(2),1:mx(3),1)
!
   END FUNCTION Source3DV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   FUNCTION SourceUpdate3DV1(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,ll,mx,h,xlower) RESULT(updateresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(0:mx(1),1:mx(2),1:mx(3),1:nfcv):: updateresult
!
      REAL(KIND=r8):: h2, tmp
!
      tmp = dt*gamma/(2.0_r8*h)
!
      updateresult(0:mx(1)  ,1:mx(2),1:mx(3),1)   &
         = 0.0!-tmp*(qold(1:mx(1)+1,1:mx(2),1:mx(3),1)   &
      !+ qold(0:mx(1)  ,1:mx(2),1:mx(3),1))  &
      !*(qnew(1:mx(1)+1,1:mx(2),1:mx(3),2)   &
      !- qnew(0:mx(1)  ,1:mx(2),1:mx(3),2))
!
   END FUNCTION SourceUpdate3DV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Operator3DV2(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,mx,h,xlower) RESULT(opresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D, ULap3DV2
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: aux
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),0:mx(2),1:mx(3),1:nfcv):: opresult
!
      REAL(KIND=r8):: h2, tmp1, tmp2
!
      h2 = h*h
      tmp1 = dt/h2
      tmp2 = epep2/h2
!
      opresult( 1:mx(1)  , 0:mx(2)  ,1:mx(3)  ,1)        &
         =   v2( 1:mx(1)  , 0:mx(2)  ,1:mx(3)  ,1)        !&
! -mu*ULap3DV2(v2( 0:mx(1)+1,-1:mx(2)+1,0:mx(3)+1,1))*dt/h2
!
   END FUNCTION Operator3DV2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Source3DV2(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,ll,mx,h,xlower) RESULT(sourceresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),0:mx(2),1:mx(3),1:nfcv):: sourceresult
!
      REAL(KIND=r8):: h2, tmp
!
      h2 = h*h
      tmp = dt/h2
!
      sourceresult(1:mx(1),0:mx(2),1:mx(3),1) &
         = v2old(1:mx(1),0:mx(2),1:mx(3),1)
!
   END FUNCTION Source3DV2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   FUNCTION SourceUpdate3DV2(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,ll,mx,h,xlower) RESULT(updateresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),0:mx(2),1:mx(3),1:nfcv):: updateresult
!
      REAL(KIND=r8):: h2, tmp
!
      h2 = h*h
      tmp = dt/h2
!
      updateresult(1:mx(1),0:mx(2),1:mx(3),1) &
         = 0.0!-tmp*(qold(1:mx(1),1:mx(2)+1,1:mx(3),1)   &
      !+ qold(1:mx(1),0:mx(2)  ,1:mx(3),1))  &
      ! *(qnew(1:mx(1),1:mx(2)+1,1:mx(3),2)   &
      !- qnew(1:mx(1),0:mx(2)  ,1:mx(3),2))
!
   END FUNCTION SourceUpdate3DV2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Operator3DV3(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,mx,h,xlower) RESULT(opresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D, ULap3DV3
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION( 0:, 0:, 0:,1:), INTENT(IN):: aux
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),0:mx(3),1:nfcv):: opresult
!
      REAL(KIND=r8):: h2, tmp1, tmp2
!
      h2 = h*h
      tmp1 = dt/h2
      tmp2 = epep2/h2
!
      opresult(1:mx(1)  , 1:mx(2)  , 0:mx(3)  ,1)        &
         =   v3(1:mx(1)  , 1:mx(2)  , 0:mx(3)  ,1)        !&
      !-mu*ULap3DV3(v3(0:mx(1)+1, 0:mx(2)+1,-1:mx(3)+1,1))*dt/h2
!
   END FUNCTION Operator3DV3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION Source3DV3(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,ll,mx,h,xlower) RESULT(sourceresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),0:mx(3),1:nfcv):: sourceresult
!
      REAL(KIND=r8):: h2, tmp
!
      h2 = h*h
      tmp = dt/h2
!
      sourceresult(1:mx(1),1:mx(2),0:mx(3),1) &
         = v3old(1:mx(1),1:mx(2),0:mx(3),1)
!
   END FUNCTION Source3DV3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   FUNCTION SourceUpdate3DV3(qnew,qold,v1,v1old,v2,v2old,v3,v3old,aux,ll,mx,h,xlower) RESULT(updateresult)
      USE NodeInfoDef
      USE GridUtilities, ONLY: ULap3D
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qnew
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: qold
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN):: v1old
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN):: v2old
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN):: v3old
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN):: aux
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      REAL(KIND=r8), INTENT(IN):: h
      REAL(KIND=r8), DIMENSION(1:3), INTENT(IN):: xlower
      REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2),0:mx(3),1:nfcv):: updateresult
!
      REAL(KIND=r8):: h2, tmp
!
      h2 = h*h
      tmp = dt/h2
!
      updateresult(1:mx(1),1:mx(2),0:mx(3),1)         &
         = 0.0!-tmp*(qold(1:mx(1)  ,1:mx(2),1:mx(3)+1,1)   &
      ! + qold(1:mx(1)  ,1:mx(2),0:mx(3)  ,1))  &
      ! *(qnew(1:mx(1)  ,1:mx(2),1:mx(3)+1,2)   &
      ! - qnew(1:mx(1)  ,1:mx(2),0:mx(3)  ,2))
!
   END FUNCTION SourceUpdate3DV3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User supplied physical boundary conditions:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE UserBC2D(q,v1,v2,ll,mx,bcnow)
      USE NodeinfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,1:), INTENT(IN OUT):: q
      REAL(KIND=r8), DIMENSION(-1:, 0:,1:), INTENT(IN OUT):: v1
      REAL(KIND=r8), DIMENSION( 0:,-1:,1:), INTENT(IN OUT):: v2
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      INTEGER, INTENT(IN):: bcnow
!
      INTEGER:: ibc
      INTEGER, DIMENSION(1:2):: ul
!
      ul(1:2) = mx(1:2)+mbc
!
      SELECT CASE(bcnow)
       CASE(1)
!
! Bndry 1:
         DO ibc = 1, mbc
            q(1-ibc,ll:ul(2),1:2) =     q(ibc,ll:ul(2),1:2)
            q(1-ibc,ll:ul(2),3)   = 2.0 - q(ibc,ll:ul(2),3)
         END DO
!
         v1(-1,0:mx(2)+1,1) = -v1(1,0:mx(2)+1,1)
         v1( 0,0:mx(2)+1,1) =  0.0_r8
!
         v2(0,-1:mx(2)+1,1) = v2(1,-1:mx(2)+1,1)
!
       CASE(2)
!
! Bndry 2:
         DO ibc = 1, mbc
            q(mx(1)+ibc,ll:ul(2),1) =       q(mx(1)-ibc+1,ll:ul(2),1)
            q(mx(1)+ibc,ll:ul(2),2) =       q(mx(1)-ibc+1,ll:ul(2),2)
            q(mx(1)+ibc,ll:ul(2),3)   =2.0 -  q(mx(1)-ibc+1,ll:ul(2),3)
         END DO
!
         v1(mx(1)+1,0:mx(2)+1,1) = -v1(mx(1)-1,0:mx(2)+1,1)
         v1(mx(1)  ,0:mx(2)+1,1) =  0.0_r8
!
         v2(mx(1)+1,-1:mx(2)+1,1) = v2(mx(1),-1:mx(2)+1,1)
!
       CASE(3)
!
! Bndry 3:
         DO ibc = 1, mbc
            q(ll:ul(1),1-ibc,1:2) =       q(ll:ul(1),ibc,1:2)
            q(ll:ul(1),1-ibc,3)   =2.0 -  q(ll:ul(1),ibc,3)
         END DO
!
         v1(-1:mx(1)+1,0,1) = 2.0_r8*shearvel-v1(-1:mx(1)+1,1,1)
!
         v2(0:mx(1)+1, 0,1) = 0.0_r8
         v2(0:mx(1)+1,-1,1) = -v2(0:mx(1)+1,1,1)
!
       CASE(4)
!
! Bndry 4:
         DO ibc = 1, mbc
            q(ll:ul(1),mx(2)+ibc,1:2) =       q(ll:ul(1),mx(2)-ibc+1,1:2)
            q(ll:ul(1),mx(2)+ibc,3)   =2.0 -  q(ll:ul(1),mx(2)-ibc+1,3)
         END DO
!
         v1(-1:mx(1)+1,mx(2)+1,1) = -2.0_r8*shearvel-v1(-1:mx(1)+1,mx(2),1)
!
         v2(0:mx(1)+1,mx(2)  ,1) = 0.0_r8
         v2(0:mx(1)+1,mx(2)+1,1) = -v2(0:mx(1)+1,mx(2)-1,1)
      END SELECT
!
   END SUBROUTINE UserBC2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE UserBC3D(q,v1,v2,v3,ll,mx,bcnow)
      USE NodeinfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(ll:,ll:,ll:,1:), INTENT(IN OUT):: q
      REAL(KIND=r8), DIMENSION(-1:, 0:, 0:,1:), INTENT(IN OUT):: v1
      REAL(KIND=r8), DIMENSION( 0:,-1:, 0:,1:), INTENT(IN OUT):: v2
      REAL(KIND=r8), DIMENSION( 0:, 0:,-1:,1:), INTENT(IN OUT):: v3
      INTEGER, INTENT(IN):: ll
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      INTEGER, INTENT(IN):: bcnow
!
      INTEGER:: ibc
      INTEGER, DIMENSION(1:3):: ul
!
      ul(1:3) = mx(1:3)+mbc
!
      SELECT CASE(bcnow)
       CASE(1)
!
! Dirichlet on bndry 1:
         DO ibc = 1, mbc
            q(1-ibc,ll:ul(2),ll:ul(3),1:2) =    q(ibc,ll:ul(2),ll:ul(3),1:2)
            q(1-ibc,ll:ul(2),ll:ul(3),3) = 2.0_r8-   q(ibc,ll:ul(2),ll:ul(3),3)
         END DO
!
         v1(-1, 0:mx(2)+1, 0:mx(3)+1,1:nfcv) = -v1( 1, 0:mx(2)+1, 0:mx(3)+1,1:nfcv)
         v1( 0, 0:mx(2)+1, 0:mx(3)+1,1:nfcv) =  0.0_r8
!
         v2( 0,-1:mx(2)+1, 0:mx(3)+1,1:nfcv) =  v2( 1,-1:mx(2)+1, 0:mx(3)+1,1:nfcv)
!
         v3( 0, 0:mx(2)+1,-1:mx(3)+1,1:nfcv) =  v3( 1, 0:mx(2)+1,-1:mx(3)+1,1:nfcv)
!
       CASE(2)
!
! Dirichlet on bndry 2:
         DO ibc = 1, mbc
            q(mx(1)+ibc,ll:ul(2),ll:ul(3),1:2)   =   q(mx(1)-ibc+1,ll:ul(2),ll:ul(3),1:2)
            q(mx(1)+ibc,ll:ul(2),ll:ul(3),3)   =   2.0_r8- q(mx(1)-ibc+1,ll:ul(2),ll:ul(3),3)
         END DO
!
         v1(mx(1)+1, 0:mx(2)+1, 0:mx(3)+1,1:nfcv) = -v1(mx(1)-1, 0:mx(2)+1, 0:mx(3)+1,1:nfcv)
         v1(mx(1)  , 0:mx(2)+1, 0:mx(3)+1,1:nfcv) =  0.0_r8
!
         v2(mx(1)+1,-1:mx(2)+1, 0:mx(3)+1,1:nfcv) =  v2(mx(1)  ,-1:mx(2)+1, 0:mx(3)+1,1:nfcv)
!
         v3(mx(1)+1, 0:mx(2)+1,-1:mx(3)+1,1:nfcv) =  v3(mx(1)  , 0:mx(2)+1,-1:mx(3)+1,1:nfcv)
!
       CASE(3)
!
! Dirichlet on bndry 3:
         DO ibc = 1, mbc
            q(ll:ul(1),1-ibc,ll:ul(3),1:2)  =    q (ll:ul(1)  ,ibc,ll:ul(3)  ,1:2)
            q(ll:ul(1),1-ibc,ll:ul(3),3)  =   2.0_r8- q (ll:ul(1)  ,ibc,ll:ul(3)  ,3)
         END DO
!
         v1(-1:mx(1)+1, 0, 0:mx(3)+1,1:nfcv) =  2.0_r8*shearvel &
            - v1(-1:mx(1)+1,  1, 0:mx(3)+1,1:nfcv)
!
         v2( 0:mx(1)+1, 0, 0:mx(3)+1,1:nfcv) =  0.0_r8
         v2( 0:mx(1)+1,-1, 0:mx(3)+1,1:nfcv) =  - v2( 0:mx(1)+1,  1, 0:mx(3)+1,1:nfcv)
!
         v3( 0:mx(1)+1, 0,-1:mx(3)+1,1:nfcv) =    v3( 0:mx(1)+1,  1,-1:mx(3)+1,1:nfcv)
!
       CASE(4)
!
! Dirichlet on bndry 4:
         DO ibc = 1, mbc
            q(ll:ul(1),mx(2)+ibc,ll:ul(3),1:2) = q (ll:ul(1),mx(2)-ibc+1,ll:ul(3),1:2)
            q(ll:ul(1),mx(2)+ibc,ll:ul(3),3) = 2.0_r8 - q (ll:ul(1),mx(2)-ibc+1,ll:ul(3),3)
         END DO
!
         v1(-1:mx(1)+1,mx(2)+1, 0:mx(3)+1,1:nfcv) = -2.0_r8*shearvel &
            -v1(-1:mx(1)+1,mx(2), 0:mx(3)+1,1:nfcv)
!
         v2( 0:mx(1)+1,mx(2)  , 0:mx(3)+1,1:nfcv) = 0.0_r8
         v2( 0:mx(1)+1,mx(2)+1, 0:mx(3)+1,1:nfcv) = -v2( 0:mx(1)+1,mx(2)-1, 0:mx(3)+1,1:nfcv)
!
         v3( 0:mx(1)+1,mx(2)+1,-1:mx(3)+1,1:nfcv) =  v3( 0:mx(1)+1,mx(2)  ,-1:mx(3)+1,1:nfcv)
!
       CASE(5)
!
! Dirichlet on bndry 5:
         DO ibc = 1, mbc
            q(ll:ul(1),ll:ul(2),1-ibc,1:2) =  q(ll:ul(1),ll:ul(2),ibc,1:2)
            q(ll:ul(1),ll:ul(2),1-ibc,3) = 2.0_r8- q(ll:ul(1),ll:ul(2),ibc,3)
         END DO
!
         v1(-1:mx(1)+1, 0:mx(2)+1, 0,1:nfcv) =   2.0_r8*shearvel &
            -   v1(-1:mx(1)+1, 0:mx(2)+1, 1,1:nfcv)
!
         v2( 0:mx(1)+1,-1:mx(2)+1, 0,1:nfcv) =   v2( 0:mx(1)+1,-1:mx(2)+1, 1,1:nfcv)
!
         v3( 0:mx(1)+1, 0:mx(2)+1, 0,1:nfcv) =   0.0
         v3( 0:mx(1)+1, 0:mx(2)+1,-1,1:nfcv) =  -v3( 0:mx(1)+1, 0:mx(2)+1, 1,1:nfcv)
!
       CASE(6)
!
! Dirichlet on bndry 6:
         DO ibc = 1, mbc
            q(ll:ul(1),ll:ul(2),mx(3)+ibc,1:2) =  q(ll:ul(1),ll:ul(2),mx(3)-ibc+1,1:2)
            q(ll:ul(1),ll:ul(2),mx(3)+ibc,3) = 2.0_r8-  q(ll:ul(1),ll:ul(2),mx(3)-ibc+1,3)
         END DO
!
         v1(-1:mx(1)+1, 0:mx(2)+1,mx(3)+1,1:nfcv) =   -2.0_r8*shearvel &
            -   v1(-1:mx(1)+1, 0:mx(2)+1,mx(3)  ,1:nfcv)
!
         v2( 0:mx(1)+1,-1:mx(2)+1,mx(3)+1,1:nfcv) =   v2( 0:mx(1)+1,-1:mx(2)+1,mx(3)  ,1:nfcv)
!
         v3( 0:mx(1)+1, 0:mx(2)+1,mx(3)  ,1:nfcv) =   0.0
         v3( 0:mx(1)+1, 0:mx(2)+1,mx(3)+1,1:nfcv) =  -v3( 0:mx(1)+1, 0:mx(2)+1,mx(3)-1,1:nfcv)
!
      END SELECT
!
   END SUBROUTINE UserBC3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Dimensionally Invariant Multigrid Routines:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ELEMENTAL FUNCTION Cube(c) RESULT(cuberesult)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), INTENT(IN):: c
      REAL(KIND=r8):: cuberesult
!
      cuberesult = c*c*c
!
   END FUNCTION Cube
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ELEMENTAL FUNCTION Sqr(c) RESULT(sqrresult)
      USE NodeInfoDef
      USE ProblemDef
      IMPLICIT NONE
!
      REAL(KIND=r8), INTENT(IN):: c
      REAL(KIND=r8):: sqrresult
!
      sqrresult = c*c
!
   END FUNCTION Sqr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ELEMENTAL FUNCTION fof1(a) RESULT(ulapdiffresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Level independent, 2D UNDIVIDED laplacian operator for non-constant
! diffusivity.
!
      REAL(KIND=r8), INTENT(IN):: a
      REAL(KIND=r8) :: ulapdiffresult
!
! Calculate the UNDIVIDED divergence of the flux:
      ulapdiffresult = 0.18_r8 *a*(a-1.0_r8)*(a-0.5_r8)
!
   END FUNCTION fof1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ELEMENTAL FUNCTION fof2(a) RESULT(ulapdiffresult)
      USE NodeInfoDef
      IMPLICIT NONE
!
! Level independent, 2D UNDIVIDED laplacian operator for non-constant
! diffusivity.
!
      REAL(KIND=r8), INTENT(IN):: a
      REAL(KIND=r8) :: ulapdiffresult
!
! Calculate the UNDIVIDED divergence of the flux:
      ulapdiffresult = 0.18_r8* ( (a-1.0_r8)*(a-0.5_r8) &
         + a*           (a-0.5_r8) &
         + a*(a-1.0_r8))
!
   END FUNCTION fof2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elemental function Mob(c) result(mobilresult)
      use ProblemDef
      implicit none
      real(KIND=r8), intent(in):: c
      real(KIND=r8):: mobilresult
      real(KIND=r8), parameter :: small = 1e-6

! This definition ensures that M stays positive

      mobilresult = c*c/epep + small

   end function Mob
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elemental function MobD(c) result(mobilresult)
      use ProblemDef
      implicit none
      real(KIND=r8), intent(in):: c
      real(KIND=r8):: mobilresult,mc
      real(KIND=r8), parameter :: small = 1e-6

! This definition ensures that M stays positive

      mc = 1.0_r8-c
!
      mobilresult = c+ DtDh*mc
!
   end function MobD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function Fof(a) result(ulapdiffresult)
      USE NodeInfoDef
      use ProblemDef
!
      implicit none
!
      REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: a
      REAL(KIND=r8), DIMENSION(1:SIZE(a,1),1:SIZE(a,2)):: d
      REAL(KIND=r8), DIMENSION(1:SIZE(a,1),1:SIZE(a,2)):: ulapdiffresult
!
      INTEGER, DIMENSION(1:2):: mx
      INTEGER:: i, j
!
      mx(1) = SIZE(a,1); mx(2) = SIZE(a,2)
!
      d(1:mx(1),1:mx(2)) =       0.18_r8*  (a(1:mx(1),1:mx(2))       )     &
         *(a(1:mx(1),1:mx(2))-1.0_r8)     &
         *(a(1:mx(1),1:mx(2))-0.5_r8)
!
      ulapdiffresult(1:mx(1),1:mx(2)) = d(1:mx(1) ,1:mx(2))
!
   end function Fof
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Problem


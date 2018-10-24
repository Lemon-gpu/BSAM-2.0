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
!
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
! File:             afasroutines.f90
! Purpose:          Adaptive FAS Multigrid module
! Contains:
! Revision History: Ver. 1.0 Oct. 2rf1006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 2.0 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------
MODULE AFASRoutines
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	SUBROUTINE MultigridIterations
USE NodeInfoDef !
USE TreeOps, ONLY: ApplyOnLevel ! 
IMPLICIT NONE
!
TYPE(funcparam):: dummy ! what is funcparam?
INTEGER:: itmg, level !this level means the level of the grid
REAL(KIND=r8), DIMENSION(1:2):: residual
!
! initializtion of the parents grid for each patch
CALL FillDown(solutionfield)
!
! If the auxiliary fields are updated in the vcycle, perform here:
DO level = finestlevel, minlevel, -1
  CALL ApplyOnLevel(level,UpdateAuxInVcycle,dummy) ! applyonlevel call the function
END DO
!
! Maybe implement this later, but it needs more work:
!CALL FillDown(auxiliaryfield)
! what is the difference between the minlevel and the rro
DO level = finestlevel, rootlevel, -1 
  CALL ApplyOnLevel(level,GetSourceFunction,dummy)
END DO
!
! If the source function is updated after each vcycle, perform here first:
DO level = finestlevel, rootlevel, -1
  CALL ApplyOnLevel(level,UpdateSourceFunction,dummy)
END DO
!
CALL FillDown(sourcefield)
!
! Perform Adaptive Full Aproximation Storage Vcycles on the grid hierarchy:
vcycleloop: DO itmg = 1, maxvcycles
!
  CALL AFASVcycle(finestlevel) ! iteration is finished by this step
!                              ! apply on level, means 
  CALL FillDown(solutionfield)
!
! If the auxiliary fields are updated in the vcycle, perform here:
  IF(MODULO(itmg-1,updateauxfreq)==0) THEN
    DO level = finestlevel, minlevel, -1
      CALL ApplyOnLevel(level,UpdateAuxInVcycle,dummy)
    END DO
  END IF
!
! Maybe implement this later, but it needs more work:
!  CALL FillDown(auxiliaryfield)
!
! If the source function is updated after each vcycle, perform here:
  DO level = finestlevel, rootlevel, -1
    CALL ApplyOnLevel(level,UpdateSourceFunction,dummy)
  END DO
!
  CALL FillDown(sourcefield)
!
! Note that it is necessary to have an up-to-date source function for 
! calculating the residual:
  residual(1:2) = ErrorAFAS(finestlevel)
  PRINT 1001, itmg, residual(1), residual(2)
  1001 FORMAT('it =', i3, 3X, 'ccv error =', ES15.7, 3X, ' fcv error =' ES15.7)
!                         
  IF(residual(1)<qerrortol) EXIT vcycleloop
!
END DO vcycleloop
!
END SUBROUTINE MultigridIterations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ErrorAFAS(level) RESULT(errorresult)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel, GetRootInfo
IMPLICIT NONE
!
TYPE(nodeinfo), POINTER:: rootinfo
!
INTEGER, INTENT(IN):: level
REAL(KIND=r8), DIMENSION(1:2):: errorresult
!
TYPE(funcparam):: dummy
INTEGER:: i, ierror
REAL(KIND=r8), DIMENSION(1:maxnccv)::  componenterror
!
SELECT CASE(errortype)
  CASE(1)
!
! Blended error calculation of cell centered variables:
    integralresult(1:2) = 0.0_r8
!
    CALL ApplyOnLevel(level,L2Error,dummy)
!
    errorresult(1) = SQRT(integralresult(1)/integralresult(2))
!
! Blended error calculation of face centered variables:
    integralresult(1:2) = 0.0_r8
!
    CALL ApplyOnLevel(level,L2ErrorFC,dummy)
!
    errorresult(2) = SQRT(integralresult(1)/integralresult(2))
  CASE(2)
!
! Component errors:
    ierror = GetRootInfo(rootinfo) ! Why? Check this!??????????????????
    integralresult(2) = 0.0_r8
    componentintegral(1:nccv) = 0.0_r8
!
    CALL ApplyOnLevel(level,L2ComponentErrors,dummy)
!
    componenterror(1:nccv) = SQRT(componentintegral(1:nccv) &
                           / integralresult(2))
!
    DO i = 1, nccv
      PRINT *, i, componenterror(i)
    END DO
!
    errorresult = MAXVAL(componenterror(1:nccv))
  CASE DEFAULT
    PRINT *, 'ErrorAFAS: Only errortype = 1,2 are supported.'
    STOP
!
END SELECT
!
END FUNCTION ErrorAFAS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION L2Error(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
! Adds the L2 error on this grid to the total.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER, DIMENSION(1:maxdims):: mx
!
L2Error = err_ok
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
!
CALL Residual(info)
!
integralresult(1) = SquareSum(info%rf( 1:mx(1),1:mx(2),1:mx(3),1:nccv)) &
                  + integralresult(1)
!
integralresult(2) = REAL(nccv*PRODUCT(mx(1:ndims)),KIND=r8) &
                  + integralresult(2)
!
END FUNCTION L2Error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION L2ErrorFC(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
! Adds the L2 error on this grid to the total.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER, DIMENSION(1:maxdims):: mx
REAL(KIND=r8), DIMENSION(:,:,:,:), ALLOCATABLE:: tmp
!
L2ErrorFC = err_ok
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
!
ALLOCATE(tmp(0:mx(1),0:mx(2),0:mx(3),1:nfcv))
!
! Not optimized:
!
! v1 residual:
         tmp(0:mx(1),1:mx(2),1:mx(3),1:nfcv) &
  = info%rf1(0:mx(1),1:mx(2),1:mx(3),1:nfcv) &
  * info%rf1(0:mx(1),1:mx(2),1:mx(3),1:nfcv)
!
            tmp(1:mx(1)  ,1:mx(2),1:mx(3),1:nfcv) &
  = 0.5_r8*(tmp(0:mx(1)-1,1:mx(2),1:mx(3),1:nfcv) &
  +         tmp(1:mx(1)  ,1:mx(2),1:mx(3),1:nfcv))
!
integralresult(1) = SUM(tmp(1:mx(1),1:mx(2),1:mx(3),1:nfcv)) &
                  + integralresult(1)
!PRINT *, 'fieldsallocated', integralresult(1)
!
! v2 residual:
         tmp(1:mx(1),0:mx(2),1:mx(3),1:nfcv) &
  = info%rf2(1:mx(1),0:mx(2),1:mx(3),1:nfcv) &
  * info%rf2(1:mx(1),0:mx(2),1:mx(3),1:nfcv)
!
            tmp(1:mx(1),1:mx(2)  ,1:mx(3),1:nfcv) &
  = 0.5_r8*(tmp(1:mx(1),0:mx(2)-1,1:mx(3),1:nfcv) &
  +         tmp(1:mx(1),1:mx(2)  ,1:mx(3),1:nfcv))
!
integralresult(1) = SUM(tmp(1:mx(1),1:mx(2),1:mx(3),1:nfcv)) &
                  + integralresult(1)
!
IF (ndims==3) THEN
!
! v3 residual:
         tmp(1:mx(1),1:mx(2),0:mx(3),1:nfcv) &
  = info%rf3(1:mx(1),1:mx(2),0:mx(3),1:nfcv) &
  * info%rf3(1:mx(1),1:mx(2),0:mx(3),1:nfcv)
!
            tmp(1:mx(1),1:mx(2),1:mx(3)  ,1:nfcv) &
  = 0.5_r8*(tmp(1:mx(1),1:mx(2),0:mx(3)-1,1:nfcv) &
  +         tmp(1:mx(1),1:mx(2),1:mx(3)  ,1:nfcv))
!
integralresult(1) = SUM(tmp(1:mx(1),1:mx(2),1:mx(3),1:nfcv)) &
                  + integralresult(1)
!
END IF
!
integralresult(2) = REAL(ndims*nfcv*PRODUCT(mx(1:ndims)),KIND=r8) &
                  + integralresult(2)
!
DEALLOCATE(tmp)
!
END FUNCTION L2ErrorFC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION L2ComponentErrors(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
! Adds the L2 component error on this grid to the component total.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: i
INTEGER, DIMENSION(1:maxdims):: mx
!
L2ComponentErrors = err_ok
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
!
CALL Residual(info)
!
DO i = 1, nccv
  componentintegral(i) = SquareSum(info%rf(1:mx(1),1:mx(2),1:mx(3),i:i)) &
                       + componentintegral(i)
END DO
!
integralresult(2) = REAL(PRODUCT(mx(1:ndims)),KIND=r8)+integralresult(2)
!
END FUNCTION L2ComponentErrors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION SquareSum(rf) RESULT(squaresumresult)
USE NodeInfoDef
IMPLICIT NONE
!
! Dimensionally invariant square `integral.`
!
REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(IN):: rf
REAL(KIND=r8):: squaresumresult
!
squaresumresult = SUM(rf*rf)
!
END FUNCTION SquareSum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE FillDown(field)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
USE Boundary, ONLY: SetGhost
IMPLICIT NONE
!
INTEGER, INTENT(IN):: field
!
TYPE(funcparam):: dummy
INTEGER:: level
!
dummy%iswitch = field
!
DO level = finestlevel, minlevel+1, -1
  CALL ApplyOnLevel(level,FillDownLevel,dummy)
  IF(field==solutionfield) CALL SetGhost(level-1,0)
END DO
!
END SUBROUTINE FillDown
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION FillDownLevel(info,dummy)!??????? nothing happens in face variables??
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY: Restriction2D, Restriction2DV1, Restriction2DV2, &
                         Restriction3D, Restriction3DV1, Restriction3DV2, Restriction3DV3
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: field, ierror, nnaxv
INTEGER, DIMENSION(1:maxdims):: cmx, mx
INTEGER, DIMENSION(1:maxdims,1:2):: mb
!
FillDownLevel = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
field = dummy%iswitch
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
!nnaxv = MAX(naxv,1)
!
ierror = GetParentInfo(parent)
!
SELECT CASE(field)
  CASE(sourcefield)
!
    SELECT CASE(ndims)
!
      CASE(2)
!
                     parent%f(  mb(1,1):mb(1,2),mb(2,1):mb(2,2),1,1:nccv) &
    = Restriction2D  ( info%f(        1:mx(1)  ,      1:mx(2)  ,1,1:nccv))
!
      CASE(3)
!
                     parent%f(  mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nccv) &
    = Restriction3D  ( info%f(        1:mx(1)  ,      1:mx(2)  ,      1:mx(3)  ,1:nccv))
!
      CASE DEFAULT
!
        PRINT *, 'FillDownLevel: Only ndims = 2,3 are supported.'
        STOP
!
    END SELECT
!
  CASE(solutionfield)
!
    SELECT CASE(ndims)
      CASE(2)
!
                    parent%q (mb(1,1)  :mb(1,2),mb(2,1)  :mb(2,2),1,1:nccv) &
    = Restriction2D  (info%q (        1:mx(1)  ,        1:mx(2)  ,1,1:nccv))
!        
                    parent%v1(mb(1,1)-1:mb(1,2),mb(2,1)  :mb(2,2),1,1:nfcv) &
    = Restriction2DV1(info%v1(        0:mx(1)  ,        1:mx(2)  ,1,1:nfcv))
!
                    parent%v2(mb(1,1)  :mb(1,2),mb(2,1)-1:mb(2,2),1,1:nfcv) &
    = Restriction2DV2(info%v2(        1:mx(1)  ,        0:mx(2)  ,1,1:nfcv))  
!                   
      CASE(3)
!
                    parent%q (mb(1,1)  :mb(1,2),mb(2,1)  :mb(2,2),mb(3,1)  :mb(3,2),1:nccv) &
    = Restriction3D  (info%q (        1:mx(1)  ,        1:mx(2)  ,        1:mx(3)  ,1:nccv))
!        
                    parent%v1(mb(1,1)-1:mb(1,2),mb(2,1)  :mb(2,2),mb(3,1)  :mb(3,2),1:nfcv) &
    = Restriction3DV1(info%v1(        0:mx(1)  ,        1:mx(2)  ,        1:mx(3)  ,1:nfcv))
!
                    parent%v2(mb(1,1)  :mb(1,2),mb(2,1)-1:mb(2,2),mb(3,1)  :mb(3,2),1:nfcv) &
    = Restriction3DV2(info%v2(        1:mx(1)  ,        0:mx(2)  ,        1:mx(3)  ,1:nfcv))  
!                   
                    parent%v3(mb(1,1)  :mb(1,2),mb(2,1)  :mb(2,2),mb(3,1)-1:mb(3,2),1:nfcv) &
    = Restriction3DV3(info%v3(        1:mx(1)  ,        1:mx(2)  ,        0:mx(3)  ,1:nfcv))  
!                   
      CASE DEFAULT
        PRINT *, 'FillDownLevel: Only ndims = 2,3 are supported.'
        STOP
    END SELECT
!  CASE(auxiliaryfield)
!    IF(naxv>0) THEN
!      SELECT CASE(ndims)
!        CASE(2)
!          parent%aux(mb(1,1):mb(1,2),mb(2,1):mb(2,2),1,1:mnaxv) &
!            = Restriction2D(info%aux(1:mx(1),1:mx(2),1,1:mnaxv))
!        CASE(3)
!          parent%aux(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:mnaxv) &
!            = Restriction3D(info%aux(1:mx(1),1:mx(2),1:mx(3),1:mnaxv))
!        CASE DEFAULT
!          PRINT *, 'FillDownLevel: Only ndims = 2,3 are supported.'
!          STOP
!      END SELECT
!    END IF
END SELECT
!
END FUNCTION FillDownLevel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION GetSourceFunction(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE Problem, ONLY: Source2D, Source2DV1, Source2DV2, &
                   Source3D, Source3DV1, Source3DV2, Source3DV3
IMPLICIT NONE
!
! Dimensionally invariant source (or partial source) function routine.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: ll, nnaxv
INTEGER, DIMENSION(1:maxdims):: amx, aul, mx, ul
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
GetSourceFunction = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
nnaxv = MAX(naxv,1)
amx = 1; IF(naxv>0) amx(1:ndims) = mx(1:ndims)
!
ll = 1-mbc
 ul = 1;  ul(1:ndims) =  mx(1:ndims)+mbc
aul = 1; aul(1:ndims) = amx(1:ndims)+mbc
!
SELECT CASE(ndims)
  CASE(2)
!
! Source for all cell centered variables:
                   info%ftmp (  1: mx(1)  , 1: mx(2)  ,1,1: nccv)  &
        = Source2D(info%q    ( ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                   info%qold ( ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                   info%v1   ( -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                   info%v1old( -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                   info%v2   (  0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                   info%v2old(  0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                   info%aux  ( ll:aul(1)  ,ll:aul(2)  ,1,1:nnaxv), &
                   ll,mx(1:2),h,xlower(1:2))
!
! Source for face centered variables v1 and v2:
                   info%f1tmp(  0: mx(1)  , 1: mx(2)  ,1,1: nfcv)  &
      = Source2DV1(info%q    ( ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                   info%qold ( ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                   info%v1   ( -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                   info%v1old( -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                   info%v2   (  0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                   info%v2old(  0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                   info%aux  ( ll:aul(1)  ,ll:aul(2)  ,1,1:nnaxv), &
                   ll,mx(1:2),h,xlower(1:2))
!
                   info%f2tmp(  1: mx(1)  , 0: mx(2)  ,1,1: nfcv)  &
      = Source2DV2(info%q    ( ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                   info%qold ( ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                   info%v1   ( -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                   info%v1old( -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                   info%v2   (  0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                   info%v2old(  0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                   info%aux  ( ll:aul(1)  ,ll:aul(2)  ,1,1:nnaxv), &
                   ll,mx(1:2),h,xlower(1:2))
  CASE(3)
!
                   info%ftmp (  1: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1: nccv)  &
        = Source3D(info%q    ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                   info%qold ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                   info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v1old( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v2   (  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v2old(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v3   (  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                   info%v3old(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                   info%aux  ( ll:aul(1)  ,ll:aul(2)  ,ll:aul(3)  ,1:nnaxv), &
                   ll,mx(1:3),h,xlower(1:3))
!
                   info%f1tmp(  0: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1: nfcv)  &
      = Source3DV1(info%q    ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                   info%qold ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                   info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v1old( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v2   (  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v2old(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v3   (  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                   info%v3old(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                   info%aux  ( ll:aul(1)  ,ll:aul(2)  ,ll:aul(3)  ,1:nnaxv), &
                   ll,mx(1:3),h,xlower(1:3))
!
                   info%f2tmp(  1: mx(1)  , 0: mx(2)  , 1: mx(3)  ,1: nfcv)  &
      = Source3DV2(info%q    ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                   info%qold ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                   info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v1old( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v2   (  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v2old(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v3   (  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                   info%v3old(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                   info%aux  ( ll:aul(1)  ,ll:aul(2)  ,ll:aul(3)  ,1:nnaxv), &
                   ll,mx(1:3),h,xlower(1:3))
!
                   info%f3tmp(  1: mx(1)  , 1: mx(2)  , 0: mx(3)  ,1: nfcv)  &
      = Source3DV3(info%q    ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                   info%qold ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                   info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v1old( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v2   (  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v2old(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v3   (  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                   info%v3old(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                   info%aux( ll:aul(1)  ,ll:aul(2)  ,ll:aul(3)  ,1:nnaxv), &
                   ll,mx(1:3),h,xlower(1:3))
!
  CASE DEFAULT
    PRINT *, 'GetSourceFunction: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION GetSourceFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION UpdateAuxInVcycle(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE Problem, ONLY: UpdateAuxVcycle2D, UpdateAuxVcycle3D
IMPLICIT NONE
!
! Used if the auxiliary variables are updated after each Vcycle iteration.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: ll, nnaxv
INTEGER, DIMENSION(1:maxdims):: amx, aul, mx, ul
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
UpdateAuxInVcycle = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
nnaxv = MAX(naxv,1)
amx = 1; IF(naxv>0) amx(1:ndims) = mx(1:ndims)
!
ll = 1-mbc
 ul = 1;  ul(1:ndims) =  mx(1:ndims)+mbc
aul = 1; aul(1:ndims) = amx(1:ndims)+mbc
!
IF(naxv<=0) RETURN
!
SELECT CASE(ndims)
  CASE(2)
    CALL UpdateAuxVcycle2D(info%q    (  ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                           info%qold (  ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                           info%v1   (  -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                           info%v1old(  -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                           info%v2   (   0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                           info%v2old(   0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                           info%aux  (  ll:aul(1)  ,ll:aul(2)  ,1,1:nnaxv), &
                           ll,mx(1:2),h,xlower(1:2))
  CASE(3)
    CALL UpdateAuxVcycle3D(info%q(   ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                           info%qold(ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                           info%aux( ll:aul(1)  ,ll:aul(2)  ,ll:aul(3)  ,1:nnaxv), &
                           ll,mx(1:3),h,xlower(1:3))
  CASE DEFAULT
    PRINT *, 'UpdateAuxInVcycle: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION UpdateAuxInVcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION UpdateSourceFunction(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE Problem, ONLY: SourceUpdate2D, SourceUpdate2DV1, SourceUpdate2DV2, &
                   SourceUpdate3D, SourceUpdate3DV1, SourceUpdate3DV2, SourceUpdate3DV3
IMPLICIT NONE
!
! Dimensionally invariant source update routine.  If no update is needed set 
! SourceUpdate# = 0.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: ll, nnaxv
INTEGER, DIMENSION(1:maxdims):: amx, aul, mx, ul
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
UpdateSourceFunction = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
nnaxv = MAX(naxv,1)
amx = 1; IF(naxv>0) amx(1:ndims) = mx(1:ndims)
!
ll = 1-mbc
 ul = 1;  ul(1:ndims) =  mx(1:ndims)+mbc
aul = 1; aul(1:ndims) = amx(1:ndims)+mbc
!
SELECT CASE(ndims)
  CASE(2)
!
! All cell centered variables:
                         info%f    (   1:mx(1)   , 1: mx(2)  ,1,1: nccv)  &
      =                  info%ftmp (   1:mx(1)   , 1: mx(2)  ,1,1: nccv)  &
      +   SourceUpdate2D(info%q    (  ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                         info%qold (  ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                         info%v1   (  -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                         info%v1old(  -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                         info%v2   (   0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                         info%v2old(   0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                         info%aux  (  ll:aul(1)  ,ll:aul(2)  ,1,1:nnaxv), &
                         ll,mx(1:2),h,xlower(1:2))
!
! East-West cell edge variables, v1:
                         info%f1   (  0: mx(1)  , 1: mx(2)  ,1,1: nfcv)  &
      =                  info%f1tmp(  0: mx(1)  , 1: mx(2)  ,1,1: nfcv)  &
      + SourceUpdate2DV1(info%q    ( ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                         info%qold ( ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                         info%v1   ( -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                         info%v1old( -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                         info%v2   (  0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                         info%v2old(  0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                         info%aux  ( ll:aul(1)  ,ll:aul(2)  ,1,1:nnaxv), &
                         ll,mx(1:2),h,xlower(1:2))
!
! North-South cell edge variables, v2:
                         info%f2   (  1: mx(1)  , 0: mx(2)  ,1,1: nfcv)  &
      =                  info%f2tmp(  1: mx(1)  , 0: mx(2)  ,1,1: nfcv)  &
      + SourceUpdate2DV2(info%q    ( ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                         info%qold ( ll: ul(1)  ,ll: ul(2)  ,1,1: nccv), &
                         info%v1   ( -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                         info%v1old( -1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                         info%v2   (  0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                         info%v2old(  0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                         info%aux  ( ll:aul(1)  ,ll:aul(2)  ,1,1:nnaxv), &
                         ll,mx(1:2),h,xlower(1:2))
  CASE(3)
!
                         info%f    (  1: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1: nccv)  &
      =                  info%ftmp (  1: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1: nccv)  &
      +   SourceUpdate3D(info%q    ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                         info%qold ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                         info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v1old( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v2   (  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v2old(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v3   (  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                         info%v3old(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                         info%aux  ( ll:aul(1)  ,ll:aul(2)  ,ll:aul(3)  ,1:nnaxv), &
                         ll,mx(1:3),h,xlower(1:3))
!
                         info%f1   (  0: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1: nfcv)  &
      =                  info%f1tmp(  0: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1: nfcv)  &
      + SourceUpdate3DV1(info%q    ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                         info%qold ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                         info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v1old( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v2   (  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v2old(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v3   (  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                         info%v3old(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                         info%aux  ( ll:aul(1)  ,ll:aul(2)  ,ll:aul(3)  ,1:nnaxv), &
                         ll,mx(1:3),h,xlower(1:3))
!
                         info%f2   (  1: mx(1)  , 0: mx(2)  , 1: mx(3)  ,1: nfcv)  &
      =                  info%f2tmp(  1: mx(1)  , 0: mx(2)  , 1: mx(3)  ,1: nfcv)  &
      + SourceUpdate3DV2(info%q    ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                         info%qold ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                         info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v1old( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v2   (  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v2old(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v3   (  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                         info%v3old(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                         info%aux  ( ll:aul(1)  ,ll:aul(2)  ,ll:aul(3)  ,1:nnaxv), &
                         ll,mx(1:3),h,xlower(1:3))
!
                         info%f3   (  1: mx(1)  , 1: mx(2)  , 0: mx(3)  ,1: nfcv)  &
      =                  info%f3tmp(  1:mx(1)   , 1: mx(2)  , 0: mx(3)  ,1: nfcv)  &
      + SourceUpdate3DV3(info%q    ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                         info%qold ( ll: ul(1)  ,ll: ul(2)  ,ll: ul(3)  ,1: nccv), &
                         info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v1old( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v2   (  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v2old(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                         info%v3   (  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                         info%v3old(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                         info%aux  ( ll:aul(1)  ,ll:aul(2)  ,ll:aul(3)  ,1:nnaxv), &
                         ll,mx(1:3),h,xlower(1:3))
!
  CASE DEFAULT
    PRINT *, 'UpdateSourceFunction: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION UpdateSourceFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
RECURSIVE SUBROUTINE AFASVcycle(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
USE Boundary, ONLY: SetGhost, GetCoarseGhostPoints
IMPLICIT NONE
!
! Dimensionally invariant recursive AFAS vcycle routine.
!
INTEGER, INTENT(IN):: level
!
TYPE(funcparam):: dummy
!
CALL LevelRelax(level)
!
IF(level>minlevel) THEN
!
  CALL RestrictSolution(level)   !1:cm(1)
!
  CALL SetGhost(level-1,0)
!
  CALL ApplyOnLevel(level,GetCoarseGhostPoints,dummy)
!
  CALL GetCoarseLevelLoading(level)
!
  CALL AFASVcycle(level-1)
!
  CALL ApplyOnLevel(level,CorrectFine,dummy) 
!
  CALL SetGhost(level,0)
!
  CALL LevelRelax(level)
!
END IF
!
END SUBROUTINE AFASVcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE LevelRelax(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
USE Boundary, ONLY: SetGhost
IMPLICIT NONE
!
! Dimensionally invariant relaxation.
!
INTEGER, INTENT(IN):: level
!
TYPE(funcparam):: dummy
INTEGER:: redblack, smoothingpass
!
! Assume Ghost points are set on entry to this routine.
DO smoothingpass = 1, nsmoothingpasses
  DO redblack = 1, 2
!
    dummy%iswitch = redblack
    CALL ApplyOnLevel(level,RelaxPatchCC,dummy)
!
    CALL ApplyOnLevel(level,RelaxPatchFC,dummy)
!
    CALL SetGhost(level,redblack)
!
  END DO
END DO
!
END SUBROUTINE LevelRelax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION RelaxPatchCC(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
USE Problem, ONLY: RelaxGrid2D, RelaxGrid3D
IMPLICIT NONE
!
! Dimensionally invariant relaxation of a patch, on RED squares for
! dummy%iswitch = 1, and BLACK squares for dummy%iswitch = 2.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: nnaxv
INTEGER, DIMENSION(1:maxdims):: amx, mx
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
RelaxPatchCC = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
nnaxv = MAX(naxv,1)
amx = 1; IF(naxv>0) amx(1:ndims) = mx(1:ndims)
!
SELECT CASE(ndims)
  CASE(2)
    CALL RelaxGrid2D(info%q    (  0: mx(1)+1, 0: mx(2)+1, 1,1: nccv), &
                     info%qold (  0: mx(1)+1, 0: mx(2)+1, 1,1: nccv), &
                     info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 1,1: nfcv), &
                     info%v1old( -1: mx(1)+1, 0: mx(2)+1, 1,1: nfcv), &
                     info%v2   (  0: mx(1)+1,-1: mx(2)+1, 1,1: nfcv), &
                     info%v2old(  0: mx(1)+1,-1: mx(2)+1, 1,1: nfcv), &
                     info%aux  (  0:amx(1)+1, 0:amx(2)+1, 1,1:nnaxv), &
                     info%f    (  1: mx(1)  , 1: mx(2)  , 1,1: nccv), &
                     mx(1:2),h,xlower(1:2),dummy%iswitch)
  CASE(3)
    CALL RelaxGrid3D(info%q    (  0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                     info%qold (  0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                     info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v1old( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v2   (  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v2old(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v3   (  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                     info%v3old(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                     info%aux  (  0:amx(1)+1, 0:amx(2)+1, 0:amx(3)+1,1:nnaxv), &
                     info%f    (  1: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1: nccv), &
                     mx(1:3),h,xlower(1:3),dummy%iswitch)
  CASE DEFAULT
    PRINT *, 'RelaxPatchCC: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION RelaxPatchCC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION RelaxPatchFC(info,dummy1)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, ApplyOnMyLevelNbrs
USE Problem, ONLY: RelaxGrid2DEdge, RelaxGrid3DEdge
USE Boundary, ONLY: TransferBCOneWayFC
IMPLICIT NONE
!
! Dimensionally invariant relaxation of a patch, on RED squares for
! dummy%iswitch = 1, and BLACK squares for dummy%iswitch = 2.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy1
!
TYPE(funcparam):: dummy2
INTEGER:: level, nnaxv
INTEGER, DIMENSION(1:maxdims):: amx, mx
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
RelaxPatchFC = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
level = info%level
mx = 1; mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
nnaxv = MAX(naxv,1)
amx = 1; IF(naxv>0) amx(1:ndims) = mx(1:ndims)
!
SELECT CASE(ndims)
CASE(2)
!
  CALL RelaxGrid2DEdge(info%q    (  0: mx(1)+1, 0: mx(2)+1, 1,1: nccv), &
                       info%qold (  0: mx(1)+1, 0: mx(2)+1, 1,1: nccv), &
                       info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 1,1: nfcv), &
                       info%v1old( -1: mx(1)+1, 0: mx(2)+1, 1,1: nfcv), &
                       info%v2   (  0: mx(1)+1,-1: mx(2)+1, 1,1: nfcv), &
                       info%v2old(  0: mx(1)+1,-1: mx(2)+1, 1,1: nfcv), &
                       info%aux  (  0:amx(1)+1, 0:amx(2)+1, 1,1:nnaxv), &
                       info%f    (  1: mx(1)  , 1: mx(2)  , 1,1: nccv), &
                       info%f1   (  0: mx(1)  , 1: mx(2)  , 1,1: nfcv), &
                       info%f2   (  1: mx(1)  , 0: mx(2)  , 1,1: nfcv), &
                       level,mx(1:2),h,xlower(1:2),dummy1%iswitch)
!
CASE(3)
!
  CALL RelaxGrid3DEdge(info%q    (  0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                       info%qold (  0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                       info%v1   ( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                       info%v1old( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                       info%v2   (  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                       info%v2old(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                       info%v3   (  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                       info%v3old(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                       info%aux  (  0:amx(1)+1, 0:amx(2)+1, 0:amx(3)+1,1:nnaxv), &
                       info%f    (  1: mx(1),   1: mx(2)  , 1: mx(3)  ,1: nccv), &
                       info%f1   (  0: mx(1),   1: mx(2)  , 1: mx(3)  ,1: nfcv), &
                       info%f2   (  1: mx(1),   0: mx(2)  , 1: mx(3)  ,1: nfcv), &
                       info%f3   (  1: mx(1),   1: mx(2)  , 0: mx(3)  ,1: nfcv), &
                       level,mx(1:3),h,xlower(1:3),dummy1%iswitch)
!
CASE DEFAULT
  PRINT *, 'BSAM 2.0 RelaxPatchFC: Only ndims = 2,3 are supported.'
  STOP
END SELECT
!
! Immediately sync-up boundary values:
dummy2%iswitch = 1 ! velocities:
IF(level > rootlevel) CALL ApplyOnMyLevelNbrs(level,TransferBCOneWayFC,dummy2)
!
END FUNCTION RelaxPatchFC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE RestrictSolution(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
IMPLICIT NONE
!
! Dimensionally invariant restriction of CC and FC variables:
!
INTEGER, INTENT(IN):: level
!
TYPE(funcparam):: dummy
!
! Cell Centered Variables: Children push down the restriction values:
CALL ApplyOnLevel(level,RestrictCCSolution,dummy)
!
! Face Centered Variables: Parents must pull up the restriction values, then
! synchronize data with other parent-level grids: 
CALL ApplyOnLevel(level-1,PullFCRestrictions,dummy)
!
END SUBROUTINE RestrictSolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION RestrictCCSolution(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY: Restriction2D, Restriction2DV1, Restriction2DV2, Restriction3D
IMPLICIT NONE
!
! Dimensionally invariant restriction.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror
INTEGER, DIMENSION(1:maxdims):: cmx, mx
INTEGER, DIMENSION(1:maxdims,1:2):: mb
!
RestrictCCSolution = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
ierror = GetParentInfo(parent)
!
SELECT CASE(ndims)
  CASE(2)
                      info%qc(1:cmx(1),1:cmx(2),1,1:nccv) &
      = Restriction2D(info%q (1: mx(1),1: mx(2),1,1:nccv))
  CASE(3)
                      info%qc(1:cmx(1),1:cmx(2),1:cmx(3),1:nccv) &
      = Restriction3D(info%q (1: mx(1),1: mx(2),1: mx(3),1:nccv))
  CASE DEFAULT
    PRINT *, 'BSAM 2.0 RestrictCCSolution: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
  parent%q (mb(1,1): mb(1,2),mb(2,1): mb(2,2),mb(3,1): mb(3,2),1:nccv) &
 =  info%qc(      1:cmx(1)  ,      1:cmx(2)  ,      1:cmx(3)  ,1:nccv)
!
END FUNCTION RestrictCCSolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION PullFCRestrictions(parent,dummy1)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, ApplyOnChildren, ApplyOnMyLevelNbrs
USE Boundary, ONLY: TransferBCOneWayFC
IMPLICIT NONE
!
! Dimensionally invariant restriction.
!
TYPE(nodeinfo):: parent
TYPE(funcparam):: dummy1
!
TYPE(funcparam):: dummy2
!
PullFCRestrictions = err_ok
!
! Pull the child-level data down to the parent:
CALL ApplyOnChildren(RestrictFCSolution,dummy2)
!
! Synchronize the parent's face centered variable data with other parent-
! level grids. iswitch = 1 indicates velocity data v#:
dummy2%iswitch = 1
IF(parent%level > rootlevel) &
      CALL ApplyOnMyLevelNbrs(parent%level,TransferBCOneWayFC,dummy2)
!
END FUNCTION PullFCRestrictions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION RestrictFCSolution(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY: Restriction2DV1, Restriction2DV2, &
                         Restriction3DV1, Restriction3DV2, Restriction3DV3
IMPLICIT NONE
!
! Dimensionally invariant restriction of face centered variable data:
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror
INTEGER, DIMENSION(1:maxdims):: cmx, mx
INTEGER, DIMENSION(1:maxdims,1:2):: mb
!
RestrictFCSolution = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
 mx = 1;  mx(1:ndims)     = info%mx     (1:ndims)
cmx = 1; cmx(1:ndims)     =           mx(1:ndims)/2 
 mb = 1;  mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
ierror = GetParentInfo(parent)
!
SELECT CASE(ndims)
CASE(2)
!
                      info%v1c(0:cmx(1),1:cmx(2),1,1:nfcv) &
    = Restriction2DV1(info%v1 (0: mx(1),1: mx(2),1,1:nfcv))
!
                      info%v2c(1:cmx(1),0:cmx(2),1,1:nfcv) &
    = Restriction2DV2(info%v2 (1: mx(1),0: mx(2),1,1:nfcv))
!    
   ! This saves the data into parent, but also keeps a 'full' approximation in 
! v#c needed for the AFAS method. We will still need to sync-up the bndry
! data of v#c:
parent%v1(mb(1,1)-1: mb(1,2),mb(2,1)  : mb(2,2),mb(3,1): mb(3,2),1:nfcv) &
  = info%v1c(     0:cmx(1)  ,        1:cmx(2)  ,      1:cmx(3)  ,1:nfcv)
parent%v2(mb(1,1)  : mb(1,2),mb(2,1)-1: mb(2,2),mb(3,1): mb(3,2),1:nfcv) &
  = info%v2c(     1:cmx(1)  ,        0:cmx(2)  ,      1:cmx(3)  ,1:nfcv) 
!    
CASE(3)
!
                      info%v1c(0:cmx(1),1:cmx(2),1:cmx(3),1:nfcv) &
    = Restriction3DV1(info%v1 (0: mx(1),1: mx(2),1: mx(3),1:nfcv))
!
                      info%v2c(1:cmx(1),0:cmx(2),1:cmx(3),1:nfcv) &
    = Restriction3DV2(info%v2 (1: mx(1),0: mx(2),1: mx(3),1:nfcv))
!
                      info%v3c(1:cmx(1),1:cmx(2),0:cmx(3),1:nfcv) &
    = Restriction3DV3(info%v3 (1: mx(1),1: mx(2),0: mx(3),1:nfcv))
!
  ! This saves the data into parent, but also keeps a 'full' approximation in 
! v#c needed for the AFAS method. We will still need to sync-up the bndry
! data of v#c:
!
parent%v1(mb(1,1)-1: mb(1,2),mb(2,1)  : mb(2,2),mb(3,1)  : mb(3,2),1:nfcv) &
  = info%v1c(     0:cmx(1)  ,        1:cmx(2)  ,        1:cmx(3)  ,1:nfcv)
parent%v2(mb(1,1)  : mb(1,2),mb(2,1)-1: mb(2,2),mb(3,1)  : mb(3,2),1:nfcv) &
  = info%v2c(     1:cmx(1)  ,        0:cmx(2)  ,        1:cmx(3)  ,1:nfcv)
parent%v3(mb(1,1)  : mb(1,2),mb(2,1)  : mb(2,2),mb(3,1)-1: mb(3,2),1:nfcv) &
  = info%v3c(     1:cmx(1)  ,        1:cmx(2)  ,        0:cmx(3)  ,1:nfcv)
!
CASE DEFAULT
  PRINT *, 'BSAM 2.0 RestrictFCSolution: Only ndims = 2,3 are supported.'
  STOP
END SELECT
!

END FUNCTION RestrictFCSolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE GetCoarseLevelLoading(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
IMPLICIT NONE
!
INTEGER, INTENT(IN):: level
!
TYPE(funcparam):: dummy
!
! 1) Calculate the standard coarse load function in AFAS:
!
! Cell centered variables:
CALL ApplyOnLevel(level,CoarseLoadingFunctionCC,dummy)
!
! Face Centered Variables: Parents must pull up the required values, then
! synchronize data with other parent-level grids: 
CALL ApplyOnLevel(level-1,PullFCLoadings,dummy)
!
! 2) Correct the CC coarse load functions to preserve CC mass:
fluxbalancing = .TRUE.
IF(fluxbalancing .AND. level > rootlevel) THEN
  ALLOCATE(zerothlayer)
  NULLIFY(zerothlayer%prevlayer)
!
  lastlayer => zerothlayer
  nmasscorrlayers = 0
!
  CALL ApplyOnLevel(level,CoarseLoadingCorrection,dummy)
!
  IF(nmasscorrlayers > 0) THEN
    CALL TransferMassCorrectionLayers(level-1)
    CALL DeleteMassCorrectionList
  END IF
!
  NULLIFY(lastlayer)
  NULLIFY(currentlayer)
  DEALLOCATE(zerothlayer)
END IF 
!
END SUBROUTINE GetCoarseLevelLoading
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION CoarseLoadingFunctionCC(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY: Restriction2D, Restriction3D
USE Problem, ONLY: Operator2D, Operator3D
IMPLICIT NONE
!
! Dimensionally invariant coarse loading function routine.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror, nnaxv
INTEGER, DIMENSION(1:maxdims):: cmx, mx
INTEGER, DIMENSION(1:maxdims,1:2):: amb, mb
REAL(KIND=r8):: ch
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
CoarseLoadingFunctionCC = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
ierror = GetParentInfo(parent)
ch = parent%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
nnaxv = MAX(naxv,1)
amb = 1; IF(naxv>0) amb(1:ndims,1:2) = mb(1:ndims,1:2)
!
CALL Residual(info)
!
SELECT CASE(ndims)
  CASE(2)
    parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),1,1:nccv) &
      = Restriction2D(info%rf(1: mx(1)  ,1: mx(2)  ,1,1:nccv)) &
      +    Operator2D(info%qc(0:cmx(1)+1,0:cmx(2)+1,1,1:nccv), &
                   parent%qold(   mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nccv), &
                   parent%v1(     mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nfcv), &
                   parent%v1old(  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nfcv), &
                   parent%v2   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1,1,1: nfcv), &
                   parent%v2old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1,1,1: nfcv), &
                   parent%aux  ( amb(1,1)-1:amb(1,2)+1, &
                                 amb(2,1)-1:amb(2,2)+1,1,1:nnaxv), &
                   cmx(1:2),ch,xlower(1:2))
  CASE(3)
    parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nccv) &
      = Restriction3D(info%rf(1: mx(1)  ,1: mx(2)  ,1: mx(3)  ,1:nccv)) &
      +    Operator3D(info%qc(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nccv), &
                   parent%qold(   mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nccv), &
                   parent%v1(     mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v1old(  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v2   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v2old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v3   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-2: mb(3,2)+1,1: nfcv), &
                   parent%v3old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-2: mb(3,2)+1,1: nfcv), &
                   parent%aux(amb(1,1)-1:amb(1,2)+1, &
                              amb(2,1)-1:amb(2,2)+1, &
                              amb(3,1)-1:amb(3,2)+1,1:nnaxv), &
                   cmx(1:3),ch,xlower(1:3))
  CASE DEFAULT
    PRINT *, 'BSAM 2.0 CoarseLoadingFunctionCC:'
    PRINT *, '         Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION CoarseLoadingFunctionCC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION CoarseLoadingCorrection(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE Problem, ONLY: Mob, MobD
USE GridUtilities, ONLY: Restriction2D, Restriction3D
USE ProblemDef
IMPLICIT NONE
!
! Dimensionally invariant coarse loading correction function routine.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
LOGICAL:: masscorrneeded
INTEGER:: ierror, i, j, k, level, nnaxv
INTEGER, DIMENSION(1:maxdims):: cmx, mx, pmx
INTEGER, DIMENSION(1:maxdims,1:2):: amb, mb, mg
REAL(KIND=r8):: ch, ch2
REAL(KIND=r8), DIMENSION(1:maxnccv):: diffcons
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
REAL(KIND=r8), DIMENSION(:,:,:,:), ALLOCATABLE:: zombiediff1, zombiediff2, zombiediff3
!:
CoarseLoadingCorrection = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
ierror = GetParentInfo(parent)
pmx = 1; pmx(1:ndims) = parent%mx(1:ndims)
ch = parent%dx(1); ch2 = ch*ch
level = info%level
xlower(1:ndims) = info%xlower(1:ndims)
!
nnaxv = MAX(naxv,1)
amb = 1; IF(naxv>0) amb(1:ndims,1:2) = mb(1:ndims,1:2)
!
! These need to come from the user. This is only for testing purposes:
diffcons(1) = 0.5_r8*dt
diffcons(2) = -epep*epep
!
SELECT CASE(ndims)
CASE(2)
  ALLOCATE(zombiediff1(1:2     ,1:cmx(2),1:1,1:3), &
           zombiediff2(1:cmx(1),1:2     ,1:1,1:3))
  zombiediff1 = 0.0_r8; zombiediff2 = 0.0_r8
!
!! Bottom correction:
  masscorrneeded = .FALSE.
  DO i = 1, cmx(1)
    IF(info%levellandscape(i,0,1) == level-1) THEN
      zombiediff2(i,1,1,2) &
        = (          (info%q(2*i  ,1,1,1)  -     info%q(2*i  ,0,1,1))   &
        +            (info%q(2*i-1,1,1,1)  -     info%q(2*i-1,0,1,1))   &    
        +  parent%q(mb(1,1)+i-1,mb(2,1)-1,1,1) &
        -  info%qcold(i,1,1,1))/ch2 
!         
      zombiediff2(i,1,1,1) &
!! 1/2 Flux^{n+1} !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        = (          (info%q(2*i  ,1,1,2)  -     info%q(2*i  ,0,1,2))   &
        * 0.5_r8*(Mob(info%q(2*i  ,1,1,1)) + Mob(info%q(2*i  ,0,1,1)))  &
!
        +            (info%q(2*i-1,1,1,2)  -     info%q(2*i-1,0,1,2))   &
        * 0.5_r8*(Mob(info%q(2*i-1,1,1,1)) + Mob(info%q(2*i-1,0,1,1)))  &
!
        - 0.5_r8*(Mob(info%qc(i,1,1,1))                                 &
                + Mob(parent%q(mb(1,1)+i-1,mb(2,1)-1,1,1)))             &
        *  (info%qc(i,1,1,2)                                            &
        -  parent%q(mb(1,1)+i-1,mb(2,1)-1,1,2)))/ch2                    &
!!! 1/2 Flux^{n}   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        + (          (info%qold(2*i  ,1,1,2)  -     info%qold(2*i  ,0,1,2))   &
        * 0.5_r8*(Mob(info%qold(2*i  ,1,1,1)) + Mob(info%qold(2*i  ,0,1,1)))  &
!
        +            (info%qold(2*i-1,1,1,2)  -     info%qold(2*i-1,0,1,2))   &
        * 0.5_r8*(Mob(info%qold(2*i-1,1,1,1)) + Mob(info%qold(2*i-1,0,1,1)))  &
!
        - 0.5_r8*(Mob(info%qcold(i,1,1,1))                                 &
                + Mob(parent%qold(mb(1,1)+i-1,mb(2,1)-1,1,1)))             &
        *  (info%qcold(i,1,1,2)                                            &
        -  parent%qold(mb(1,1)+i-1,mb(2,1)-1,1,2)))/ch2
!


      zombiediff2(i,1,1,3) &
        = (          (info%q(2*i  ,1,1,3)  -     info%q(2*i  ,0,1,3))   &
        *0.5_r8*(MobD(info%q(2*i  ,1,1,1)) +MobD(info%q(2*i  ,0,1,1)))  &
!
        +            (info%q(2*i-1,1,1,3)  -     info%q(2*i-1,0,1,3))   &
        *0.5_r8*(MobD(info%q(2*i-1,1,1,1)) +MobD(info%q(2*i-1,0,1,1)))  &
!      
        - 0.5_r8*(MobD(info%qc(i,1,1,1)) + MobD(parent%q(mb(1,1)+i-1,mb(2,1)-1,1,1)))     &
                *(     info%qc(i,1,1,3)  -      parent%q(mb(1,1)+i-1,mb(2,1)-1,1,3)))/ch2        
!
      zombiediff2(i,1,1,1) = zombiediff2(i,1,1,1)*diffcons(1)
!
      zombiediff2(i,1,1,2) = zombiediff2(i,1,1,2)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  IF(masscorrneeded) THEN
    IF(mb(2,1)-1 >= 1) THEN
!    
        parent%f(mb(1,1):mb(1,2),mb(2,1)-1,1,1:3) &
      = parent%f(mb(1,1):mb(1,2),mb(2,1)-1,1,1:3) &
      + zombiediff2(1:cmx(1),1,1,1:3)
!
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)+mb(1,1:2)-1
      mg(2,1:2) = parent%mglobal(2,1)-1
      CALL AddStripToList(mg,zombiediff2(1:cmx(1),1:1,1:1,1:3))
    END IF
  END IF
!
!! Top correction:
  masscorrneeded = .FALSE.
  DO i = 1, cmx(1)
    IF(info%levellandscape(i,cmx(2)+1,1) == level-1) THEN
!
          zombiediff2(i,2,1,2) &
        = (       (info%q(2*i  , mx(2)  ,1,1)-info%q(2*i  ,mx(2)+1  ,1,1)) &
        +         (info%q(2*i-1, mx(2)  ,1,1)-info%q(2*i-1,mx(2)+1  ,1,1)) &    
        +  parent%q(mb(1,1)+i-1,mb(2,2)+1,1,1) &
        -  info%qc(i,cmx(2),1,1))/ch2
!        
      zombiediff2(i,2,1,1) &
       = (       (     info%q(2*i  ,mx(2)  ,1,2)  -     info%q(2*i  ,mx(2)+1,1,2))  &
        *  0.5_r8*(Mob(info%q(2*i  ,mx(2)  ,1,1)) + Mob(info%q(2*i  ,mx(2)+1,1,1))) &
!
        +        (     info%q(2*i-1,mx(2)  ,1,2)  -     info%q(2*i-1,mx(2)+1,1,2))  &
        *  0.5_r8*(Mob(info%q(2*i-1,mx(2)  ,1,1)) + Mob(info%q(2*i-1,mx(2)+1,1,1))) &
!      
        -  0.5_r8*(Mob(info%qc(i,cmx(2),1,1)) + Mob(parent%q(mb(1,1)+i-1,mb(2,2)+1,1,1)))          &
        *             (info%qc(i,cmx(2),1,2)  -     parent%q(mb(1,1)+i-1,mb(2,2)+1,1,2)))/ch2 &
!!  
        + (       (    info%qold(2*i  ,mx(2)  ,1,2)  -     info%qold(2*i  ,mx(2)+1,1,2))  &
        *  0.5_r8*(Mob(info%qold(2*i  ,mx(2)  ,1,1)) + Mob(info%qold(2*i  ,mx(2)+1,1,1))) &
!
        +        (     info%qold(2*i-1,mx(2)  ,1,2)  -     info%qold(2*i-1,mx(2)+1,1,2))  &
        *  0.5_r8*(Mob(info%qold(2*i-1,mx(2)  ,1,1)) + Mob(info%qold(2*i-1,mx(2)+1,1,1))) &
!      
        -  0.5_r8*(Mob(info%qcold(i,cmx(2),1,1)) + Mob(parent%qold(mb(1,1)+i-1,mb(2,2)+1,1,1)))          &
        *             (info%qcold(i,cmx(2),1,2)  -     parent%qold(mb(1,1)+i-1,mb(2,2)+1,1,2)))/ch2
!  
		zombiediff2(i,2,1,3) &
       = (        (     info%q(2*i  ,mx(2)  ,1,3)  -      info%q(2*i ,mx(2) +1 ,1,3))  &
        *  0.5_r8*(MobD(info%q(2*i  ,mx(2)  ,1,1)) + MobD(info%q(2*i ,mx(2) +1 ,1,1))) &
!
        +         (     info%q(2*i-1,mx(2)  ,1,3)  -      info%q(2*i-1,mx(2) +1  ,1,3))  &
        *  0.5_r8*(MobD(info%q(2*i-1,mx(2)  ,1,1)) + MobD(info%q(2*i-1,mx(2) +1  ,1,1))) &
!
        -  0.5_r8*(MobD(info%qc(i,cmx(2),1,1)) + MobD(parent%q(mb(1,1)+i-1,mb(2,2)+1,1,1)))          &
        *  (            info%qc(i,cmx(2),1,3)  -      parent%q(mb(1,1)+i-1,mb(2,2)+1,1,3)))/ch2
!   
      zombiediff2(i,2,1,1) = zombiediff2(i,2,1,1)*diffcons(1)
!
      zombiediff2(i,2,1,2) = zombiediff2(i,2,1,2)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  IF(masscorrneeded) THEN
    IF(mb(2,2)+1 <= pmx(2)) THEN
        parent%f(mb(1,1):mb(1,2),mb(2,2)+1,1,1:3) &
      = parent%f(mb(1,1):mb(1,2),mb(2,2)+1,1,1:3) &
      + zombiediff2(1:cmx(1),2,1,1:3)
!            
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)+mb(1,1:2)-1
      mg(2,1:2) = parent%mglobal(2,1)+pmx(2)
      CALL AddStripToList(mg,zombiediff2(1:cmx(1),2:2,1:1,1:3))
    END IF
  END IF
!!
!! Left Correction:
  masscorrneeded = .FALSE.
  DO j = 1, cmx(2)
    IF(info%levellandscape(0,j,1) == level-1) THEN
   
      zombiediff1(1,j,1,2) &
        =  (       (info%q(1,2*j  ,1,1)-info%q(0,2*j  ,1,1))        &      
        +          (info%q(1,2*j-1,1,1)-info%q(0,2*j-1,1,1))        &    
        +   parent%q(mb(1,1)-1,mb(2,1)+j-1,1,1)             &
        -    info%qc(1,j,1,1))  /ch2                             
!
       zombiediff1(1,j,1,1) &
        =  (          (info%q(1,2*j  ,1,2)  -     info%q(0,2*j  ,1,2))  &
		*  0.5_r8*(Mob(info%q(1,2*j  ,1,1)) + Mob(info%q(0,2*j  ,1,1))) &
!
        +             (info%q(1,2*j-1,1,2)	-     info%q(0,2*j-1,1,2))  &
		*  0.5_r8*(Mob(info%q(1,2*j-1,1,1)) + Mob(info%q(0,2*j-1,1,1))) &
!             
        -  0.5_r8*(Mob(info%qc(1,j,1,1)) + Mob(parent%q(mb(1,1)-1,mb(2,1)+j-1,1,1)))         	 &
        *             (info%qc(1,j,1,2)  -     parent%q(mb(1,1)-1,mb(2,1)+j-1,1,2)))/ch2		 &
        +  (          (info%qold(1,2*j  ,1,2)  -     info%qold(0,2*j  ,1,2))  &
		*  0.5_r8*(Mob(info%qold(1,2*j  ,1,1)) + Mob(info%qold(0,2*j  ,1,1))) &
!
        +             (info%qold(1,2*j-1,1,2)	-     info%qold(0,2*j-1,1,2))  &
		*  0.5_r8*(Mob(info%qold(1,2*j-1,1,1)) + Mob(info%qold(0,2*j-1,1,1))) &
!             
        -  0.5_r8*(Mob(info%qcold(1,j,1,1)) + Mob(parent%qold(mb(1,1)-1,mb(2,1)+j-1,1,1)))         	 &
        *             (info%qcold(1,j,1,2)  -     parent%qold(mb(1,1)-1,mb(2,1)+j-1,1,2)))/ch2							
       
 
	zombiediff1(1,j,1,3) &
        =  (           (info%q(1,2*j  ,1,3)       - info%q(0,2*j ,1,3))  &
		*  0.5_r8*(MobD(info%q(1,2*j  ,1,1)) + MobD(info%q(0,2*j ,1,1))) &
!
        +              (info%q(1,2*j-1,1,3)		  - info%q(0,2*j-1,1,3))        &
		*  0.5_r8*(MobD(info%q(1,2*j-1,1,1)) + MobD(info%q(0,2*j-1,1,1))) &
!             
        -  0.5_r8*(MobD(info%qc(1,j,1,1)) + MobD(parent%q(mb(1,1)-1,mb(2,1)+j-1,1,1)))          &
        *  (            info%qc(1,j,1,3)  -      parent%q(mb(1,1)-1,mb(2,1)+j-1,1,3)))/ch2									
!        
      zombiediff1(1,j,1,1) = zombiediff1(1,j,1,1)*diffcons(1)
!
      zombiediff1(1,j,1,2) = zombiediff1(1,j,1,2)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  IF(masscorrneeded) THEN
    IF(mb(1,1)-1 >= 1) THEN
        parent%f(mb(1,1)-1,mb(2,1):mb(2,2),1,1:3) &
      = parent%f(mb(1,1)-1,mb(2,1):mb(2,2),1,1:3) &
      + zombiediff1(1,1:cmx(2),1,1:3)
!
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)-1
      mg(2,1:2) = parent%mglobal(2,1)+mb(2,1:2)-1
      CALL AddStripToList(mg,zombiediff1(1:1,1:cmx(2),1:1,1:3))
    END IF
  END IF
!!!!
!!! Right Correction:
  masscorrneeded = .FALSE.
  DO j = 1, cmx(2)
    IF(info%levellandscape(cmx(1)+1,j,1) == level-1) THEN

      zombiediff1(2,j,1,2) &
        = ( (info%q(mx(1),2*j  ,1,1)-info%q(mx(1)+1,2*j  ,1,1))        &      
        +   (info%q(mx(1),2*j-1,1,1)-info%q(mx(1)+1,2*j-1,1,1))        &    
        +   parent%q(mb(1,2)+1,mb(2,1)+j-1,1,1)             				&
        -    info%qc(cmx(1),j,1,1))  /ch2

	zombiediff1(2,j,1,1) &
        =  (          (info%q(mx(1),2*j  ,1,2)     - info%q(mx(1)+1,2*j ,1,2))        &
		*  0.5_r8*(Mob(info%q(mx(1),2*j  ,1,1))+ Mob(info%q(mx(1)+1,2*j ,1,1))) &
! 
        +             (info%q(mx(1),2*j-1,1,2)	   - info%q(mx(1)+1,2*j-1,1,2))        &
		*  0.5_r8*(Mob(info%q(mx(1),2*j-1,1,1))+ Mob(info%q(mx(1)+1,2*j-1,1,1))) &
!             
        -  0.5_r8*(Mob(info%qc(cmx(1),j,1,1)) + Mob(parent%q(mb(1,2)+1,mb(2,1)+j-1,1,1)))          &
        *             (info%qc(cmx(1),j,1,2)  -     parent%q(mb(1,2)+1,mb(2,1)+j-1,1,2)))/ch2      &
! 
        +  (          (info%qold(mx(1),2*j  ,1,2)     - info%qold(mx(1)+1,2*j ,1,2))        &
		*  0.5_r8*(Mob(info%qold(mx(1),2*j  ,1,1))+ Mob(info%qold(mx(1)+1,2*j ,1,1))) &
! 
        +             (info%qold(mx(1),2*j-1,1,2)	   - info%qold(mx(1)+1,2*j-1,1,2))        &
		*  0.5_r8*(Mob(info%qold(mx(1),2*j-1,1,1))+ Mob(info%qold(mx(1)+1,2*j-1,1,1))) &
!             
        -  0.5_r8*(Mob(info%qcold(cmx(1),j,1,1)) + Mob(parent%qold(mb(1,2)+1,mb(2,1)+j-1,1,1)))          &
        *             (info%qcold(cmx(1),j,1,2)  -     parent%qold(mb(1,2)+1,mb(2,1)+j-1,1,2)))/ch2
! 
         zombiediff1(2,j,1,3) &
        =  (           (info%q(mx(1),2*j  ,1,3)       - info%q(mx(1)+1,2*j  ,1,3))        &
		*  0.5_r8*(MobD(info%q(mx(1),2*j  ,1,1)) + MobD(info%q(mx(1)+1,2*j  ,1,1))) &
 
        +              (info%q(mx(1),2*j-1,1,3)		  - info%q(mx(1)+1,2*j-1,1,3))        &
		*  0.5_r8*(MobD(info%q(mx(1)  ,2*j-1   ,1,1)) + MobD(info%q(mx(1)+1  ,2*j-1 ,1,1))) &
!             
        -  0.5_r8*(MobD(info%qc(cmx(1),j,1,1)) + MobD(parent%q(mb(1,2)+1,mb(2,1)+j-1,1,1)))          &
        *              (info%qc(cmx(1),j,1,3) -       parent%q(mb(1,2)+1,mb(2,1)+j-1,1,3)))/ch2  
!
      zombiediff1(2,j,1,1) = zombiediff1(2,j,1,1)*diffcons(1)
!
      zombiediff1(2,j,1,2) = zombiediff1(2,j,1,2)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  IF(masscorrneeded) THEN
    IF(mb(1,2)+1 <= pmx(1)) THEN
        parent%f(mb(1,2)+1,mb(2,1):mb(2,2),1,1:3) &
      = parent%f(mb(1,2)+1,mb(2,1):mb(2,2),1,1:3) &
      + zombiediff1(2,1:cmx(2),1,1:3)
!
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)+pmx(1)
      mg(2,1:2) = parent%mglobal(2,1)+mb(2,1:2)-1
      CALL AddStripToList(mg,zombiediff1(2:2,1:cmx(2),1:1,1:3))
    END IF
  END IF
!  
  DEALLOCATE(zombiediff1,zombiediff2)
!
CASE(3)
!
  ALLOCATE(zombiediff1(1:2     ,1:cmx(2),1:cmx(3),1:nccv), &
           zombiediff2(1:cmx(1),1:cmx(2),1:2     ,1:nccv), &
           zombiediff3(1:cmx(1),1:2     ,1:cmx(3),1:nccv))
!           
  zombiediff1 = 0.0_r8; zombiediff2 = 0.0_r8; zombiediff3 = 0.0_r8
!
! Bottom correction:
  masscorrneeded = .FALSE.
  DO i = 1, cmx(1)
  DO j = 1, cmx(2)
    IF(info%levellandscape(i,j,0) == level-1) THEN
      zombiediff2(i,j,1,2) &
        = (0.5_r8*(info%q(2*i  ,2*j  ,1,1)-info%q(2*i  ,2*j  ,0,1)) &
        +  0.5_r8*(info%q(2*i-1,2*j  ,1,1)-info%q(2*i-1,2*j  ,0,1)) &
        +  0.5_r8*(info%q(2*i  ,2*j-1,1,1)-info%q(2*i  ,2*j-1,0,1)) &
        +  0.5_r8*(info%q(2*i-1,2*j-1,1,1)-info%q(2*i-1,2*j-1,0,1)) &
        +  parent%q(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,1)-1,1) &
        -  info%qc(i,j,1,1))/ch2
!
      zombiediff2(i,j,1,1) &
        = (0.5_r8*(    info%q(2*i  ,2*j  ,1,2)  -     info%q(2*i  ,2*j  ,0,2))  &
        *  0.5_r8*(Mob(info%q(2*i  ,2*j  ,1,1)) + Mob(info%q(2*i  ,2*j  ,0,1))) &
!
        +  0.5_r8*(    info%q(2*i-1,2*j  ,1,2)  -     info%q(2*i-1,2*j  ,0,2))  &
        *  0.5_r8*(Mob(info%q(2*i-1,2*j  ,1,1)) + Mob(info%q(2*i-1,2*j  ,0,1))) &
!
        +  0.5_r8*(    info%q(2*i  ,2*j-1,1,2)  -     info%q(2*i  ,2*j-1,0,2))  &
        *  0.5_r8*(Mob(info%q(2*i  ,2*j-1,1,1)) + Mob(info%q(2*i  ,2*j-1,0,1))) &
!
        +  0.5_r8*(    info%q(2*i-1,2*j-1,1,2)  -     info%q(2*i-1,2*j-1,0,2))  &
        *  0.5_r8*(Mob(info%q(2*i-1,2*j-1,1,1)) + Mob(info%q(2*i-1,2*j-1,0,1))) &
        -  0.5_r8*(Mob(info%qc(i,j,1,1))                                        &
                 + Mob(parent%q(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,1)-1,1)))          &
        *  (info%qc(i,j,1,2)                                                    &
        -  parent%q(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,1)-1,2)))/ch2                  &
!
        + (0.5_r8*(    info%qold(2*i  ,2*j  ,1,2)  -     info%qold(2*i  ,2*j  ,0,2))  &
        *  0.5_r8*(Mob(info%qold(2*i  ,2*j  ,1,1)) + Mob(info%qold(2*i  ,2*j  ,0,1))) &
!
        +  0.5_r8*(    info%qold(2*i-1,2*j  ,1,2)  -     info%qold(2*i-1,2*j  ,0,2))  &
        *  0.5_r8*(Mob(info%qold(2*i-1,2*j  ,1,1)) + Mob(info%qold(2*i-1,2*j  ,0,1))) &
!
        +  0.5_r8*(    info%qold(2*i  ,2*j-1,1,2)  -     info%qold(2*i  ,2*j-1,0,2))  &
        *  0.5_r8*(Mob(info%qold(2*i  ,2*j-1,1,1)) + Mob(info%qold(2*i  ,2*j-1,0,1))) &
!
        +  0.5_r8*(    info%qold(2*i-1,2*j-1,1,2)  -     info%qold(2*i-1,2*j-1,0,2))  &
        *  0.5_r8*(Mob(info%qold(2*i-1,2*j-1,1,1)) + Mob(info%qold(2*i-1,2*j-1,0,1))) &
        -  0.5_r8*(Mob(info%qcold(i,j,1,1))                                        &
                 + Mob(parent%qold(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,1)-1,1)))          &
        *  (info%qcold(i,j,1,2)                                                    &
        -  parent%qold(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,1)-1,2)))/ch2
!


      zombiediff2(i,j,1,3) &
        = (0.5_r8*(     info%q(2*i  ,2*j  ,1,3)  -      info%q(2*i  ,2*j  ,0,3))  &
        *  0.5_r8*(MobD(info%q(2*i  ,2*j  ,1,1)) + MobD(info%q(2*i  ,2*j  ,0,1))) &
!
        +  0.5_r8*(     info%q(2*i-1,2*j  ,1,3)  -      info%q(2*i-1,2*j  ,0,3))  &
        *  0.5_r8*(MobD(info%q(2*i-1,2*j  ,1,1)) + MobD(info%q(2*i-1,2*j  ,0,1))) &
!
        +  0.5_r8*(     info%q(2*i  ,2*j-1,1,3)  -      info%q(2*i  ,2*j-1,0,3))  &
        *  0.5_r8*(MobD(info%q(2*i  ,2*j-1,1,1)) + MobD(info%q(2*i  ,2*j-1,0,1))) &
!
        +  0.5_r8*(     info%q(2*i-1,2*j-1,1,3)  -      info%q(2*i-1,2*j-1,0,3))  &
        *  0.5_r8*(MobD(info%q(2*i-1,2*j-1,1,1)) + MobD(info%q(2*i-1,2*j-1,0,1))) &
        -  0.5_r8*(MobD(info%qc(i,j,1,1))                                         &
                 + MobD(parent%q(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,1)-1,1)))           &
        *  (info%qc(i,j,1,3)                                                    &
        -  parent%q(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,1)-1,3)))/ch2
!
      zombiediff2(i,j,1,1) = zombiediff2(i,j,1,1)*diffcons(1)
      zombiediff2(i,j,1,2) = zombiediff2(i,j,1,2)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  END DO
  IF(masscorrneeded) THEN
    IF(mb(3,1)-1 >= 1) THEN
        parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1)-1,1:nccv) &
      = parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1)-1,1:nccv) &
      + zombiediff2(1:cmx(1),1:cmx(2),1,1:nccv)
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)+mb(1,1:2)-1
      mg(2,1:2) = parent%mglobal(2,1)+mb(2,1:2)-1
      mg(3,1:2) = parent%mglobal(3,1)-1
      CALL AddStripToList(mg,zombiediff2(1:cmx(1),1:cmx(2),1:1,1:nccv))
    END IF
  END IF
!!
! Top correction:
  masscorrneeded = .FALSE.
  DO i = 1, cmx(1)
  DO j = 1, cmx(2)
    IF(info%levellandscape(i,j,cmx(3)+1) == level-1) THEN
      zombiediff2(i,j,2,2) &
        = (0.5_r8*(info%q(2*i  ,2*j  ,mx(3),1)-info%q(2*i  ,2*j  ,mx(3)+1,1)) &
        +  0.5_r8*(info%q(2*i-1,2*j  ,mx(3),1)-info%q(2*i-1,2*j  ,mx(3)+1,1)) &
        +  0.5_r8*(info%q(2*i  ,2*j-1,mx(3),1)-info%q(2*i  ,2*j-1,mx(3)+1,1)) &
        +  0.5_r8*(info%q(2*i-1,2*j-1,mx(3),1)-info%q(2*i-1,2*j-1,mx(3)+1,1)) &
        +  parent%q(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,2)+1,1) &
        -  info%qc(i,j,cmx(3),1))/ch2
!
      zombiediff2(i,j,2,1) &
        = (0.5_r8*(    info%q(2*i  ,2*j  ,mx(3),2)  -     info%q(2*i  ,2*j  ,mx(3)+1,2))   &
        *  0.5_r8*(Mob(info%q(2*i  ,2*j  ,mx(3),1)) + Mob(info%q(2*i  ,2*j  ,mx(3)+1,1)))  &
        +  0.5_r8*(    info%q(2*i-1,2*j  ,mx(3),2)  -     info%q(2*i-1,2*j  ,mx(3)+1,2))   &
        *  0.5_r8*(Mob(info%q(2*i-1,2*j  ,mx(3),1)) + Mob(info%q(2*i-1,2*j  ,mx(3)+1,1)))  &
        +  0.5_r8*(    info%q(2*i  ,2*j-1,mx(3),2)  -     info%q(2*i  ,2*j-1,mx(3)+1,2))   &
        *  0.5_r8*(Mob(info%q(2*i  ,2*j-1,mx(3),1)) + Mob(info%q(2*i  ,2*j-1,mx(3)+1,1)))  &
        +  0.5_r8*(    info%q(2*i-1,2*j-1,mx(3),2)  -     info%q(2*i-1,2*j-1,mx(3)+1,2))   &
        *  0.5_r8*(Mob(info%q(2*i-1,2*j-1,mx(3),1)) + Mob(info%q(2*i-1,2*j-1,mx(3)+1,1)))  &
        +  0.5_r8*(Mob(parent%q(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,2)+1,1))                      &
        +          Mob(info%qc(i,j,cmx(3),1)))                                             &
        * (parent%q(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,2)+1,2)                                   &
        -              info%qc(i,j,cmx(3),2)))/ch2     &
        + (0.5_r8*(    info%qold(2*i  ,2*j  ,mx(3),2)  -     info%qold(2*i  ,2*j  ,mx(3)+1,2))   &
        *  0.5_r8*(Mob(info%qold(2*i  ,2*j  ,mx(3),1)) + Mob(info%qold(2*i  ,2*j  ,mx(3)+1,1)))  &
        +  0.5_r8*(    info%qold(2*i-1,2*j  ,mx(3),2)  -     info%qold(2*i-1,2*j  ,mx(3)+1,2))   &
        *  0.5_r8*(Mob(info%qold(2*i-1,2*j  ,mx(3),1)) + Mob(info%qold(2*i-1,2*j  ,mx(3)+1,1)))  &
        +  0.5_r8*(    info%qold(2*i  ,2*j-1,mx(3),2)  -     info%qold(2*i  ,2*j-1,mx(3)+1,2))   &
        *  0.5_r8*(Mob(info%qold(2*i  ,2*j-1,mx(3),1)) + Mob(info%qold(2*i  ,2*j-1,mx(3)+1,1)))  &
        +  0.5_r8*(    info%qold(2*i-1,2*j-1,mx(3),2)  -     info%qold(2*i-1,2*j-1,mx(3)+1,2))   &
        *  0.5_r8*(Mob(info%qold(2*i-1,2*j-1,mx(3),1)) + Mob(info%qold(2*i-1,2*j-1,mx(3)+1,1)))  &
        +  0.5_r8*(Mob(parent%qold(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,2)+1,1))                      &
        +          Mob(info%qcold(i,j,cmx(3),1)))                                             &
        * (parent%qold(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,2)+1,2)                                   &
        -              info%qcold(i,j,cmx(3),2)))/ch2

!
      zombiediff2(i,j,2,3) &
        = (0.5_r8*(     info%q(2*i  ,2*j  ,mx(3),3)  -      info%q(2*i  ,2*j  ,mx(3)+1,3))   &
        *  0.5_r8*(MobD(info%q(2*i  ,2*j  ,mx(3),1)) + MobD(info%q(2*i  ,2*j  ,mx(3)+1,1)))  &
        +  0.5_r8*(     info%q(2*i-1,2*j  ,mx(3),3)  -      info%q(2*i-1,2*j  ,mx(3)+1,3))   &
        *  0.5_r8*(MobD(info%q(2*i-1,2*j  ,mx(3),1)) + MobD(info%q(2*i-1,2*j  ,mx(3)+1,1)))  &
        +  0.5_r8*(     info%q(2*i  ,2*j-1,mx(3),3)  -      info%q(2*i  ,2*j-1,mx(3)+1,3))   &
        *  0.5_r8*(MobD(info%q(2*i  ,2*j-1,mx(3),1)) + MobD(info%q(2*i  ,2*j-1,mx(3)+1,1)))  &
        +  0.5_r8*(     info%q(2*i-1,2*j-1,mx(3),3)  -      info%q(2*i-1,2*j-1,mx(3)+1,3))   &
        *  0.5_r8*(MobD(info%q(2*i-1,2*j-1,mx(3),1)) + MobD(info%q(2*i-1,2*j-1,mx(3)+1,1)))  &
        +  0.5_r8*(MobD(parent%q(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,2)+1,1))                       &
        +          MobD(info%qc(i,j,cmx(3),1)))                                              &
        * (parent%q(mb(1,1)+i-1,mb(2,1)+j-1,mb(3,2)+1,3)                                     &
        -              info%qc(i,j,cmx(3),3)))/ch2
!
      zombiediff2(i,j,2,1) = zombiediff2(i,j,2,1)*diffcons(1)
      zombiediff2(i,j,2,2) = zombiediff2(i,j,2,2)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  END DO
  IF(masscorrneeded) THEN
    IF(mb(3,2)+1 <= pmx(3)) THEN
        parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,2)+1,1:nccv) &
      = parent%f(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,2)+1,1:nccv) &
      + zombiediff2(1:cmx(1),1:cmx(2),2,1:nccv)
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)+mb(1,1:2)-1
      mg(2,1:2) = parent%mglobal(2,1)+mb(2,1:2)-1
      mg(3,1:2) = parent%mglobal(3,1)+pmx(3)
      CALL AddStripToList(mg,zombiediff2(1:cmx(1),1:cmx(2),2:2,1:nccv))
    END IF
  END IF
!!
! Left Correction:
  masscorrneeded = .FALSE.
  DO j = 1, cmx(2)
  DO k = 1, cmx(3)
    IF(info%levellandscape(0,j,k) == level-1) THEN
      zombiediff1(1,j,k,2) &
        = (0.5_r8*(info%q(1,2*j  ,2*k  ,1)-info%q(0,2*j  ,2*k  ,1)) &
        +  0.5_r8*(info%q(1,2*j-1,2*k  ,1)-info%q(0,2*j-1,2*k  ,1)) &
        +  0.5_r8*(info%q(1,2*j  ,2*k-1,1)-info%q(0,2*j  ,2*k-1,1)) &
        +  0.5_r8*(info%q(1,2*j-1,2*k-1,1)-info%q(0,2*j-1,2*k-1,1)) &
        +  parent%q(mb(1,1)-1,mb(2,1)+j-1,mb(3,1)+k-1,1) &
        -  info%qc(1,j,k,1))/ch2
!
      zombiediff1(1,j,k,1) &
        = (0.5_r8*(    info%q(1,2*j  ,2*k  ,2)  -     info%q(0,2*j  ,2*k  ,2))   &
        *  0.5_r8*(Mob(info%q(1,2*j  ,2*k  ,1)) + Mob(info%q(0,2*j  ,2*k  ,1)))  &
        +  0.5_r8*(    info%q(1,2*j-1,2*k  ,2)  -     info%q(0,2*j-1,2*k  ,2))   &
        *  0.5_r8*(Mob(info%q(1,2*j-1,2*k  ,1)) + Mob(info%q(0,2*j-1,2*k  ,1)))  &
        +  0.5_r8*(    info%q(1,2*j  ,2*k-1,2)  -     info%q(0,2*j  ,2*k-1,2))   &
        *  0.5_r8*(Mob(info%q(1,2*j  ,2*k-1,1)) + Mob(info%q(0,2*j  ,2*k-1,1)))  &
        +  0.5_r8*(    info%q(1,2*j-1,2*k-1,2)  -     info%q(0,2*j-1,2*k-1,2))   &
        *  0.5_r8*(Mob(info%q(1,2*j-1,2*k-1,1)) + Mob(info%q(0,2*j-1,2*k-1,1)))  &
        + (parent%q(mb(1,1)-1,mb(2,1)+j-1,mb(3,1)+k-1,2)                         &
        -  info%qc(1,j,k,2))                                                     &
        *0.5_r8*(Mob(parent%q(mb(1,1)-1,mb(2,1)+j-1,mb(3,1)+k-1,1))              &
        +  Mob(info%qc(1,j,k,1))))/ch2                                           &
        + (0.5_r8*(    info%qold(1,2*j  ,2*k  ,2)  -     info%qold(0,2*j  ,2*k  ,2))   &
        *  0.5_r8*(Mob(info%qold(1,2*j  ,2*k  ,1)) + Mob(info%qold(0,2*j  ,2*k  ,1)))  &
        +  0.5_r8*(    info%qold(1,2*j-1,2*k  ,2)  -     info%qold(0,2*j-1,2*k  ,2))   &
        *  0.5_r8*(Mob(info%qold(1,2*j-1,2*k  ,1)) + Mob(info%qold(0,2*j-1,2*k  ,1)))  &
        +  0.5_r8*(    info%qold(1,2*j  ,2*k-1,2)  -     info%qold(0,2*j  ,2*k-1,2))   &
        *  0.5_r8*(Mob(info%qold(1,2*j  ,2*k-1,1)) + Mob(info%qold(0,2*j  ,2*k-1,1)))  &
        +  0.5_r8*(    info%qold(1,2*j-1,2*k-1,2)  -     info%qold(0,2*j-1,2*k-1,2))   &
        *  0.5_r8*(Mob(info%qold(1,2*j-1,2*k-1,1)) + Mob(info%qold(0,2*j-1,2*k-1,1)))  &
        + (parent%qold(mb(1,1)-1,mb(2,1)+j-1,mb(3,1)+k-1,2)                         &
        -  info%qcold(1,j,k,2))                                                     &
        *0.5_r8*(Mob(parent%qold(mb(1,1)-1,mb(2,1)+j-1,mb(3,1)+k-1,1))              &
        +  Mob(info%qcold(1,j,k,1))))/ch2

!
      zombiediff1(1,j,k,3) &
        = (0.5_r8*(     info%q(1,2*j  ,2*k  ,3)  -      info%q(0,2*j  ,2*k  ,3))   &
        *  0.5_r8*(MobD(info%q(1,2*j  ,2*k  ,1)) + MobD(info%q(0,2*j  ,2*k  ,1)))  &
        +  0.5_r8*(     info%q(1,2*j-1,2*k  ,3)  -      info%q(0,2*j-1,2*k  ,3))   &
        *  0.5_r8*(MobD(info%q(1,2*j-1,2*k  ,1)) + MobD(info%q(0,2*j-1,2*k  ,1)))  &
        +  0.5_r8*(     info%q(1,2*j  ,2*k-1,3)  -      info%q(0,2*j  ,2*k-1,3))   &
        *  0.5_r8*(MobD(info%q(1,2*j  ,2*k-1,1)) + MobD(info%q(0,2*j  ,2*k-1,1)))  &
        +  0.5_r8*(     info%q(1,2*j-1,2*k-1,3)  -      info%q(0,2*j-1,2*k-1,3))   &
        *  0.5_r8*(MobD(info%q(1,2*j-1,2*k-1,1)) + MobD(info%q(0,2*j-1,2*k-1,1)))  &
        + (parent%q(mb(1,1)-1,mb(2,1)+j-1,mb(3,1)+k-1,3)                         &
        -  info%qc(1,j,k,3))                                                     &
        *0.5_r8*(MobD(parent%q(mb(1,1)-1,mb(2,1)+j-1,mb(3,1)+k-1,1))              &
        +  MobD(info%qc(1,j,k,1))))/ch2
!
      zombiediff1(1,j,k,1) = zombiediff1(1,j,k,1)*diffcons(1)
      zombiediff1(1,j,k,2) = zombiediff1(1,j,k,2)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  END DO
  IF(masscorrneeded) THEN
    IF(mb(1,1)-1 >= 1) THEN
        parent%f(mb(1,1)-1,mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nccv) &
      = parent%f(mb(1,1)-1,mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nccv) &
      + zombiediff1(1,1:cmx(2),1:cmx(3),1:nccv)
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)-1
      mg(2,1:2) = parent%mglobal(2,1)+mb(2,1:2)-1
      mg(3,1:2) = parent%mglobal(3,1)+mb(3,1:2)-1
      CALL AddStripToList(mg,zombiediff1(1:1,1:cmx(2),1:cmx(3),1:nccv))
    END IF
  END IF
!
! Right Correction:
  masscorrneeded = .FALSE.
  DO j = 1, cmx(2)
  DO k = 1, cmx(3)
    IF(info%levellandscape(cmx(1)+1,j,k) == level-1) THEN
      zombiediff1(2,j,k,2) &
        = (0.5_r8*(info%q(mx(1),2*j  ,2*k  ,1)-info%q(mx(1)+1,2*j  ,2*k  ,1)) &
        +  0.5_r8*(info%q(mx(1),2*j-1,2*k  ,1)-info%q(mx(1)+1,2*j-1,2*k  ,1)) &
        +  0.5_r8*(info%q(mx(1),2*j  ,2*k-1,1)-info%q(mx(1)+1,2*j  ,2*k-1,1)) &
        +  0.5_r8*(info%q(mx(1),2*j-1,2*k-1,1)-info%q(mx(1)+1,2*j-1,2*k-1,1)) &
        +  parent%q(mb(1,2)+1,mb(2,1)+j-1,mb(3,1)+k-1,1) &
        -  info%qc(cmx(1),j,k,1))/ch2
!
      zombiediff1(2,j,k,1) &
        = (0.5_r8*(    info%q(mx(1),2*j  ,2*k  ,2)  -     info%q(mx(1)+1,2*j  ,2*k  ,2))   &
        *  0.5_r8*(Mob(info%q(mx(1),2*j  ,2*k  ,1)) + Mob(info%q(mx(1)+1,2*j  ,2*k  ,1)))  &
        +  0.5_r8*(    info%q(mx(1),2*j-1,2*k  ,2)  -     info%q(mx(1)+1,2*j-1,2*k  ,2))   &
        *  0.5_r8*(Mob(info%q(mx(1),2*j-1,2*k  ,1)) + Mob(info%q(mx(1)+1,2*j-1,2*k  ,1)))  &
        +  0.5_r8*(    info%q(mx(1),2*j  ,2*k-1,2)  -     info%q(mx(1)+1,2*j  ,2*k-1,2))   &
        *  0.5_r8*(Mob(info%q(mx(1),2*j  ,2*k-1,1)) + Mob(info%q(mx(1)+1,2*j  ,2*k-1,1)))  &
        +  0.5_r8*(    info%q(mx(1),2*j-1,2*k-1,2)  -     info%q(mx(1)+1,2*j-1,2*k-1,2))   &
        *  0.5_r8*(Mob(info%q(mx(1),2*j-1,2*k-1,1)) + Mob(info%q(mx(1)+1,2*j-1,2*k-1,1)))  &
        +  (parent%q(mb(1,2)+1,mb(2,1)+j-1,mb(3,1)+k-1,2)                    &
        -              info%qc(cmx(1),j,k,2))                                &
        *  0.5_r8*(Mob(parent%q(mb(1,2)+1,mb(2,1)+j-1,mb(3,1)+k-1,1))                    &
        +          Mob(info%qc(cmx(1),j,k,1))))/ch2                          &
        + (0.5_r8*(    info%qold(mx(1),2*j  ,2*k  ,2)  -     info%qold(mx(1)+1,2*j  ,2*k  ,2))   &
        *  0.5_r8*(Mob(info%qold(mx(1),2*j  ,2*k  ,1)) + Mob(info%qold(mx(1)+1,2*j  ,2*k  ,1)))  &
        +  0.5_r8*(    info%qold(mx(1),2*j-1,2*k  ,2)  -     info%qold(mx(1)+1,2*j-1,2*k  ,2))   &
        *  0.5_r8*(Mob(info%qold(mx(1),2*j-1,2*k  ,1)) + Mob(info%qold(mx(1)+1,2*j-1,2*k  ,1)))  &
        +  0.5_r8*(    info%qold(mx(1),2*j  ,2*k-1,2)  -     info%qold(mx(1)+1,2*j  ,2*k-1,2))   &
        *  0.5_r8*(Mob(info%qold(mx(1),2*j  ,2*k-1,1)) + Mob(info%qold(mx(1)+1,2*j  ,2*k-1,1)))  &
        +  0.5_r8*(    info%qold(mx(1),2*j-1,2*k-1,2)  -     info%qold(mx(1)+1,2*j-1,2*k-1,2))   &
        *  0.5_r8*(Mob(info%qold(mx(1),2*j-1,2*k-1,1)) + Mob(info%qold(mx(1)+1,2*j-1,2*k-1,1)))  &
        +  (parent%qold(mb(1,2)+1,mb(2,1)+j-1,mb(3,1)+k-1,2)                    &
        -              info%qcold(cmx(1),j,k,2))                                &
        *  0.5_r8*(Mob(parent%qold(mb(1,2)+1,mb(2,1)+j-1,mb(3,1)+k-1,1))                    &
        +          Mob(info%qcold(cmx(1),j,k,1))))/ch2
!
      zombiediff1(2,j,k,3) &
        = (0.5_r8*(     info%q(mx(1),2*j  ,2*k  ,3)  -      info%q(mx(1)+1,2*j  ,2*k  ,3))   &
        *  0.5_r8*(MobD(info%q(mx(1),2*j  ,2*k  ,1)) + MobD(info%q(mx(1)+1,2*j  ,2*k  ,1)))  &
        +  0.5_r8*(     info%q(mx(1),2*j-1,2*k  ,3)  -      info%q(mx(1)+1,2*j-1,2*k  ,3))   &
        *  0.5_r8*(MobD(info%q(mx(1),2*j-1,2*k  ,1)) + MobD(info%q(mx(1)+1,2*j-1,2*k  ,1)))  &
        +  0.5_r8*(     info%q(mx(1),2*j  ,2*k-1,3)  -      info%q(mx(1)+1,2*j  ,2*k-1,3))   &
        *  0.5_r8*(MobD(info%q(mx(1),2*j  ,2*k-1,1)) + MobD(info%q(mx(1)+1,2*j  ,2*k-1,1)))  &
        +  0.5_r8*(     info%q(mx(1),2*j-1,2*k-1,3)  -      info%q(mx(1)+1,2*j-1,2*k-1,3))   &
        *  0.5_r8*(MobD(info%q(mx(1),2*j-1,2*k-1,1)) + MobD(info%q(mx(1)+1,2*j-1,2*k-1,1)))  &
        +  (parent%q(mb(1,2)+1,mb(2,1)+j-1,mb(3,1)+k-1,3)                    &
        -              info%qc(cmx(1),j,k,3))                                &
        *  0.5_r8*(MobD(parent%q(mb(1,2)+1,mb(2,1)+j-1,mb(3,1)+k-1,1))                    &
        +          MobD(info%qc(cmx(1),j,k,1))))/ch2
!
      zombiediff1(2,j,k,1) = zombiediff1(2,j,k,1)*diffcons(1)
      zombiediff1(2,j,k,2) = zombiediff1(2,j,k,2)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  END DO
  IF(masscorrneeded) THEN
    IF(mb(1,2)+1 <= pmx(1)) THEN
        parent%f(mb(1,2)+1,mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nccv) &
      = parent%f(mb(1,2)+1,mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nccv) &
      + zombiediff1(2,1:cmx(2),1:cmx(3),1:nccv)
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)+pmx(1)
      mg(2,1:2) = parent%mglobal(2,1)+mb(2,1:2)-1
      mg(3,1:2) = parent%mglobal(3,1)+mb(3,1:2)-1
      CALL AddStripToList(mg,zombiediff1(2:2,1:cmx(2),1:cmx(3),1:nccv))
    END IF
  END IF
!
! Front Correction:
  masscorrneeded = .FALSE.
  DO i = 1, cmx(1)
  DO k = 1, cmx(3)
    IF(info%levellandscape(i,0,k) == level-1) THEN
      zombiediff3(i,1,k,2) &
        = (0.5_r8*(info%q(2*i  ,1,2*k  ,1)-info%q(2*i  ,0,2*k  ,1)) &
        +  0.5_r8*(info%q(2*i-1,1,2*k  ,1)-info%q(2*i-1,0,2*k  ,1)) &
        +  0.5_r8*(info%q(2*i  ,1,2*k-1,1)-info%q(2*i  ,0,2*k-1,1)) &
        +  0.5_r8*(info%q(2*i-1,1,2*k-1,1)-info%q(2*i-1,0,2*k-1,1)) &
        +  parent%q(mb(1,1)+i-1,mb(2,1)-1,mb(3,1)+k-1,1) &
        -  info%qc(i,1,k,1))/ch2
!
      zombiediff3(i,1,k,1) &
        = (0.5_r8*(    info%q(2*i  ,1,2*k  ,2)  -     info%q(2*i  ,0,2*k  ,2))   &
        *  0.5_r8*(Mob(info%q(2*i  ,1,2*k  ,1)) + Mob(info%q(2*i  ,0,2*k  ,1)))  &
        +  0.5_r8*(    info%q(2*i-1,1,2*k  ,2)  -     info%q(2*i-1,0,2*k  ,2))   &
        *  0.5_r8*(Mob(info%q(2*i-1,1,2*k  ,1)) + Mob(info%q(2*i-1,0,2*k  ,1)))  &
        +  0.5_r8*(    info%q(2*i  ,1,2*k-1,2)  -     info%q(2*i  ,0,2*k-1,2))   &
        *  0.5_r8*(Mob(info%q(2*i  ,1,2*k-1,1)) + Mob(info%q(2*i  ,0,2*k-1,1)))  &
        +  0.5_r8*(    info%q(2*i-1,1,2*k-1,2)  -     info%q(2*i-1,0,2*k-1,2))   &
        *  0.5_r8*(Mob(info%q(2*i-1,1,2*k-1,1)) + Mob(info%q(2*i-1,0,2*k-1,1)))  &
        + (parent%q(mb(1,1)+i-1,mb(2,1)-1,mb(3,1)+k-1,2)          &
        -              info%qc(i,1,k,2))                          &
        *  0.5_r8*(Mob(parent%q(mb(1,1)+i-1,mb(2,1)-1,mb(3,1)+k-1,1))          &
        +          Mob(info%qc(i,1,k,1))))/ch2                     &
        + (0.5_r8*(    info%qold(2*i  ,1,2*k  ,2)  -     info%qold(2*i  ,0,2*k  ,2))   &
        *  0.5_r8*(Mob(info%qold(2*i  ,1,2*k  ,1)) + Mob(info%qold(2*i  ,0,2*k  ,1)))  &
        +  0.5_r8*(    info%qold(2*i-1,1,2*k  ,2)  -     info%qold(2*i-1,0,2*k  ,2))   &
        *  0.5_r8*(Mob(info%qold(2*i-1,1,2*k  ,1)) + Mob(info%qold(2*i-1,0,2*k  ,1)))  &
        +  0.5_r8*(    info%qold(2*i  ,1,2*k-1,2)  -     info%qold(2*i  ,0,2*k-1,2))   &
        *  0.5_r8*(Mob(info%qold(2*i  ,1,2*k-1,1)) + Mob(info%qold(2*i  ,0,2*k-1,1)))  &
        +  0.5_r8*(    info%qold(2*i-1,1,2*k-1,2)  -     info%qold(2*i-1,0,2*k-1,2))   &
        *  0.5_r8*(Mob(info%qold(2*i-1,1,2*k-1,1)) + Mob(info%qold(2*i-1,0,2*k-1,1)))  &
        + (parent%qold(mb(1,1)+i-1,mb(2,1)-1,mb(3,1)+k-1,2)          &
        -              info%qcold(i,1,k,2))                          &
        *  0.5_r8*(Mob(parent%qold(mb(1,1)+i-1,mb(2,1)-1,mb(3,1)+k-1,1))          &
        +          Mob(info%qcold(i,1,k,1))))/ch2
!
      zombiediff3(i,1,k,3) &
        = (0.5_r8*(     info%q(2*i  ,1,2*k  ,3)  -      info%q(2*i  ,0,2*k  ,3))   &
        *  0.5_r8*(MobD(info%q(2*i  ,1,2*k  ,1)) + MobD(info%q(2*i  ,0,2*k  ,1)))  &
        +  0.5_r8*(     info%q(2*i-1,1,2*k  ,3)  -      info%q(2*i-1,0,2*k  ,3))   &
        *  0.5_r8*(MobD(info%q(2*i-1,1,2*k  ,1)) + MobD(info%q(2*i-1,0,2*k  ,1)))  &
        +  0.5_r8*(     info%q(2*i  ,1,2*k-1,3)  -      info%q(2*i  ,0,2*k-1,3))   &
        *  0.5_r8*(MobD(info%q(2*i  ,1,2*k-1,1)) + MobD(info%q(2*i  ,0,2*k-1,1)))  &
        +  0.5_r8*(     info%q(2*i-1,1,2*k-1,3)  -      info%q(2*i-1,0,2*k-1,3))   &
        *  0.5_r8*(MobD(info%q(2*i-1,1,2*k-1,1)) + MobD(info%q(2*i-1,0,2*k-1,1)))  &
        + (parent%q(mb(1,1)+i-1,mb(2,1)-1,mb(3,1)+k-1,3)          &
        -              info%qc(i,1,k,3))                          &
        *  0.5_r8*(MobD(parent%q(mb(1,1)+i-1,mb(2,1)-1,mb(3,1)+k-1,1))          &
        +          MobD(info%qc(i,1,k,1))))/ch2
!
      zombiediff3(i,1,k,1) = zombiediff3(i,1,k,1)*diffcons(1)
      zombiediff3(i,1,k,2) = zombiediff3(i,1,k,2)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  END DO
  IF(masscorrneeded) THEN
    IF(mb(2,1)-1 >= 1) THEN
        parent%f(mb(1,1):mb(1,2),mb(2,1)-1,mb(3,1):mb(3,2),1:nccv) &
      = parent%f(mb(1,1):mb(1,2),mb(2,1)-1,mb(3,1):mb(3,2),1:nccv) &
      + zombiediff3(1:cmx(1),1,1:cmx(3),1:nccv)
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)+mb(1,1:2)-1
      mg(2,1:2) = parent%mglobal(2,1)-1
      mg(3,1:2) = parent%mglobal(3,1)+mb(3,1:2)-1
      CALL AddStripToList(mg,zombiediff3(1:cmx(1),1:1,1:cmx(3),1:nccv))
    END IF
  END IF
!
!
! Back Correction:
  masscorrneeded = .FALSE.
  DO i = 1, cmx(1)
  DO k = 1, cmx(3)
    IF(info%levellandscape(i,cmx(2)+1,k) == level-1) THEN
      zombiediff3(i,2,k,2) &
        = (0.5_r8*(info%q(2*i  ,mx(2),2*k  ,1)-info%q(2*i  ,mx(2)+1,2*k  ,1)) &
        +  0.5_r8*(info%q(2*i-1,mx(2),2*k  ,1)-info%q(2*i-1,mx(2)+1,2*k  ,1)) &
        +  0.5_r8*(info%q(2*i  ,mx(2),2*k-1,1)-info%q(2*i  ,mx(2)+1,2*k-1,1)) &
        +  0.5_r8*(info%q(2*i-1,mx(2),2*k-1,1)-info%q(2*i-1,mx(2)+1,2*k-1,1)) &
        +  parent%q(mb(1,1)+i-1,mb(2,2)+1,mb(3,1)+k-1,1) &
        -  info%qc(i,cmx(2),k,1))/ch2
!
      zombiediff3(i,2,k,1) &
        = (0.5_r8*(    info%q(2*i  ,mx(2),2*k  ,2)  -     info%q(2*i  ,mx(2)+1,2*k  ,2))   &
        *  0.5_r8*(Mob(info%q(2*i  ,mx(2),2*k  ,1)) + Mob(info%q(2*i  ,mx(2)+1,2*k  ,1)))  &
        +  0.5_r8*(    info%q(2*i-1,mx(2),2*k  ,2)  -     info%q(2*i-1,mx(2)+1,2*k  ,2))   &
        *  0.5_r8*(Mob(info%q(2*i-1,mx(2),2*k  ,1)) + Mob(info%q(2*i-1,mx(2)+1,2*k  ,1)))  &
        +  0.5_r8*(    info%q(2*i  ,mx(2),2*k-1,2)  -     info%q(2*i  ,mx(2)+1,2*k-1,2))   &
        *  0.5_r8*(Mob(info%q(2*i  ,mx(2),2*k-1,1)) + Mob(info%q(2*i  ,mx(2)+1,2*k-1,1)))  &
        +  0.5_r8*(    info%q(2*i-1,mx(2),2*k-1,2)  -     info%q(2*i-1,mx(2)+1,2*k-1,2))   &
        *  0.5_r8*(Mob(info%q(2*i-1,mx(2),2*k-1,1)) + Mob(info%q(2*i-1,mx(2)+1,2*k-1,1)))  &
        +  (parent%q(mb(1,1)+i-1,mb(2,2)+1,mb(3,1)+k-1,2)                    &
        -              info%qc(i,cmx(2),k,2))                                &
        *  0.5_r8*(Mob(parent%q(mb(1,1)+i-1,mb(2,2)+1,mb(3,1)+k-1,1))                    &
        +          Mob(info%qc(i,cmx(2),k,1))))/ch2                          &
        + (0.5_r8*(    info%qold(2*i  ,mx(2),2*k  ,2)  -     info%qold(2*i  ,mx(2)+1,2*k  ,2))   &
        *  0.5_r8*(Mob(info%qold(2*i  ,mx(2),2*k  ,1)) + Mob(info%qold(2*i  ,mx(2)+1,2*k  ,1)))  &
        +  0.5_r8*(    info%qold(2*i-1,mx(2),2*k  ,2)  -     info%qold(2*i-1,mx(2)+1,2*k  ,2))   &
        *  0.5_r8*(Mob(info%qold(2*i-1,mx(2),2*k  ,1)) + Mob(info%qold(2*i-1,mx(2)+1,2*k  ,1)))  &
        +  0.5_r8*(    info%qold(2*i  ,mx(2),2*k-1,2)  -     info%qold(2*i  ,mx(2)+1,2*k-1,2))   &
        *  0.5_r8*(Mob(info%qold(2*i  ,mx(2),2*k-1,1)) + Mob(info%qold(2*i  ,mx(2)+1,2*k-1,1)))  &
        +  0.5_r8*(    info%qold(2*i-1,mx(2),2*k-1,2)  -     info%qold(2*i-1,mx(2)+1,2*k-1,2))   &
        *  0.5_r8*(Mob(info%qold(2*i-1,mx(2),2*k-1,1)) + Mob(info%qold(2*i-1,mx(2)+1,2*k-1,1)))  &
        +  (parent%qold(mb(1,1)+i-1,mb(2,2)+1,mb(3,1)+k-1,2)                    &
        -              info%qcold(i,cmx(2),k,2))                                &
        *  0.5_r8*(Mob(parent%qold(mb(1,1)+i-1,mb(2,2)+1,mb(3,1)+k-1,1))                    &
        +          Mob(info%qcold(i,cmx(2),k,1))))/ch2
!
      zombiediff3(i,2,k,3) &
        = (0.5_r8*(     info%q(2*i  ,mx(2),2*k  ,3)  -      info%q(2*i  ,mx(2)+1,2*k  ,3))   &
        *  0.5_r8*(MobD(info%q(2*i  ,mx(2),2*k  ,1)) + MobD(info%q(2*i  ,mx(2)+1,2*k  ,1)))  &
        +  0.5_r8*(     info%q(2*i-1,mx(2),2*k  ,3)  -      info%q(2*i-1,mx(2)+1,2*k  ,3))   &
        *  0.5_r8*(MobD(info%q(2*i-1,mx(2),2*k  ,1)) + MobD(info%q(2*i-1,mx(2)+1,2*k  ,1)))  &
        +  0.5_r8*(     info%q(2*i  ,mx(2),2*k-1,3)  -      info%q(2*i  ,mx(2)+1,2*k-1,3))   &
        *  0.5_r8*(MobD(info%q(2*i  ,mx(2),2*k-1,1)) + MobD(info%q(2*i  ,mx(2)+1,2*k-1,1)))  &
        +  0.5_r8*(     info%q(2*i-1,mx(2),2*k-1,3)  -      info%q(2*i-1,mx(2)+1,2*k-1,3))   &
        *  0.5_r8*(MobD(info%q(2*i-1,mx(2),2*k-1,1)) + MobD(info%q(2*i-1,mx(2)+1,2*k-1,1)))  &
        +  (parent%q(mb(1,1)+i-1,mb(2,2)+1,mb(3,1)+k-1,3)                    &
        -              info%qc(i,cmx(2),k,3))                                &
        *  0.5_r8*(MobD(parent%q(mb(1,1)+i-1,mb(2,2)+1,mb(3,1)+k-1,1))                    &
        +          MobD(info%qc(i,cmx(2),k,1))))/ch2
!
      zombiediff3(i,2,k,1) = zombiediff3(i,2,k,1)*diffcons(1)
      zombiediff3(i,2,k,2) = zombiediff3(i,2,k,2)*diffcons(2)
!
      masscorrneeded = .TRUE.
    END IF
  END DO
  END DO
  IF(masscorrneeded) THEN
    IF(mb(2,2)+1 <= pmx(2)) THEN
        parent%f(mb(1,1):mb(1,2),mb(2,2)+1,mb(3,1):mb(3,2),1:nccv) &
      = parent%f(mb(1,1):mb(1,2),mb(2,2)+1,mb(3,1):mb(3,2),1:nccv) &
      + zombiediff3(1:cmx(1),2,1:cmx(3),1:nccv)
    ELSE
      mg = 1
      mg(1,1:2) = parent%mglobal(1,1)+mb(1,1:2)-1
      mg(2,1:2) = parent%mglobal(2,1)+pmx(2)
      mg(3,1:2) = parent%mglobal(3,1)+mb(3,1:2)-1
      CALL AddStripToList(mg,zombiediff3(1:cmx(1),2:2,1:cmx(3),1:nccv))
    END IF
  END IF
!!

  DEALLOCATE(zombiediff1,zombiediff2,zombiediff3)

CASE DEFAULT
  PRINT *, 'BSAM 2.0 CoarseLoadingCorrection:'
  PRINT *, '         Only ndims = 2,3 are supported.'
  STOP
END SELECT
!
END FUNCTION CoarseLoadingCorrection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AddStripToList(mg,masscorr)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN):: mg
REAL(KIND=r8), DIMENSION(mg(1,1):mg(1,2), &
                         mg(2,1):mg(2,2), &
                         mg(3,1):mg(3,2),1:nccv), INTENT(IN):: masscorr
!
! This routine adds the mass correction strip to a global linked list. Each
! strip is indexed with global coordinates.
!
nmasscorrlayers =  nmasscorrlayers+1
ALLOCATE(currentlayer)
currentlayer%id = nmasscorrlayers
currentlayer%mg(1:maxdims,1:2) = mg(1:maxdims,1:2)
!
ALLOCATE(currentlayer%masscorr(mg(1,1):mg(1,2), &
                               mg(2,1):mg(2,2), &
                               mg(3,1):mg(3,2),1:nccv))
!
currentlayer%masscorr(mg(1,1):mg(1,2),mg(2,1):mg(2,2),mg(3,1):mg(3,2),1:nccv) &
           = masscorr(mg(1,1):mg(1,2),mg(2,1):mg(2,2),mg(3,1):mg(3,2),1:nccv)
!
currentlayer%prevlayer => lastlayer
lastlayer => currentlayer
!
END SUBROUTINE AddStripToList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE TransferMassCorrectionLayers(level)
USE NodeInfoDef
USE TreeOps, ONLY: ApplyOnLevel
USE Boundary, ONLY: GetPeriodicOffsets
IMPLICIT NONE
!
INTEGER, INTENT(IN):: level
!
TYPE(funcparam):: dummy
LOGICAL:: periodicbuffer
INTEGER:: n, offset, polarity
INTEGER, DIMENSION(1:maxdims,1:2):: mgsave
!
mgsave = 1
currentlayer => lastlayer
searchloop: DO
  IF(.NOT. ASSOCIATED(currentlayer%prevlayer)) EXIT searchloop
!
  dummy%offset = 0
  CALL ApplyOnLevel(level,TransferLoads,dummy)
!
! Buffering of periodic edge tags. Check to see if the buffer area cuts across 
! a periodic boundary.  If so, add offset and apply buffer:
  mgsave(1:ndims,1:2) = currentlayer%mg(1:ndims,1:2)
  IF(periodicboundaryconditions) THEN
    CALL GetPeriodicOffsets(level)
    DO polarity = -1, 1, 2
      DO offset = 1, nperiodicoffsets
        dummy%offset = 0
        DO n = 1, ndims
          dummy%offset(n) = polarity*poffset(n,offset)
          currentlayer%mg(n,1:2) = mgsave(n,1:2)+dummy%offset(n)
        END DO
        CALL ApplyOnLevel(level,TransferLoads,dummy)
      END DO
    END DO
  END IF
  currentlayer%mg(1:ndims,1:2) = mgsave(1:ndims,1:2)
!
  currentlayer => currentlayer%prevlayer
END DO searchloop
!
END SUBROUTINE TransferMassCorrectionLayers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION TransferLoads(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
INTEGER:: n
INTEGER, DIMENSION(1:maxdims,1:2):: mglobal, mlf, mls, movr
!
TransferLoads = err_ok
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mglobal = 1; mglobal(1:ndims,1:2) = info%mglobal(1:ndims,1:2)
!
! 1. Find overlap region in global index space:
movr = 1
movr(1:ndims,1) = MAX(currentlayer%mg(1:ndims,1),mglobal(1:ndims,1))
movr(1:ndims,2) = MIN(currentlayer%mg(1:ndims,2),mglobal(1:ndims,2))
!
! 2. Check for nonempty intersection:
IF(ANY(movr(1:maxdims,2)-movr(1:maxdims,1)<0)) RETURN
!
! 3. Transform common index space to grid index spaces:
mlf = 1; mls = 1
DO n = 1, ndims
  mlf(n,1:2) = movr(n,1:2)-mglobal(n,1)+1
  mls(n,1:2) = movr(n,1:2)-dummy%offset(n)
END DO
!
                 info%f(mlf(1,1):mlf(1,2),  &
                        mlf(2,1):mlf(2,2),  &
                        mlf(3,1):mlf(3,2),1:nccv)  &
=                info%f(mlf(1,1):mlf(1,2),  &
                        mlf(2,1):mlf(2,2),  &
                        mlf(3,1):mlf(3,2),1:nccv)  & 
+ currentlayer%masscorr(mls(1,1):mls(1,2),  &
                        mls(2,1):mls(2,2),  &
                        mls(3,1):mls(3,2),1:nccv)
!
END FUNCTION TransferLoads
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE DeleteMassCorrectionList
USE NodeInfoDef
IMPLICIT NONE
!
searchloop: DO
  currentlayer => lastlayer%prevlayer
  DEALLOCATE(lastlayer%masscorr)
  DEALLOCATE(lastlayer)
  lastlayer => currentlayer
  IF(.NOT. ASSOCIATED(lastlayer%prevlayer)) EXIT searchloop
END DO searchloop
!
END SUBROUTINE DeleteMassCorrectionList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION PullFCLoadings(parent,dummy1)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, ApplyOnChildren, ApplyOnMyLevelNbrs
USE Boundary, ONLY: TransferBCOneWayFC
IMPLICIT NONE
!
! Dimensionally invariant restriction.
!
TYPE(nodeinfo):: parent
TYPE(funcparam):: dummy1
!
TYPE(funcparam):: dummy2
!
PullFCLoadings = err_ok
!
! Pull the child-level data down to the parent:
CALL ApplyOnChildren(CoarseLoadingFunctionFC,dummy2)
!
! Synchronize the parent's face centered variable data with other parent-
! level grids. iswitch = 2 indicates loading data f#:
dummy2%iswitch = 2
IF(parent%level > rootlevel) &
      CALL ApplyOnMyLevelNbrs(parent%level,TransferBCOneWayFC,dummy2)  
!
END FUNCTION PullFCLoadings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION CoarseLoadingFunctionFC(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY: Restriction2DV1, Restriction2DV2, &
                         Restriction3DV1, Restriction3DV2, Restriction3DV3
USE Problem, ONLY: Operator2DV1, Operator2DV2, Operator3DV1, Operator3DV2, Operator3DV3
IMPLICIT NONE
!
! Dimensionally invariant coarse loading function routine for face variables.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror, nnaxv
INTEGER, DIMENSION(1:maxdims):: cmx, mx
INTEGER, DIMENSION(1:maxdims,1:2):: amb, mb
REAL(KIND=r8):: ch
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
CoarseLoadingFunctionFC = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
ierror = GetParentInfo(parent)
ch = parent%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
nnaxv = MAX(naxv,1)
amb = 1; IF(naxv>0) amb(1:ndims,1:2) = mb(1:ndims,1:2)
!
CALL Residual(info)
!
SELECT CASE(ndims)
CASE(2)
  parent%f1(mb(1,1)-1:mb(1,2),mb(2,1):mb(2,2),1,1:nfcv) &
    = Restriction2DV1(info%rf1(0: mx(1)  ,1: mx(2)  ,1,1:nfcv)) &
    +    Operator2DV1(info%qc (0:cmx(1)+1,0:cmx(2)+1,1,1:nccv), &
                   parent%qold (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nccv), &
                   parent%v1   (  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nfcv), &
                   parent%v1old(  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nfcv), &
                   parent%v2   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1,1,1: nfcv), &
                   parent%v2old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1,1,1: nfcv), &
                   parent%aux  ( amb(1,1)-1:amb(1,2)+1, &
                                 amb(2,1)-1:amb(2,2)+1,1,1:nnaxv), &
                   cmx(1:2),ch,xlower(1:2))
  parent%f2(mb(1,1):mb(1,2),mb(2,1)-1:mb(2,2),1,1:nfcv) &
    = Restriction2DV2(info%rf2(1: mx(1)  ,0: mx(2)  ,1,1:nfcv)) &
    +    Operator2DV2(info%qc (0:cmx(1)+1,0:cmx(2)+1,1,1:nccv), &
                   parent%qold (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nccv), &
                   parent%v1   (  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nfcv), &
                   parent%v1old(  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nfcv), &
                   parent%v2   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1,1,1: nfcv), &
                   parent%v2old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1,1,1: nfcv), &
                   parent%aux  ( amb(1,1)-1:amb(1,2)+1, &
                                 amb(2,1)-1:amb(2,2)+1,1,1:nnaxv), &
                   cmx(1:2),ch,xlower(1:2))
  CASE(3)
!
      parent%f1(mb(1,1)-1:mb(1,2),mb(2,1):mb(2,2),mb(3,1):mb(3,2),1:nfcv) &
      = Restriction3DV1(info%rf1(0: mx(1)  ,1: mx(2)  ,1: mx(3)  ,1:nfcv)) &
      +    Operator3DV1(info%qc (0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nccv), &
                   parent%qold (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nccv), &
                   parent%v1   (  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v1old(  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v2   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v2old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v3   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-2: mb(3,2)+1,1: nfcv), &
                   parent%v3old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-2: mb(3,2)+1,1: nfcv), &
                   parent%aux  ( amb(1,1)-1:amb(1,2)+1, &
                                 amb(2,1)-1:amb(2,2)+1, &
                                 amb(3,1)-1:amb(3,2)+1,1:nnaxv), &
                   cmx(1:3),ch,xlower(1:3))
!
      parent%f2(mb(1,1):mb(1,2),mb(2,1)-1:mb(2,2),mb(3,1):mb(3,2),1:nfcv) &
      = Restriction3DV2(info%rf2(1: mx(1)  ,0: mx(2)  ,1: mx(3)  ,1:nfcv)) &
      +    Operator3DV2(info%qc (0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nccv), &
                   parent%qold (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nccv), &
                   parent%v1   (  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v1old(  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v2   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v2old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v3   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-2: mb(3,2)+1,1: nfcv), &
                   parent%v3old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-2: mb(3,2)+1,1: nfcv), &
                   parent%aux  ( amb(1,1)-1:amb(1,2)+1, &
                                 amb(2,1)-1:amb(2,2)+1, &
                                 amb(3,1)-1:amb(3,2)+1,1:nnaxv), &
                   cmx(1:3),ch,xlower(1:3))
!
      parent%f3(mb(1,1):mb(1,2),mb(2,1):mb(2,2),mb(3,1)-1:mb(3,2),1:nfcv) &
      = Restriction3DV3(info%rf3(1: mx(1)  ,1: mx(2)  ,0: mx(3)  ,1:nfcv)) &
      +    Operator3DV3(info%qc (0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nccv), &
                   parent%qold (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nccv), &
                   parent%v1   (  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v1old(  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v2   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v2old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v3   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-2: mb(3,2)+1,1: nfcv), &
                   parent%v3old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-2: mb(3,2)+1,1: nfcv), &
                   parent%aux  ( amb(1,1)-1:amb(1,2)+1, &
                                 amb(2,1)-1:amb(2,2)+1, &
                                 amb(3,1)-1:amb(3,2)+1,1:nnaxv), &
                   cmx(1:3),ch,xlower(1:3))
!
  CASE DEFAULT
    PRINT *, 'BSAM 2.0 CoarseLoadingFunctionFC:'
    PRINT *, '         Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION CoarseLoadingFunctionFC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION RelativeTruncationError(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY: Restriction2D, Restriction3D
USE Problem, ONLY: Operator2D, Operator3D
IMPLICIT NONE
!
! Dimensionally invariant relative trunctation error:
!
!   L_{2h}(I_h^{2h} q_h)-I_h^{2h}(L_h q_h).
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror, nnaxv
INTEGER, DIMENSION(1:maxdims):: amx, cmx, mx
INTEGER, DIMENSION(1:maxdims,1:2):: amb, mb
REAL(KIND=r8):: ch, h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
RelativeTruncationError = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
nnaxv = MAX(naxv,1)
amx = 1; IF(naxv>0) amx(1:ndims) = mx(1:ndims)
!
ierror = GetParentInfo(parent)
ch = parent%dx(1)
!
amb = 1; IF(naxv>0) amb(1:ndims,1:2) = mb(1:ndims,1:2)
!
SELECT CASE(ndims)
  CASE(2)
    info%qrte(1:cmx(1),1:cmx(2),1,1:nccv) &
      = Restriction2D( &
        Operator2D(info%q    (  0: mx(1)+1, 0: mx(2)+1,1,1: nccv), &
                   info%qold (  0: mx(1)+1, 0: mx(2)+1,1,1: nccv), &
                   info%v1   ( -1: mx(1)+1, 0: mx(2)+1,1,1: nccv), &
                   info%v1old( -1: mx(1)+1, 0: mx(2)+1,1,1: nccv), &
                   info%v2   (  0: mx(1)+1,-1: mx(2)+1,1,1: nccv), &
                   info%v2old(  0: mx(1)+1,-1: mx(2)+1,1,1: nccv), &
                   info%aux  (  0:amx(1)+1, 0:amx(2)+1,1,1:nnaxv), &
                   mx(1:2),h,xlower(1:2))) &
      - Operator2D(info%qc(0:cmx(1)+1,0:cmx(2)+1,1,1:nccv), &
                   parent%qold (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nccv), &
                   parent%v1   (  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nccv), &
                   parent%v1old(  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1,1,1: nccv), &
                   parent%v2   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1,1,1: nccv), &
                   parent%v2old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1,1,1: nccv), &
                   parent%aux  ( amb(1,1)-1:amb(1,2)+1, &
                                 amb(2,1)-1:amb(2,2)+1,1,1:nnaxv), &
                   cmx(1:2),ch,xlower(1:2))
  CASE(3)
    info%qrte(1:cmx(1),1:cmx(2),1:cmx(3),1:nccv) &
      = Restriction3D( &
        Operator3D(info%q    ( 0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                   info%qold ( 0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                   info%v1   (-1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v1old(-1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v2   ( 0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v2old( 0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                   info%v3   ( 0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                   info%v3old( 0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                   info%aux  ( 0:amx(1)+1, 0:amx(2)+1, 0:amx(3)+1,1:nnaxv), &
                   mx(1:3),h,xlower(1:3))) &
      - Operator3D(info%qc(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nccv), &
                   parent%qold(   mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nccv), &
                   parent%v1(     mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v1old(  mb(1,1)-2: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v2   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v2old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-2: mb(2,2)+1, &
                                  mb(3,1)-1: mb(3,2)+1,1: nfcv), &
                   parent%v3   (  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-2: mb(3,2)+1,1: nfcv), &
                   parent%v3old(  mb(1,1)-1: mb(1,2)+1, &
                                  mb(2,1)-1: mb(2,2)+1, &
                                  mb(3,1)-2: mb(3,2)+1,1: nfcv), &
                   parent%aux(   amb(1,1)-1:amb(1,2)+1, &
                                 amb(2,1)-1:amb(2,2)+1, &
                                 amb(3,1)-1:amb(3,2)+1,1:nnaxv), &
                   cmx(1:3),ch,xlower(1:3))
  CASE DEFAULT
    PRINT *, 'RelativeTruncationError: Only ndims = 2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION RelativeTruncationError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Residual(info)
USE NodeInfoDef
USE Problem, ONLY: Operator2D, Operator2DV1, Operator2DV2, &
                   Operator3D, Operator3DV1, Operator3DV2, Operator3DV3
IMPLICIT NONE
!
TYPE(nodeinfo):: info
!
INTEGER:: nnaxv
INTEGER, DIMENSION(1:maxdims):: amx, mx
REAL(KIND=r8):: h
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
h = info%dx(1)
xlower(1:ndims) = info%xlower(1:ndims)
!
nnaxv = MAX(naxv,1)
amx = 1; IF(naxv>0) amx(1:ndims) = mx(1:ndims)
!
SELECT CASE(ndims)
  CASE(2)
                     info%rf   ( 1: mx(1)  , 1: mx(2)  ,1,1: nccv)  &
        =            info%f    ( 1: mx(1)  , 1: mx(2)  ,1,1: nccv)  &
        - Operator2D(info%q    ( 0: mx(1)+1, 0: mx(2)+1,1,1: nccv), &
                     info%qold ( 0: mx(1)+1, 0: mx(2)+1,1,1: nccv), &
                     info%v1   (-1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                     info%v1old(-1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                     info%v2   ( 0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                     info%v2old( 0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                     info%aux  ( 0:amx(1)+1, 0:amx(2)+1,1,1:nnaxv), &
                     mx(1:2),h,xlower(1:2))
!
                     info%rf1  ( 0: mx(1)  , 1: mx(2)  ,1,1: nfcv)  &
      =              info%f1   ( 0: mx(1)  , 1: mx(2)  ,1,1: nfcv)  &
      - Operator2DV1(info%q    ( 0: mx(1)+1, 0: mx(2)+1,1,1: nccv), &
                     info%qold ( 0: mx(1)+1, 0: mx(2)+1,1,1: nccv), &
                     info%v1   (-1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                     info%v1old(-1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                     info%v2   ( 0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                     info%v2old( 0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                     info%aux  ( 0:amx(1)+1, 0:amx(2)+1,1,1:nnaxv), &
                     mx(1:2),h,xlower(1:2))
!
                     info%rf2  ( 1: mx(1)  , 0: mx(2)  ,1,1: nfcv)  &
      =              info%f2   ( 1: mx(1)  , 0: mx(2)  ,1,1: nfcv)  &
      - Operator2DV2(info%q    ( 0: mx(1)+1, 0: mx(2)+1,1,1: nccv), &
                     info%qold ( 0: mx(1)+1, 0: mx(2)+1,1,1: nccv), &
                     info%v1   (-1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                     info%v1old(-1: mx(1)+1, 0: mx(2)+1,1,1: nfcv), &
                     info%v2   ( 0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                     info%v2old( 0: mx(1)+1,-1: mx(2)+1,1,1: nfcv), &
                     info%aux  ( 0:amx(1)+1, 0:amx(2)+1,1,1:nnaxv), &
                     mx(1:2),h,xlower(1:2))
  CASE(3)
                     info%rf   ( 1: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1: nccv)  &
        =            info%f    ( 1: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1: nccv)  &
        - Operator3D(info%q    ( 0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                     info%qold ( 0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                     info%v1   (-1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v1old(-1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v2   ( 0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v2old( 0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v3   ( 0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                     info%v3old( 0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                     info%aux  ( 0:amx(1)+1, 0:amx(2)+1, 0:amx(3)+1,1:nnaxv), &
                     mx(1:3),h,xlower(1:3))
!
                     info%rf1  ( 0: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1: nfcv)  &
        =            info%f1   ( 0: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1: nfcv)  &
      - Operator3DV1(info%q    ( 0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                     info%qold ( 0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                     info%v1   (-1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v1old(-1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v2   ( 0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v2old( 0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v3   ( 0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                     info%v3old( 0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                     info%aux  ( 0:amx(1)+1, 0:amx(2)+1, 0:amx(3)+1,1:nnaxv), &
                     mx(1:3),h,xlower(1:3))
!
                     info%rf2  ( 1: mx(1)  , 0: mx(2)  , 1: mx(3)  ,1: nfcv)  &
      =              info%f2   ( 1: mx(1)  , 0: mx(2)  , 1: mx(3)  ,1: nfcv)  &
      - Operator3DV2(info%q    ( 0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                     info%qold ( 0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                     info%v1   (-1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v1old(-1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v2   ( 0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v2old( 0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v3   ( 0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                     info%v3old( 0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                     info%aux  ( 0:amx(1)+1, 0:amx(2)+1, 0:amx(3)+1,1:nnaxv), &
                     mx(1:3),h,xlower(1:3))
!
                     info%rf3  ( 1: mx(1)  , 1: mx(2)  , 0: mx(3)  ,1: nfcv)  &
      =              info%f3   ( 1: mx(1)  , 1: mx(2)  , 0: mx(3)  ,1: nfcv)  &
      - Operator3DV3(info%q    ( 0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                     info%qold ( 0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nccv), &
                     info%v1   (-1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v1old(-1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v2   ( 0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v2old( 0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1: nfcv), &
                     info%v3   ( 0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                     info%v3old( 0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1: nfcv), &
                     info%aux  ( 0:amx(1)+1, 0:amx(2)+1, 0:amx(3)+1,1:nnaxv), &
                     mx(1:3),h,xlower(1:3))
!
  CASE DEFAULT
    PRINT *, 'Residua! fields at all above root levels, saving the data to uniformgrid(level)%q'
    STOP
END SELECT
!
END SUBROUTINE Residual
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION CorrectFine(info,dummy)
USE NodeInfoDef
USE TreeOps, ONLY: err_ok, GetParentInfo
USE GridUtilities, ONLY:       Prolongation2D    ,       Prolongation3D,     &
                               Prolongation2DV1Ex,       Prolongation2DV2Ex, &
                          BiLinProlongationP1    ,  BiLinProlongationP2,     &
                         TriLinProlongationP1    , TriLinProlongationP2,     &
                               Prolongation3DV1Ex,       Prolongation3DV2Ex, & 
                               Prolongation3DV3Ex
IMPLICIT NONE
!
! Dimensionally invariant coarse-grid-correction routine.
!
TYPE(nodeinfo):: info
TYPE(funcparam):: dummy
!
TYPE(nodeinfo), POINTER:: parent
INTEGER:: ierror, ll
INTEGER, DIMENSION(1:maxdims):: cmx, mx, ul
INTEGER, DIMENSION(1:maxdims,1:2):: mb
!
CorrectFine = err_ok
!
IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
ierror = GetParentInfo(parent)
!
mx = 1; mx(1:ndims) = info%mx(1:ndims)
cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
SELECT CASE(ndims)
  CASE(2)
       info%qc(       1-mbc:cmx(1)+mbc ,          &
                      1-mbc:cmx(2)+mbc ,1,1:nccv) &
    = parent%q( mb(1,1)-mbc:mb(1,2)+mbc,          &
                mb(2,1)-mbc:mb(2,2)+mbc,1,1:nccv) &
    -  info%qc(       1-mbc:cmx(1)+mbc ,          &
                      1-mbc:cmx(2)+mbc ,1,1:nccv)
!
      info%v1c(          -1:cmx(1)+1   ,          &
                          0:cmx(2)+1   ,1,1:nfcv) &
    = parent%v1(  mb(1,1)-2:mb(1,2)+1  ,          &
                  mb(2,1)-1:mb(2,2)+1  ,1,1:nfcv) &
    -  info%v1c(         -1:cmx(1)+1   ,          &
                          0:cmx(2)+1,1,1:nfcv)
!
       info%v2c(0:cmx(1)+1,-1:cmx(2)+1,1,1:nfcv) &
    = parent%v2(mb(1,1)-1:mb(1,2)+1, &
                mb(2,1)-2:mb(2,2)+1,1,1:nfcv) &
    -  info%v2c(0:cmx(1)+1,-1:cmx(2)+1,1,1:nfcv)
!
	SELECT CASE(mbc)
		CASE(1)
!
! Bilinear prolongation:
          info%q( 0: mx(1)+1,0: mx(2)+1,1,1:nccv) &
        = info%q( 0: mx(1)+1,0: mx(2)+1,1,1:nccv) &
        + BiLinProlongationP1( &
          info%qc(0:cmx(1)+1,0:cmx(2)+1,1,1:nccv))
!
! Simple injection:
!          info%q( 1: mx(1),1: mx(2),1,1:nccv) &
!        = info%q( 1: mx(1),1: mx(2),1,1:nccv) &
!        + Prolongation2D( &
!          info%qc(1:cmx(1),1:cmx(2),1,1:nccv))
!
		CASE(2)
!
! Bilinear prolongation:
          info%q( -1: mx(1)+2,-1: mx(2)+2,1,1:nccv) &
        = info%q( -1: mx(1)+2,-1: mx(2)+2,1,1:nccv) &
        + BiLinProlongationP2( &
          info%qc(-1:cmx(1)+2,-1:cmx(2)+2,1,1:nccv))
	END SELECT
!
          info%v1( -1: mx(1)+1, 0: mx(2)+1,1,1:nfcv) &
        = info%v1( -1: mx(1)+1, 0: mx(2)+1,1,1:nfcv) &
        + Prolongation2DV1Ex( &
          info%v1c(-1:cmx(1)+1, 0:cmx(2)+1,1,1:nfcv))
!
          info%v2(  0: mx(1)+1,-1: mx(2)+1,1,1:nfcv) &
        = info%v2(  0: mx(1)+1,-1: mx(2)+1,1,1:nfcv) &
        + Prolongation2DV2Ex( &
          info%v2c( 0:cmx(1)+1,-1:cmx(2)+1,1,1:nfcv))
!
  CASE(3)
!
       info%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1-mbc:cmx(3)+mbc,1:nccv) &
    = parent%q(mb(1,1)-mbc:mb(1,2)+mbc, &
               mb(2,1)-mbc:mb(2,2)+mbc, &
               mb(3,1)-mbc:mb(3,2)+mbc,1:nccv) &
    - info%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1-mbc:cmx(3)+mbc,1:nccv)
!
       info%v1c(-1:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nfcv) &
    = parent%v1(mb(1,1)-2:mb(1,2)+1, &
                mb(2,1)-1:mb(2,2)+1, &
                mb(3,1)-1:mb(3,2)+1,1:nfcv) &
    -  info%v1c(-1:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nfcv)
!
       info%v2c(0:cmx(1)+1,-1:cmx(2)+1,0:cmx(3)+1,1:nfcv) &
    = parent%v2(mb(1,1)-1:mb(1,2)+1, &
                mb(2,1)-2:mb(2,2)+1, &
                mb(3,1)-1:mb(3,2)+1,1:nfcv) &
    -  info%v2c(0:cmx(1)+1,-1:cmx(2)+1,0:cmx(3)+1,1:nfcv)
!
       info%v3c(0:cmx(1)+1,0:cmx(2)+1,-1:cmx(3)+1,1:nfcv) &
    = parent%v3(mb(1,1)-1:mb(1,2)+1, &
                mb(2,1)-1:mb(2,2)+1, &
                mb(3,1)-2:mb(3,2)+1,1:nfcv) &
    -  info%v3c(0:cmx(1)+1,0:cmx(2)+1,-1:cmx(3)+1,1:nfcv)
!
	SELECT CASE(mbc)
		CASE(1)
!
! Trilinear prolongation:
        info%q( 0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nccv) &
      = info%q( 0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nccv) &
      + TriLinProlongationP1( &
        info%qc(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1,1:nccv))
!
! Simple injection:
!        info%q( 1: mx(1),1: mx(2),1: mx(3),1:nccv) &
!      = info%q( 1: mx(1),1: mx(2),1: mx(3),1:nccv) &
!      + Prolongation3D( &
!        info%qc(1:cmx(1),1:cmx(2),1:cmx(3),1:nccv))
		CASE(2)
!
! Trilinear prolongation:
        info%q( -1: mx(1)+2,-1: mx(2)+2,-1: mx(3)+2,1:nccv) &
      = info%q( -1: mx(1)+2,-1: mx(2)+2,-1: mx(3)+2,1:nccv) &
      + TriLinProlongationP2( &
        info%qc(-1:cmx(1)+2,-1:cmx(2)+2,-1:cmx(3)+2,1:nccv))
	END SELECT
!
       info%v1( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1:nfcv) &
     = info%v1( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1:nfcv) &
     + Prolongation3DV1Ex( &
      info%v1c( -1:cmx(1)+1, 0:cmx(2)+1, 0:cmx(3)+1,1:nfcv))
!
       info%v2(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1:nfcv) &
     = info%v2(  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1:nfcv) &
     + Prolongation3DV2Ex( &
       info%v2c( 0:cmx(1)+1,-1:cmx(2)+1, 0:cmx(3)+1,1:nfcv))
!
       info%v3(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1:nfcv) &
     = info%v3(  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1:nfcv) &
     + Prolongation3DV3Ex( &
       info%v3c( 0:cmx(1)+1, 0:cmx(2)+1,-1:cmx(3)+1,1:nfcv))
!
  CASE DEFAULT
    PRINT *, 'CorrectFine: Only ndims=2,3 are supported.'
    STOP
END SELECT
!
END FUNCTION CorrectFine
END MODULE AFASRoutines

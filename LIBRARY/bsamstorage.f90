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
! Portions of the code
!
! (c) Copyright Sorin Mitran, 2002
! Department of Mathematics
! University of North Carolina at Chapel Hill
! mitran@amath.unc.edu
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
! File:             bsamstorage.f90
! Purpose:          BSAM memory allocation module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 2.0 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------
MODULE BSAMStorage
!
CONTAINS
!
SUBROUTINE AllocFields(info,parent)
USE NodeInfoDef
IMPLICIT NONE
!
TYPE(nodeinfo):: info
TYPE(nodeinfo), INTENT(IN), OPTIONAL:: parent
!
INTEGER:: ierror, nnaxv, n
INTEGER, DIMENSION(1:maxdims):: amx, cmx, mx
!
mx = info%mx
cmx = 1
DO n = 1, ndims
  IF(MODULO(mx(n),2)==1) THEN
    cmx(n) = 1
  ELSE
    cmx(n) = mx(n)/2
  END IF
END DO
!
nnaxv = MAX(naxv,1)
amx = 1; IF(naxv>0) amx(1:ndims) = mx(1:ndims)
!
! Allocate space for all POINTER components from nodeinfodef:
SELECT CASE(ndims)
  CASE(2)
    ALLOCATE( &
      info%errorflags(1    : mx(1)    ,1    : mx(2)    ,1:1        ), &
      info%q(         1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1:1,1: nccv), &
      info%qold(      1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1:1,1: nccv), &
      info%qc(        1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1:1,1: nccv), &
      info%qcold(     1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1:1,1: nccv), &
      info%qrte(      1    :cmx(1)    ,1    :cmx(2)    ,1:1,1: nccv), &
      info%aux(       1-mbc:amx(1)+mbc,1-mbc:amx(2)+mbc,1:1,1:nnaxv), &
      info%f(         1    : mx(1)    ,1    : mx(2)    ,1:1,1: nccv), &
      info%rf(        1    : mx(1)    ,1    : mx(2)    ,1:1,1: nccv), &
      info%ftmp(      1    : mx(1)    ,1    : mx(2)    ,1:1,1: nccv), &
        STAT=ierror)
!
    IF(ierror/=0) THEN
      PRINT *, 'BSAM 2.0: Error in allocation of cc variables in AllocFields'
      STOP
    END IF
!
! These variables will not be present if there is no velocity problem:
    ALLOCATE(info%v1   (-1: mx(1)+1, 0: mx(2)+1,1:1,1:nfcv), &
             info%v1old(-1: mx(1)+1, 0: mx(2)+1,1:1,1:nfcv), &
             info%v1c  (-1:cmx(1)+1, 0:cmx(2)+1,1:1,1:nfcv), &
             info%f1   ( 0: mx(1)  , 1: mx(2)  ,1:1,1:nfcv), &
             info%f1tmp( 0: mx(1)  , 1: mx(2)  ,1:1,1:nfcv), &
             info%rf1  ( 0: mx(1)  , 1: mx(2)  ,1:1,1:nfcv), &
             info%v2   ( 0: mx(1)+1,-1: mx(2)+1,1:1,1:nfcv), &
             info%v2old( 0: mx(1)+1,-1: mx(2)+1,1:1,1:nfcv), &
             info%v2c  ( 0:cmx(1)+1,-1:cmx(2)+1,1:1,1:nfcv), &
             info%f2   ( 1: mx(1)  , 0: mx(2)  ,1:1,1:nfcv), &
             info%f2tmp( 1: mx(1)  , 0: mx(2)  ,1:1,1:nfcv), &
             info%rf2  ( 1: mx(1)  , 0: mx(2)  ,1:1,1:nfcv),STAT=ierror)
!
    IF(ierror/=0) THEN
      PRINT *, 'BSAM 2.0: Error in allocation of face variables in AllocFields'
      STOP
    END IF
!
	ALLOCATE(info%levellandscape(0:cmx(1)+1,0:cmx(2)+1,1:1),STAT=ierror)
!
    IF(ierror/=0) THEN
      PRINT *, 'BSAM 2.0: Error in allocation of levellandscape in AllocFields'
      STOP
    END IF
!
! Initialize solution arrays:
info%q  = 0.0_r8; info%qold  = 0.0_r8
!
info%v1 = 0.0_r8; info%v1old = 0.0_r8
!
info%v2 = 0.0_r8; info%v2old = 0.0_r8
!
  CASE(3)
    ALLOCATE( &
      info%errorflags( 1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)            ), &
      info%q         ( 1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1-mbc: mx(3)+mbc,1: nccv), &
      info%qold      ( 1-mbc: mx(1)+mbc,1-mbc: mx(2)+mbc,1-mbc: mx(3)+mbc,1: nccv), &
      info%qc        ( 1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1-mbc:cmx(3)+mbc,1: nccv), &
      info%qcold     ( 1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1-mbc:cmx(3)+mbc,1: nccv), &
      info%qrte      ( 1    :cmx(1)    ,1    :cmx(2)    ,1    :cmx(3)    ,1: nccv), &
      info%aux       ( 1-mbc:amx(1)+mbc,1-mbc:amx(2)+mbc,1-mbc:amx(3)+mbc,1:nnaxv), &
      info%f         ( 1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)    ,1: nccv), &
      info%rf        ( 1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)    ,1: nccv), &
      info%ftmp      ( 1    : mx(1)    ,1    : mx(2)    ,1    : mx(3)    ,1: nccv),STAT=ierror)
!
    IF(ierror/=0) THEN
      PRINT *, 'BSAM 2.0: Error in allocation of cc variables in AllocFields'
      STOP
    END IF
!
! These variables will not be present if there is no velocity problem:
    ALLOCATE(info%v1   (  -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1:nfcv), &
             info%v1old(  -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1:nfcv), &
             info%v1c  (  -1:cmx(1)+1, 0:cmx(2)+1, 0:cmx(3)+1,1:nfcv), &
             info%f1   (   0: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1:nfcv), &
             info%f1tmp(   0: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1:nfcv), &
             info%rf1  (   0: mx(1)  , 1: mx(2)  , 1: mx(3)  ,1:nfcv), &
!           
             info%v2   (   0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1:nfcv), &
             info%v2old(   0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1:nfcv), &
             info%v2c  (   0:cmx(1)+1,-1:cmx(2)+1, 0:cmx(3)+1,1:nfcv), &
             info%f2   (   1: mx(1)  , 0: mx(2)  , 1: mx(3)  ,1:nfcv), &
             info%f2tmp(   1: mx(1)  , 0: mx(2)  , 1: mx(3)  ,1:nfcv), &
             info%rf2  (   1: mx(1)  , 0: mx(2)  , 1: mx(3)  ,1:nfcv), &
!
             info%v3   (   0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1:nfcv), &
             info%v3old(   0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1:nfcv), &
             info%v3c  (   0:cmx(1)+1, 0:cmx(2)+1,-1:cmx(3)+1,1:nfcv), &
             info%f3   (   1: mx(1)  , 1: mx(2)  , 0: mx(3)  ,1:nfcv), &
             info%f3tmp(   1: mx(1)  , 1: mx(2)  , 0: mx(3)  ,1:nfcv), &
             info%rf3  (   1: mx(1)  , 1: mx(2)  , 0: mx(3)  ,1:nfcv),STAT=ierror)
!
    IF(ierror/=0) THEN
      PRINT *, 'BSAM 2.0: Error in allocation of face variables in AllocFields'
      STOP
    END IF
!
	ALLOCATE(info%levellandscape(0:cmx(1)+1,0:cmx(2)+1,0:cmx(3)+1),STAT=ierror)
!
    IF(ierror/=0) THEN
      PRINT *, 'BSAM 2.0: Error in allocation of levellandscape in AllocFields'
      STOP
    END IF
!
! Initialize solution arrays:
info%q  = 0.0_r8; info%qold  = 0.0_r8
!
info%v1 = 0.0_r8; info%v1old = 0.0_r8
!
info%v2 = 0.0_r8; info%v2old = 0.0_r8
!
info%v3 = 0.0_r8; info%v3old = 0.0_r8
!
  CASE DEFAULT
    PRINT *, 'AllocFields: Only ndims = 2,3 are supported.'
    STOP
END SELECT
! Initialize error flags:
info%errorflags = 0
!
info%fieldsallocated = .TRUE.
!
END SUBROUTINE AllocFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE DeAllocFields(info)
USE NodeInfoDef
IMPLICIT NONE
!
TYPE(nodeinfo):: info
!
IF(.NOT. info%fieldsallocated) RETURN
!
IF(ASSOCIATED(info%errorflags))     DEALLOCATE(info%errorflags)
!
IF(ASSOCIATED(info%q))              DEALLOCATE(info%q)
IF(ASSOCIATED(info%qold))           DEALLOCATE(info%qold)
IF(ASSOCIATED(info%qc))             DEALLOCATE(info%qc)
IF(ASSOCIATED(info%qcold))          DEALLOCATE(info%qcold)
IF(ASSOCIATED(info%qrte))           DEALLOCATE(info%qrte)
!
IF(ASSOCIATED(info%aux))            DEALLOCATE(info%aux)
!
IF(ASSOCIATED(info%f))              DEALLOCATE(info%f)
IF(ASSOCIATED(info%rf))             DEALLOCATE(info%rf)
IF(ASSOCIATED(info%ftmp))           DEALLOCATE(info%ftmp)
!
IF(ASSOCIATED(info%v1))             DEALLOCATE(info%v1)
IF(ASSOCIATED(info%v1old))          DEALLOCATE(info%v1old)
IF(ASSOCIATED(info%v1c))            DEALLOCATE(info%v1c)
IF(ASSOCIATED(info%f1))             DEALLOCATE(info%f1)
IF(ASSOCIATED(info%f1tmp))          DEALLOCATE(info%f1tmp)
IF(ASSOCIATED(info%rf1))            DEALLOCATE(info%rf1)
!
IF(ASSOCIATED(info%v2))             DEALLOCATE(info%v2)
IF(ASSOCIATED(info%v2old))          DEALLOCATE(info%v2old)
IF(ASSOCIATED(info%v2c))            DEALLOCATE(info%v2c)
IF(ASSOCIATED(info%f2))             DEALLOCATE(info%f2)
IF(ASSOCIATED(info%f2tmp))          DEALLOCATE(info%f2tmp)
IF(ASSOCIATED(info%rf2))            DEALLOCATE(info%rf2)
!
IF(ASSOCIATED(info%levellandscape)) DEALLOCATE(info%levellandscape)
!
! Nullify the pointers:
NULLIFY(info%errorflags,                                                 &
        info%q, info%qold, info%qc, info%qcold, info%f, info%ftmp, info%rf,          &
        info%v1,info%v1old,info%v1c,info%f1,info%f1tmp,info%rf1,         &
        info%v2,info%v2old,info%v2c,info%f2,info%f2tmp,info%rf2,         &
        info%qrte,info%aux,info%levellandscape)
!
If (ndims == 3) THEN
IF(ASSOCIATED(info%v3))             DEALLOCATE(info%v3)
IF(ASSOCIATED(info%v3old))          DEALLOCATE(info%v3old)
IF(ASSOCIATED(info%v3c))            DEALLOCATE(info%v3c)
IF(ASSOCIATED(info%f3))             DEALLOCATE(info%f3)
IF(ASSOCIATED(info%f3tmp))          DEALLOCATE(info%f3tmp)
IF(ASSOCIATED(info%rf3))            DEALLOCATE(info%rf3)

!
! Nullify the pointers:
NULLIFY(info%v3,info%v3old,info%v3c,info%f3,info%f3tmp,info%rf3)
!
END IF
info%fieldsallocated = .FALSE.  ! Don't compute on this node:
info%tobedeleted     = .TRUE.   ! Mark for garbage collection:
!
END SUBROUTINE DeAllocFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AllocPeriodicBCStorage(np)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, INTENT(IN):: np
!
INTEGER:: ierror
!
ALLOCATE(periodicoffsetindex(1:np),poffset(1:maxdims,1:np),STAT=ierror)
!
IF(ierror/=0) THEN
  PRINT *,'BSAM 2.0 Error in AllocPeriodicBCStorage.'
  STOP
END IF
!
END SUBROUTINE AllocPeriodicBCStorage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE DeallocPeriodicBCStorage
USE NodeInfoDef
IMPLICIT NONE
!
IF(ALLOCATED(periodicoffsetindex)) DEALLOCATE(periodicoffsetindex)
IF(ALLOCATED(poffset)) DEALLOCATE(poffset)
!
END SUBROUTINE DeallocPeriodicBCStorage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE AllocUniformGrids(mxroot)
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER, DIMENSION(1:maxdims), INTENT(IN):: mxroot ! root-level grid size.
!
INTEGER:: ierror, level, n, r
INTEGER, DIMENSION(1:maxdims):: high, low, mx
!
mx = 1; mx(1:ndims) = mxroot(1:ndims)
low = 1; low(1:ndims) = 1-mbc
high = 1
!
DO level = 0, maxlevel
!
  IF(level>0) mx(1:ndims) = mx(1:ndims)*2
  uniformgrid(level)%mx = 1; uniformgrid(level)%mx(1:ndims) = mx(1:ndims)
!
  high(1:ndims) = mx(1:ndims)+mbc
!
  ALLOCATE( &
    uniformgrid(level)%q(low(1):high(1),low(2):high(2),low(3):high(3),1:nccv), &
    STAT=ierror)
!
  IF(ierror/=0) THEN
    PRINT *, 'AllocUniformGrids: Error allocating uniformgrid on level', &
             level
    STOP
  END IF
!  
END DO
!
END SUBROUTINE AllocUniformGrids
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE DeallocUniformGrids
USE NodeInfoDef
IMPLICIT NONE
!
INTEGER:: level
!
DO level = rootlevel, maxlevel
  uniformgrid(level)%mx = 0
  IF(ASSOCIATED(uniformgrid(level)%q)) DEALLOCATE(uniformgrid(level)%q)
END DO
!
END SUBROUTINE DeallocUniformGrids
!
END MODULE BSAMStorage

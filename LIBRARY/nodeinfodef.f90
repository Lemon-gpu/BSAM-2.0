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
! (c) Copyright Steven M. Wise, 2015
! Department of Mathematics
! University of Tennessee
! swise@math.utk.edu
!
! Portions of the code
!
! (c) Copyright Sorin Mitran, 2002
! Department of Mathematics
! University of North Carolina at Chapel Hill
! mitran@amath.unc.edu
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
! File:             nodeinfodef.f90
! Purpose:          BSAM node data structures and global allocation module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 2.0 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------
MODULE NodeInfoDef
IMPLICIT NONE
!
SAVE
PRIVATE
!
LOGICAL, PUBLIC:: fluxbalancing, getafterstepstats, makeconformingmesh, &
                  meshbuildcomplete, outputinitialdata, outputuniformmesh, &
                  periodicboundaryconditions, restart, syncelliptic
!
! Double (r8) and extended (rext) precision if available:
INTEGER, PARAMETER, PUBLIC:: r8 = SELECTED_REAL_KIND(15,307), &
                             int2 = SELECTED_INT_KIND(2)
!
INTEGER, PARAMETER, PUBLIC:: errflagdefault = 1, &
                             errflaguser = 10, &
                             internalbc = 999, &
                             maxsubgrids = 1024, &
                             maxdims = 3, &
                             maxdepth = 10, &
                             maxnccv = 10, &
                             rootlevel = 0, &
                             sourcefield = 1, &
                             solutionfield = 2, &
                             auxiliaryfield = 3
LOGICAL, DIMENSION(0:maxdepth), PUBLIC:: defectivegridlevel
INTEGER, PUBLIC:: amrrestarts, errortype, finestlevel, gridnumber, maxvcycles, &
                  maxlevel, mbc, minlevel, ndims, nmasscorrlayers, naxv, &
                  nccv, nfcv, nperiodicoffsets, nrootgrids, nsmoothingpasses, &
                  ntaggedcells, outframes, restartframe, timeiterations, &
                  totalmeshsize, updateauxfreq
INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC:: periodicoffsetindex
INTEGER, DIMENSION(0:maxdepth), PUBLIC:: errflagopt, ibuffer, minimumgridpoints
INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC:: poffset
INTEGER, DIMENSION(0:maxdepth,1:maxdims), PUBLIC:: mxmax
!
REAL(KIND=r8), PUBLIC:: currenttime, dt, finaltime, omega, qerrortol, &
                        restarttime
REAL(KIND=r8), DIMENSION(1:2), PUBLIC:: integralresult
REAL(KIND=r8), DIMENSION(1:maxnccv), PUBLIC:: componentintegral
REAL(KIND=r8), DIMENSION(0:maxdepth), PUBLIC:: desiredfillratios, qtolerance
!
TYPE, PUBLIC:: taggedcell
  INTEGER:: id
  INTEGER, DIMENSION(1:maxdims):: coordinate
  TYPE(taggedcell), POINTER:: prevcell
END TYPE taggedcell
!
TYPE(taggedcell), POINTER, PUBLIC:: zerothtaggedcell
TYPE(taggedcell), POINTER, PUBLIC:: currenttaggedcell
TYPE(taggedcell), POINTER, PUBLIC:: lasttaggedcell
!
TYPE, PUBLIC:: masscorrlayer
  INTEGER:: id
  INTEGER, DIMENSION(1:maxdims,1:2):: mg
  REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: masscorr
  TYPE(masscorrlayer), POINTER:: prevlayer
END TYPE masscorrlayer
!
TYPE(masscorrlayer), POINTER, PUBLIC:: zerothlayer
TYPE(masscorrlayer), POINTER, PUBLIC:: currentlayer
TYPE(masscorrlayer), POINTER, PUBLIC:: lastlayer
!
! Uniform grids used for restarting and output:
TYPE, PUBLIC:: uniformgridtype
  INTEGER, DIMENSION(1:maxdims):: mx
  REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: q
END TYPE uniformgridtype
TYPE(uniformgridtype), DIMENSION(0:maxdepth), PUBLIC:: uniformgrid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE, PUBLIC:: nodeinfo
!
! This must be the first component to ensure proper parallel communication:
INTEGER:: nodeinfostart
!
! LOGICAL variables !!!!!!!!!
!
! A necessary component. Mark for garbage collection:
LOGICAL:: tobedeleted
!
! Flag to show whether grid may accept values from elder siblings on this 
!   level:
LOGICAL:: activegrid
!
! All patches are probationary until the mesh hierarchy is fully constructed:
LOGICAL:: defective
!
! Flag to show whether grid has been created during a restart from checkpoint 
!   file:
LOGICAL:: restartgrid
!
! Flag showing allocation status of fields within this node:
LOGICAL:: fieldsallocated
!
! Flag to show whether this is an initial grid, i.e. created during start-up:
LOGICAL:: initialgrid
!
! INTEGER variables !!!!!!!!!
!
INTEGER:: maxlevel               ! Maximum level to which this grid may be 
                                 !   refined:
INTEGER:: ngrid                  ! Number of this grid:
INTEGER:: nsubgrids              ! Number of child grids:
INTEGER:: level                  ! Level on which this node lives:
!
! Number of grid cells in q along each dimension:
INTEGER, DIMENSION(1:maxdims):: mx
!
! Index bounds within parent where this child was created:
INTEGER, DIMENSION(1:maxdims,1:2):: mbounds
!
! Index bounds of this grid in global indexing system:
INTEGER, DIMENSION(1:maxdims,1:2):: mglobal
!       
! Boundary condition codes:
!   1 - left, 2 - right, 3 - bottom, 4 - top, 5 - back, 6 - front:
INTEGER, DIMENSION(1:2*maxdims):: mthbc
!
! Array of error flags:
INTEGER, DIMENSION(:,:,:), POINTER:: errorflags
!
! Pointers to the grid's nearest neighbors' levels:
INTEGER(KIND=int2), DIMENSION(:,:,:), POINTER::  levellandscape
!
! REAL variables !!!!!!!!!!!!!!!!!!!!!!
!
! Lower coordinates for this grid:
REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
! Upper coordinates for this grid:
REAL(KIND=r8), DIMENSION(1:maxdims):: xupper
!
! The current time at which this grid exists:
REAL(KIND=r8):: gridtime
!
! Grid spacings:
REAL(KIND=r8), DIMENSION(1:maxdims):: dx
!
! Pointer to field variable arrays:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: q
!
! Pointer to the velocity variable on the 1-2 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v1
!
! Pointer to the velocity variable on the 3-4 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v2
!
! Pointer to the velocity variable on the 5-6 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v3
!
! Pointer to field variable arrays at previous time:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: qold
!
! Pointer to the velocity variable on the 1-2 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v1old
!
! Pointer to the velocity variable on the 3-4 faces at previous time:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v2old
!
! Pointer to the velocity variable on the 5-6 faces at previous time:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v3old
!
! Pointer to the coarse-under-fine field variable arrays:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: qc
!
! Pointer to the coarse-under-fine field variable arrays:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: qcold
!
! Pointer to the coarse-under-fine v1 variable array:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v1c
!
! Pointer to the coarse-under-fine v2 variable array:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v2c
!
! Pointer to the coarse-under-fine v2 variable array:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v3c
!
! Pointer to the coarse-under-fine relative truncation error:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: qrte
!
! Pointer to the load function:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f
!
! Pointer to the temporary function.  Same size as f and rf:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: ftmp
!
! Pointer to the load function on the 1-2 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f1
!
! Pointer to the load function on the 1-2 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f1tmp 
!
! Pointer to the load function on the 3-4 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f2
!
! Pointer to the load function on the 3-4 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f2tmp
!
! Pointer to the load function on the 3-4 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f3
!
! Pointer to the load function on the 3-4 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f3tmp
!
! Pointer to the cell-centered residual:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: rf
!
! Pointer to the residual on the 1-2 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: rf1
!
! Pointer to the residual on the 3-4 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: rf2
!
! Pointer to the residual on the 3-4 faces:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: rf3
!
! Pointer to auxiliary arrays:
REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: aux
!    
! This must be the last component to ensure proper parallel communication:
INTEGER:: nodeinfoend
!
END TYPE nodeinfo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE, PUBLIC:: funcparam
!
INTEGER:: iswitch
INTEGER, DIMENSION(1:maxdims):: offset
!
TYPE(nodeinfo), POINTER:: info
!
END TYPE funcparam
!
END MODULE NodeInfoDef

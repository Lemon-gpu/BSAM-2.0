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
!END IF
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
! File:             bsamroutines.f90
! Purpose:          BSAM control module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2006 Steven Wise
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 2.0 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------
MODULE BSAMRoutines
! ============================= 中文概览 =============================
! 模块：BSAMRoutines（BSAM 控制/调度模块）
! 作用：
! - 统一驱动 BSAM 2.0 的一次完整求解流程：读取输入、构建网格森林、
!   初始化、时间推进、AMR 细化/回收、输出与清理。
! - 封装核心调度子程序：BSAMSolver、TakeTimeSteps、AMR 等。
!
! 重要外部依赖（来自其它模块，按功能分组）：
! - 网格/树与全局参数：NodeInfoDef（定义 nodeinfo、全局开关与尺寸等）
! - 树操作：TreeOps（InitForest/KillForest、ApplyOnLevel(…)/Pairs、
!   CreateBelowSeedLevels、CreateChild、DeleteMarkedNode 等）
! - 存储与 IO：BSAMStorage（Alloc/DeAllocFields、周期边界存储）、
!   BSAMInputOutput（ReadQ/WriteQ/WriteUniformMeshQ）
! - 边界条件：Boundary（SetGhost、周期偏移等）
! - 多重网格：AFASRoutines（MultigridIterations、误差评估/限制/加密）
! - 物理问题：Problem（SetupProblem、Initialize2D/3D、SetAux/SetSrc、AfterRun）
!
! 文件内主要过程（简要说明）：
! - BSAMSolver：主入口；读入运行参数，建树和根网格，初始化字段，
!   执行时间推进（TakeTimeSteps），最终清理资源。
! - RootInit：读取 griddata.dat（Namelist）并构造根网格的几何与边界设置，
!   进行一致性检查并分配内存。
! - InitSeed：基于细网格回推生成多重网格低层“种子”网格的信息。
! - Initialize/InitializeFields：根据问题维度以用户回调方式初始化状态场。
! - SetAuxFields/SetSrcFields：设置辅助变量与源项。
! - CopyQToQold：保存本地 q 到 qold（包含面心速度的快照）。
! - TakeTimeSteps：时间推进的高层循环；每帧中做 AMR、FAS V-cycle、统计、输出。
! - AMR：递归的自适应网格细化/修复流程；标记误差、生成新子网格、转移场数据。
! - 误差标记相关：EstimateLevelErrors、ErrFlag、SetErrFlags2D/3D、
!   BufferAndList/InflateEdgeTags/BufferTaggedCells 等。
! - 网格生成/传输：NewSubGrids（改 Berger–Rigoutsos 切分）、
!   MakeNewGrid/InitFields、Transferq/TransferOverlap。
! - 其它工具：FindCoarseLevelNeighbors、FindMeshDefects 等。
!
! 约定：
! - 本文件所有新增为纯注释（以 ! 开头），不改变任何现有语句与逻辑。
! - 变量名/接口与原始代码保持一致，注释尽量解释“含义与使用时机”。
! ==================================================================
   IMPLICIT NONE
!
   SAVE
   PRIVATE
   PUBLIC BSAMSolver
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE BSAMSolver
! 子程序：BSAMSolver（全局驱动）
! 作用：
! - 读取运行参数（rundata.dat），初始化网格森林与根网格，
! - 初始化物理问题、边界与场，按时间推进执行 AFAS 多重网格迭代，
! - 周期性写出结果，最终回收内存并销毁森林。
      USE NodeInfoDef
      USE TreeOps, ONLY: AddRootLevelNode, ApplyOnForest, ApplyOnLevel, &
         CreateBelowSeedLevels, DeleteMarkedNode, InitForest, &
         KillForest
      USE BSAMInputOutput, ONLY: ReadQ
      USE BSAMStorage, ONLY: DeallocPeriodicBCStorage
      USE Problem, ONLY: SetupProblem, AfterRun
      USE Boundary, ONLY: SetGhost
      IMPLICIT NONE
!
      TYPE(funcparam):: dummy
      CHARACTER(LEN=5):: zone
      CHARACTER(LEN=8):: date
      CHARACTER(LEN=10):: time
      INTEGER:: i, ierror, level, ilevel
      INTEGER, DIMENSION(1:8):: values
!
! 局部变量说明：
! - dummy            ：在 ApplyOnLevel/… 调用中承载简单开关/参数的通用结构。
! - zone/date/time   ：DATE_AND_TIME 返回的时区/日期/时间字符串。
! - values(1:8)      ：DATE_AND_TIME 返回的年月日时分秒等数值部件。
! - i,ierror,level… ：通用整型临时变量/错误码/层级游标。
!
      NAMELIST/rundata/dt, errortype, getafterstepstats, maxvcycles, &
         nsmoothingpasses, omega, outframes, outputuniformmesh, &
         qerrortol, restart, restartframe, syncelliptic, &
         timeiterations, updateauxfreq

! Namelist [rundata] 关键参数释义：
! - dt               ：时间步长（若 syncelliptic 为真，初始会临时置 0）。
! - errortype        ：误差类型选择（用于 AMR/多重网格中的误差评估策略）。
! - getafterstepstats：是否统计每步后的积分量/网格规模等数据。
! - maxvcycles       ：每个时间步执行的 AFAS V-cycle 最大次数。
! - nsmoothingpasses ：每次松弛迭代的遍数（多重网格平滑器参数）。
! - omega            ：松弛参数（如加权 Jacobi/Gauss–Seidel 的权重）。
! - outframes        ：输出帧的总数（时间推进将分配等量步数到这些帧）。
! - outputuniformmesh：是否同时输出统一网格（便于可视化/重启）。
! - qerrortol        ：解向量 L2 残差阈值（控制 V-cycle 收敛/提前退出）。
! - restart          ：是否从已有输出（重启帧）继续计算。
! - restartframe     ：重启时读取的帧编号。
! - syncelliptic     ：是否在起步阶段对椭圆变量做一次同步（dt=0 的步骤）。
! - timeiterations   ：总时间步数（会按 outframes 等分到每帧）。
! - updateauxfreq    ：辅助场更新的频率（每多少 V-cycle 或步数更新一次）。
!
! Initializations:
      errortype = 1
! dt = 0.0001_r8
      getafterstepstats = .FALSE.
      maxvcycles = 20
      nsmoothingpasses = 2
      omega = 1.0_r8
      outframes = 1
      outputuniformmesh = .FALSE.
      qerrortol = 1.0E-06_r8
      restart = .FALSE.
      restartframe = 0
      syncelliptic = .FALSE.
      timeiterations = 1
      updateauxfreq = 1
!
! Read general input data:
      OPEN(UNIT=75,FILE='rundata.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
      ! OPEN() 这里的每一个参数的意思是：
         ! UNIT: 指定文件单元号（75），文件单元号的意思是：
            ! 在 Fortran 中，文件单元号是一个整数，用于标识和管理打开的文件。
            ! 每个打开的文件都需要一个唯一的文件单元号，以便程序能够正确地读写数据。
            ! 通过指定文件单元号，程序可以区分不同的文件，并执行相应的输入输出操作。
         ! FILE: 指定要打开的文件名（'rundata.dat'），该文件包含程序运行所需的参数和配置。
         ! STATUS: 指定文件的状态（'OLD'），表示该文件必须已经存在，否则会报错。
         ! ACTION: 指定对文件的操作类型（'READ'），表示该文件将被读取。
         ! IOSTAT: 指定一个整数变量（ierror），用于捕获文件操作的错误状态。
            ! 如果文件操作成功，ierror 将被设置为 0；
            ! 如果发生错误，ierror 将被设置为一个非零值，表示具体的错误类型。
            ! 通过检查 ierror 的值，程序可以判断文件操作是否成功，并采取相应的措施（如报错或继续执行）。
      IF(ierror/=0) THEN
         PRINT *,'Error opening input file rundata.dat. Program stop.'
         STOP
      END IF
      READ(UNIT=75,NML=rundata)
      CLOSE(75)
      PRINT *,'dt', dt
!
      CALL DATE_AND_TIME(date,time,zone,values)
      OPEN(UNIT=76,FILE='output.dat',STATUS='UNKNOWN',ACTION='WRITE', &
         FORM='FORMATTED',POSITION='APPEND')
      ! STATUS='UNKNOWN'：如果文件不存在则创建，存在则打开并追加内容。
      ! FORM='FORMATTED'：指定文件为文本格式（而非二进制）。
      ! POSITION='APPEND'：将新内容追加到文件末尾，而不是覆盖原有内容。
      WRITE(76,1001) date, time
      ! Write() 这里的每一个参数的意思是：
         ! 76: 指定文件单元号（76），表示将数据写入该单元号对应的文件。
         ! 1001: 指定格式标签，表示使用标签为 1001 的格式进行数据输出。
            ! 该格式标签在后续的 FORMAT 语句中定义，控制输出数据的布局和格式。
         ! date, time: 指定要写入的数据变量，这里是两个字符串变量，分别表示当前的日期和时间。
1001  FORMAT(' '/'New run at date ',A8,' and time ',A10/' ')
! 格式标签 1001 定义了输出的具体格式：
! - ' '：输出一个空格。
! - 'New run at date '：输出字符串 "New run at date "。
! - A8：输出一个长度为 8 的字符串（对应 date 变量）。
! - ' and time '：输出字符串 " and time "。
! - A10：输出一个长度为 10 的字符串（对应 time 变量）。
! - '/'：换行符，表示输出结束后换行。
! - ' '：输出一个空格。
      WRITE(76,NML=rundata)
      CLOSE(76)
! 注：以上将本次运行的关键参数追记到 output.dat，便于复现实验。
!
! By default only one rootlevel grid.  This may change in the future:
      nrootgrids = 1
! 目前默认只有一个根网格（rootlevel）。如需多根网格，需扩展相关逻辑。
!
! Initialize a forest of trees.  One root level grid generated in this call:
      CALL InitForest
! 初始化“森林”数据结构（层级树）。该调用内会创建至少一个根层节点。
!
! Read data file to initialize root level grids:
      CALL ApplyOnLevel(rootlevel,RootInit,dummy)
! RootInit：读取 griddata.dat 并构造根网格（尺寸/边界/坐标与存储分配）。
!
! Create the levels below the seed needed for multigrid:
      CALL CreateBelowSeedLevels(minlevel)
! 为多重网格准备更粗的“种子层”（minlevel 通常为负或 0）。
!
! Initialize the levels below the forest seed:
      DO ilevel = rootlevel-1, minlevel, -1
         CALL ApplyOnLevel(ilevel,InitSeed,dummy)
      END DO
! InitSeed：由更细层向下推导 coarse 层几何/索引信息并分配存储。
!
! Set user problem parameters:
      CALL SetupProblem
      PRINT *,'dt', dt
! 调用户问题模块设置方程/物理常数等（可能调整 dt）。
!
! Initialize the rootlevel data fields.  In the case of a restart, we read the
! fields at all above root levels, saving the data to uniformgrid(level)%q':
      IF(restart) THEN
         PRINT *, 'Restart of computation from plot frame ', restartframe
         CALL ReadQ
         outputinitialdata = .FALSE.
         syncelliptic = .FALSE.
      ELSE
         CALL ApplyOnLevel(rootlevel,InitializeFields,dummy)
         outputinitialdata = .TRUE.
      END IF
! 重启：从输出文件装载多层数据；新开：调用 Initialize… 以用户自定义方式初始化。
!
      CALL SetGhost(rootlevel,0)
      CALL ApplyOnLevel(rootlevel,SetAuxFields,dummy)
      CALL ApplyOnLevel(rootlevel,SetSrcFields,dummy)
      CALL ApplyOnLevel(rootlevel,CopyQToQold,dummy)
! 设置根层边界、辅助变量、源项，并保存初始快照（q -> qold）。
!
! Initialization complete. Start run:
      PRINT *, 'BSAM 2.0 is running ... '
      PRINT *, ' '
!
      PRINT *,'dt', dt
! (Core) Carry out the time steps:
      CALL TakeTimeSteps
! 时间推进主循环（见 TakeTimeSteps）。
!
      PRINT *,'dt', dt
! User-specified actions before program ends:
      CALL AfterRun
! 用户自定义的收尾动作（统计/导出等）。
!
! Delete the below-seed-level grids:
      DO level = minlevel, maxlevel
         CALL ApplyOnLevel(level,MarkNodeInactive,dummy)
         CALL ApplyOnLevel(level,ReleaseInactiveFields,dummy)
         PRINT *,' Level ',level,' fields have been released.'
      END DO
! 回收所有层的场内存（将仍存在的非活动网格释放其字段）。
!
! Delete the forest of trees:
      CALL KillForest
! 销毁森林（释放树结构节点）。
!
! Delete forest seed and below-seed levels:
      CALL ApplyOnLevel(minlevel,MarkNodeToBeDeleted,dummy)
      CALL ApplyOnLevel(minlevel,DeleteMarkedNode,dummy)
! 确保最粗层节点被标记并真正删除。
!
      CALL DeallocPeriodicBCStorage
! 释放周期边界辅助存储。
!
   END SUBROUTINE BSAMSolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION MarkNodeToBeDeleted(info,dummy)
! 功能：标记当前网格节点为“待删除”。
! 参数：
! - info   ：节点信息（值引用传递，此处修改其 tobedeleted 标志）。
! - dummy  ：通用参数（未使用）。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      MarkNodeToBeDeleted = err_ok
!
      info%tobedeleted = .TRUE.
!
   END FUNCTION MarkNodeToBeDeleted
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION MarkNodeInactive(info,dummy)
! 功能：将网格标记为非活动（activegrid = .FALSE.），以便稍后释放字段或删除。
! 说明：AMR 重建新网格后，旧网格通常会被标成非活动，等待回收。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      MarkNodeInactive = err_ok
!
      info%activegrid = .FALSE.
!
   END FUNCTION MarkNodeInactive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION MarkNodeNonInitial(info,dummy)
! 功能：将网格的 initialgrid 标志设为 .FALSE.，表示已不再是“初始网格”。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      MarkNodeNonInitial = err_ok
!
      info%initialgrid = .FALSE.
!
   END FUNCTION MarkNodeNonInitial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION RootInit(rootinfo,dummy)
! 功能：读取 griddata.dat 的 Namelist，初始化根层网格的尺寸/边界/坐标等基本信息，
!       并执行一致性检查与存储分配。
! 关键 Namelist [griddata] 字段说明（与 NodeInfoDef 全局变量关联）：
! - ndims              ：维数（2 或 3）。
! - mx(1:ndims)        ：根网格各向单元数（需为 2 的倍数）。
! - mglobal(:,1:2)     ：根网格在全局索引系中的左右/上下(/前后)范围。
! - xlower/xupper      ：物理域边界坐标（各向）。
! - mbc                ：幽灵层层数（仅支持 1 或 2）。
! - nccv/nfcv/naxv     ：单元中心/面心/辅助变量数量。
! - mthbc(1:2*ndims)   ：边界条件编码（2 表示周期，其它为物理/内部）。
! - desiredfillratios  ：期望填充率（AMR 子网格生成阈值参考）。
! - minimumgridpoints  ：子网格的最小单元数（偶数，不小于 4）。
! - maxlevel/minlevel  ：多重网格层级上限/下限。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      USE Boundary, ONLY: PeriodicSetup
      USE BSAMStorage, ONLY: AllocFields
      IMPLICIT NONE
!
      TYPE(nodeinfo):: rootinfo
      TYPE(funcparam):: dummy
!
      CHARACTER(LEN=12):: filename
      INTEGER:: i, ierror, n
      INTEGER, DIMENSION(1:maxdims):: mx
      INTEGER, DIMENSION(1:2*maxdims):: mthbc
      INTEGER, DIMENSION(1:maxdims,1:2):: mglobal
      REAL(KIND=r8), DIMENSION(1:maxdims):: dx, xlower, xupper
!
      NAMELIST/griddata/desiredfillratios, errflagopt, ibuffer, maxlevel, mbc, &
         mglobal, minimumgridpoints, minlevel, mthbc, mx, ndims, &
         naxv, nccv, nfcv, qtolerance, xlower, xupper
!
      PRINT *, 'Reading grid data for root-level grid.'
      ! 这里的`*`表示标准输出设备，通常是控制台或终端。
!
! Set default values !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      desiredfillratios = 8.0E-01_r8
      errflagopt = errflagdefault
      ibuffer = 0
      maxlevel = 0
      mbc = 1
      mglobal = 1
      minimumgridpoints = 2
      minlevel = 0
      mthbc = 10
      mx = 1
      ndims = 2
      naxv = 0
      nccv = 1
      nfcv = 0
      qtolerance = 1.0E-06_r8
      xlower = 0.0_r8
      xupper = 1.0_r8
!
! Read from namelist !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      OPEN(UNIT=75,FILE='griddata.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
      IF(ierror/=0) THEN
         PRINT *,'Error opening input file griddata.dat. Program stop.'
         STOP
      END IF
      READ(75,NML=griddata)
      CLOSE(75)
      OPEN(UNIT=76,FILE='output.dat',STATUS='OLD',ACTION='WRITE',FORM='FORMATTED', &
         POSITION='APPEND')
      WRITE(76,NML=griddata)
      CLOSE(76)
!
! Check input variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      IF(ANY(desiredfillratios<0.0) .OR. ANY(desiredfillratios>1.0)) THEN
         PRINT *, 'BSAM 2.0 Error: Code only supports desiredfilratios between 0 and 1.'
         STOP
      END IF
!
      IF(ANY(ibuffer<0)) THEN
         PRINT *, 'BSAM 2.0 Error: Code only supports ibuffers>=0.'
         STOP
      END IF
!
      IF(maxlevel<0) THEN
         PRINT *, 'BSAM 2.0 Error: Code only supports maxlevel>=0.'
         STOP
      END IF
!
      IF(naxv<0) THEN
         PRINT *, 'BSAM 2.0 Error: Code only supports naxv>=0.'
         STOP
      END IF
!
      IF(mbc<1 .OR. mbc>2) THEN
         PRINT *, 'BSAM 2.0 Error: Code only supports mbc=1,2.'
         STOP
      END IF
!
      IF(minlevel>0) THEN
         PRINT *, 'BSAM 2.0 Error: Code only supports minlevel<=0.'
         STOP
      END IF
!
      IF(ndims<2 .OR. ndims>3) THEN
         PRINT *, 'BSAM 2.0 Error: Code only supports ndims=2,3.'
         STOP
      END IF
!
      DO i = 1, 2*ndims-1, 2
         IF(nrootgrids==1 .AND. &
            mthbc(i)==2 .AND. mthbc(i+1)/=2) THEN
            PRINT *, 'BSAM 2.0 Error: Incorrect periodic boundary conditions.'
            STOP
         END IF
         IF(nrootgrids==1 .AND. &
            mthbc(i+1)==2 .AND. mthbc(i)/=2) THEN
            PRINT *, 'BSAM 2.0 Error: Incorrect periodic boundary conditions.'
            STOP
         END IF
      END DO
!
      DO i = 1, ndims
         IF(mx(i)<=0) THEN
            PRINT *, 'BSAM 2.0 Error: mx<=0 along dim', i, '.'
            STOP
         END IF
      END DO
!
      DO i = 1, ndims
         IF(MODULO(mx(i),2)/=0) THEN
            PRINT *, 'BSAM 2.0 Error: Initial grid must have dimensions which are a multiple'
            PRINT *, '       of the refinement ratio 2.'
            STOP
         END IF
      END DO
!
      DO i = 1, ndims
         IF(mglobal(i,2)-mglobal(i,1)+1/=mx(i)) THEN
            PRINT *, 'BSAM 2.0 Error: mglobal(i,2)-mglobal(i,1)+1/=mx(i), i=', i, '.'
            STOP
         END IF
      END DO
!
      IF(ANY(minimumgridpoints<1)) THEN
         PRINT *, 'BSAM 2.0 Error: Code only supports minimumgridpoints>=1.'
         STOP
      END IF
!
      IF(nccv>maxnccv .OR. nccv<1) THEN
         PRINT *, 'BSAM 2.0 Error: Code only supports nccv<=maxnccv and nccv>=1.'
         STOP
      END IF
!
      IF(ANY(qtolerance<=0.0)) THEN
         PRINT *, 'BSAM 2.0 Error: Code only supports qtolerance>0.0.'
         STOP
      END IF
!
      DO i = 1, ndims
         IF(xupper(i)<=xlower(i)) THEN
            PRINT *, 'BSAM 2.0 Error: Code only supports xupper > xlower.'
            STOP
         END IF
      END DO
!
      dx = 0.0_r8
      dx(1:ndims) = (xupper(1:ndims)-xlower(1:ndims))/REAL(mx(1:ndims),KIND=r8)
!
      DO i = 2, ndims
         IF(ABS(dx(1)-dx(i))>1.0E-10_r8) THEN
            PRINT *, 'BSAM 2.0 Error: dx(1)\=dx(i) along dim', i, '.'
            STOP
         END IF
      END DO
!
      mxmax = 1; mxmax(0,1:ndims) = mx(1:ndims)
      DO i = rootlevel+1, maxlevel
         mxmax(i,1:ndims) = 2*mxmax(i-1,1:ndims)
      END DO
!
! Copy data to global storage !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      rootinfo%ngrid = 0
      rootinfo%level = rootlevel
      rootinfo%maxlevel = maxlevel
!
      rootinfo%tobedeleted = .FALSE.
      rootinfo%initialgrid = .TRUE.
      rootinfo%activegrid = .TRUE.
      rootinfo%defective = .FALSE.
!
      rootinfo%mx = 1
      rootinfo%mx(1:ndims) = mx(1:ndims)
      rootinfo%mglobal = 1
      rootinfo%mglobal(1:ndims,1:2) = mglobal(1:ndims,1:2)
!
      rootinfo%gridtime = 0.0_r8
!
      rootinfo%xlower = 0.0_r8
      rootinfo%xlower(1:ndims) = xlower(1:ndims)
      rootinfo%xupper = 0.0_r8
      rootinfo%xupper(1:ndims) = xupper(1:ndims)
      rootinfo%dx = 0.0_r8
      rootinfo%dx(1:ndims) = dx(1:ndims)
      rootinfo%mthbc = 10
      rootinfo%mthbc(1:2*ndims) = mthbc(1:2*ndims)
!
      rootinfo%mbounds = 1; rootinfo%mbounds(1:ndims,2) = rootinfo%mx(1:ndims)
!
! Allocate storage !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Set up the periodic offsets used in transferring periodic bc's:
      CALL PeriodicSetup(rootinfo)
!
      CALL AllocFields(rootinfo)
      rootinfo%levellandscape = 0
!
! Finished initialization of root node info structure
      RootInit = err_ok
!
      PRINT *, 'Finished reading root-level grid data.'
!
   END FUNCTION RootInit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION InitSeed(seedinfo,dummy)
! 功能：基于较细层子网格的信息构造上一层（更粗层，level-1）的“种子”网格信息，
!       用于多重网格的层级构建（仅保留多重网格所需的信息）。
! 说明：
! - 本过程会检查细网格尺寸是否能被 2 整除（保证可粗化）。
! - 设置 mx、dx、mglobal、mbounds 等，并分配基本存储。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok, GetChildInfo
      USE BSAMStorage, ONLY: AllocFields
      IMPLICIT NONE
!
      TYPE(nodeinfo):: seedinfo
      TYPE(funcparam) :: dummy
!
      TYPE(nodeinfo), POINTER:: child
      INTEGER:: ierror, level
!
! Only data needed for multigrid need be copied from the child grid.  Parent
! (seedinfo) is born from child:
!
      ierror = GetChildInfo(child)
!
      level = child%level-1
!
      seedinfo%maxlevel = maxlevel
      seedinfo%level = level
!
      PRINT *, 'level=', level, 'ndims=', ndims
!
      IF(ANY(MODULO(child%mx(1:ndims),2)/=0)) THEN
         PRINT *, 'BSAM 2.0 Error in InitSeed: grid on level', child%level, 'will not coarsen.'
         STOP
      END IF
!
      seedinfo%mx(1:ndims) = child%mx(1:ndims)/2
      seedinfo%dx(1:ndims) = child%dx(1:ndims)*REAL(2,KIND=r8)
      seedinfo%mglobal(1:ndims,1) = 1
      seedinfo%mglobal(1:ndims,2) = child%mglobal(1:ndims,2)/2
      child%mbounds(1:ndims,1) = 1; child%mbounds(1:ndims,2) = seedinfo%mx(1:ndims)
!
      seedinfo%mbounds(1:ndims,1) = 1; seedinfo%mbounds(1:ndims,2) = seedinfo%mx(1:ndims)
!
      IF(ANY(seedinfo%mx(1:ndims)/=seedinfo%mglobal(1:ndims,2))) THEN
         PRINT *, 'BSAM 2.0 Error in InitSeed: on level-1', seedinfo%level
         PRINT *, 'mx(:)/=mglobal(:,2)'
         STOP
      END IF
!
      seedinfo%ngrid = 0
!
      seedinfo%tobedeleted = .FALSE.; seedinfo%level = child%level-1
      seedinfo%initialgrid = .FALSE.
      seedinfo%activegrid = .TRUE.
      seedinfo%defective = .FALSE.
!
      seedinfo%gridtime = child%gridtime
!
      seedinfo%xlower = child%xlower; seedinfo%xupper = child%xupper
      seedinfo%mthbc = child%mthbc
!
      CALL AllocFields(seedinfo)
!
! Finished initialization of seed node info structure:
      InitSeed = err_ok
!
   END FUNCTION InitSeed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION SetAuxFields(info,dummy)
! 功能：若 naxv>0，则调用问题模块 SetAux(info) 设置辅助（aux）字段。
! 条件：跳过待删除或非活动网格。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      USE Problem, ONLY: SetAux
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      SetAuxFields = err_ok
      IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
      IF(naxv>0) CALL SetAux(info)
!
   END FUNCTION SetAuxFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION SetSrcFields(info,dummy)
! 功能：调用问题模块 SetSrc(info) 设置源项（source）字段。
! 条件：跳过待删除或非活动网格。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      USE Problem, ONLY: SetSrc
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      SetSrcFields = err_ok
      IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
      CALL SetSrc(info)
!
   END FUNCTION SetSrcFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION InitializeFields(info,dummy)
! 功能：调用 Initialize 子程序，按维度分派到 Initialize2D/3D 来生成初始场。
! 条件：跳过待删除或非活动网格。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      InitializeFields = err_ok
      IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
      CALL Initialize(info)
!
   END FUNCTION InitializeFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE Initialize(info)
! 功能：根据 ndims 分派到问题模块的 Initialize2D/Initialize3D，
!       在局部（含幽灵层或不含，视接口而定）数组片段上填充初值。
! 变量：
! - mx    ：本网格各向单元数（来自 info%mx）。
! - h     ：网格步长（info%dx(1)）。
! - xlower：左下(后)角坐标（info%xlower）。
      USE NodeInfoDef
      USE PROBLEM, ONLY: Initialize2D, Initialize3D
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
!
      INTEGER, DIMENSION(1:maxdims):: mx
      REAL(KIND=r8):: h
      REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
      mx(1:ndims) = info%mx(1:ndims)
      h = info%dx(1)
      xlower(1:ndims) = info%xlower(1:ndims)
!
      SELECT CASE(ndims)
       CASE(2)
         CALL Initialize2D(info%q (1:mx(1),1:mx(2),1      ,1:nccv), &
            info%v1(0:mx(1),1:mx(2),1      ,1:nfcv), &
            info%v2(1:mx(1),0:mx(2),1      ,1:nfcv), &
            mx(1:2),h,xlower(1:2))
       CASE(3)
         CALL Initialize3D(info%q (1:mx(1),1:mx(2),1:mx(3),1:nccv), &
            info%v1(0:mx(1),1:mx(2),1:mx(3),1:nfcv), &
            info%v2(1:mx(1),0:mx(2),1:mx(3),1:nfcv), &
            info%v3(1:mx(1),1:mx(2),0:mx(3),1:nfcv), &
            mx(1:3),h,xlower(1:3))
      END SELECT
!
   END SUBROUTINE Initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION CopyQToQold(info,dummy)
! 功能：保存当前解/速度到对应的 *old 缓存（q->qold、v*->v*old，及 qc->qcold）。
! 场景：时间推进或网格转移前的备份，以便计算截断误差或回滚。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      CopyQToQold = err_ok
      IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
      SELECT CASE(ndims)
       CASE(2)
         info% qold = info%q
         info%qcold = info%qc
!
         info%v1old = info%v1
         info%v2old = info%v2
       CASE(3)
         info% qold = info%q
         info%qcold = info%qc
!
         info%v1old = info%v1
         info%v2old = info%v2
         info%v3old = info%v3
      END SELECT
!
   END FUNCTION CopyQToQold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE TakeTimeSteps
! 功能：时间推进主循环（按 outframes 帧输出，timeiterations 总步数）。
! 流程概要：
! 1) 重启/新运行的起始时间与帧区间设置；必要时 syncelliptic（dt=0 一步）。
! 2) 对每个输出帧：
!    - 将旧层网格标为非活动，执行 AMR 构建新层级；
!    - 删除非活动网格字段并真正删除节点；
!    - 首次 AMR 后自上而下 FillDown 以补齐粗层的 qold；
!    - 保存当前层级数据到 qold；
!    - 可选：输出初始帧数据与统计；
!    - 执行 MultigridIterations（AFAS V-cycle）；
!    - 更新时间 currenttime，并可选统计；
!    - 若完成 syncelliptic，则恢复 dt；
!    - 为各层设置 info%gridtime=currenttime；
!    - 达到帧间步数后输出本帧数据（WriteQ/WriteUniformMeshQ）。
      USE NodeInfoDef
      USE TreeOps, ONLY: ApplyOnLevel, DeleteMarkedNode
      USE BSAMInputOutput, ONLY: WriteQ, WriteUniformMeshQ
      USE AFASRoutines, ONLY: FillDown, MultigridIterations
      IMPLICIT NONE
!
! currenttime is the time of the rootlevel grid.
! restarttime is the time at restart.
! finaltime is the calculated final time.
!
      TYPE(funcparam):: dummy
      LOGICAL:: firstamr
      INTEGER:: firstframe, it, itperprint, lastframe, level, n
      REAL(KIND=r8):: starttime, dtsave
!
      dtsave = dt
!
      IF(restart) THEN
         IF(restartframe>=outframes) THEN
            PRINT *, 'BSAM 2.0 Error in rundata.dat: outframes=', outframes, ' restartframe=', &
               restartframe
            PRINT *, 'outframes must be greater than restartframe.'
            STOP
         END IF
!
         currenttime = restarttime
         firstframe = restartframe+1; lastframe = outframes
!
         PRINT 1002, restartframe, outframes
1002     FORMAT('================================================'/ &
            ' BSAM 2.0: Restart frame ', i4,' out of ', i4,' requested frames'/ &
            '================================================')
      ELSE
         currenttime = 0.0_r8
         firstframe = 1; lastframe = outframes
         IF(syncelliptic) THEN
            dt = 0.0_r8
         END IF
      END IF
!
      itperprint = timeiterations/outframes
!
      starttime = currenttime
      finaltime = starttime &
         + REAL((lastframe-firstframe+1)*itperprint,KIND=r8)*dtsave
!
      firstamr = .TRUE.
!
      printloop: DO n = firstframe, lastframe
!
         it = 1
         timesteploop: DO
!
! For testing purposes only:
!
! Print out grid information for memory checking:
!    DO level = rootlevel+1, maxlevel
!      PRINT *, 'Grid information, level =', level
!      gridnumber = 0
!      CALL ApplyOnLevel(level,PrintGridInfo,dummy)
!      READ(*,*)
!    END DO
!
! Mark all above-root-level grids from the old mesh inactive:
            DO level = rootlevel+1, maxlevel
               CALL ApplyOnLevel(level,MarkNodeInactive,dummy)
            END DO
!
! Make a new mesh and move the solution to the new mesh:
            finestlevel = rootlevel
            defectivegridlevel(0:maxlevel) = .FALSE.
            amrrestarts = 0; meshbuildcomplete = .FALSE.
!
            CALL AMR(rootlevel)
!
! After exit from AMR delete all inactive grids (the previous mesh):
            meshbuildcomplete = .TRUE.
            DO level = maxlevel, rootlevel+1, -1
               CALL ApplyOnLevel(level,ReleaseInactiveFields,dummy)
               CALL ApplyOnLevel(level,DeleteMarkedNode,dummy)
            END DO
!
! This ensures that the below rootlevel values of qold are filled after a
! clean start or restart:
            IF(firstamr) THEN
               CALL FillDown(solutionfield)
               firstamr = .FALSE.
               CALL ApplyOnLevel(rootlevel,MarkNodeNonInitial,dummy)
            END IF
!
! Save a copy of data at the last time step:
            DO level = minlevel, finestlevel
               CALL ApplyOnLevel(level,CopyQToQold,dummy)
            END DO
!
! If initial grid, write data:
            IF(outputinitialdata .AND. (.NOT. syncelliptic)) THEN
!
               outputinitialdata = .FALSE.
               PRINT 799
799            FORMAT('==============================='/ &
                  ' BSAM 2.0: Writing initial data'/ &
                  '===============================')
               IF(getafterstepstats) THEN
                  integralresult = 0.0_r8
                  totalmeshsize = 0
                  DO level = finestlevel, rootlevel, -1
                     CALL ApplyOnLevel(level,AfterStepStatistics,dummy)
                     CALL ApplyOnLevel(level,GetMeshSize,dummy)
                  END DO
                  OPEN(UNIT=65,FILE='OUT/stats.dat',STATUS='UNKNOWN',ACTION='WRITE', &
                     FORM='FORMATTED',POSITION='APPEND')
                  WRITE(65,'(3(F25.12),1x,I8)') currenttime, integralresult(1:2), &
                     totalmeshsize
                  CLOSE(65)
               END IF
               totalmeshsize = 0
               DO level = finestlevel, rootlevel, -1
                  CALL ApplyOnLevel(level,GetMeshSize,dummy)
               END DO
               PRINT *, 'Initial composite mesh size =', totalmeshsize
               CALL WriteQ(0,currenttime)
               IF(outputuniformmesh) CALL WriteUniformMeshQ(0,currenttime)
!
            END IF
!
            IF(syncelliptic) THEN
               PRINT *, ' '
               PRINT *, ' '
               PRINT *, 'BSAM 2.0: Synchronizing elliptic fields at start.'
               PRINT *, ' '
            ELSE
               PRINT *, ' '
               PRINT *, ' '
               PRINT 8001, currenttime+dt, finaltime
8001           FORMAT('BSAM 2.0: Advancing to time =', ES15.7, 1X, 'of', ES15.7)
               PRINT *, ' '
            END IF
!
! For testing purposes only:
!
! Print out grid information for memory checking:
!    DO level = rootlevel+1, maxlevel
!      PRINT *, 'Grid information, level =', level
!      gridnumber = 0
!      CALL ApplyOnLevel(level,PrintGridInfo,dummy)
!      READ(*,*)
!    END DO
!
! Perform Multigrid on the multilevel mesh:
            CALL MultigridIterations
!
            currenttime = starttime+REAL((n-firstframe)*itperprint+it,KIND=r8)*dt
!
! Calculate after step statistics:
            IF(getafterstepstats .AND. (.NOT. syncelliptic)) THEN
!
               integralresult = 0.0_r8
               totalmeshsize = 0
               DO level = finestlevel, rootlevel, -1
                  CALL ApplyOnLevel(level,AfterStepStatistics,dummy)
                  CALL ApplyOnLevel(level,GetMeshSize,dummy)
               END DO
               OPEN(UNIT=65,FILE='OUT/stats.dat',STATUS='UNKNOWN',ACTION='WRITE', &
                  FORM='FORMATTED',POSITION='APPEND')
               WRITE(65,'(3(F25.12),1x,I8)') currenttime, integralresult(1:2), &
                  totalmeshsize
               CLOSE(65)
            END IF
!
            IF(syncelliptic) THEN
               syncelliptic = .false.
               dt = dtsave
               totalmeshsize = 0
               DO level = finestlevel, rootlevel, -1
                  CALL ApplyOnLevel(level,GetMeshSize,dummy)
               END DO
               PRINT *, 'BSAM 2.0: Synchronization composite mesh size =', totalmeshsize
            ELSE
               it = it+1
            END IF
!
! Set current time on the grid:
            DO level = finestlevel, rootlevel, -1
               CALL ApplyOnLevel(level,SetCurrentTime,dummy)
            END DO
!
            IF(it>itperprint) EXIT timesteploop
!
         END DO timesteploop
!
         totalmeshsize = 0
         DO level = finestlevel, rootlevel, -1
            CALL ApplyOnLevel(level,GetMeshSize,dummy)
         END DO
         PRINT *, 'Composite mesh size =', totalmeshsize
!
         PRINT 1001, n, outframes, currenttime
1001     FORMAT( &
            '===================================================================='/ &
            ' Writing frame ',I3,' out of ',I3,' requested frames at t=',ES15.7,   / &
            '====================================================================')
!
         CALL WriteQ(n,currenttime)
         IF(outputuniformmesh) CALL WriteUniformMeshQ(n,currenttime)
!
      END DO printloop
!
   END SUBROUTINE TakeTimeSteps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION PrintGridInfo(info,dummy)
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      PrintGridInfo = err_ok
!
      gridnumber = gridnumber+1
      PRINT *, ' '
      PRINT *, 'Grid number =', gridnumber, 'level=', info%level
      PRINT *, 'mx =', info%mx
      PRINT *, 'mb =', info%mbounds
      PRINT *, 'ngrid', info%ngrid
      PRINT *, 'activegrid', info%activegrid
      PRINT *, 'initialgrid', info%initialgrid
      PRINT *, 'tobedeleted', info%tobedeleted
      PRINT *, 'fieldsallocated', info%fieldsallocated
      PRINT *, ' '
!
   END FUNCTION PrintGridInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION GetMeshSize(info,dummy)
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      INTEGER, DIMENSION(1:maxdims):: mx, cmx
!
      GetMeshSize = err_ok
      IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
      mx(1:ndims) = info%mx(1:ndims)
      cmx(1:ndims) = mx(1:ndims)/2
!
      IF(info%level==0) THEN
         totalmeshsize = totalmeshsize+PRODUCT(mx(1:ndims))
      ELSE
         totalmeshsize = totalmeshsize+PRODUCT(mx(1:ndims))-PRODUCT(cmx(1:ndims))
      END IF
!
   END FUNCTION GetMeshSize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION ReleaseInactiveFields(info,dummy)
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      USE BSAMStorage, ONLY: DeAllocFields
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      ReleaseInactiveFields = err_ok
!
      IF(.NOT. info%activegrid) CALL DeAllocFields(info)
!
   END FUNCTION ReleaseInactiveFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION ReleaseActiveFields(info,dummy)
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      USE BSAMStorage, ONLY: DeAllocFields
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      ReleaseActiveFields = err_ok
!
      IF(info%activegrid) CALL DeAllocFields(info)
!
   END FUNCTION ReleaseActiveFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION SetCurrentTime(info,dummy)
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      SetCurrentTime = err_ok
!
      info%gridtime = currenttime
!
   END FUNCTION SetCurrentTime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION AfterStepStatistics(info,dummy)
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok, GetParentInfo
      USE PROBLEM, ONLY: AfterStep2D, AfterStep3D
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      TYPE(nodeinfo), POINTER:: parent
      INTEGER:: ierror, level
      INTEGER, DIMENSION(1:maxdims):: mx
      INTEGER, DIMENSION(1:maxdims,1:2):: mb
      REAL(KIND=r8):: h
      REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
      AfterStepStatistics = err_ok
      IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
      mx(1:ndims) = info%mx(1:ndims)
      h = info%dx(1)
      xlower(1:ndims) = info%xlower(1:ndims)
      level = info%level
      mb = 1; mb(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
!
      ierror = GetParentInfo(parent)
!
      SELECT CASE(ndims)
       CASE(2)
         CALL AfterStep2D(info%q(0:mx(1)+1,0:mx(2)+1,1,1:nccv), &
            parent%q(mb(1,1)-1:mb(1,2)+1, &
            mb(2,1)-1:mb(2,2)+1,1,1:nccv), &
            mx(1:2),h,xlower(1:2),level)
       CASE(3)
         CALL AfterStep3D(info%q(0:mx(1)+1,0:mx(2)+1,0:mx(3)+1,1:nccv), &
            parent%q(mb(1,1)-1:mb(1,2)+1, &
            mb(2,1)-1:mb(2,2)+1, &
            mb(3,1)-1:mb(3,2)+1,1:nccv), &
            mx(1:3),h,xlower(1:3),level)
       CASE DEFAULT
         PRINT *, 'BSAM 2.0: AfterStepStatistics: Only ndims = 2,3 are supported.'
         STOP
      END SELECT
!
   END FUNCTION AfterStepStatistics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   RECURSIVE SUBROUTINE AMR(level,estimateswitch)
! 功能：递归执行自适应网格细化/修复。
! 输入：
! - level            ：当前递归处理的层级。
! - estimateswitch?  ：可选开关，存在时跳过误差估计（用于回退重建）。
! 流程：
! - 调 SetGhost 同步边界；若需要，检测并修复“非相容网格”（hanging node）
!   缺陷（FindMeshDefects）。
! - 若修复触发，回退删除最近层级，并从更粗层重启 AMR 构建。
! - 若需误差估计：EstimateLevelErrors（计算截断误差/用户误差并标记）。
! - GridAdapt：依据标记生成子网格并转移场数据；更新 finestlevel；
! - 若还有更细层需处理，递归调用 AMR(level+1)。
      USE NodeInfoDef
      USE TreeOps, ONLY: ApplyOnLevel, DeleteMarkedNode, SetLevelNodeNumbers
      USE Boundary, ONLY: SetGhost
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: level
      INTEGER, OPTIONAL, INTENT(IN):: estimateswitch
      ! OPTIONAL means that the caller may choose to not provide this argument
!
      INTEGER:: lvl, noestimates
      TYPE(funcparam):: dummy
!
!PRINT *, 'AMR level', level
!
! Hard-coded for now, this will be a user option in the near future.
! to make the mesh having at most one hanging node.
      makeconformingmesh = .TRUE.
!
      IF(.NOT. makeconformingmesh) meshbuildcomplete = .TRUE.
!
! Detect coarse-level neighbors of the fine-level patches:
      IF(level >= 1 .AND. makeconformingmesh) &
         CALL ApplyOnLevel(level,FindCoarseLevelNeighbors,dummy)
!
! Fill in the ghost points:
      CALL SetGhost(level,0)
!
! Before we adapt the mesh further, check to see that the mesh is conforming:
       IF(level >= 2 .AND. makeconformingmesh) CALL FindMeshDefects(level)
       ! 这里：
         ! level: 当前处理层级
         ! makeconformingmesh: 是否强制一致网格（目前硬编码为真，未来可做成用户选项）
         ! FindMeshDefects: 检测并标记本层的“非相容网格缺陷”（defectivegridlevel(level)=.TRUE.）
!
! If the mesh is defective, then we need to fix it before moving on:
      IF(defectivegridlevel(level) .AND. amrrestarts < 100) THEN
!
! Delete the current and last active levels. Note, need to keep inactive grids:
         DO lvl = level, level-1, -1
            CALL ApplyOnLevel(lvl,ReleaseActiveFields,dummy)
! 设定：当前强制构建“至多一个悬挂点”的一致网格（未来可做成用户选项）。
            CALL ApplyOnLevel(lvl,DeleteMarkedNode,dummy)
         END DO
! 若不强制一致网格，则直接标记“网格构建完成”。
         defectivegridlevel(level-1:level) = .FALSE.
!
!  PRINT *, 'Restarting AMR at level =', level, ',    amrrestarts =', amrrestarts
! 目的：为细层 patch 查找对应的粗层邻居关系（用于后续幽灵填充/跨层插值等）。
! 条件：level>=1 且启用一致网格修复。
         amrrestarts = amrrestarts+1
         finestlevel = level-2
         CALL AMR(level-2,noestimates)
         ! here AMR will call himself so many times untill the top level is reached
! 为当前层所有活动网格填充幽灵点；第二实参 0 为模式/标志（项目中常用 0 表示标准幽灵填充）。
         RETURN
!
      END IF
! 在进一步细化前检查网格一致性：当层数≥2 且启用一致网格，检测是否存在“非相容网格缺陷”。
!
! Tag cells for refinement:
      IF(level < maxlevel .AND. (.NOT. PRESENT(estimateswitch))) THEN
! 若本层被标记为“有缺陷”，且回退次数未超过 100，则回退并修复。
         CALL EstimateLevelErrors(level)
      END IF
!
! 释放当前层与上一层的“活动网格”字段，并删除已标记节点（注意：保留非活动网格以便后续转移/对齐）。
! Create new subgrids if necessary:
      CALL GridAdapt(level)
      CALL SetLevelNodeNumbers(level+1)
!
! 清除此两层的缺陷标记。
      IF(level < finestlevel) CALL AMR(level+1)
!
   END SUBROUTINE AMR
! 递增回退计数，回退两层重置最细层，并从更粗的 level-2 重新进入 AMR；
! 递归调用时传入 noestimates（仅为“存在性”），跳过误差估计以加速恢复。
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE EstimateLevelErrors(level)
! 功能：在指定层做误差估计与标记（构建 tagged cell 链表并进行缓冲/扩张）。
! 返回上一层调用点（本分支修复流程结束）。
! 步骤：
! 1) 若采用默认误差策略：限制到粗层、补齐粗层幽灵、复制 q->qold，
!    计算相对截断误差（RelativeTruncationError）。
! 2) 遍历网格调用 EstimateError（ErrFlag），生成标记链表。
! 3) 对标记进行边缘膨胀（InflateEdgeTags），处理周期偏移；最后销毁链表。
! 若未达最大层，且本次未请求“跳过估计”，则进行误差估计并标记需要细化的单元。
      USE NodeInfoDef
      USE TreeOps, ONLY: ApplyOnLevel
      USE AFASRoutines, ONLY: RestrictCCSolution, RelativeTruncationError
      USE Boundary, ONLY: SetGhost, GetCoarseGhostPoints
      IMPLICIT NONE
! 根据误差标记与缓冲结果生成新子网格、转移/插值数据，并更新 finestlevel。
!
! 新层生成后，为 level+1 层重新分配并设置节点编号，保持层内一致性。
      INTEGER, INTENT(IN):: level
!
! 若当前层尚不是最细层，则递归处理更细一层。
      TYPE(funcparam):: dummy
!
! 1) Calculate the relative truncation error if needed:
      IF(errflagopt(level) == errflagdefault) THEN
         CALL ApplyOnLevel(level,RestrictCCSolution,dummy)
!
         CALL SetGhost(level-1,0)
!
         CALL ApplyOnLevel(level,GetCoarseGhostPoints,dummy)
!
         CALL ApplyOnLevel(level,CopyQToQold,dummy)
!
         CALL ApplyOnLevel(level-1,CopyQToQold,dummy)
!
         CALL ApplyOnLevel(level,RelativeTruncationError,dummy)
      END IF
!
! 2) Refine based on the size of the error estimator.  Make a linked list of
!    tagged cells:
      ALLOCATE(zerothtaggedcell)
      NULLIFY(zerothtaggedcell%prevcell)
!
      lasttaggedcell => zerothtaggedcell
      ntaggedcells = 0
!
      CALL ApplyOnLevel(level,EstimateError,dummy)
!
! 3) Inflate tagged regions and destroy the list:
      IF(ntaggedcells > 0) THEN
         CALL InflateEdgeTags(level)
         CALL DeleteTaggedCellsList
      END IF
!
      NULLIFY(lasttaggedcell)
      NULLIFY(currenttaggedcell)
      DEALLOCATE(zerothtaggedcell)
!
   END SUBROUTINE EstimateLevelErrors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION EstimateError(info,dummy)
! 功能：对单个网格执行误差标记流程 ErrFlag（默认或用户策略）。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      EstimateError = err_ok
      IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
      CALL ErrFlag(info)
!
   END FUNCTION EstimateError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE ErrFlag(info)
! 功能：根据 errflagopt(level) 选择默认或用户自定义的误差标记算法，
!       写入 info%errorflags。若配置了 ibuffer(level)>0，则触发 BufferAndList
!       将“跨 patch 的缓冲层”映射到全局坐标并加入链表。
! 关键变量：
! - qrte ：相对截断误差在粗层（coarse cell）的度量（维度依 ndims）。
! - errorflags：与细层单元对齐的二值标记。
! - mglobal：当前网格在全局索引系的范围（用于跨格缓冲）。
      USE NodeInfoDef
      USE Problem, ONLY: SetErrFlagsUser2D, SetErrFlagsUser3D
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
!
      INTEGER:: level
      INTEGER, DIMENSION(1:maxdims):: mx, cmx
      INTEGER, DIMENSION(1:maxdims,1:2):: mglobal
      REAL(KIND=r8):: h
      REAL(KIND=r8), DIMENSION(1:maxdims):: xlower
!
      level = info%level
      h = info%dx(1)
      xlower(1:ndims) = info%xlower(1:ndims)
      mx = 1; mx(1:ndims) = info%mx(1:ndims)
      cmx = 1; cmx(1:ndims) = mx(1:ndims)/2
      mglobal = 1; mglobal(1:ndims,1:2) = info%mglobal(1:ndims,1:2)
!
      SELECT CASE(errflagopt(level))
       CASE(errflagdefault)
!
! Default error tagging based on the relative truncation error:
         SELECT CASE(ndims)
          CASE(2)
            CALL SetErrFlags2D(info%qrte(1:cmx(1),1:cmx(2),1,1:nccv), &
               info%errorflags(1:mx(1),1:mx(2),1), &
               mx(1:2),cmx(1:2),h,level)
          CASE(3)
            CALL SetErrFlags3D(info%qrte(1:cmx(1),1:cmx(2),1:cmx(3),1:nccv), &
               info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
               mx(1:3),cmx(1:3),h,level)
!
          CASE DEFAULT
            PRINT *, 'ErrFlag: only ndims=2d,3d supported'
            STOP
         END SELECT
!
       CASE(errflaguser)
!
! User chooses the error tagging proceedure:
         SELECT CASE(ndims)
          CASE(2)
            CALL SetErrFlagsUser2D(info%qrte(1:cmx(1)  ,1:cmx(2)  ,1,1:nccv), &
               info%q   (0: mx(1)+1,0: mx(2)+1,1,1:nccv), &
               info%errorflags(1:mx(1),1:mx(2),1), &
               mx(1:2),cmx(1:2),h,xlower(1:2),level)
          CASE(3)
            CALL SetErrFlagsUser3D(info%qrte(1:cmx(1)  ,1:cmx(2)  ,1:cmx(3)  ,1:nccv), &
               info%q   (0: mx(1)+1,0: mx(2)+1,0: mx(3)+1,1:nccv), &
               info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
               mx(1:3),cmx(1:3),h,xlower(1:3),level)
!
          CASE DEFAULT
            PRINT *, 'ErrFlag: only ndims=2d,3d supported'
            STOP
         END SELECT
!
       CASE DEFAULT
         PRINT *, 'ErrFlag: No error flagging algorithm selected'
         STOP
      END SELECT
!
! Add the tagged cells to a linked list. Buffering requires a global approach,
! since buffer layers might go into neighboring grids:
      IF(ibuffer(level) > 0) THEN
         CALL BufferAndList(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
            mglobal(1:3,1:2),mx(1:3),level)
      END IF
!
   END SUBROUTINE ErrFlag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE SetErrFlags2D(qrte,errorflags,mx,cmx,h,level)
! 功能：2D 情况下，将粗单元 qrte 超阈值的 2x2 细单元块标记为 1。
! 阈值：tol = qtolerance(level) / h^2。
      USE NodeInfoDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN):: qrte
      INTEGER, DIMENSION(1:,1:), INTENT(OUT):: errorflags
      INTEGER, DIMENSION(1:2), INTENT(IN):: mx
      INTEGER, DIMENSION(1:2), INTENT(IN):: cmx
      REAL(KIND=r8), INTENT(IN):: h
      INTEGER, INTENT(IN):: level
!
      INTEGER:: i, j
      REAL(KIND=r8):: tol
!
      tol = qtolerance(level)/h/h
!
      errorflags(1:mx(1),1:mx(2)) = 0
!
      DO i = 1, cmx(1)
         DO j = 1, cmx(2)
            IF(MAXVAL(ABS(qrte(i,j,1:nccv)))>tol) THEN
               errorflags(2*i  ,2*j  ) = 1
               errorflags(2*i-1,2*j  ) = 1
               errorflags(2*i  ,2*j-1) = 1
               errorflags(2*i-1,2*j-1) = 1
            END IF
         END DO
      END DO
!
   END SUBROUTINE SetErrFlags2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE SetErrFlags3D(qrte,errorflags,mx,cmx,h,level)
! 功能：3D 情况下，将粗单元 qrte 超阈值的 2x2x2 细单元块标记为 1。
! 阈值：tol = qtolerance(level) / h^3。
      USE NodeInfoDef
      IMPLICIT NONE
!
      REAL(KIND=r8), DIMENSION(1:,1:,1:,1:), INTENT(IN):: qrte
      INTEGER, DIMENSION(1:,1:,1:), INTENT(OUT):: errorflags
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      INTEGER, DIMENSION(1:3), INTENT(IN):: cmx
      REAL(KIND=r8), INTENT(IN):: h
      INTEGER, INTENT(IN):: level
!
      INTEGER:: i, j, k
      REAL(KIND=r8):: tol
!
      tol = qtolerance(level)/h/h/h
!
      errorflags(1:mx(1),1:mx(2),1:mx(3)) = 0
!
      DO i = 1, cmx(1)
         DO j = 1, cmx(2)
            DO k = 1, cmx(3)
               IF(MAXVAL(ABS(qrte(i,j,k,1:nccv)))>tol) THEN
                  errorflags(2*i  ,2*j  ,2*k  ) = 1
                  errorflags(2*i-1,2*j  ,2*k  ) = 1
                  errorflags(2*i  ,2*j-1,2*k  ) = 1
                  errorflags(2*i-1,2*j-1,2*k  ) = 1
                  errorflags(2*i  ,2*j  ,2*k-1) = 1
                  errorflags(2*i-1,2*j  ,2*k-1) = 1
                  errorflags(2*i  ,2*j-1,2*k-1) = 1
                  errorflags(2*i-1,2*j-1,2*k-1) = 1
               END IF
            END DO
         END DO
      END DO
!
   END SUBROUTINE SetErrFlags3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE BufferAndList(errorflags,mglobal,mx,level)
! 功能：在当前 patch 内对已标记单元进行 ibuffer(level) 层的局部扩张。
!       若缓冲区跨越 patch 边界，则将相应的“中心单元”全局坐标加入链表，
!       供后续在其它 patch 上进行缓冲处理（跨格一致）。
      USE NodeInfoDef
      IMPLICIT NONE
!
      INTEGER, DIMENSION(1:,1:,1:), INTENT(IN OUT):: errorflags
      INTEGER, DIMENSION(1:3,1:2), INTENT(IN):: mglobal
      INTEGER, DIMENSION(1:3), INTENT(IN):: mx
      INTEGER, INTENT(IN):: level
!
      INTEGER:: i, j, k
      INTEGER, DIMENSION(1:maxdims):: index
      INTEGER, DIMENSION(1:maxdims,1:2):: mtg
      INTEGER, DIMENSION(1:mx(1),1:mx(2),1:mx(3)):: errorflagstmp
!
      errorflagstmp = errorflags
      mtg = 1
!
      DO k = 1, mx(3)
         index(3) = k
         DO j = 1, mx(2)
            index(2) = j
            DO i = 1, mx(1)
               index(1) = i
               IF(errorflagstmp(i,j,k)==1) THEN
                  mtg(1:ndims,1) = index(1:ndims)-ibuffer(level)
                  mtg(1:ndims,2) = index(1:ndims)+ibuffer(level)
                  IF(ANY(mtg(1:3,1)<1) .OR. ANY(mtg(1:3,2)>mx(1:3))) THEN
!
! If the buffer area overlaps with any edge of the patch, then record the
! global coordinates of the cell:
                     Call AddTaggedCellToList(index(1:maxdims)+mglobal(1:maxdims,1)-1)
                  ELSE
!
! If the buffer area lies within the patch, apply the buffer:
                     errorflags(mtg(1,1):mtg(1,2),mtg(2,1):mtg(2,2),mtg(3,1):mtg(3,2)) = 1
                  END IF
               END IF
            END DO
         END DO
      END DO
!
   END SUBROUTINE BufferAndList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE AddTaggedCellToList(globalindex)
! 功能：将一个全局单元坐标加入“标记单元”的链表（栈式，尾插为新头）。
! 参数：globalindex(1:maxdims) 为全局索引系中的 cell 坐标。
      USE NodeInfoDef
      IMPLICIT NONE
!
      INTEGER, DIMENSION(1:maxdims), INTENT(IN):: globalindex
!
      ntaggedcells = ntaggedcells+1
      ALLOCATE(currenttaggedcell)
      currenttaggedcell%id = ntaggedcells
      currenttaggedcell%coordinate(1:maxdims) = 1
      currenttaggedcell%coordinate(1:ndims) = globalindex(1:ndims)
      currenttaggedcell%prevcell => lasttaggedcell
      lasttaggedcell => currenttaggedcell
!
   END SUBROUTINE AddTaggedCellToList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE InflateEdgeTags(level)
! 功能：对链表中的全局标记单元执行跨 patch 的缓冲膨胀。
! 步骤：
! - 在当前 level 上，通过 ApplyOnLevel(BufferTaggedCells) 应用 ibuffer 层膨胀；
! - 若开启周期边界，则对每个周期偏移（正负方向）重复该缓冲；
! - 逐个回溯链表直至头结点。
      USE NodeInfoDef
      USE TreeOps, ONLY: ApplyOnLevel
      USE Boundary, ONLY: GetPeriodicOffsets
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: level
!
      TYPE(funcparam):: dummy
      LOGICAL:: periodicbuffer
      INTEGER:: offset, polarity
      INTEGER, DIMENSION(1:maxdims):: coordinatesave
!
      coordinatesave = 1
      currenttaggedcell => lasttaggedcell
      searchloop: DO
         IF(.NOT. ASSOCIATED(currenttaggedcell%prevcell)) EXIT searchloop
!
! Ordinary buffering of an edge tag:
         dummy%iswitch = ibuffer(level)
         CALL ApplyOnLevel(level,BufferTaggedCells,dummy)
!
! Buffering of periodic edge tags. Check to see if the buffer area cuts across
! a periodic boundary.  If so, add offset and apply buffer:
         coordinatesave(1:ndims) = currenttaggedcell%coordinate(1:ndims)
         IF(periodicboundaryconditions) THEN
            CALL GetPeriodicOffsets(level)
            DO polarity = -1, 1, 2
               DO offset = 1, nperiodicoffsets
                  currenttaggedcell%coordinate(1:ndims) &
                     = coordinatesave(1:ndims)+polarity*poffset(1:ndims,offset)
                  dummy%iswitch = ibuffer(level)
                  CALL ApplyOnLevel(level,BufferTaggedCells,dummy)
               END DO
            END DO
         END IF
         currenttaggedcell%coordinate(1:ndims) = coordinatesave(1:ndims)
!
         currenttaggedcell => currenttaggedcell%prevcell
      END DO searchloop
!
   END SUBROUTINE InflateEdgeTags
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION BufferTaggedCells(info,dummy)
! 功能：将 currenttaggedcell 指定的全局“缓冲区”映射到本地 patch 索引范围，
!       并把落在本 patch 内的单元在 info%errorflags 中置 1。
! 说明：dummy%iswitch = ibuff 表示缓冲半径（单位：细网格单元）。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      INTEGER:: ibuff, n
      INTEGER, DIMENSION(1:maxdims,1:2):: mglobal, mglobaltag, mlocal, moverlap
!
      BufferTaggedCells = err_ok
      IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
! ibuff stands for the buffer layer of the tagged cell, the default value is 1
      ibuff = dummy%iswitch
! info%mglobal(1:ndims,1:2) is the global index of the current patch
      mglobal = 1; mglobal(1:ndims,1:2) = info%mglobal(1:ndims,1:2)
      mglobaltag = 1
! mglobaltag is the index of the buffer layer
      mglobaltag(1:ndims,1) = currenttaggedcell%coordinate(1:ndims)-ibuff
      mglobaltag(1:ndims,2) = currenttaggedcell%coordinate(1:ndims)+ibuff
!
! 1. Find overlap region in global index space:
      moverlap = 1
      moverlap(1:ndims,1) = MAX(mglobaltag(1:ndims,1),mglobal(1:ndims,1))
      moverlap(1:ndims,2) = MIN(mglobaltag(1:ndims,2),mglobal(1:ndims,2))
!
! 2. Check for nonempty intersection:
      IF(ANY(moverlap(:,2)-moverlap(:,1)<0)) RETURN
!
! 3. Transform gobal index space to local grid index spaces:
      mlocal = 1
      DO n = 1, ndims
         mlocal(n,1:2) = moverlap(n,1:2)-mglobal(n,1)+1
      END DO
!
      info%errorflags(mlocal(1,1):mlocal(1,2), &
         mlocal(2,1):mlocal(2,2), &
         mlocal(3,1):mlocal(3,2)) = 1
!
   END FUNCTION BufferTaggedCells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE DeleteTaggedCellsList
! 功能：销毁“标记单元”链表（自尾到头逐个释放）。
      USE NodeInfoDef
      IMPLICIT NONE
!
      searchloop: DO
         currenttaggedcell => lasttaggedcell%prevcell
         DEALLOCATE(lasttaggedcell)
         lasttaggedcell => currenttaggedcell
         IF(.NOT. ASSOCIATED(lasttaggedcell%prevcell)) EXIT searchloop
      END DO searchloop
!
   END SUBROUTINE DeleteTaggedCellsList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE GridAdapt(level)
! 功能：依据误差标记生成 level+1 层的新子网格，并在新旧网格间转移场数据。
! 步骤：
! - RefineGrid：在每个活动网格上根据 errorflags 计算子网格范围并创建；
! - ApplyOnLevelPairs：在同层新旧网格对之间转移重叠区域的场值。
      USE NodeInfoDef
      USE TreeOps, ONLY: ApplyOnLevel, ApplyOnLevelPairs, DeleteMarkedNode
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: level
!
! Generate new subgrids of level:
!
      TYPE(funcparam):: dummy
      INTEGER:: ilevel
!
! Generate new grids in accordance with error flags:
      CALL ApplyOnLevel(level,RefineGrid,dummy)
!
! Transfer field values from previous grids on this level to the newly created
! grids:
      CALL ApplyOnLevelPairs(level+1,TransferValues,dummy)
!
   END SUBROUTINE GridAdapt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION RefineGrid(info,dummy)
! 功能：对单个网格，根据 errorflags 调用 NewSubGrids 生成子网格。
! 变量：
! - mbounds：初始为整个网格（1..mx），随后在 NewSubGrids 中分裂裁剪。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      INTEGER, DIMENSION(1:maxdims,1:2):: mbounds
!
      RefineGrid = err_ok
      IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
! Generate new subgrids of the current grid:
      mbounds(1:ndims,1) = 1; mbounds(1:ndims,2) = info%mx(1:ndims)
      CALL NewSubGrids(info,mbounds)
!
   END FUNCTION RefineGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE NewSubGrids(info,mbounds)
! 功能：改造的 Berger–Rigoutsos 算法实现。
! 过程：
! - 计算签名（沿各坐标方向对误差标记求和），裁剪外层零带；
! - 若填充率不足，优先在“空洞”处分裂，否则在“拐点”处分裂；
! - 保证偶数长度、最小长度 mgp 等约束；直到不可再分或达到上限；
! - 对得到的每个子块调用 MakeNewGrid 生成子网格。
      USE NodeInfoDef
      IMPLICIT NONE
!
! Modified implementation of Berger-Rigoutsos algorithm (IEEE Trans. Systems,
! Man. & Cyber., 21(5):1278-1286, 1991):
!
      TYPE(nodeinfo):: info
      INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN OUT):: mbounds
!
      LOGICAL havesplit, cansplitgrid
      LOGICAL, DIMENSION(1:maxsubgrids):: cansplit
      INTEGER, PARAMETER:: maxsplitpasses = 15
      INTEGER:: del, dist0, dist1, i, i1, i2, ierror, igrid, inflect, level, maxm, &
         mgp, minm, n, nn, ngrid, npass
      INTEGER, DIMENSION(1:maxdims):: isplit, mx
      INTEGER, DIMENSION(1:maxdims,1:2,1:maxsubgrids):: msubbounds
      INTEGER, DIMENSION(:,:), ALLOCATABLE:: signature ,ddsignature
      REAL(KIND=r8):: fillratio, desfillratio
!
      mx = info%mx
      level = info%level; info%nsubgrids = 0
      desfillratio = desiredfillratios(level)
!
      mgp = minimumgridpoints(level)
      IF(mgp<4) THEN
         PRINT *,'BSAM 2.0: Error on level', level, 'minimumgridpoints cannot be less than 4.'
         STOP
      END IF
      IF(MODULO(mgp,2)==1) THEN
         PRINT *,'BSAM 2.0: Error on level', level, 'minimumgridpoints must be even.'
         STOP
      END IF
!
! Compute fill ratio for this grid:
      fillratio = GridFlagRatio(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
         mbounds)
!
! Don't generate a new grid if we don't have flagged points:
      IF(fillratio < 1.0E-10_r8) RETURN
!
! Allocate space for signatures:
      maxm = MAXVAL(mx(1:ndims))
      ALLOCATE(signature(1:maxm,1:ndims),ddsignature(1:maxm,1:ndims),STAT=ierror)
      IF(ierror /= 0) THEN
         PRINT *,'BSAM 2.0: Error allocating signatures arrays in NewSubGrids'
         STOP
      END IF
      signature=0; ddsignature=0
!
! Initialize list of subgrids:
      ngrid = 1; cansplit(:) = .TRUE.
      msubbounds(1:ndims,1:2,ngrid) = mbounds(1:ndims,1:2)
!
! Loop until no better grid splitting can be found:
      igrid = 1
      DO WHILE (ngrid<maxsubgrids .AND. igrid<=ngrid)
         npass = 0
         DO WHILE (cansplit(igrid) .AND. npass<maxsplitpasses)
            npass = npass+1
            signature = GetSignatures(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
               msubbounds(:,:,igrid),maxm)
!
! Trim unflagged points on the edges of this grid:
            DO n = 1, ndims
               i1 = msubbounds(n,1,igrid); i2=msubbounds(n,2,igrid)
!
               DO WHILE(signature(i1,n)==0 .AND. i1<msubbounds(n,2,igrid) .AND. &
                  i2-i1+1>mgp)
                  i1 = i1+1
               END DO
               DO WHILE(signature(i2,n)==0 .AND. i2>msubbounds(n,1,igrid) .AND. &
                  i2-i1+1>mgp)
                  i2 = i2-1
               END DO
!
! We only make grid adjustments in increments of 2 (coarse) grid points:
               msubbounds(n,1,igrid) = i1-MODULO(i1+1,2)
               msubbounds(n,2,igrid) = i2+MODULO(i2  ,2)
            END DO
!
            fillratio = GridFlagRatio(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
               msubbounds(:,:,igrid))
            minm = MINVAL(msubbounds(1:ndims,2,igrid)-msubbounds(1:ndims,1,igrid))+1
            IF(fillratio<desfillratio) THEN
               signature = GetSignatures(info%errorflags(1:mx(1),1:mx(2),1:mx(3)), &
                  msubbounds(:,:,igrid),maxm)
!
! Look for holes along which to split grid:
               isplit = 0; havesplit = .FALSE.
               DO n = 1, ndims
                  i1 = msubbounds(n,1,igrid); i2 = msubbounds(n,2,igrid)
                  DO i = i1+mgp-1, i2-mgp, 2
!
! i is the second (terminal) index of the first half of the split grid.
                     IF(signature(i,n)==0 .AND. MIN(i-i1+1,i2-i)>=mgp) THEN
                        isplit(n) = i
                        havesplit = .TRUE.
                        EXIT
                     END IF
                  END DO
                  IF(havesplit) EXIT
               END DO
!
               IF(.NOT. havesplit) THEN
!
! No split along a hole. Try split along inflection point:
                  DO n = 1, ndims
                     i1 = msubbounds(n,1,igrid); i2 = msubbounds(n,2,igrid)
                     DO i = i1+1, i2-1
                        ddsignature(i,n) = signature(i-1,n)-2*signature(i,n) &
                           + signature(i+1,n)
                     END DO
                  END DO
!
                  inflect = 0; dist0 = 0
                  DO n = 1, ndims
                     i1 = msubbounds(n,1,igrid); i2 = msubbounds(n,2,igrid)
                     DO i = i1+mgp, i2-mgp+1, 2
!
! Here i is the first index of the second half of the split grid.
                        del = ABS(ddsignature(i,n)-ddsignature(i-1,n))
                        IF(del>inflect) THEN
                           inflect = del; isplit = 0; isplit(n) = i-1; havesplit = .TRUE.
                           dist0 = MIN(i-i1,i2-i+1)
                        ELSE IF(del==inflect .AND. inflect>0) THEN
                           dist1 = MIN(i-i1,i2-i+1)
                           IF(dist1>dist0) THEN
                              isplit = 0; isplit(n) = i-1; havesplit = .TRUE.
                              dist0 = dist1
                           END IF
                        END IF
                     END DO
                  END DO
               END IF
!
               IF(havesplit) THEN
!
! Split the grid along a determined line:
                  DO n = 1, ndims
                     IF(isplit(n)>0 .AND. &
                        MIN(msubbounds(n,2,igrid)-isplit(n), &
                        isplit(n)-msubbounds(n,1,igrid)+1)>=mgp) THEN
!
! Add a new subgrid to the end of the grid list:
                        ngrid = ngrid+1
                        cansplit(ngrid) = .TRUE.
                        msubbounds(1:ndims,1:2,ngrid) = msubbounds(1:ndims,1:2,igrid)
                        msubbounds(n,1,ngrid) = isplit(n)+1
!
! Replace current grid with a subgrid:
                        msubbounds(n,2,igrid) = isplit(n)
                        EXIT
                     END IF
                  END DO
               ELSE
!
! Mark grid if no split is possible:
                  cansplit(igrid) = .FALSE.
               END IF
            ELSE
               cansplit(igrid) = .FALSE.
            END IF
         END DO
         igrid = igrid+1
      END DO
!
! Generate the newly determined subgrids:
      info%nsubgrids = ngrid
      DO i = 1, ngrid
         CALL MakeNewGrid(info,msubbounds(:,:,i))
         IF(MINVAL(msubbounds(1:ndims,2,i)-msubbounds(1:ndims,1,i)+1)<mgp) THEN
            PRINT *, 'BSAM 2.0: Error in NewSubGrids, grid smaller than minimumgridpoints'
            STOP
         END IF
         DO n = 1, ndims
            IF(MODULO(msubbounds(n,2,i)-msubbounds(n,1,i)+1,2)==1) THEN
               PRINT *, 'BSAM 2.0: Error in NewSubGrids, igrid=', i, ', grid length not even number'
               STOP
            END IF
            IF(MODULO(msubbounds(n,1,i),2)==0) THEN
               PRINT *, 'BSAM 2.0: Error in NewSubGrids, igrid=', i, ', index starts on even number'
               STOP
            END IF
            IF(MODULO(msubbounds(n,2,i),2)==1) THEN
               PRINT *, 'BSAM 2.0: Error in NewSubGrids, igrid=', i, ', index end on odd number'
               STOP
            END IF
         END DO
      END DO
!
      DEALLOCATE(signature,ddsignature,STAT=ierror)
      IF(ierror/=0) THEN
         PRINT *,'BSAM 2.0: Error deallocating signatures arrays in NewSubGrids'
         STOP
      END IF
!
   END SUBROUTINE NewSubGrids
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION GetSignatures(errorflags,msubbounds,maxm) RESULT(gsresult)
! 功能：计算给定子区域在各方向上的“签名”（将另两维压缩求和后的 1D 序列）。
      USE NodeInfoDef
      IMPLICIT NONE
!
      INTEGER, DIMENSION(1:,1:,1:), INTENT(IN):: errorflags
      INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN):: msubbounds
      INTEGER, INTENT(IN):: maxm
      INTEGER, DIMENSION(1:maxm,1:ndims):: gsresult
!
      INTEGER:: i, n
      INTEGER, DIMENSION(1:maxdims):: i1, i2
!
      i1 = 1; i2 = 1
      DO n = 1, ndims
         i1(1:ndims) = msubbounds(1:ndims,1)
         i2(1:ndims) = msubbounds(1:ndims,2)
         DO i = msubbounds(n,1), msubbounds(n,2)
            i1(n) = i; i2(n) = i
!
            gsresult(i,n) = SUM(errorflags(i1(1):i2(1),i1(2):i2(2),i1(3):i2(3)))
!
         END DO
      END DO
!
   END FUNCTION GetSignatures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION GridFlagRatio(errorflags,mbounds) RESULT(gfrresult)
! 功能：统计子区域内被标记单元占比（用于决定是否继续分裂或接受当前子块）。
      USE NodeInfoDef
      IMPLICIT NONE
!
      INTEGER, DIMENSION(1:,1:,1:), INTENT(IN):: errorflags
      INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN OUT):: mbounds
      REAL(KIND=r8):: gfrresult
!
      REAL(KIND=r8):: flagged, total
!
      total = REAL(PRODUCT(mbounds(1:ndims,2)-mbounds(1:ndims,1)+1),KIND=r8)
!
      mbounds(ndims+1:maxdims,1:2) = 1
!
      flagged = REAL(SUM(errorflags(mbounds(1,1):mbounds(1,2), &
         mbounds(2,1):mbounds(2,2), &
         mbounds(3,1):mbounds(3,2))),KIND=r8)
!
      IF(flagged<-1.0E-08_r8) THEN
         PRINT *, 'BSAM 2.0: Error in GridFlagRatio: flagged < 0.'
         STOP
      END IF
!
      gfrresult = flagged/total
!
   END FUNCTION GridFlagRatio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE MakeNewGrid(parent,mbounds)
! 功能：在父网格的 mbounds 粗单元范围内创建一个新的细化子网格节点，
!       设置几何/边界/索引信息，分配存储，并通过 InitFields 从父网格插值初始化。
! 要点：
! - 子网格的 mglobal 依据父网格 mglobal 与 mbounds 映射到全局索引（×2）。
! - mthbc 先设为 internal，再继承靠近父物理边界的一侧物理边界条件。
! - field 分配后，levellandscape 初始化为本层级，便于缺陷检测。
      USE NodeInfoDef
      USE TreeOps, ONLY: CreateChild, GetChildInfo
      USE BSAMStorage, ONLY: AllocFields
      USE Problem, ONLY: SetAux, SetSrc
      IMPLICIT NONE
!
      TYPE(nodeinfo):: parent
      INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN):: mbounds
!
! Generate a new, finer grid within mbounds of the grid info:
!
      TYPE(nodeinfo), POINTER:: child
      INTEGER:: i, ierror, n, nb
      INTEGER, DIMENSION(1:maxdims):: cmx
      INTEGER, DIMENSION(1:maxdims,1:2):: mglobalbounds
      REAL(KIND=r8):: rand
!
! Create a child of the currentnode:
      CALL CreateChild
!
! Initialize the grid information for the new node:
! Get pointer to youngest child's info:
      ierror = GetChildInfo(child)
!
! Start filling in the fields of child's info:
      child%tobedeleted = .FALSE.
      child%activegrid = .TRUE.
      child%defective = .FALSE.
      child%fieldsallocated = .FALSE.
      child%initialgrid = parent%initialgrid
!
      child%maxlevel = parent%maxlevel; child%nsubgrids = 0
      child%level = parent%level+1
      IF(child%level>finestlevel) finestlevel = child%level
!
      child%mx = 1; child%mx(1:ndims) = (mbounds(1:ndims,2)-mbounds(1:ndims,1)+1)*2
      child%mbounds = 1; child%mbounds(1:ndims,:) = mbounds(1:ndims,:)
      child%mglobal = 1
      mglobalbounds(1:ndims,1) = parent%mglobal(1:ndims,1)+mbounds(1:ndims,1)-1
      mglobalbounds(1:ndims,2) = mglobalbounds(1:ndims,1)+mbounds(1:ndims,2) &
         - mbounds(1:ndims,1)
      child%mglobal(1:ndims,1) = (mglobalbounds(1:ndims,1)-1)*2+1
      child%mglobal(1:ndims,2) = mglobalbounds(1:ndims,2)*2
!
! First assume all boundaries are internal:
      child%mthbc = internalbc
!
! Now check if we have any physical boundaries:
      DO n = 1, ndims
         nb = 2*n-1
!
! If parent boundary condition is physical and child left is same as parent's
! left:
         IF((parent%mthbc(nb)<internalbc) .AND. &
            (child%mbounds(n,1)==1)) THEN
!
! Then child left boundary is physical and inherited from parent:
            child%mthbc(nb)=parent%mthbc(nb)
         END IF
         nb = nb+1
!
! If parent boundary condition is physical and child right is same as parent's
! right:
         IF((parent%mthbc(nb)<internalbc) .AND. &
            (child%mbounds(n,2)==parent%mx(n))) THEN
!
! Then child right boundary is physical and inherited from parent
            child%mthbc(nb) = parent%mthbc(nb)
         END IF
      END DO
!
! ID the grids randomly:
      CALL RANDOM_NUMBER(rand)
      child%ngrid = NINT(10000000*rand)
!
      child%gridtime = parent%gridtime
!
      child%xlower = 0.0_r8; child%xupper = 0.0_r8
      child%xlower(1:ndims) = REAL(mbounds(1:ndims,1)-1,KIND=r8)*parent%dx(1:ndims) &
         + parent%xlower(1:ndims)
      child%xupper(1:ndims) = REAL(mbounds(1:ndims,2)  ,KIND=r8)*parent%dx(1:ndims) &
         + parent%xlower(1:ndims)
      child%dx = 0.0_r8; child%dx(1:ndims) = parent%dx(1:ndims)/REAL(2,KIND=r8)
!
! Allocate dynamic space:
      CALL AllocFields(child,parent)
!
      cmx = 1
      cmx(1:ndims) = child%mx(1:ndims)/2
      child%levellandscape = parent%level
      child%levellandscape(1:cmx(1),1:cmx(2),1:cmx(3)) = child%level
!
! Initialize field variables with values from parent or from initial data:
      CALL InitFields(parent,child)
!
! Call user routine to initialize auxiliary and source-term array values:
      IF(naxv>0) CALL SetAux(child)
      CALL SetSrc(child)
!
! Delete later:
      child%rf = 0.0
      child%rf1 = 0.0
!
      child%qold = child%q
      child%qcold = child%qc
!
   END SUBROUTINE MakeNewGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE InitFields(parent,child)
! 功能：从父网格到子网格的场值初始化：
! - CC 变量：先将父层 CC 值拷贝到子层 coarse 存储（qc），再按 mbc 调用
!   双/三线性带质量修正的 ProlongationP1MC/P2MC 或其 3D 版本生成细网格 q。
! - FC 变量：调用 Prolongation2D/3D V1/V2(/V3) 的扩展版本在边/面上插值。
! - 若 child%initialgrid 为真，则直接调用 Initialize(child) 由用户初始化。
      USE NodeInfoDef
      USE GridUtilities, ONLY:       Prolongation2DV1Ex,       Prolongation2DV2Ex, &
         BiLinProlongationP1MC  ,  BiLinProlongationP2MC,   &
         TriLinProlongationP1MC  , TriLinProlongationP2MC,   &
         Prolongation3DV1Ex,       Prolongation3DV2Ex, &
         Prolongation3DV3Ex
      IMPLICIT NONE
!
! Redone to support only bilinear interpolation and 1st layer updating.
!
      TYPE(nodeinfo):: parent, child
!
      INTEGER:: ierror, interpopt
      INTEGER, DIMENSION(1:maxdims):: cmx, mx
      INTEGER, DIMENSION(1:maxdims,1:2):: mb
!
! Check whether we can apply user-provided Initialize for initialization.
      IF(child%initialgrid) THEN
         CALL Initialize(child)
         RETURN
      END IF
!
      mx = 1;  mx(1:ndims)     = child%mx(1:ndims)
      cmx = 1; cmx(1:ndims)     =       mx(1:ndims)/2
      mb = 1;  mb(1:ndims,1:2) = child%mbounds(1:ndims,1:2)
!
      SELECT CASE(ndims)
       CASE(2)
         child%qc(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1,1:nccv) &
            = parent%q(mb(1,1)-mbc:mb(1,2)+mbc, &
            mb(2,1)-mbc:mb(2,2)+mbc,1,1:nccv)
!
         child%qcold(1-mbc:cmx(1)+mbc,1-mbc:cmx(2)+mbc,1,1:nccv) &
            = parent%qold(mb(1,1)-mbc:mb(1,2)+mbc, &
            mb(2,1)-mbc:mb(2,2)+mbc,1,1:nccv)
!
         child%v1c(-1:cmx(1)+1,0:cmx(2)+1,1,1:nfcv) &
            = parent%v1(mb(1,1)-2:mb(1,2)+1, &
            mb(2,1)-1:mb(2,2)+1,1,1:nfcv)
!
         child%v2c(0:cmx(1)+1,-1:cmx(2)+1,1,1:nfcv) &
            = parent%v2(mb(1,1)-1:mb(1,2)+1, &
            mb(2,1)-2:mb(2,2)+1,1,1:nfcv)
!
         SELECT CASE(mbc)
          CASE(1)
            child%q(  0: mx(1)+1, 0: mx(2)+1,1,1:nccv) &
               = BiLinProlongationP1MC( &
               child%qc( 0:cmx(1)+1, 0:cmx(2)+1,1,1:nccv))
          CASE(2)
            child%q( -1: mx(1)+2,-1: mx(2)+2,1,1:nccv) &
               = BiLinProlongationP2MC( &
               child%qc(-1:cmx(1)+2,-1:cmx(2)+2,1,1:nccv))
          CASE DEFAULT
            PRINT *, 'InitFields: only supports mbc=1,2.'
         END SELECT
!
! Face centered variables:
         child%v1( -1: mx(1)+1, 0: mx(2)+1,1,1:nfcv) &
            = Prolongation2DV1Ex( &
            child%v1c(-1:cmx(1)+1, 0:cmx(2)+1,1,1:nfcv))
!
         child%v2(  0: mx(1)+1,-1: mx(2)+1,1,1:nfcv) &
            = Prolongation2DV2Ex( &
            child%v2c( 0:cmx(1)+1,-1:cmx(2)+1,1,1:nfcv))
!
       CASE(3)
         child%qc (      1-mbc:cmx(1)  +mbc,        &
            1-mbc:cmx(2)  +mbc,        &
            1-mbc:cmx(3)  +mbc,1:nccv) &
            = parent%q (mb(1,1)-mbc: mb(1,2)+mbc,        &
            mb(2,1)-mbc: mb(2,2)+mbc,        &
            mb(3,1)-mbc: mb(3,2)+mbc,1:nccv)
!
         child%qcold(      1-mbc:cmx(1)  +mbc,        &
            1-mbc:cmx(2)  +mbc,        &
            1-mbc:cmx(3)  +mbc,1:nccv) &
            = parent%qold(mb(1,1)-mbc: mb(1,2)+mbc,        &
            mb(2,1)-mbc: mb(2,2)+mbc,        &
            mb(3,1)-mbc: mb(3,2)+mbc,1:nccv)
!
         child%v1c(         -1:cmx(1)  +1  ,        &
            0:cmx(2)  +1  ,        &
            0:cmx(3)  +1  ,1:nfcv) &
            = parent%v1(mb(1,1)-2  : mb(1,2)+1  ,        &
            mb(2,1)-1  : mb(2,2)+1  ,        &
            mb(3,1)-1  : mb(3,2)+1  ,1:nfcv)
!
         child%v2c(          0:cmx(1)  +1  ,        &
            -1:cmx(2)  +1  ,        &
            0:cmx(3)  +1  ,1:nfcv) &
            = parent%v2(mb(1,1)-1  : mb(1,2)+1  ,        &
            mb(2,1)-2  : mb(2,2)+1  ,        &
            mb(3,1)-1  : mb(3,2)+1  ,1:nfcv)
!
         child%v3c(          0:cmx(1)  +1  ,        &
            0:cmx(2)  +1  ,        &
            -1:cmx(3)  +1  ,1:nfcv) &
            = parent%v3(mb(1,1)-1  : mb(1,2)+1  ,        &
            mb(2,1)-1  : mb(2,2)+1  ,        &
            mb(3,1)-2  : mb(3,2)+1  ,1:nfcv)
!
         SELECT CASE(mbc)
          CASE(1)
            child%q(  0: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1:nccv) &
               = TriLinProlongationP1MC( &
               child%qc( 0:cmx(1)+1, 0:cmx(2)+1, 0:cmx(3)+1,1:nccv))
          CASE(2)
            child%q( -1: mx(1)+2,-1: mx(2)+2,-1: mx(3)+2,1:nccv) &
               = TriLinProlongationP2MC( &
               child%qc(-1:cmx(1)+2,-1:cmx(2)+2,-1:cmx(3)+2,1:nccv))
          CASE DEFAULT
            PRINT *, 'InitFields: only supports mbc=1,2.'
         END SELECT
!
! Face centered variables:
!
         child%v1 ( -1: mx(1)+1, 0: mx(2)+1, 0: mx(3)+1,1:nfcv) &
            = Prolongation3DV1Ex( &
            child%v1c( -1:cmx(1)+1, 0:cmx(2)+1, 0:cmx(3)+1,1:nfcv))
!
         child%v2 (  0: mx(1)+1,-1: mx(2)+1, 0: mx(3)+1,1:nfcv) &
            = Prolongation3DV2Ex( &
            child%v2c(  0:cmx(1)+1,-1:cmx(2)+1, 0:cmx(3)+1,1:nfcv))
!
         child%v3 (  0: mx(1)+1, 0: mx(2)+1,-1: mx(3)+1,1:nfcv) &
            = Prolongation3DV3Ex( &
            child%v3c(  0:cmx(1)+1, 0:cmx(2)+1,-1:cmx(3)+1,1:nfcv))
!
       CASE DEFAULT
         PRINT *, 'InitFields: only supports ndims=2,3.'
      END SELECT
!
   END SUBROUTINE InitFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION TransferValues(grid1,grid2,dummy)
! 功能：在同层新旧网格对中，查找重叠并将“旧网格”的解转移到“新网格”。
! 规则：当 grid1 活动而 grid2 非活动时，从 grid2→grid1；反之亦然。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
      TYPE(nodeinfo):: grid1, grid2
      TYPE(funcparam):: dummy
!
! Transfer values from previous grids on this level to newly created grids:
!
      TransferValues = err_ok
!
      IF(grid1%activegrid .AND. (.NOT. grid2%activegrid)) THEN
!
! Look for overlap and transfer grid values from Grid2 to Grid1
         CALL Transferq(grid2,grid1)
      END IF
!
      IF(grid2%activegrid .AND. (.NOT. grid1%activegrid)) THEN
!
! Look for overlap and transfer grid values from Grid1 to Grid2
         CALL Transferq(grid1,grid2)
      END IF
!
! We either have two new grids or two old grids, nothing to be done:
      RETURN
   END FUNCTION TransferValues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE Transferq(sourceinfo,targetinfo)
! 功能：对两网格求全局索引的重叠区域，在重叠内将 sourceinfo%q 赋给 targetinfo%q。
! 注意：此处仅处理单元中心变量 q（可按需扩展到 v* 等）。
      USE NodeInfoDef
      IMPLICIT NONE
!
      TYPE(nodeinfo):: sourceinfo
      TYPE(nodeinfo):: targetinfo
!
! Look for overlap and transfer grid values from source to target
!
      INTEGER, DIMENSION(1:maxdims,1:2):: mtarget, msource
!
      IF(.NOT. sourceinfo%fieldsallocated) THEN
         PRINT *, 'BSAM 2.0 Error: Trying to transfer values from unallocated'
         PRINT *, '                sourceinfo in Transferq.'
         STOP
      END IF
!
      IF(.NOT. targetinfo%fieldsallocated) THEN
         PRINT *, 'BSAM 2.0 Error: Trying to transfer values to unallocated'
         PRINT *, '                targetinfo in Transferq.'
         STOP
      END IF
!
! Look for overlap (Later, maybe copy over ghost cells as well):
      msource = sourceinfo%mglobal; mtarget = targetinfo%mglobal
!
      CALL TransferOverlap(msource,mtarget,sourceinfo,targetinfo)
!
   END SUBROUTINE Transferq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE TransferOverlap(msource,mtarget,sourceinfo,targetinfo)
! 功能：给定两网格的全局索引范围 msource/mtarget，求其交集 moverlap，
!       再映射到各自局部索引 ms/mt，并复制相应子数组。
      USE NodeInfoDef
      IMPLICIT NONE
!
      INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN):: msource
      INTEGER, DIMENSION(1:maxdims,1:2), INTENT(IN):: mtarget
      TYPE(nodeinfo):: sourceinfo
      TYPE(nodeinfo):: targetinfo
!
      INTEGER:: n
      INTEGER, DIMENSION(1:maxdims,1:2):: moverlap, ms, mt
!
! Transfer values from source to target in overlap region:
! 1. Find overlap region in global index space:
      moverlap = 1
      moverlap(1:ndims,1) = MAX(msource(1:ndims,1),mtarget(1:ndims,1))
      moverlap(1:ndims,2) = MIN(msource(1:ndims,2),mtarget(1:ndims,2))
!
! 2. Check for nonempty intersection:
      IF(ANY(moverlap(:,2)-moverlap(:,1)<0)) RETURN
!
! 3. Transform common index space to grid index spaces:
      ms = 1; mt = 1
      DO n = 1, ndims
         ms(n,:) = moverlap(n,:)-sourceinfo%mglobal(n,1)+1
         mt(n,:) = moverlap(n,:)-targetinfo%mglobal(n,1)+1
      END DO
!
! 4. Carry out the transfer:
      targetinfo%q(mt(1,1):mt(1,2),mt(2,1):mt(2,2),mt(3,1):mt(3,2),:) &
         = sourceinfo%q(ms(1,1):ms(1,2),ms(2,1):ms(2,2),ms(3,1):ms(3,2),:)
!
   END SUBROUTINE TransferOverlap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION FindCoarseLevelNeighbors(info,dummy)
! 功能：在粗层（父层）上查询与当前细网格相邻的“对面块”的层级数，
!       将其写入当前网格的 levellandscape（用于后续缺陷检测）。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok, GetParentInfo
      IMPLICIT NONE
!
! This routine determines the level values of neighbor grids touching the
! current grid on the coarse level. Fine level neighbors are determined in
! MODULE Boundary:
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      TYPE(nodeinfo), POINTER:: parent
      INTEGER:: ierror
      INTEGER, DIMENSION(1:maxdims):: cmx
      INTEGER, DIMENSION(1:maxdims,1:2):: mboc, mbounds
!
      cmx = 1
      cmx(1:ndims) = info%mx(1:ndims)/2
      mbounds = 1
      mbounds(1:ndims,1:2) = info%mbounds(1:ndims,1:2)
      mboc = 1
      mboc(1:ndims,1) = (mbounds(1:ndims,1)+1)/2
      mboc(1:ndims,2) =  mbounds(1:ndims,2)   /2
!
      FindCoarseLevelNeighbors = err_ok
      IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
      ierror = GetParentInfo(parent)
!
      SELECT CASE(ndims)
       CASE(2)
!
! Right (x = cmx+1):
         info%levellandscape( cmx(1)  +1,          1: cmx(2)+1:2,1) &
            = parent%levellandscape(mboc(1,2)+1,mboc(2,1)  :mboc(2,2)+1,1)
         info%levellandscape( cmx(1)  +1,          0: cmx(2)  :2,1) &
            = parent%levellandscape(mboc(1,2)+1,mboc(2,1)-1:mboc(2,2)  ,1)
!
! Left  (x = 0):
         info%levellandscape(          0,          1: cmx(2)+1:2,1) &
            = parent%levellandscape(mboc(1,1)-1,mboc(2,1)  :mboc(2,2)+1,1)
         info%levellandscape(          0,          0: cmx(2)  :2,1) &
            = parent%levellandscape(mboc(1,1)-1,mboc(2,1)-1:mboc(2,2)  ,1)
!
! Top   (y = cmx+1):
         info%levellandscape(          1: cmx(1)+1:2, cmx(2)  +1,1) &
            = parent%levellandscape(mboc(1,1)  :mboc(1,2)+1,mboc(2,2)+1,1)
         info%levellandscape(          0: cmx(1)  :2, cmx(2)  +1,1) &
            = parent%levellandscape(mboc(1,1)-1:mboc(1,2)  ,mboc(2,2)+1,1)
!
! Bottom(y = 0):
         info%levellandscape(          1: cmx(1)+1:2,          0,1) &
            = parent%levellandscape(mboc(1,1)  :mboc(1,2)+1,mboc(2,1)-1,1)
         info%levellandscape(          0: cmx(1)  :2,          0,1) &
            = parent%levellandscape(mboc(1,1)-1:mboc(1,2)  ,mboc(2,1)-1,1)
!
       CASE(3)
!
! Right (x = cmx+1):
         info%levellandscape( cmx(1)  +1,          1: cmx(2)+1:2,          1: cmx(3)+1:2) &
            = parent%levellandscape(mboc(1,2)+1,mboc(2,1)  :mboc(2,2)+1,mboc(3,1)  :mboc(3,2)+1)
         info%levellandscape( cmx(1)  +1,          0: cmx(2)  :2,          1: cmx(3)+1:2) &
            = parent%levellandscape(mboc(1,2)+1,mboc(2,1)-1:mboc(2,2)  ,mboc(3,1)  :mboc(3,2)+1)
         info%levellandscape( cmx(1)  +1,          1: cmx(2)+1:2,          0: cmx(3)  :2) &
            = parent%levellandscape(mboc(1,2)+1,mboc(2,1)  :mboc(2,2)+1,mboc(3,1)-1:mboc(3,2)  )
         info%levellandscape( cmx(1)  +1,          0: cmx(2)  :2,          0: cmx(3)  :2) &
            = parent%levellandscape(mboc(1,2)+1,mboc(2,1)-1:mboc(2,2)  ,mboc(3,1)-1:mboc(3,2)  )
!
! Left  (x = 0):
         info%levellandscape(          0,          1: cmx(2)+1:2,          1: cmx(3)+1:2) &
            = parent%levellandscape(mboc(1,1)-1,mboc(2,1)  :mboc(2,2)+1,mboc(3,1)  :mboc(3,2)+1)
         info%levellandscape(          0,          0: cmx(2)  :2,          1: cmx(3)+1:2) &
            = parent%levellandscape(mboc(1,1)-1,mboc(2,1)-1:mboc(2,2)  ,mboc(3,1)  :mboc(3,2)+1)
         info%levellandscape(          0,          1: cmx(2)+1:2,          0: cmx(3)  :2) &
            = parent%levellandscape(mboc(1,1)-1,mboc(2,1)  :mboc(2,2)+1,mboc(3,1)-1:mboc(3,2)  )
         info%levellandscape(          0,          0: cmx(2)  :2,          0: cmx(3)  :2) &
            = parent%levellandscape(mboc(1,1)-1,mboc(2,1)-1:mboc(2,2)  ,mboc(3,1)-1:mboc(3,2)  )
!
! Top   (y = cmx+1):
         info%levellandscape(          1: cmx(1)+1:2, cmx(2)  +1,          1: cmx(3)+1:2) &
            = parent%levellandscape(mboc(1,1)  :mboc(1,2)+1,mboc(2,2)+1,mboc(3,1)  :mboc(3,2)+1)
         info%levellandscape(          0: cmx(1)  :2, cmx(2)  +1,          1: cmx(3)+1:2) &
            = parent%levellandscape(mboc(1,1)-1:mboc(1,2)  ,mboc(2,2)+1,mboc(3,1)  :mboc(3,2)+1)
         info%levellandscape(          1: cmx(1)+1:2, cmx(2)  +1,          0: cmx(3)  :2) &
            = parent%levellandscape(mboc(1,1)  :mboc(1,2)+1,mboc(2,2)+1,mboc(3,1)-1:mboc(3,2)  )
         info%levellandscape(          0: cmx(1)  :2, cmx(2)  +1,          0: cmx(3)  :2) &
            = parent%levellandscape(mboc(1,1)-1:mboc(1,2)  ,mboc(2,2)+1,mboc(3,1)-1:mboc(3,2)  )
!
! Bottom(y = 0):
         info%levellandscape(          1: cmx(1)+1:2,          0,          1: cmx(3)+1:2) &
            = parent%levellandscape(mboc(1,1)  :mboc(1,2)+1,mboc(2,1)-1,mboc(3,1)  :mboc(3,2)+1)
         info%levellandscape(          0: cmx(1)  :2,          0,          1: cmx(3)+1:2) &
            = parent%levellandscape(mboc(1,1)-1:mboc(1,2)  ,mboc(2,1)-1,mboc(3,1)  :mboc(3,2)+1)
         info%levellandscape(          1: cmx(1)+1:2,          0,          0: cmx(3)  :2) &
            = parent%levellandscape(mboc(1,1)  :mboc(1,2)+1,mboc(2,1)-1,mboc(3,1)-1:mboc(3,2)  )
         info%levellandscape(          0: cmx(1)  :2,          0,          0: cmx(3)  :2) &
            = parent%levellandscape(mboc(1,1)-1:mboc(1,2)  ,mboc(2,1)-1,mboc(3,1)-1:mboc(3,2)  )
!
! Front (z = 0):
         info%levellandscape(          1: cmx(1)+1:2,          1: cmx(2)+1:2,          0) &
            = parent%levellandscape(mboc(1,1)  :mboc(1,2)+1,mboc(2,1)  :mboc(2,2)+1,mboc(3,1)-1)
         info%levellandscape(          0: cmx(1)  :2,          1: cmx(2)+1:2,          0) &
            = parent%levellandscape(mboc(1,1)-1:mboc(1,2)  ,mboc(2,1)  :mboc(2,2)+1,mboc(3,1)-1)
         info%levellandscape(          1: cmx(1)+1:2,          0: cmx(2)  :2,          0) &
            = parent%levellandscape(mboc(1,1)  :mboc(1,2)+1,mboc(2,1)-1:mboc(2,2)  ,mboc(3,1)-1)
         info%levellandscape(          0: cmx(1)  :2,          0: cmx(2)  :2,          0) &
            = parent%levellandscape(mboc(1,1)-1:mboc(1,2)  ,mboc(2,1)-1:mboc(2,2)  ,mboc(3,1)-1)
!
! Back  (z = cmx+1):
         info%levellandscape(          1: cmx(1)+1:2,          1: cmx(2)+1:2, cmx(3)  +1) &
            = parent%levellandscape(mboc(1,1)  :mboc(1,2)+1,mboc(2,1)  :mboc(2,2)+1,mboc(3,2)+1)
         info%levellandscape(          0: cmx(1)  :2,          1: cmx(2)+1:2, cmx(3)  +1) &
            = parent%levellandscape(mboc(1,1)-1:mboc(1,2)  ,mboc(2,1)  :mboc(2,2)+1,mboc(3,2)+1)
         info%levellandscape(          1: cmx(1)+1:2,          0: cmx(2)  :2, cmx(3)  +1) &
            = parent%levellandscape(mboc(1,1)  :mboc(1,2)+1,mboc(2,1)-1:mboc(2,2)  ,mboc(3,2)+1)
         info%levellandscape(          0: cmx(1)  :2,          0: cmx(2)  :2, cmx(3)  +1) &
            = parent%levellandscape(mboc(1,1)-1:mboc(1,2)  ,mboc(2,1)-1:mboc(2,2)  ,mboc(3,2)+1)
!
       CASE DEFAULT
         PRINT *, 'BSAM 2.0, FindCoarseLevelNeighbors: only supports ndims=2,3.'
      END SELECT
!
   END FUNCTION FindCoarseLevelNeighbors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE FindMeshDefects(level)
! 功能：基于 levellandscape 在 level 层检测“非相容网格”（例如跨越两层落差的连接），
!       并在 level-2 层标记需要再细化的单元以修复缺陷。
! 步骤：收集缺陷→在 level-2 扩张一层缓冲（含周期偏移）→销毁链表。
      USE NodeInfoDef
      USE Boundary, ONLY: GetPeriodicOffsets
      USE TreeOps, ONLY: ApplyOnLevel
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: level
!
      TYPE(funcparam):: dummy
      LOGICAL:: periodicbuffer
      INTEGER:: offset, polarity
      INTEGER, DIMENSION(1:maxdims):: coordinatesave
!
      ALLOCATE(zerothtaggedcell)
      NULLIFY(zerothtaggedcell%prevcell)
      lasttaggedcell => zerothtaggedcell
      ntaggedcells = 0
!
! Get a list of defects. Tagged cells, if any will be located on level-2:
      CALL ApplyOnLevel(level,FindAndListPatchDefects,dummy)
!
      IF(ntaggedcells > 0) THEN
         coordinatesave = 1
         currenttaggedcell => lasttaggedcell
         searchloop: DO
            IF(.NOT. ASSOCIATED(currenttaggedcell%prevcell)) EXIT searchloop
!
! Ordinary buffering of an edge tag:
            dummy%iswitch = 1 ! Buffer by one cell only:
            CALL ApplyOnLevel(level-2,BufferTaggedCells,dummy)
!
! Buffering of periodic edge tags. Check to see if the buffer area cuts across
! a periodic boundary.  If so, add offset and apply buffer:
            coordinatesave(1:ndims) = currenttaggedcell%coordinate(1:ndims)
            IF(periodicboundaryconditions) THEN
               CALL GetPeriodicOffsets(level-2) ! Tags are on level-2:
               DO polarity = -1, 1, 2
                  DO offset = 1, nperiodicoffsets
                     currenttaggedcell%coordinate(1:ndims) &
                        = coordinatesave(1:ndims)+polarity*poffset(1:ndims,offset)
                     dummy%iswitch = 1 ! Buffer by one cell only:
                     CALL ApplyOnLevel(level-2,BufferTaggedCells,dummy)
                  END DO
               END DO
            END IF
            currenttaggedcell%coordinate(1:ndims) = coordinatesave(1:ndims)
!
            currenttaggedcell => currenttaggedcell%prevcell
         END DO searchloop
!
         CALL DeleteTaggedCellsList
      END IF
!
      NULLIFY(lasttaggedcell)
      NULLIFY(currenttaggedcell)
      DEALLOCATE(zerothtaggedcell)
!
   END SUBROUTINE FindMeshDefects
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   INTEGER FUNCTION FindAndListPatchDefects(info,dummy)
! 功能：在单个网格上，根据 levellandscape 的周边值是否等于 mylevel-2 来判定是否
!       存在缺陷；若有，则将相应 coarse-coarse 单元的全局坐标加入链表。
! 说明：2D 检查四边与四角；3D 则检查六面、十二条棱与八个角。
      USE NodeInfoDef
      USE TreeOps, ONLY: err_ok
      IMPLICIT NONE
!
! This routine finds mesh defects and tags for re-refinement if necessary:
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      INTEGER:: i, ii, j, jj, k, kk, mylevel
      INTEGER, DIMENSION(1:maxdims):: cmx, ind
      INTEGER, DIMENSION(1:maxdims,1:2):: mg
!
      FindAndListPatchDefects = err_ok
!
! Check for inactive grid awaiting garbage collection
      IF(info%tobedeleted .OR. (.NOT. info%activegrid)) RETURN
!
      cmx = 1
      cmx(1:ndims) = info%mx(1:ndims)/2
!
      mg = 1
      mg(1:ndims,1:2) = info%mglobal(1:ndims,1:2)
!
      mylevel = info%level
!
      ind = 1
!
      SELECT CASE(ndims)
       CASE(2)
!
! Faces:
! Right and Left:
         DO j = 1, cmx(2)
            DO i = 1, 2
               ii = cmx(1)+1-(2-i)*(cmx(1)+1) ! 0 and cmx(1)+1
               IF(info%levellandscape(ii,j,1) == mylevel-2) THEN
                  defectivegridlevel(mylevel) = .TRUE.
                  info%defective = .TRUE.
                  ind(1:maxdims) = CoarseCoarseGlobalIndex(mg(1,i),mg(2,1)+2*j-1,1)
                  CALL AddTaggedCellToList(ind(1:maxdims))
               END IF
            END DO
         END DO
!
! Top and Bottom:
         DO i = 1, cmx(1)
            DO j = 1, 2
               jj = cmx(2)+1-(2-j)*(cmx(2)+1)
               IF(info%levellandscape(i,jj,1) == mylevel-2) THEN
                  defectivegridlevel(mylevel) = .TRUE.
                  info%defective = .TRUE.
                  ind(1:maxdims) = CoarseCoarseGlobalIndex(mg(1,1)+2*i-1,mg(2,j),1)
                  CALL AddTaggedCellToList(ind(1:maxdims))
               END IF
            END DO
         END DO
!
! Four corners points:
         DO i = 1, 2
            DO j = 1, 2
               ii = cmx(1)+1-(2-i)*(cmx(1)+1)
               jj = cmx(2)+1-(2-j)*(cmx(2)+1)
               IF(info%levellandscape(ii,jj,1) == mylevel-2) THEN
                  defectivegridlevel(mylevel) = .TRUE.
                  info%defective = .TRUE.
                  ind(1:maxdims) = CoarseCoarseGlobalIndex(mg(1,i),mg(2,j),1)
                  CALL AddTaggedCellToList(ind(1:maxdims))
               END IF
            END DO
         END DO
!
       CASE(3)
!
! Faces:
! Right and Left
         DO j = 1, cmx(2)
            DO k = 1, cmx(3)
               DO i = 1, 2
                  ii = cmx(1)+1-(2-i)*(cmx(1)+1) ! 0 and cmx(1)+1
                  IF(info%levellandscape(ii,j,k) == mylevel-2) THEN
                     defectivegridlevel(mylevel) = .TRUE.
                     info%defective = .TRUE.
                     ind(1:maxdims) = CoarseCoarseGlobalIndex(mg(1,i),mg(2,1)+2*j-1,mg(3,1)+2*k-1)
                     CALL AddTaggedCellToList(ind(1:maxdims))
                  END IF
               END DO
            END DO
         END DO
!
! Top and Bottom
         DO i = 1, cmx(1)
            DO k = 1, cmx(3)
               DO j = 1, 2
                  jj = cmx(2)+1-(2-j)*(cmx(2)+1)
                  IF(info%levellandscape(i,jj,k) == mylevel-2) THEN
                     defectivegridlevel(mylevel) = .TRUE.
                     info%defective = .TRUE.
                     ind(1:maxdims) = CoarseCoarseGlobalIndex(mg(1,1)+2*i-1,mg(2,j),mg(3,1)+2*k-1)
                     CALL AddTaggedCellToList(ind(1:maxdims))
                  END IF
               END DO
            END DO
         END DO
!
! Front and Back
         DO i = 1, cmx(1)
            DO j = 1, cmx(2)
               DO k = 1, 2
                  kk = cmx(3)+1-(2-k)*(cmx(3)+1)
                  IF(info%levellandscape(i,j,kk) == mylevel-2) THEN
                     defectivegridlevel(mylevel) = .TRUE.
                     info%defective = .TRUE.
                     ind(1:maxdims) = CoarseCoarseGlobalIndex(mg(1,1)+2*i-1,mg(2,1)+2*j-1,mg(3,k))
                     CALL AddTaggedCellToList(ind(1:maxdims))
                  END IF
               END DO
            END DO
         END DO
!
! Twelve edges
! x and y fixed
!
         DO k = 1, cmx(3)
            DO i = 1, 2
               DO j = 1, 2
                  ii = cmx(1)+1-(2-i)*(cmx(1)+1)
                  jj = cmx(2)+1-(2-j)*(cmx(2)+1)
                  IF(info%levellandscape(ii,jj,k) == mylevel-2) THEN
                     defectivegridlevel(mylevel) = .TRUE.
                     info%defective = .TRUE.
                     ind(1:maxdims) = CoarseCoarseGlobalIndex(mg(1,i),mg(2,j),mg(3,1)+2*k-1)
                     CALL AddTaggedCellToList(ind(1:maxdims))
                  END IF
               END DO
            END DO
         END DO
!
! x and z fixed
!
         DO j = 1, cmx(2)
            DO i = 1, 2
               DO k = 1, 2
                  ii = cmx(1)+1-(2-i)*(cmx(1)+1)
                  kk = cmx(3)+1-(2-k)*(cmx(3)+1)
                  IF(info%levellandscape(ii,j,kk) == mylevel-2) THEN
                     defectivegridlevel(mylevel) = .TRUE.
                     info%defective = .TRUE.
                     ind(1:maxdims) = CoarseCoarseGlobalIndex(mg(1,i),mg(2,1)+2*j-1,mg(3,k))
                     CALL AddTaggedCellToList(ind(1:maxdims))
                  END IF
               END DO
            END DO
         END DO
!
! y and z fixed
!
         DO i = 1, cmx(1)
            DO j = 1, 2
               DO k = 1, 2
                  jj = cmx(2)+1-(2-j)*(cmx(2)+1)
                  kk = cmx(3)+1-(2-k)*(cmx(3)+1)
                  IF(info%levellandscape(i,jj,kk) == mylevel-2) THEN
                     defectivegridlevel(mylevel) = .TRUE.
                     info%defective = .TRUE.
                     ind(1:maxdims) = CoarseCoarseGlobalIndex(mg(1,1)+2*i-1,mg(2,j),mg(3,k))
                     CALL AddTaggedCellToList(ind(1:maxdims))
                  END IF
               END DO
            END DO
         END DO
!
! Eight corners points:
!
         DO i = 1, 2
            DO j = 1, 2
               DO k = 1, 2
                  ii = cmx(1)+1-(2-i)*(cmx(1)+1)
                  jj = cmx(2)+1-(2-j)*(cmx(2)+1)
                  kk = cmx(3)+1-(2-k)*(cmx(3)+1)
                  IF(info%levellandscape(ii,jj,kk) == mylevel-2) THEN
                     defectivegridlevel(mylevel) = .TRUE.
                     info%defective = .TRUE.
                     ind(1:maxdims) = CoarseCoarseGlobalIndex(mg(1,i),mg(2,j),mg(3,k))
                     CALL AddTaggedCellToList(ind(1:maxdims))
                  END IF
               END DO
            END DO
         END DO
!
!
      END SELECT
!
   END FUNCTION FindAndListPatchDefects
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   FUNCTION CoarseCoarseGlobalIndex(i,j,k) RESULT(indexres)
! 功能：给定细网格全局点 (i,j,k)，返回“向下两层”粗网格中包含该点的 cell 的全局索引。
! 实现：两次除以 2（带奇偶修正），在各向上落到 coarse-coarse 的索引格点上。
      USE NodeInfoDef
      IMPLICIT NONE
!
! Given global indices of point p=(i,j,k), find the global indices of a cell
! two levels below (coarse-coarse grid) that contains p:
!
      INTEGER:: i, j,k
      INTEGER, DIMENSION(1:maxdims):: indexres
!
      indexres = 1
!
      indexres(1) = (i+MODULO(i,2))/2
      indexres(2) = (j+MODULO(j,2))/2
      indexres(3) = (k+MODULO(k,2))/2
!
      indexres(1:ndims) = (indexres(1:ndims)+MODULO(indexres(1:ndims),2))/2
!
   END FUNCTION CoarseCoarseGlobalIndex
!
END MODULE BSAMRoutines

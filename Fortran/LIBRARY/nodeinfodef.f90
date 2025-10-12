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
      PRIVATE   ! 本模块默认私有（仅显式标注为 PUBLIC 的符号对外可见）
!
      LOGICAL, PUBLIC:: fluxbalancing, getafterstepstats, makeconformingmesh, &
            meshbuildcomplete, outputinitialdata, outputuniformmesh, &
            periodicboundaryconditions, restart, syncelliptic
      ! fluxbalancing                    : （可选）是否进行通量守恒/质量平衡修正
      ! getafterstepstats                : 是否在每个时间步后收集统计量（积分、网格规模等）
      ! makeconformingmesh               : 构网时是否强制网格“符合性”（最多一层悬挂点）
      ! meshbuildcomplete                : 当前 AMR 网格层级是否完成构建（用于控制删除旧网格）
      ! outputinitialdata                : 是否输出初始场数据（第一帧）
      ! outputuniformmesh                : 是否输出各层统一网格（Uniform grid）数据
      ! periodicboundaryconditions       : 是否启用周期性边界条件
      ! restart                          : 是否从重启文件恢复
      ! syncelliptic                     : 是否执行椭圆场同步（例如首次将 dt 置零的同步阶段）
!
! Double (r8) and extended (rext) precision if available:
      INTEGER, PARAMETER, PUBLIC:: r8 = SELECTED_REAL_KIND(15,307), &
            int2 = SELECTED_INT_KIND(2)
      ! r8   : 双精度实数 kind（约 15 位十进制有效数字，指数范围至 10^307）
      ! int2 : 短整型 kind（可存储小整数，常用于紧凑存储层级景观等）
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
      ! errflagdefault : 默认误差标记策略（相对截断误差）
      ! errflaguser    : 用户自定义误差标记策略的代号
      ! internalbc     : 内部边界的占位编码（非物理边界），用于块间接口
      ! maxsubgrids    : 单次分裂生成子网格的最大数量上限
      ! maxdims        : 最大空间维度（支持 2D/3D，此处为 3 的上限）
      ! maxdepth       : 最大层级深度（0..maxdepth）
      ! maxnccv        : 最大 CC 分量数上限（如 q 的通道数）
      ! rootlevel      : 根层级编号（约定为 0）
      ! sourcefield    : 泛型“字段选择器”常量：源项字段编号
      ! solutionfield  : 泛型“字段选择器”常量：解场（solution）字段编号
      ! auxiliaryfield : 泛型“字段选择器”常量：辅助场（aux）字段编号
      LOGICAL, DIMENSION(0:maxdepth), PUBLIC:: defectivegridlevel   ! 标记各层是否存在不合规网格（需先修复）
      INTEGER, PUBLIC:: amrrestarts, errortype, finestlevel, gridnumber, maxvcycles, &
            maxlevel, mbc, minlevel, ndims, nmasscorrlayers, naxv, &
            nccv, nfcv, nperiodicoffsets, nrootgrids, nsmoothingpasses, &
            ntaggedcells, outframes, restartframe, timeiterations, &
            totalmeshsize, updateauxfreq
      ! amrrestarts        : AMR 重建重试计数（用于修复网格缺陷时回退重构）
      ! errortype          : 误差度量类型（1=合成/混合 L2；2=分量级 L2 等）
      ! finestlevel        : 当前构网后最细层级编号
      ! gridnumber         : 调试/打印中使用的网格计数器
      ! maxvcycles         : 每个时间步内最多执行的 V-cycle 次数
      ! maxlevel           : 允许细化到的最大层级
      ! mbc                : 幽灵层格点厚度（支持 1 或 2）
      ! minlevel           : 多重网格的最粗层级（≤ 0，通常为 0 或负数作为 seed）
      ! ndims              : 维度（2 或 3）
      ! nmasscorrlayers    : 质量修正层（边带）链表的层数计数
      ! naxv               : 辅助变量通道数（aux 的分量数）
      ! nccv               : CC 变量通道数（q 的分量数）
      ! nfcv               : FC 变量通道数（v1/v2/v3 的分量数）
      ! nperiodicoffsets   : 周期性偏移向量的数量（每层）
      ! nrootgrids         : 根层块数（当前实现通常为 1）
      ! nsmoothingpasses   : 每层每次 LevelRelax 的平滑步（红黑各一次为 2）
      ! ntaggedcells       : 被误差估计标记的单元数（用于 BR 分裂）
      ! outframes          : 计划输出帧数（时间段内总帧数）
      ! restartframe       : 重启时的起始帧编号
      ! timeiterations     : 时间步总数（整体迭代次数）
      ! totalmeshsize      : 合成网格尺寸（总单元数）
      ! updateauxfreq      : V-cycle 内调用 UpdateAux 的频率（步数间隔）
      INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC:: periodicoffsetindex   ! 周期偏移索引表（索引到 poffset 列表）
      INTEGER, DIMENSION(0:maxdepth), PUBLIC:: errflagopt, ibuffer, minimumgridpoints
      ! errflagopt(level)     : 每层误差标记策略选择（默认/用户）
      ! ibuffer(level)        : 对被标记单元的缓冲层厚度（用于膨胀标记区域）
      ! minimumgridpoints(L)  : 子块的最小尺寸（每维）约束
      INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC:: poffset              ! 周期性偏移向量集合（ndims × nperiodicoffsets）
      INTEGER, DIMENSION(0:maxdepth,1:maxdims), PUBLIC:: mxmax            ! 各层最大可能的网格尺寸（按 r=2 递推）
!
      REAL(KIND=r8), PUBLIC:: currenttime, dt, finaltime, omega, qerrortol, &
            restarttime
      ! currenttime  : 根层当前时间（等级间共享参考时间）
      ! dt           : 时间步长（若 syncelliptic 初始同步阶段为 0.0）
      ! finaltime    : 目标结束时间（基于帧数与 itperprint 推导）
      ! omega        : 松弛因子（如 SOR/红黑 GS 的加权系数）
      ! qerrortol    : 收敛阈值（误差度量的容差）
      ! restarttime  : 重启时刻（读档恢复的时间）
      REAL(KIND=r8), DIMENSION(1:2), PUBLIC:: integralresult   ! 统计量积累（例如质量/能量等 2 个指标的积分值）
      REAL(KIND=r8), DIMENSION(1:maxnccv), PUBLIC:: componentintegral   ! 各 CC 分量的分量级积分统计
      REAL(KIND=r8), DIMENSION(0:maxdepth), PUBLIC:: desiredfillratios, qtolerance
      ! desiredfillratios(L) : 目标填充率（子块面积/标记区域）阈值，用于 BR 分裂停止准则
      ! qtolerance(L)        : 误差标记/截断误差阈值（维度一致性见 2D/3D 处理）
!
   TYPE, PUBLIC:: taggedcell
            INTEGER:: id                                           ! 节点在链表中的递增 ID
            INTEGER, DIMENSION(1:maxdims):: coordinate             ! 标记单元的全局坐标（按 coarse/coarse 对齐）
            TYPE(taggedcell), POINTER:: prevcell                   ! 指向前一个标记单元（后进先出栈式链表）
   END TYPE taggedcell
!
      TYPE(taggedcell), POINTER, PUBLIC:: zerothtaggedcell      ! 标记单元链表的哨兵/零节点
      TYPE(taggedcell), POINTER, PUBLIC:: currenttaggedcell     ! 遍历/构造过程中的当前节点指针
      TYPE(taggedcell), POINTER, PUBLIC:: lasttaggedcell        ! 链表尾（最新加入）的指针
!
   TYPE, PUBLIC:: masscorrlayer
            INTEGER:: id                                           ! 质量修正层的编号
            INTEGER, DIMENSION(1:maxdims,1:2):: mg                 ! 修正条带在全局/局部的索引范围 [min,max]
            REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: masscorr  ! 质量修正的通量/累积数组（与 f/rf 同型）
            TYPE(masscorrlayer), POINTER:: prevlayer               ! 指向上一层修正条带（链表）
   END TYPE masscorrlayer
!
      TYPE(masscorrlayer), POINTER, PUBLIC:: zerothlayer        ! 质量修正链表的哨兵节点
      TYPE(masscorrlayer), POINTER, PUBLIC:: currentlayer       ! 当前修正层指针
      TYPE(masscorrlayer), POINTER, PUBLIC:: lastlayer          ! 最新的修正层指针
!
! Uniform grids used for restarting and output:
   TYPE, PUBLIC:: uniformgridtype
            INTEGER, DIMENSION(1:maxdims):: mx                     ! 该层统一网格的尺寸（各维）
            REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: q         ! 统一网格上的 CC 场缓存（用于重启/统一输出）
   END TYPE uniformgridtype
      TYPE(uniformgridtype), DIMENSION(0:maxdepth), PUBLIC:: uniformgrid  ! 各层统一网格缓冲区
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   TYPE, PUBLIC:: nodeinfo
!
! This must be the first component to ensure proper parallel communication:
      INTEGER:: nodeinfostart                                ! 结构体起始标记（通信/序列化时用于边界定位）
!
! LOGICAL variables !!!!!!!!!
!
! A necessary component. Mark for garbage collection:
      LOGICAL:: tobedeleted                                  ! 标记此网格是否待删除（垃圾回收）
!
! Flag to show whether grid may accept values from elder siblings on this
!   level:
      LOGICAL:: activegrid                                   ! 是否为当前活动网格（可参与计算/收发）
!
! All patches are probationary until the mesh hierarchy is fully constructed:
      LOGICAL:: defective                                    ! 是否为不合规网格（构网阶段发现挂点需修复）
!
! Flag to show whether grid has been created during a restart from checkpoint
!   file:
      LOGICAL:: restartgrid                                  ! 是否为重启时新建的网格
!
! Flag showing allocation status of fields within this node:
      LOGICAL:: fieldsallocated                              ! 本节点内场变量是否已分配内存
!
! Flag to show whether this is an initial grid, i.e. created during start-up:
      LOGICAL:: initialgrid                                  ! 是否为初始化阶段生成的网格（非重启）
!
! INTEGER variables !!!!!!!!!
!
      INTEGER:: maxlevel               ! Maximum level to which this grid may be  ! 该网格允许细化到的最大层级
      !   refined:
      INTEGER:: ngrid                  ! Number of this grid:                      ! 网格（块）的随机/唯一编号（调试/追踪用）
      INTEGER:: nsubgrids              ! Number of child grids:                    ! 子网格数量（由 BR 分裂生成）
      INTEGER:: level                  ! Level on which this node lives:           ! 所在层级编号
!
! Number of grid cells in q along each dimension:
      INTEGER, DIMENSION(1:maxdims):: mx                   ! q 的各向尺寸（单元数），决定 CC 网格大小
!
! Index bounds within parent where this child was created:
      INTEGER, DIMENSION(1:maxdims,1:2):: mbounds          ! 在父网格中的创建范围（父网格局部索引）
!
! Index bounds of this grid in global indexing system:
      INTEGER, DIMENSION(1:maxdims,1:2):: mglobal          ! 全局索引范围（用于同层重叠/转移/统计）
!
! Boundary condition codes:
!   1 - left, 2 - right, 3 - bottom, 4 - top, 5 - back, 6 - front:
      INTEGER, DIMENSION(1:2*maxdims):: mthbc              ! 物理/内部边界编码（每维左右各一：1..6）
!
! Array of error flags:
      INTEGER, DIMENSION(:,:,:), POINTER:: errorflags       ! 误差标记掩码（1=需要细化，0=不需要）
!
! Pointers to the grid's nearest neighbors' levels:
      INTEGER(KIND=int2), DIMENSION(:,:,:), POINTER::  levellandscape  ! 邻域触及的“对面块”的层级景观（用于缺陷检测）
!
! REAL variables !!!!!!!!!!!!!!!!!!!!!!
!
! Lower coordinates for this grid:
      REAL(KIND=r8), DIMENSION(1:maxdims):: xlower         ! 网格物理域左下（或左下后）坐标
!
! Upper coordinates for this grid:
      REAL(KIND=r8), DIMENSION(1:maxdims):: xupper         ! 网格物理域右上（或右上前）坐标
!
! The current time at which this grid exists:
      REAL(KIND=r8):: gridtime                             ! 网格当前时间戳（与根层 currenttime 同步）
!
! Grid spacings:
      REAL(KIND=r8), DIMENSION(1:maxdims):: dx             ! 网格步长（各向，通常等距且各向相等）
!
! Pointer to field variable arrays:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: q       ! CC 变量（解场）数组 [i,j,(k),nccv]，含幽灵层
!
! Pointer to the velocity variable on the 1-2 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v1      ! FC 变量 x-向（1-2 面） [i-face,j,(k),nfcv]
!
! Pointer to the velocity variable on the 3-4 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v2      ! FC 变量 y-向（3-4 面） [i,j-face,(k),nfcv]
!
! Pointer to the velocity variable on the 5-6 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v3      ! FC 变量 z-向（5-6 面） [i,j,k-face,nfcv]
!
! Pointer to field variable arrays at previous time:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: qold    ! 上一时间步的 CC 变量快照
!
! Pointer to the velocity variable on the 1-2 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v1old   ! 上一时间步的 v1 快照
!
! Pointer to the velocity variable on the 3-4 faces at previous time:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v2old   ! 上一时间步的 v2 快照
!
! Pointer to the velocity variable on the 5-6 faces at previous time:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v3old   ! 上一时间步的 v3 快照
!
! Pointer to the coarse-under-fine field variable arrays:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: qc      ! 粗网-覆盖-细网的 CC 拷贝（CUF，用于插值/限制）
!
! Pointer to the coarse-under-fine field variable arrays:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: qcold   ! qc 的旧值快照（用于误差估计/比较）
!
! Pointer to the coarse-under-fine v1 variable array:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v1c     ! v1 的粗网-覆盖-细网拷贝（CUF）
!
! Pointer to the coarse-under-fine v2 variable array:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v2c     ! v2 的粗网-覆盖-细网拷贝（CUF）
!
! Pointer to the coarse-under-fine v2 variable array:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: v3c     ! v3 的粗网-覆盖-细网拷贝（CUF）
!
! Pointer to the coarse-under-fine relative truncation error:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: qrte    ! 相对截断误差（coarse 下的估计，维度按 cmx）
!
! Pointer to the load function:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f       ! CC 负载/右端项（FAS/粗层加载）
!
! Pointer to the temporary function.  Same size as f and rf:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: ftmp    ! 临时负载缓冲，与 f/rf 同形
!
! Pointer to the load function on the 1-2 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f1      ! x-向 FC 负载（对应 v1 的一维面量）
!
! Pointer to the load function on the 1-2 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f1tmp   ! f1 的临时缓冲
!
! Pointer to the load function on the 3-4 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f2      ! y-向 FC 负载（对应 v2 的一维面量）
!
! Pointer to the load function on the 3-4 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f2tmp   ! f2 的临时缓冲
!
! Pointer to the load function on the 3-4 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f3      ! z-向 FC 负载（对应 v3 的一维面量）
!
! Pointer to the load function on the 3-4 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: f3tmp   ! f3 的临时缓冲
!
! Pointer to the cell-centered residual:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: rf      ! CC 残差数组（Operator(q)-f 等）
!
! Pointer to the residual on the 1-2 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: rf1     ! x-向 FC 残差
!
! Pointer to the residual on the 3-4 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: rf2     ! y-向 FC 残差
!
! Pointer to the residual on the 3-4 faces:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: rf3     ! z-向 FC 残差
!
! Pointer to auxiliary arrays:
      REAL(KIND=r8), DIMENSION(:,:,:,:), POINTER:: aux     ! 辅助场（用户问题相关的额外变量）
!
! This must be the last component to ensure proper parallel communication:
      INTEGER:: nodeinfoend                                  ! 结构体结束标记（通信/序列化边界）
!
   END TYPE nodeinfo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   TYPE, PUBLIC:: funcparam
!
      INTEGER:: iswitch                                     ! 通用整型开关（在回调/遍历中传参）
      INTEGER, DIMENSION(1:maxdims):: offset               ! 通用偏移量（如周期偏移/方向选择）
!
      TYPE(nodeinfo), POINTER:: info                        ! 可选：关联的 nodeinfo 指针（用于特定回调）
!
   END TYPE funcparam
!
END MODULE NodeInfoDef

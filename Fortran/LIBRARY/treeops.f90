!==========================================================================
! BSAM 2.0: Block-Structured Adaptive Multigrid Solver
!==========================================================================
!
! WiseSoft: Innovators, Brothers.
!
! (c) Copyright Sorin Mitran, 2000
! Department of Applied Mathematics
! University of Washington
! mitran@amath.washington.edu
!
! (c) Copyright Sorin Mitran, 2002
! Department of Mathematics
! University of North Carolina at Chapel Hill
! mitran@amath.unc.edu
!
! Portions of the code
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
! This software is made available "as is" without any assurance that it
! will work for your purposes. The software may in fact have defects,
! so use the software at your own risk.
!
! -----------------------------------------------------------------------
! File:             treeops.f90
! Purpose:          Tree operations module.
! Contains:
! Revision History: Ver. 1.0 Oct. 2000 Sorin Mitran
! Revision History: Ver. 1.1 May. 2007 Steven Wise
! Revision History: Ver. 2.0 Jul. 2015 Steven Wise
! -----------------------------------------------------------------------

!*************
MODULE TreeOps
!*************
!
! ============================= 中文总览 / 使用说明 =============================
! 模块名称: TreeOps（树操作）
! 适用范围: BSAM 2.0 自适应多重网格（Block-Structured Adaptive Multigrid）框架
!
! 本模块实现了一个“森林/树/节点”的层次数据结构及其遍历与维护操作：
! - ForestSeed(森林种子) 是所有树的虚拟根，level = FOREST_SEED（负值）；
! - Root(根结点) 是实际问题域的根，level = rootlevel（通常为 0，由 NodeInfoDef 提供）；
! - 每个节点 Node 之间通过 4 类指针建立关系：
!     Parent（父亲，指向上一层）、Child（孩子，指向最年轻的孩子）、
!     Sibling（兄长，指向比自己年长的兄弟）、Neighbor（同层相邻节点，横向链表）。
! - Level 语义：
!     0 = 根层；正整数 = 向下更细的层级；负整数 = 位于根之上的“种子层”(仅用于包装/对齐)；
!     MaxDepth/RootLevel 在 NodeInfoDef 中定义。
! - 同层节点通过 Neighbor 形成从“最年轻(Youngest) -> … -> 最年长(Eldest)”的单向链；
!   Sibling 是“同父兄弟”之间按年龄的单向链：NewChild%Sibling => Parent%Child（即指向更年长者）。
!
! 功能要点：
! 1) 初始化/销毁：InitForest, KillForest, CreateBelowSeedLevels 等；
! 2) 结构修改：AddRootLevelNode, CreateChild, KillNode, DeleteMarkedNode；
! 3) 遍历应用：ApplyOnForest/Levels/Leaves/Level/LevelPairs/Children/MyLevelNbrs；
! 4) 查询与工具：Get* 家族、ExistLevel、PushForest/PopForest、SetLevelNodeNumbers；
! 5) 遍历策略：PreEval(先对子节点/同级准备再评估) 或 PostEval(先评估当前再下探)。
!
! 重要约定：
! - 当前遍历上下文由 CurrentNode/CurrentLevel 表示，并通过 Stack 保存递归状态；
! - YoungestOnLevel/LdestOnLevel/NrOfNodes 维护每层的横向邻接链首尾与计数；
! - globalnodenum 提供全局递增节点编号；levelnodenum 提供同层顺序号（通过 SetLevelNodeNumbers 设置）。
!
! 开发者提示：
! - 本补充为中文注释，仅新增注释行，不改变任何原有可执行代码；
! - NodeInfoDef 模块定义了 nodeinfo 类型、根层级 rootlevel、最大深度 MaxDepth、
!   以及用户回调所需的 funcparam 等；
! - 错误码统一放在本模块顶部，ErrCode 保存最近一次错误，使用 err_OK 等常量判断。
! ==========================================================================
   ! User definitions of data structures associated with each node
   ! 用户定义的节点信息（nodeinfo）、根层级 rootlevel、最大深度 MaxDepth 等均在 NodeInfoDef 中声明。
   USE NodeInfoDef
   IMPLICIT NONE
   PRIVATE  ! Everything is implicitly PRIVATE
   PUBLIC   InitForest,KillForest,AddRootLevelNode,          &
      CreateChild,KillNode,DeleteMarkedNode,           &
      ApplyOnLeaves,ApplyOnLevels,ApplyOnForest,       &
      ApplyOnLevel,ApplyOnLevelPairs,ApplyOnChildren,  &
      ApplyOnMyLevelNbrs,GetTreeOpsErrCode,            &
      GetNodeInfo,SetNodeInfo,GetCurrentNode,          &
      Getglobalnodenum,GetChildInfo,GetParentInfo,GetParent,  &
      GetSibling,GetChild,GetLevel,GetSiblingIndex,    &
      GetNrOfChildren,GetRootInfo,SetRootInfo,         &
      ExistLevel, PushForest,PopForest,                &
      CreateBelowSeedLevels,CurrentNodeToYoungest,     &
      SetLevelNodeNumbers
! InitForest: 初始化森林结构，创建 ForestSeed 与 Root，并设置全局遍历上下文
! KillForest: 递归释放整棵树（从 Root 起）及其节点信息内存
! AddRootLevelNode: 在根层新增最年轻根节点并更新 Root 与同层邻接链
! CreateChild: 为 CurrentNode 新建一个子节点并维护层级、兄弟链与同层邻接链
! KillNode: 递归删除指定节点及其全部子孙并修复父/兄弟/同层链接
! DeleteMarkedNode: 若 CurrentNode 的 info 被标记删除则调用 KillNode 删除之
! ApplyOnLeaves: 对所有叶子节点调用用户回调函数
! ApplyOnLevels: 对所有层（满足条件的节点）调用用户回调函数
! ApplyOnForest: 遍历整个森林并对所有节点调用用户回调函数
! ApplyOnLevel: 对指定层的每个节点调用用户回调函数
! ApplyOnLevelPairs: 对指定层内所有有序节点对依次调用双参回调
! ApplyOnChildren: 对 CurrentNode 的每个子节点调用用户回调函数
! ApplyOnMyLevelNbrs: 让 CurrentNode 与同层所有其他节点逐一交互调用双参回调
! GetTreeOpsErrCode: 返回最近一次树操作的错误码
! GetNodeInfo: 读取给定节点的 info（值拷贝）
! SetNodeInfo: 写入给定节点的 info（值拷贝）
! GetCurrentNode: 返回当前遍历节点的指针
! Getglobalnodenum: 返回 CurrentNode 的全局节点编号
! GetChildInfo: 返回 CurrentNode 最年轻子节点的 info（指针）
! GetParentInfo: 返回 CurrentNode 父节点的 info（指针）
! GetParent: 返回给定节点的父节点指针
! GetSibling: 返回给定节点的兄长（更年长兄弟）指针
! GetChild: 返回给定节点的最年轻子节点指针
! GetLevel: 返回当前层级（并自检与 CurrentNode%level 一致性）
! GetSiblingIndex: 计算 CurrentNode 在兄弟链中的序号并校验 ChildNo 一致性
! GetNrOfChildren: 计算并返回 CurrentNode 的子节点数量（并与存储值核对）
! GetRootInfo: 返回 Root 的 info（指针）
! SetRootInfo: 写入 Root 的 info（值拷贝）
! ExistLevel: 判断某层是否存在任何节点
! PushForest: 保存当前遍历上下文以便执行子遍历
! PopForest: 恢复最近保存的遍历上下文
! CreateBelowSeedLevels: 在 Root 之上创建到给定最小层的“种子层”包装节点
! CurrentNodeToYoungest: 将当前节点切换到指定层的最年轻节点
! SetLevelNodeNumbers: 为指定层节点按邻接顺序设置 levelnodenum

!
! Error handling
! ----------------------------- 错误处理与错误码 ------------------------------
! ErrCode 保存最近一次错误码（模块内全局）。常用错误码：
! - err_OK:              正常
! - err_UndefinedNode:   节点指针未关联/不存在
! - err_NoParent:        无父节点（已到最上层或 forest seed）
! - err_NoSibling:       无兄弟节点（是家族中最年长）
! - err_NoChild:         无子节点（叶子）
! - err_BadLevel:        请求的层级越界/非法
! 其中 ErrPrint 用于控制错误打印的阈值。
   INTEGER  ErrCode      ! The most recent error code
   ! Error codes
   PUBLIC   err_OK,                                          &
      err_UndefinedNode,                               &
      err_NoParent,err_NoSibling,err_NoChild
   INTEGER, PARAMETER :: ErrPrint          = 1000
   INTEGER, PARAMETER :: err_OK            =    0
   INTEGER, PARAMETER :: err_UndefinedNode =  100
   INTEGER, PARAMETER :: err_NoParent      =  201
   INTEGER, PARAMETER :: err_NoSibling     =  202
   INTEGER, PARAMETER :: err_NoChild       =  203
   INTEGER, PARAMETER :: err_BadLevel      =  301

   INTEGER, PARAMETER :: IntTrue=0, IntFalse=-1
!
   ! ------------------------------ 核心数据结构 ------------------------------
   ! 类型: node
   ! 说明: 抽象出树中的一个网格块/结点，维护家族关系与层级信息。
   TYPE, PUBLIC:: node
      PRIVATE
      TYPE(nodeinfo), POINTER:: info      ! 用户数据载体，由 NodeInfoDef 定义
      TYPE(node), POINTER:: parent        ! 父节点（上一层）。在 forest seed 上可能为空
      TYPE(node), POINTER:: sibling       ! 兄长（更年长的兄弟）。单向链到最年长
      TYPE(node), POINTER:: child         ! 最年轻的孩子（下一层）。孩子链头
      TYPE(node), POINTER:: neighbor      ! 同层相邻节点（横向遍历的单向链）
      INTEGER:: level                     ! 所在层级。rootlevel 为 0，向下递增；森林种子为负数
      INTEGER:: childno                   ! 家族内的“出生序”。Eldest=1，Next=2，…（与 sibling 链一致）
      INTEGER:: nrofchildren              ! 子女数量（与实际 sibling 遍历一致性在某些检查中核验）
      INTEGER:: leafdist                  ! 距离叶子的距离（本实现中常作标记/保留）
      INTEGER:: globalnodenum             ! 全局节点唯一编号（创建时自增）
      INTEGER:: levelnodenum              ! 同层的顺序号（通过 SetLevelNodeNumbers 计算填充）
   END TYPE node
!
   ! 辅助指针包装类型：用于栈/数组中存放指向 node 的指针
   TYPE:: nodepointer
      TYPE(node), POINTER:: nodeptr
   END TYPE nodepointer
!
   ! 指向 nodeinfo 的辅助指针封装
   TYPE, PUBLIC:: infopointer
      TYPE (nodeinfo), POINTER :: infoptr
   END TYPE infopointer
!
   ! ------------------------------- 常量与约定 -------------------------------
   INTEGER, PARAMETER, PUBLIC :: FOREST_SEED = -999   ! 特殊标记：森林种子/虚拟层
   INTEGER, PARAMETER :: NO_CHILDREN = 0
   INTEGER, PARAMETER :: FIRST_CHILD = 1
   INTEGER, PARAMETER :: NOT_A_CHILD = 0
   INTEGER, PARAMETER :: BAD_CHILD_NO = 0
   INTEGER, PARAMETER :: ZERO_LEAF_NODE_DIST = 0

   LOGICAL, PARAMETER :: NoInfoInit = .FALSE.
   LOGICAL, PARAMETER :: PreEvalNext = .TRUE.
   LOGICAL, PARAMETER :: PostEvalNext = .FALSE.

   ! -------------------------- 遍历上下文（全局） ----------------------------
   ! ForestSeed: 森林虚拟根（其 child 指向第一个 Root）。
   ! Root:       第一棵树的根节点（rootlevel）。
   ! CurrentNode:当前遍历焦点结点（遍历/回调时使用）。
   ! Stack:      遍历栈，按层记录当前节点（便于在回调里查询祖先等）。
   ! YoungestOnLevel/EldestOnLevel: 每层的横向链首/尾；NrOfNodes: 每层节点数。
   TYPE (Node), POINTER, SAVE :: ForestSeed
   TYPE (Node), POINTER, SAVE :: Root
   TYPE (Node), POINTER, SAVE :: CurrentNode
   TYPE (NodePointer), DIMENSION(-MaxDepth:MaxDepth) :: Stack
   TYPE (NodePointer), DIMENSION(-MaxDepth:MaxDepth) :: YoungestOnLevel
   TYPE (NodePointer), DIMENSION(-MaxDepth:MaxDepth) :: EldestOnLevel
   INTEGER, DIMENSION(-MaxDepth:MaxDepth) :: NrOfNodes

   INTEGER:: currentlevel, globalnodeindex, levelnodeindex   ! 当前层级/全局节点计数/层内顺序

   ! 嵌套遍历支持：允许在一次遍历中暂存上下文，发起子遍历，再恢复
   INTEGER, PARAMETER :: MaxForestStacks=32
   INTEGER, SAVE :: StackLevel
   TYPE (NodePointer), SAVE, DIMENSION(0:MaxForestStacks) :: ForestStack

   ! 全局 I/O（供内部调试打印使用）
   INTEGER InputUnit,OutputUnit

   ! 内存使用（占位/历史遗留）
   INTEGER, PARAMETER :: iPrec = SELECTED_INT_KIND(9)
   !INTEGER (KIND=iPrec) :: MemFree, MemAllocated

   !-------
CONTAINS
   ! 这个关键词`CONTAINS`表示模块内包含子程序和函数的实现。上面是模块的声明部分，声明了函数和用到的变量
   !-------
   ! Interface routines are defined first.
   ! ============================= 过程与遍历设计 =============================
   ! 遍历核心由 ForestTraversal + TreeTraversal 组成：
   ! - ForestTraversal: 对根层（root 或多棵树的根）进行横向遍历，逐棵树调用 TreeTraversal；
   ! - TreeTraversal: 深度优先遍历一棵树，根据 cond() 判断是否对当前节点调用回调 f(Info,Param)。
   !   借助 EvalNextBefore_f 控制“先下子树再评估”或“先评估再下子树”（前序/后序的变体）。
   !   遍历时 CurrentNode/CurrentLevel 始终保持一致，用 Stack 记录以支持在回调中追溯。
   ! 常用高阶接口：
   ! - ApplyOnForest: 对整片森林的所有节点应用回调（通常用于全量初始化/收尾）；
   ! - ApplyOnLevels: 对所有层/特定层按 cond 逻辑应用回调；
   ! - ApplyOnLeaves: 仅对叶子层应用回调；
   ! - ApplyOnLevel/ApplyOnLevelPairs/ApplyOnMyLevelNbrs/ApplyOnChildren：
   !   分别对“同一层的每个节点 / 同层所有成对组合 / 当前节点与同层邻居 / 当前节点的子女”应用回调。
   ! ========================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: InitForest
   ! 作用: 初始化森林/树的基本骨架，创建 ForestSeed 与第一棵树的 Root，
   !       建立根层的 Youngest/Eldest/NrOfNodes，并为 Root 分配 nodeinfo。
   ! 参数: InfoInit [逻辑, 可选] 若提供且为 .TRUE.，表示调用者希望随后进行 info 的初始化；
   !       本过程不负责具体填充 info，仅分配占位。
   ! 前置: 需已 USE NodeInfoDef，以便获得 rootlevel、MaxDepth 等。
   ! 后置: CurrentNode 指向 Root，CurrentLevel = rootlevel，globalnodeindex = 1。
   ! 错误: 内存分配失败会 STOP；
   ! 注意: 本函数只搭建结构，不对用户数据 info 做实际初始化。
   SUBROUTINE InitForest(InfoInit)
      LOGICAL, INTENT(IN), OPTIONAL :: InfoInit
      INTEGER :: iError
      INTEGER :: L
      !
      ALLOCATE (Root,ForestSeed,STAT=iError)
      IF (iError /= 0) THEN
         PRINT *,"Error allocating tree/forest roots in InitForest."
         STOP
      END IF
      StackLevel=0

      ! Seed of all trees
      NULLIFY(ForestSeed%Parent) ! NULLIFY的意思是将指针设为空
      NULLIFY(ForestSeed%Sibling)
      NULLIFY(ForestSeed%Neighbor)

      ForestSeed%Child=>Root
      ForestSeed%level=FOREST_SEED
      ForestSeed%ChildNo=NOT_A_CHILD
      ForestSeed%NrOfChildren=1
      ForestSeed%LeafDist=ZERO_LEAF_NODE_DIST
      ForestSeed%globalnodenum=FOREST_SEED

      ! First tree root
      NULLIFY(Root%Sibling)
      NULLIFY(Root%Child)
      NULLIFY(Root%Neighbor)

      Root%Parent=>ForestSeed
      Root%level=rootlevel
      Root%ChildNo=FIRST_CHILD
      Root%NrOfChildren=NO_CHILDREN
      Root%LeafDist=ZERO_LEAF_NODE_DIST
      globalnodeindex=1
      Root%globalnodenum=globalnodeindex

      DO L=-MaxDepth,MaxDepth
         NULLIFY(YoungestOnLevel(L)%NodePtr)
         NULLIFY(EldestOnLevel(L)%NodePtr)
         NrOfNodes(L)=0
      END DO
      YoungestOnLevel(rootlevel)%NodePtr=>Root
      EldestOnLevel(rootlevel)%NodePtr=>Root
      NrOfNodes(rootlevel)=1
      ALLOCATE(Root%Info,STAT=iError)  ! Allocate space for NodeInfo
      IF (iError /= 0) THEN
         PRINT *,"Error allocating tree root Info in InitForest."
         STOP
      END IF
      ! Set current node pointer
      CurrentNode => Root; CurrentLevel=rootlevel
      !
      IF (.NOT. PRESENT(InfoInit)) RETURN
      IF (.NOT. InfoInit) RETURN
   END SUBROUTINE InitForest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: CreateBelowSeedLevels
   ! 作用: 在 ForestSeed 与 Root 之间创建额外的“种子层”（level = -1, -2, …），
   !       使得森林可以拥有位于根之上的层级包装。
   ! 参数: MinLevel [整型] 最小层级（负值，绝对值不得超过 MaxDepth）；
   !       InfoInit [逻辑, 可选] 同 InitForest 的注释；仅分配 info，不填充。
   ! 语义: 创建的每层仅有一个“父包装节点”，其 Child 指向下层（更靠近 Root）的单一孩子节点。
   ! 注意: 调用完成后会将 CurrentNode/CurrentLevel 复位回 Root/rootlevel。
   SUBROUTINE CreateBelowSeedLevels(MinLevel,InfoInit)
      INTEGER, INTENT(IN):: MinLevel
      LOGICAL, INTENT(IN), OPTIONAL :: InfoInit
      INTEGER :: iError
      !
      IF(MinLevel<-MaxDepth) THEN
         PRINT *,"Error MinLevel less than MaxDepth in CreateBelowSeedLevels."
         STOP
      END IF
      !
      ALLOCATE(ForestSeed%Info,STAT=iError)  ! Allocate space for NodeInfo
      IF (iError /= 0) THEN
         PRINT *,"Error allocating forest seed Info in CreateBelowSeedLevels."
         STOP
      END IF
      !
      ForestSeed%Level = -1
      !
      YoungestOnLevel(-1)%NodePtr => ForestSeed
      EldestOnLevel(-1)%NodePtr => ForestSeed
      NrOfNodes(-1) = 1

      IF(-1>MinLevel) THEN
         ForestSeed%ChildNo = 1
         CurrentNode => ForestSeed; CurrentLevel=-1
         CALL CreateBelowSeedLevelNode(MinLevel,InfoInit)
         ! Set current node pointer back to root.
         CurrentNode => Root; CurrentLevel=rootlevel
      END IF
      !
      IF (.NOT. PRESENT(InfoInit)) RETURN
      IF (.NOT. InfoInit) RETURN
   END SUBROUTINE CreateBelowSeedLevels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: CreateBelowSeedLevelNode (递归)
   ! 作用: 实际执行“在当前节点之上再包一层”的递归创建逻辑；每层仅生成一个父节点，
   !       并把当前节点挂到该父节点的 Child。
   ! 参数: MinLevel [整型] 递归到此层后停止（负值）；InfoInit [逻辑, 可选] 同上。
   RECURSIVE SUBROUTINE CreateBelowSeedLevelNode(MinLevel,InfoInit)
      ! Interface declarations
      INTEGER, INTENT(IN):: MinLevel
      LOGICAL, INTENT(IN), OPTIONAL:: InfoInit
      ! Internal declarations
      TYPE (Node), POINTER:: BelowSeedNode
      INTEGER:: iError, BelowSeedLevel
      !
      ! There are only one of these grids at each level below the root.
      ! Current node is assumed to be root node or lower.
      !
      ALLOCATE(BelowSeedNode,STAT=iError)
      IF (iError /= err_OK) THEN
         PRINT *,"Error allocating BelowSeedNode in treeops::CreateBelowSeedLevelNode."
         STOP
      END IF
      !
      NULLIFY(BelowSeedNode%Parent); NULLIFY(BelowSeedNode%Sibling)
      NULLIFY(BelowSeedNode%Child);  NULLIFY(BelowSeedNode%Neighbor)
      !
      BelowSeedNode%globalnodenum = FOREST_SEED
      BelowSeedLevel = CurrentNode%level-1       ! BelowSeedNode is one level up.
      !
      NrOfNodes(BelowSeedLevel) = 1
      BelowSeedNode%Child => CurrentNode         ! Parent is created for CurrentNode
      BelowSeedNode%NrOfChildren = 1             ! BelowSeedNodes have only one child.
      CurrentNode%Parent => BelowSeedNode
      !
      BelowSeedNode%ChildNo = NOT_A_CHILD
      !
      YoungestOnLevel(BelowSeedLevel)%NodePtr => BelowSeedNode  ! Start of level stack
      EldestOnLevel(BelowSeedLevel)%NodePtr => BelowSeedNode    ! End of level stack
      !
      BelowSeedNode%level = BelowSeedLevel       ! Set the below-root level counter
      BelowSeedNode%LeafDist = FOREST_SEED
      !
      ! Finished with internal tree structure maintenance
      ! Allocate space for NodeInfo
      ALLOCATE(BelowSeedNode%Info,STAT=iError)
      IF (iError /= err_OK) THEN
         PRINT *,"Error allocating BelowSeedNode Info."
         STOP
      END IF
      !
      IF(BelowSeedLevel>MinLevel) THEN
         BelowSeedNode%ChildNo = 1
         CurrentNode => BelowSeedNode; CurrentLevel = BelowSeedLevel
         CALL CreateBelowSeedLevelNode(MinLevel,InfoInit)
      END IF
      !
      IF (.NOT. PRESENT(InfoInit)) RETURN
      IF (.NOT. InfoInit) RETURN
   END SUBROUTINE CreateBelowSeedLevelNode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: KillForest
   ! 作用: 释放整棵树（从 Root 起）及其所有子孙节点，并释放其 info。
   ! 注意: 调用后 Root 及其结构被释放，ForestSeed 自身不在此处释放。
   SUBROUTINE KillForest
      CALL KillNode(Root)
   END SUBROUTINE KillForest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: AddRootLevelNode
   ! 作用: 在根层（rootlevel）添加一个新的“最年轻”根节点，并更新 Root 指向该新节点；
   !       同时维护该层的 Neighbor 链与 Youngest/Eldest/NrOfNodes 统计。
   ! 参数: InfoInit [逻辑, 可选] 若为真，仅分配新节点的 info。
   ! 注意: 新节点将成为新的 Root，原 Root 成为其 Sibling/Neighbor。
   SUBROUTINE AddRootLevelNode(InfoInit)
      LOGICAL, INTENT(IN), OPTIONAL :: InfoInit
      ! Internal declarations
      TYPE (Node), POINTER :: NewRootLevelNode
      INTEGER :: iError
      !
      ALLOCATE (NewRootLevelNode,STAT=iError)
      IF (iError /= 0) THEN
         PRINT *,"Error allocating node in AddRootLevelNode."
         STOP
      END IF
      NewRootLevelNode%level=rootlevel
      NrOfNodes(rootlevel)=NrOfNodes(rootlevel)+1
      globalnodeindex=globalnodeindex+1; NewRootLevelNode%globalnodenum=globalnodeindex
      NULLIFY(NewRootLevelNode%Parent); NULLIFY(NewRootLevelNode%Child)
      NewRootLevelNode%NrOfChildren=NO_CHILDREN
      NewRootLevelNode%ChildNo=Root%ChildNo+1
      NewRootLevelNode%LeafDist=ZERO_LEAF_NODE_DIST
      ! The newly created node becomes the first node on this level and
      ! the start of the tree
      !  1. Set the Sibling and Neighbor of the newly created node to
      !     point to the old Root node
      NewRootLevelNode%Sibling => Root; NewRootLevelNode%Neighbor => Root
      !  2. Set the global Root to point to the newly created root level node
      Root => NewRootLevelNode; ForestSeed%Child=>Root; CurrentNode=>Root
      !  3. The newly created node starts this level's neighbor list
      YoungestOnLevel(rootlevel)%NodePtr=>NewRootLevelNode
      ALLOCATE(NewRootLevelNode%Info,STAT=iError)  ! Allocate space for NodeInfo
      IF (iError /= 0) THEN
         PRINT *,"Error allocating Info in AddRootLevelNode."
         STOP
      END IF
      IF (.NOT. PRESENT(InfoInit)) RETURN
      IF (.NOT. InfoInit) RETURN
   END SUBROUTINE AddRootLevelNode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: CreateChild
   ! 作用: 为 CurrentNode 创建一个子节点（作为该父亲最年轻的孩子），
   !       并维护同层的 Youngest/Eldest/NrOfNodes 与 Neighbor 链。
   ! 参数: InfoInit [逻辑, 可选] 若为真，仅分配 info；
   !       ReadMode [逻辑, 可选] 若为真，表示“从文件读入还原树”模式，此时新节点追加到该层尾部
   !       （EldestOnLevel%Neighbor 指向新节点）；否则作为该层新的最年轻节点挂在链首。
   ! 错误: 超过最大深度 MaxDepth 会 STOP；分配失败也会 STOP。
   SUBROUTINE CreateChild(InfoInit,ReadMode)
      ! Interface declarations
      LOGICAL, INTENT(IN), OPTIONAL :: InfoInit
      LOGICAL, INTENT(IN), OPTIONAL :: ReadMode
      ! Internal declarations
      LOGICAL UpdateYoungest
      TYPE (Node), POINTER :: NewNode
      INTEGER :: iError,ThisLevel
      !
      ALLOCATE(NewNode,STAT=iError)
      IF (iError /= err_OK) THEN
         PRINT *,"Error allocating NewNode in treeops::CreateChild."
         STOP
      END IF
      NULLIFY(NewNode%Parent); NULLIFY(NewNode%Sibling)
      NULLIFY(NewNode%Child);  NULLIFY(NewNode%Neighbor)
      globalnodeindex=globalnodeindex+1; NewNode%globalnodenum=globalnodeindex
      ThisLevel = CurrentNode%level+1   ! Child is one level down
      IF (CurrentLevel /= CurrentNode%level) THEN
         WRITE(1,*)'Internal inconsistency in level counters in treeops::CreateChild'
         WRITE(1,1000)CurrentLevel,CurrentNode%level
1000     FORMAT('Global level=',i2,' CurrentNode%level=',i2)
         STOP
      END IF
      IF (ThisLevel > MaxDepth) THEN
         PRINT *,'Error in treeops::CreateChild. Maximum tree depth exceeded'
         STOP
      END IF
      NrOfNodes(ThisLevel)=NrOfNodes(ThisLevel)+1
      NewNode%Parent => CurrentNode         ! Child is born of Current node
      CurrentNode%NrOfChildren=CurrentNode%NrOfChildren+1
      NewNode%Sibling => CurrentNode%Child  ! Point to next elder sibling
      ! First node of this parent?
      IF (.NOT. ASSOCIATED(NewNode%Sibling)) THEN
         NewNode%ChildNo=FIRST_CHILD
      ELSE
         NewNode%ChildNo=NewNode%Sibling%ChildNo+1
      END IF
      ! Determine if we're reading nodes from a file, in which case
      ! we'll be adding elements to the end of the level list
      IF (.NOT. PRESENT(ReadMode)) THEN
         UpdateYoungest = .TRUE.
      ELSE
         IF (ReadMode) THEN
            UpdateYoungest = .FALSE.
         ELSE
            UpdateYoungest = .TRUE.
         END IF
      END IF
      ! Is this the first child created on this level?
      IF (NrOfNodes(ThisLevel)==1) THEN
         ! Yes. Start stack spanning this level
         YoungestOnLevel(ThisLevel)%NodePtr => NewNode  ! Start of level stack
         EldestOnLevel(ThisLevel)%NodePtr => NewNode    ! End of level stack
         NULLIFY(NewNode%Neighbor)   ! ApplyOnLevel will end on this node
      ELSE
         ! No. Update the Neighbor list spanning this level
         IF (UpdateYoungest) THEN
            ! We're creating a new node during program execution
            ! Set the new node's Neighbor pointer to the previous YoungestOnLevel
            NewNode%Neighbor => YoungestOnLevel(thisLevel)%NodePtr
         ELSE
            ! We're creating a new node while reading from file
            EldestOnLevel(ThisLevel)%NodePtr%Neighbor => NewNode
            NULLIFY(NewNode%Neighbor)
            EldestOnLevel(ThisLevel)%NodePtr => NewNode
         END IF
      END IF
      CurrentNode%Child => NewNode           ! Parent points to youngest child
      NULLIFY(NewNode%Child)                 ! Just born doesn't have children
      NewNode%NrOfChildren=0
      NewNode%level = ThisLevel              ! Set the child's level counter
      NewNode%LeafDist = ZERO_LEAF_NODE_DIST ! New node is a leaf
      ! Save the most recently created child on this level. Used in maintaining
      ! the Neighbor list spanning a level
      IF (UpdateYoungest) YoungestOnLevel(ThisLevel)%NodePtr => NewNode
      ! Finished with internal tree structure maintenance
      ! Allocate space for NodeInfo
      ALLOCATE(NewNode%Info,STAT=iError)
      IF (iError /= err_OK) THEN
         PRINT *,"Error allocating NewNode Info."
         STOP
      END IF
      IF (.NOT. PRESENT(InfoInit)) RETURN
      IF (.NOT. InfoInit) RETURN
   END SUBROUTINE CreateChild
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! Kill aNode and all of its children
   ! ------------------------------------------------------------------------
   ! 子程序: KillNode (递归)
   ! 作用: 递归删除指定节点及其全部子孙：
   !   1) 先深度优先删除所有孩子；
   !   2) 从父亲子链与同层 Neighbor 链中移除该节点；
   !   3) 释放其 info 与节点本身内存。
   ! 参数: aNode [指针] 要删除的节点（可为 Root 或任意子树根）。
   ! 注意: 会同步维护 NrOfNodes/YoungestOnLevel/EldestOnLevel 与链表一致性。
   RECURSIVE SUBROUTINE KillNode(aNode)
      ! Interface declarations
      TYPE (Node), POINTER :: aNode
      ! Internal declarations
      TYPE (Node), POINTER :: Child,Sibling,Prev
      LOGICAL NeighborUpdated
      INTEGER :: iError,ThisLevel
      ! If aNode has children kill those first
      Child => aNode%Child
      DO
         IF (.NOT. ASSOCIATED(Child)) EXIT
         Sibling => Child%Sibling
         Call KillNode(Child)
         Child => Sibling
      END DO
      ! Remove leaf node
      IF (ASSOCIATED(aNode%Parent)) THEN
         aNode%Parent%NrOfChildren = aNode%Parent%NrOfChildren - 1
         IF (ASSOCIATED(aNode%Parent%Child,aNode)) THEN
            ! Update parent child pointer to next eldest
            aNode%Parent%Child => aNode%Sibling
         ELSE
            ! Search for position of this node withing sibling list
            Sibling => aNode%Parent%Child
            DO
               IF (ASSOCIATED(Sibling%Sibling,aNode)) EXIT
               Sibling => Sibling%Sibling
            END DO
            ! Remove this node from the sibling list
            Sibling%Sibling => aNode%Sibling
         END IF
      END IF
      ThisLevel = aNode%level
      NrOfNodes(ThisLevel)=NrOfNodes(ThisLevel)-1
      ! Update the Neighbor list spanning this level
      NeighborUpdated=.FALSE.
      ! Was this node the last node on this level?
      IF (NrOfNodes(ThisLevel)==0) THEN
         NULLIFY(YoungestOnLevel(ThisLevel)%NodePtr)
         NULLIFY(EldestOnLevel(ThisLevel)%NodePtr)
         NeighborUpdated=.TRUE.
      END IF
      ! Was it the start of the Neighbor list?
      IF ((.NOT. NeighborUpdated) .AND. &
         (ASSOCIATED(aNode,YoungestOnLevel(ThisLevel)%NodePtr))) THEN
         ! Set the start of the list to the next node on this level
         YoungestOnLevel(ThisLevel)%NodePtr => &
            YoungestOnLevel(ThisLevel)%NodePtr%Neighbor
         NeighborUpdated=.TRUE.
      END IF
      IF (.NOT. NeighborUpdated) THEN
         ! aNode definitely has a previous element. Find it.
         Prev => YoungestOnLevel(ThisLevel)%NodePtr
         DO WHILE (.NOT. ASSOCIATED(Prev%Neighbor,aNode))
            Prev => Prev%Neighbor
         END DO
      END IF
      ! Was it the end of the Neighbor stack?
      IF ((.NOT. NeighborUpdated) .AND. &
         (ASSOCIATED(aNode,EldestOnLevel(ThisLevel)%NodePtr))) THEN
         NULLIFY(Prev%Neighbor)
         EldestOnLevel(ThisLevel)%NodePtr => Prev
         NeighborUpdated = .TRUE.
      END IF
      IF (.NOT. NeighborUpdated) THEN
         ! If we're here, aNode is neither the last nor the first in Neighbor list
         Prev%Neighbor => aNode%Neighbor ! Skip aNode in Neighbor list
      END IF
      DEALLOCATE(aNode%Info,STAT=iError)
      IF (iError /= 0) THEN
         PRINT *,"Error deallocating aNode Info."
         STOP
      END IF
      DEALLOCATE(aNode,STAT=iError)
      IF (iError /= 0) THEN
         PRINT *,"Error deallocating aNode."
         STOP
      END IF
   END SUBROUTINE KillNode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 函数: DeleteMarkedNode
   ! 作用: 若当前节点的 info.tobedeleted = .TRUE.，则删除 CurrentNode。
   ! 签名: INTEGER FUNCTION(info, param)
   ! 参数: info [nodeinfo, 值传]，dummy [funcparam]
   ! 返回: err_ok
   INTEGER FUNCTION DeleteMarkedNode(info,dummy)
      IMPLICIT NONE
!
      TYPE(nodeinfo):: info
      TYPE(funcparam):: dummy
!
      DeleteMarkedNode = err_ok
!
      IF(.NOT. info%tobedeleted) RETURN
!
      CALL KillNode(currentnode)
!
   END FUNCTION DeleteMarkedNode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: CurrentNodeToYoungest
   ! 作用: 将 CurrentNode/CurrentLevel 设置为给定层 level 上的最年轻节点。
   SUBROUTINE CurrentNodeToYoungest(level)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: level
!
      currentnode => youngestonlevel(level)%nodeptr
      currentlevel = level
!
   END SUBROUTINE CurrentNodeToYoungest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: SetLevelNodeNumbers
   ! 作用: 为指定层的所有节点设置按 Neighbor 链顺序递增的 levelnodenum。
   ! 注意: 要求 level >= rootlevel；不会创建节点，仅写入该层所有节点的 levelnodenum。
   SUBROUTINE SetLevelNodeNumbers(level)
      IMPLICIT NONE
!
      INTEGER:: level
!
      TYPE(node), POINTER:: nextnode
!
      IF(level < rootlevel) THEN
         PRINT *,'Error in SetLevelNodeNumbers'
         STOP
      END IF
!
      levelnodeindex = 0
!
      currentnode => youngestonlevel(level)%nodeptr
!
      DO
         IF(.NOT. ASSOCIATED(currentnode)) EXIT
         nextnode => currentnode%neighbor
         levelnodeindex = levelnodeindex+1
         currentnode%levelnodenum = levelnodeindex
         currentnode => nextnode
      END DO
!
   END SUBROUTINE SetLevelNodeNumbers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: ApplyOnLevel
   ! 作用: 对给定层 level 上的所有节点逐一调用回调 FncN(info, param)。
   ! 参数: level [整型] 目标层（|level| 不得超过 MaxDepth）；
   !       FncN  [接口] 用户提供的过程，签名: INTEGER FUNCTION(info, param)
   !       fparam[funcparam] 传递的用户参数。
   ! 错误: FncN 返回非 err_ok 时将写入 ErrCode 并返回；最大深度越界会 STOP。
   SUBROUTINE ApplyOnLevel(level,FncN,fparam)
      IMPLICIT NONE
!
      INTEGER:: level
      INTERFACE
         INTEGER FUNCTION FncN(info,param)
            USE NodeInfoDef
            TYPE(nodeinfo):: info
            TYPE(funcparam):: param
         END FUNCTION FncN
      END INTERFACE
      TYPE(funcparam):: fparam
!
      TYPE(node), POINTER:: nextnode
      INTEGER:: ferrcode
!
      IF(ABS(level)>maxdepth) THEN
         PRINT *,'Error in ApplyOnLevel. Maximum tree depth exceeded.'
         STOP
      END IF
!
      currentnode => youngestonlevel(level)%nodeptr
      currentlevel = level
!
      DO
         IF(.NOT. ASSOCIATED(currentnode)) EXIT
         nextnode => currentnode%neighbor
!
         ferrcode = FncN(currentnode%info,fparam)
!
         IF(ferrcode /= err_ok) THEN
            errcode = ferrcode
            IF(ferrcode<errprint) PRINT *, 'Error in ApplyOnLevel. ferrcode=', &
               ferrcode
            RETURN
         END IF
         currentnode => nextnode
      END DO
!
   END SUBROUTINE ApplyOnLevel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: ApplyOnLevelPairs
   ! 作用: 对同一层 level 上的所有“有序成对组合 (node1, node2, node2 在 node1 之后)”
   !       调用回调 FncN(info1, info2, param)。常用于同层相互作用。
   ! 注:   自反对（node1 与自身）被跳过，成对数量为 n(n-1)/2。
   SUBROUTINE ApplyOnLevelPairs(level,FncN,fparam)
      IMPLICIT NONE
!
      INTEGER:: level
      INTERFACE
         INTEGER FUNCTION FncN(info1,info2,param)
            USE NodeInfoDef
            TYPE(nodeinfo):: info1, info2
            TYPE(funcparam):: param
         END FUNCTION FncN
      END INTERFACE
      TYPE(funcparam):: fparam
!
      TYPE(node), POINTER:: node1, node2
      INTEGER:: ferrcode
!
      IF(ABS(level) > MaxDepth) THEN
         PRINT *, 'Error in ApplyOnLevelPairs. Maximum tree depth exceeded.'
         STOP
      END IF
!
      node1 => youngestonlevel(level)%nodeptr
!
      DO
         IF(.NOT. ASSOCIATED(node1)) EXIT
!
         node2 => node1%neighbor
         DO
            IF(.NOT. ASSOCIATED(node2)) EXIT
!
            ferrcode = FncN(node1%info,node2%info,fparam)
!
            IF(ferrcode /= err_ok) THEN
               errcode = ferrcode
               IF(ferrcode < errprint) THEN
                  PRINT *, 'Error in ApplyOnLevelPairs. ferrcode=', ferrcode
               END IF
               RETURN
            END IF
!
            node2 => node2%neighbor
         END DO
         node1 => node1%neighbor
      END DO
!
   END SUBROUTINE ApplyOnLevelPairs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: ApplyOnMyLevelNbrs
   ! 作用: 以 CurrentNode 为“左操作数”，与同层所有其他节点（排除自身）组合，
   !       调用回调 FncN(currentnode%info, nbr%info, param)。
   SUBROUTINE ApplyOnMyLevelNbrs(level,FncN,fparam)
      IMPLICIT NONE
!
! Enables the preset currentnode to interact with all the other nodes on level:
!
      INTEGER:: level
      INTERFACE
         INTEGER FUNCTION FncN(info1,info2,param)
            USE NodeInfoDef
            TYPE(nodeinfo):: info1, info2
            TYPE(funcparam):: param
         END FUNCTION FncN
      END INTERFACE
      TYPE(funcparam):: fparam
!
      TYPE(node), POINTER:: nbrnode
      INTEGER:: ferrcode
!
      IF(ABS(level) > MaxDepth) THEN
         PRINT *, 'Error in ApplyOnLevelPairs. Maximum tree depth exceeded.'
         STOP
      END IF
!
      nbrnode => youngestonlevel(level)%nodeptr
!
      DO
         IF(.NOT. ASSOCIATED(nbrnode)) EXIT
!
! Avoid self interactions:
         IF(currentnode%levelnodenum /= nbrnode%levelnodenum) THEN
!
            ferrcode = FncN(currentnode%info,nbrnode%info,fparam)
!
            IF(ferrcode /= err_ok) THEN
               errcode = ferrcode
               IF(ferrcode < errprint) THEN
                  PRINT *, 'Error in ApplyOnLevelPairs. ferrcode=', ferrcode
               END IF
               RETURN
            END IF
         END IF
!
         nbrnode => nbrnode%neighbor
      END DO
!
   END SUBROUTINE ApplyOnMyLevelNbrs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: ApplyOnChildren
   ! 作用: 针对 CurrentNode 的所有孩子，逐一调用 FncN(child%info, param)。
   SUBROUTINE ApplyOnChildren(FncN,fparam)
      IMPLICIT NONE
!
      INTERFACE
         INTEGER FUNCTION FncN(info,param)
            USE NodeInfoDef
            TYPE(nodeinfo):: info
            TYPE(funcparam):: param
         END FUNCTION FncN
      END INTERFACE
      TYPE(funcparam):: fparam
!
      INTEGER ferrcode
      TYPE(node), POINTER :: savecurrentnode
!
      savecurrentnode => currentnode
      currentnode => currentnode%child
!
      DO
         IF(.NOT. ASSOCIATED(currentnode)) EXIT
!
         ferrcode = FncN(currentnode%info,fparam)
!
         IF(ferrcode /= err_ok) THEN
            errcode = ferrcode
            IF(ferrcode < errprint) PRINT *, 'Error in ApplyOnChildren. ferrcode=', &
               ferrcode
            RETURN
         END IF
!
         currentnode => currentnode%sibling
      END DO
!
      currentnode => savecurrentnode
!
   END SUBROUTINE ApplyOnChildren
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: ApplyOnForest
   ! 作用: 对森林中的所有节点（逐棵树，从 Root 开始）执行树遍历并在满足 TrueCond() 时
   !       调用用户回调 FncN(info, param)。
   ! 参数: PreEval [逻辑, 可选] 若 .TRUE.（默认），表示先确定 Next 再评估当前（子先/前序变体）。
   SUBROUTINE ApplyOnForest(FncN,fparam,PreEval)
      IMPLICIT NONE
!
      INTERFACE
         INTEGER FUNCTION FncN(info,param)
            USE NodeInfoDef
            TYPE(nodeinfo):: info
            TYPE(funcparam):: param
         END FUNCTION FncN
      END INTERFACE
      TYPE(funcparam):: fparam
      LOGICAL, OPTIONAL:: preeval
!
      LOGICAL evalnextbefore_f
!
      currentnode => root
      currentlevel = rootlevel
!
      IF(PRESENT(preeval)) THEN
         evalnextbefore_f = preeval
      ELSE
         evalnextbefore_f = .TRUE.
      END IF
!
      CALL ForestTraversal(root,FncN,truecond,fparam,evalnextbefore_f)
!
   END SUBROUTINE ApplyOnForest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! ------------------------------------------------------------------------
   ! 子程序: ApplyOnLevels
   ! 作用: 对所有可遍历层执行树遍历，并在 LevelCond() 恒为真时调用回调 f(Info, Param)。
   ! 参数: PreEval 同 ApplyOnForest。
   SUBROUTINE ApplyOnLevels(f,fparam,PreEval)
      ! Level  0: Root
      ! Level  1: 1st generation children
      ! ....
      ! Level  n: n-th generation children
      ! Level -1: Parents of leaves
      ! Level -2: Grandparents of leaves
      INTERFACE
         INTEGER FUNCTION f(Info,Param)
            USE NodeInfoDef
            ! Interface declarations
            TYPE (NodeInfo) :: Info
            TYPE (FuncParam) :: Param
         END FUNCTION f
      END INTERFACE
      TYPE (FuncParam) :: fparam
      LOGICAL, OPTIONAL :: PreEval
      LOGICAL EvalNextBefore_f
      CurrentNode => Root; CurrentLevel=rootlevel
      IF (PRESENT(PreEval)) THEN
         EvalNextBefore_f=PreEval
      ELSE
         EvalNextBefore_f=.TRUE.
      END IF
      CALL ForestTraversal(Root,f,LevelCond,fparam,EvalNextBefore_f)
   END SUBROUTINE ApplyOnLevels

   ! ------------------------------------------------------------------------
   ! 子程序: ApplyOnLeaves
   ! 作用: 仅对叶子结点（无子节点）执行回调 f(Info, Param)。
   ! 参数: PreEval 同 ApplyOnForest。
   SUBROUTINE ApplyOnLeaves(f,fparam,PreEval)
      ! Leaves = youngest generation in existence
      INTERFACE
         INTEGER FUNCTION f(Info,Param)
            USE NodeInfoDef
            ! Interface declarations
            TYPE (NodeInfo) :: Info
            TYPE (FuncParam) :: Param
         END FUNCTION f
      END INTERFACE
      TYPE (FuncParam) :: fparam
      LOGICAL, OPTIONAL :: PreEval
      LOGICAL EvalNextBefore_f
      ! First executable statement
      CurrentNode => Root; CurrentLevel=rootlevel
      IF (PRESENT(PreEval)) THEN
         EvalNextBefore_f=PreEval
      ELSE
         EvalNextBefore_f=.TRUE.
      END IF
      CALL ForestTraversal(Root,f,LeafCond,fparam,EvalNextBefore_f)
   END SUBROUTINE ApplyOnLeaves

   ! Routines interior to the module

   ! ------------------------------------------------------------------------
   ! 子程序: ForestTraversal
   ! 作用: 遍历同层的根节点（可能有多棵树），对每棵树调用 TreeTraversal。
   ! 参数: aNode = 根层起点；f/cond = 回调与条件；EvalNextBefore_f 控制先/后评估策略。
   SUBROUTINE ForestTraversal(aNode,f,cond,fparam,EvalNextBefore_f)
      ! Traverse the root level nodes
      TYPE (Node), POINTER :: aNode
      INTERFACE
         INTEGER FUNCTION f(Info,Param)
            USE NodeInfoDef
            ! Interface declarations
            TYPE (NodeInfo) :: Info
            TYPE (FuncParam) :: Param
         END FUNCTION f
         LOGICAL FUNCTION cond()
         END FUNCTION cond
      END INTERFACE
      TYPE (FuncParam) :: fparam
      LOGICAL EvalNextBefore_f
      ! Internal declarations
      TYPE (Node), POINTER :: Next
      Next => aNode
      DO
         IF (.NOT. ASSOCIATED(Next)) EXIT
         CALL TreeTraversal(Next,f,cond,fparam,EvalNextBefore_f)
         Next=>Next%Sibling
      END DO
   END SUBROUTINE ForestTraversal

   ! ------------------------------------------------------------------------
   ! 子程序: TreeTraversal (递归)
   ! 作用: 深度优先遍历：
   !   - 将 CurrentNode 设为 aNode，并根据 EvalNextBefore_f 决定先取 Next(=Child) 还是先回调；
   !   - 若 cond() 为真，调用 f(CurrentNode%Info, fparam)；
   !   - 递归下探到子节点，之后回到兄弟节点；
   !   - 返回时维护 CurrentLevel/CurrentNode 的一致性。
   ! 错误: 若 f 返回非 err_OK，ErrCode 被设置并立刻回退整个递归。
   RECURSIVE SUBROUTINE TreeTraversal(aNode,f,cond,fparam,EvalNextBefore_f)
      ! The core tree function in terms of which all others are defined:
      !   1) Traverse tree
      !   2) If condition *Cond* is satisfied, apply function *Func* to node
      !      *aNode*
      !   3) Update *Level* counter
      TYPE (Node), TARGET, INTENT(IN) :: aNode
      INTERFACE
         INTEGER FUNCTION f(Info,Param)
            USE NodeInfoDef
            ! Interface declarations
            TYPE (NodeInfo) :: Info
            TYPE (FuncParam) :: Param
         END FUNCTION f
         LOGICAL FUNCTION cond()
         END FUNCTION cond
      END INTERFACE
      TYPE (FuncParam) :: fparam
      LOGICAL :: EvalNextBefore_f
      ! Internal declarations
      TYPE (Node), POINTER :: Next
      TYPE (Node), POINTER, SAVE :: aNodeSave
      INTEGER fErrCode,Applyf,level
      ! First executable statement
      IF (ErrCode /= err_OK) RETURN   ! Go up recursion upon error
      CurrentNode => aNode            ! Set current node and maintain stack
      Stack(CurrentLevel)%NodePtr => CurrentNode ! so other functions know what to work on
      IF (EvalNextBefore_f) Next => CurrentNode%Child
      ! Do work on this node if condition cond() is satisfied
      IF (cond()) THEN
         fErrCode=f(CurrentNode%Info,fparam)
         IF (fErrCode /= err_OK) THEN
            ErrCode = fErrCode
            IF (fErrCode<ErrPrint ) PRINT *,'Error on applying f. fErrCode=',fErrCode
            RETURN
         END IF
      END IF
      IF (.NOT. EvalNextBefore_f) Next => CurrentNode%Child
      ! Find next node to work on
      DO
         IF (.NOT. ASSOCIATED(Next)) THEN
            ! Reached a leaf
            CurrentLevel=MAX(CurrentLevel-1,rootlevel); level=CurrentLevel
            ! Go up tree, exit clause from DO loop
            EXIT
         ELSE
            CurrentLevel=CurrentLevel+1; level=CurrentLevel  ! Node has children; go down
            CALL TreeTraversal(Next,f,cond,fparam,EvalNextBefore_f)
            ! Reset global context from local context when returning from recursion
            CurrentNode => aNode
            Next => Next%Sibling       ! After children have been exhausted
            ! work on next sibling
         END IF
      END DO
   END SUBROUTINE TreeTraversal

   ! 简单条件：恒为真（用于 ApplyOnForest）
   LOGICAL FUNCTION TrueCond()
      TrueCond = .TRUE.
   END FUNCTION TrueCond

   ! 简单条件：恒为真（用于 ApplyOnLevels）
   LOGICAL FUNCTION LevelCond()
      LevelCond = .TRUE.
   END FUNCTION LevelCond

   ! 叶子条件：仅在无子节点时返回真（用于 ApplyOnLeaves）
   LOGICAL FUNCTION LeafCond()
      IF (.NOT. ASSOCIATED(CurrentNode%Child)) THEN
         LeafCond = .TRUE.
      ELSE
         LeafCond = .FALSE.
      END IF
   END FUNCTION LeafCond

   ! ------------------------------------------------------------------------
   ! 函数: SetChildNo
   ! 作用: 若当前节点 ChildNo 未就绪（= BAD_CHILD_NO），则为其所在家庭按“从最年长到最年轻”
   !       进行编号一致化设置（保证 ChildNo 与 sibling 链位置一致）。
   INTEGER FUNCTION SetChildNo(Info,Param)
      TYPE (NodeInfo) :: Info
      TYPE (FuncParam) :: Param
      !
      TYPE (Node), POINTER :: Youngest
      INTEGER ListNo,NrOfChildren
      IF (CurrentNode%Level==rootlevel) THEN
         CurrentNode%ChildNo=FIRST_CHILD
         SetChildNo=err_OK
         RETURN
      END IF
      IF (CurrentNode%ChildNo==BAD_CHILD_NO) THEN
         ! Number all children within this family
         Youngest => CurrentNode%Parent%Child
         ListNo=1
         NrOfChildren=CurrentNode%Parent%NrOfChildren
         Youngest%ChildNo=NrOfChildren
         DO
            IF (.NOT. ASSOCIATED(Youngest%Sibling)) EXIT
            Youngest => Youngest%Sibling
            ListNo=ListNo+1
            Youngest%ChildNo=NrOfChildren-ListNo+1
         END DO
      END IF
      SetChildNo=err_OK
   END FUNCTION SetChildNo

   ! 返回最近一次的模块错误码 ErrCode
   INTEGER FUNCTION GetTreeOpsErrCode()
      GetTreeOpsErrCode = ErrCode
   END FUNCTION GetTreeOpsErrCode

   ! 读取指定节点的 info（值拷贝），若节点无效返回 err_UndefinedNode
   INTEGER FUNCTION GetNodeInfo(aNode,aNodeInfo)
      TYPE (Node), POINTER :: aNode
      TYPE (NodeInfo) :: aNodeInfo
      IF (ASSOCIATED(aNode)) THEN
         aNodeInfo = aNode%Info
         GetNodeInfo = err_OK
      ELSE
         GetNodeInfo = err_UndefinedNode
      END IF
   END FUNCTION GetNodeInfo

   ! 设置指定节点的 info（值拷贝），若节点无效返回 err_UndefinedNode
   INTEGER FUNCTION SetNodeInfo(aNode,aNodeInfo)
      TYPE (Node), POINTER :: aNode
      TYPE (NodeInfo) :: aNodeInfo
      IF (ASSOCIATED(aNode)) THEN
         aNode%Info = aNodeInfo
         SetNodeInfo = err_OK
      ELSE
         SetNodeInfo = err_UndefinedNode
      END IF
   END FUNCTION SetNodeInfo

   ! 读取 CurrentNode 的 info（值拷贝）；无当前节点则 err_UndefinedNode
   INTEGER FUNCTION GetCurrentNodeInfo(aNodeInfo)
      TYPE (NodeInfo) :: aNodeInfo
      IF (ASSOCIATED(CurrentNode)) THEN
         aNodeInfo = CurrentNode%Info
         GetCurrentNodeInfo = err_OK
      ELSE
         GetCurrentNodeInfo = err_UndefinedNode
      END IF
   END FUNCTION GetCurrentNodeInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! 获取 CurrentNode 的全局编号 globalnodenum
   INTEGER FUNCTION Getglobalnodenum( )
      IMPLICIT NONE
!
      Getglobalnodenum = currentnode%globalnodenum
!
   END FUNCTION Getglobalnodenum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! 设置 CurrentNode 的 info（值拷贝）
   INTEGER FUNCTION SetCurrentNodeInfo(aNodeInfo)
      IMPLICIT NONE
!
      TYPE(nodeinfo):: anodeinfo
!
      IF(ASSOCIATED(currentnode)) THEN
         currentnode%info = anodeinfo
         SetCurrentNodeInfo = err_ok
      ELSE
         SetCurrentNodeInfo = err_undefinednode
      END IF
!
   END FUNCTION SetCurrentNodeInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! 获取 Root 的 info（指针引用返回）
   INTEGER FUNCTION GetRootInfo(anodeinfo)
      IMPLICIT NONE
!
      TYPE(nodeinfo), POINTER:: anodeinfo
!
      IF(ASSOCIATED(root)) THEN
         anodeinfo => root%info
         GetRootInfo = err_ok
      ELSE
         GetRootInfo = err_undefinednode
      END IF
!
   END FUNCTION GetRootInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! 设置 Root 的 info（值拷贝）
   INTEGER FUNCTION SetRootInfo(aNodeInfo)
      TYPE (NodeInfo) :: aNodeInfo
      IF (ASSOCIATED(Root)) THEN
         Root%Info = aNodeInfo
         SetRootInfo = err_OK
      ELSE
         SetRootInfo = err_UndefinedNode
      END IF
   END FUNCTION SetRootInfo

   ! 获取 CurrentNode 父节点的 info（指针引用返回）。若到达 ForestSeed 或无父节点，返回 err_NoParent。
   INTEGER FUNCTION GetParentInfo(aNodeInfo)
      TYPE (NodeInfo), POINTER :: aNodeInfo
      IF (ASSOCIATED(CurrentNode) .AND. ASSOCIATED(CurrentNode%Parent) .AND.  &
         CurrentNode%Parent%level > FOREST_SEED) THEN
         aNodeInfo => CurrentNode%Parent%Info
         GetParentInfo = err_OK
      ELSE
         GetParentInfo = err_NoParent
      END IF
   END FUNCTION GetParentInfo

   ! 获取当前节点 CurrentNode（指针引用返回）
   INTEGER FUNCTION GetCurrentNode(aNode)
      TYPE (Node), POINTER :: aNode
      GetCurrentNode = err_OK
      aNode => CurrentNode
   END FUNCTION GetCurrentNode

   ! 获取 aNode 的父节点 Parent（指针引用）。无父节点则返回 err_NoParent。
   INTEGER FUNCTION GetParent(aNode,Parent)
      TYPE (Node), POINTER :: aNode
      TYPE (Node), POINTER :: Parent
      IF (ASSOCIATED(aNode)) THEN
         IF (ASSOCIATED(aNode%Parent)) THEN
            Parent => aNode%Parent
            GetParent = err_OK
         ELSE
            NULLIFY(Parent)
            GetParent = err_NoParent
         END IF
      ELSE
         GetParent = err_UndefinedNode
      END IF
   END FUNCTION GetParent

   ! 获取 aNode 的兄长 Sibling（指针引用）。无 Sibling 返回 err_NoSibling。
   INTEGER FUNCTION GetSibling(aNode,Sibling)
      TYPE (Node), POINTER :: aNode
      TYPE (Node), POINTER :: Sibling
      IF (ASSOCIATED(aNode)) THEN
         IF (ASSOCIATED(aNode%Sibling)) THEN
            Sibling => aNode%Sibling
            GetSibling = err_OK
         ELSE
            NULLIFY(Sibling)
            GetSibling = err_NoSibling
         END IF
      ELSE
         GetSibling = err_UndefinedNode
      END IF
   END FUNCTION GetSibling

   ! 获取 aNode 的“最年轻孩子” Child（指针引用）。无子返回 err_NoChild。
   INTEGER FUNCTION GetChild(aNode,Child)
      TYPE (Node), POINTER :: aNode
      TYPE (Node), POINTER :: Child
      IF (ASSOCIATED(aNode)) THEN
         IF (ASSOCIATED(aNode%Child)) THEN
            Child => aNode%Child
            GetChild = err_OK
         ELSE
            NULLIFY(Child)
            GetChild = err_NoChild
         END IF
      ELSE
         GetChild = err_UndefinedNode
      END IF
   END FUNCTION GetChild

   ! 获取 CurrentNode 最年轻孩子的 info（指针引用）。无子返回 err_NoChild。
   INTEGER FUNCTION GetChildInfo(aNodeInfo)
      TYPE (NodeInfo), POINTER :: aNodeInfo
      IF (ASSOCIATED(CurrentNode) .AND. ASSOCIATED(CurrentNode%Child)) THEN
         aNodeInfo => CurrentNode%Child%Info
         GetChildInfo = err_OK
      ELSE
         GetChildInfo = err_NoChild
      END IF
   END FUNCTION GetChildInfo

   ! 返回当前层级（额外自检 CurrentLevel 与 CurrentNode%level 的一致性）
   INTEGER FUNCTION GetLevel(ThisLevel)
      INTEGER ThisLevel
      ThisLevel = CurrentLevel
      IF (ThisLevel /= CurrentNode%level) THEN
         PRINT *,'Internal inconsistency in level counter'
         STOP
      END IF
      GetLevel = err_OK
   END FUNCTION GetLevel

   ! 计算 CurrentNode 在“同父兄弟链”中的位置（1=最年长），并交叉验证 ChildNo 一致性。
   INTEGER FUNCTION GetSiblingIndex(SiblingIndex)
      INTEGER SiblingIndex
      INTEGER NrOfSiblings
      TYPE (Node), POINTER :: NextSibling,ThisSibling
      ThisSibling => CurrentNode
      SiblingIndex=1; NrOfSiblings=1
      IF (.NOT. ASSOCIATED(CurrentNode%Parent)) THEN
         ! We are on the forest seed level
         GetSiblingIndex = err_OK
         RETURN
      ELSE
         ! We are below the forest seed level
         NextSibling => CurrentNode%Parent%Child
         DO
            IF (ASSOCIATED(NextSibling,ThisSibling)) EXIT
            SiblingIndex=SiblingIndex+1
            NextSibling => NextSibling%Sibling
         END DO
         GetSiblingIndex = err_OK
         NextSibling => CurrentNode%Parent%Child
         DO
            IF (.NOT. ASSOCIATED(nextSibling%Sibling)) EXIT
            NrOfSiblings=NrOfSiblings+1
            NextSibling => NextSibling%Sibling
         END DO
      END IF
      IF ((NrOfSiblings+1-SiblingIndex) /= CurrentNode%ChildNo) THEN
         PRINT *,'Internal inconsistency in ChildNo counter'
         PRINT *,'Level=',CurrentNode%level
         PRINT *,'ChildNo=',CurrentNode%ChildNo
         PRINT *,'SiblingIndex=',SiblingIndex
         PRINT *,'NrOfSiblings=',NrOfSiblings
         STOP
      END IF
   END FUNCTION GetSiblingIndex

   ! 计算 CurrentNode 的孩子数量（遍历 sibling 链），并与存储的 nrofchildren 一致性检查。
   INTEGER FUNCTION GetNrOfChildren(NrOfChildren)
      INTEGER NrOfChildren
      TYPE (Node), POINTER :: NextChild
      NrOfChildren=0
      NextChild => CurrentNode%Child
      DO
         IF (.NOT. ASSOCIATED(NextChild)) EXIT
         NrOfChildren=NrOfChildren+1
         NextChild => NextChild%Sibling
      END DO
      IF (NrOfChildren /= CurrentNode%NrOfChildren) THEN
         WRITE(1,1001)
1001     FORMAT('GetNrOfChildren: Internal inconsistency in NrOfChildren')
         WRITE(1,1002)NrOfChildren,CurrentNode%NrOfChildren
1002     FORMAT('Computed from links:',i3,' stored in Node:',i3)
      END IF
      GetNrOfChildren = err_OK
   END FUNCTION GetNrOfChildren
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! 查询某层是否存在任何节点（通过 EldestOnLevel(level) 是否关联判断）
   LOGICAL FUNCTION ExistLevel(level)
      IMPLICIT NONE
!
      INTEGER:: level
!
      IF(ABS(level) > maxdepth) THEN
         errcode = err_badlevel
         PRINT *, 'Bad level in call to ExistLevel.'
         STOP
      END IF
!
      IF(ASSOCIATED(eldestonlevel(level)%nodeptr)) THEN
         ExistLevel = .TRUE.
      ELSE
         ExistLevel = .FALSE.
      END IF
!
   END FUNCTION ExistLevel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: PushForest
   ! 作用: 将当前遍历的 CurrentNode 保存到 foreststack，以便在子遍历完成后恢复。
   ! 注意: stacklevel 计数上限 MaxForestStacks。
   SUBROUTINE PushForest
      IMPLICIT NONE
!
! Save current tree traversal state to allow a subsidiary tree traversal to
! take place:
!
      foreststack(stacklevel)%nodeptr => currentnode
      stacklevel = stacklevel+1
!
      IF(stacklevel > maxforeststacks) THEN
         PRINT *, 'Too many forest traversal recursions.'
         PRINT *, 'Increase MaxForestStacks in treeops.f90'
         STOP
      END IF
!
      currentnode => root
!
   END SUBROUTINE PushForest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   ! ------------------------------------------------------------------------
   ! 子程序: PopForest
   ! 作用: 从 foreststack 中恢复先前保存的 CurrentNode。
   ! 错误: 若无对应的 PushForest 即调用，则 STOP。
   SUBROUTINE PopForest
      IMPLICIT NONE
!
! Restore tree traversal state on return from a subsidiary tree traversal:
!
      stacklevel = stacklevel-1
      IF(stacklevel < 0) THEN
         PRINT *,'PopForest requested without prior PushForest'
         STOP
      END IF
!
      currentnode => foreststack(stacklevel)%nodeptr
!
   END SUBROUTINE PopForest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
END MODULE TreeOps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

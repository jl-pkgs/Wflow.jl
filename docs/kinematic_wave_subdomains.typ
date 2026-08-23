#set document(title: "汇流过程与运动波子流域划分原理", author: "Wflow / RiverNetwork")
#set page(
  paper: "a4",
  margin: (x: 24mm, y: 22mm),
  numbering: "1",
  number-align: center,
)
#set text(font: ("Noto Serif CJK SC", "Noto Serif CJK JP"), size: 10.5pt, lang: "zh")
#set par(justify: true, leading: 0.72em, first-line-indent: 2em)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)")
#show heading: set text(font: "Noto Sans CJK SC", weight: "bold")
#show raw: set text(font: "DejaVu Sans Mono", size: 8.5pt)

#align(center)[
  #text(font: "Noto Sans CJK SC", size: 19pt, weight: "bold")[
    汇流过程与运动波子流域划分原理
  ]

  #v(0.5em)
  #text(size: 11pt)[——以 Wflow / RiverNetwork 为例]
  #v(1.2em)
]

#block(
  fill: rgb("f2f6f8"),
  stroke: 0.6pt + rgb("9aadb8"),
  radius: 3pt,
  inset: 10pt,
)[
  *摘要*　运动波汇流的核心是：按河网拓扑顺序，将每个网格上游出流之和作为其入流，再求解质量守恒方程与曼宁型流量—过水断面关系。所谓“运动波子流域划分”并不改变水动力学方程，而是将有向无环河网切分为若干内部串行、彼此分层并行的计算单元。RiverNetwork 以 Strahler 河流级别跃迁点为分界，以子流域图到总出口的距离构造并行层，从而同时满足汇流依赖关系与多线程并行要求。
]

= 问题本质

栅格汇流同时包含两个层次：

- *物理层*：水量沿坡面或河道向下游传播，并在汇流点叠加；
- *计算层*：下游单元依赖上游单元的本时段出流，因此必须满足“上游先算、下游后算”。

RiverNetwork 负责第二层中的河网拓扑及任务划分，Wflow 的运动波模块负责第一层中的流量更新。二者通过拓扑序连接：河网结构决定计算顺序，但不改变连续方程、曼宁关系和时间离散格式。

需要特别指出：这里的“子流域”首先是*并行调度单元*。它由河流级别阈值和拓扑依赖确定，不等同于以任意水文站或行政断面划定的传统子流域。

= 从流向栅格到汇流图

== 有向无环图

将每个有效栅格记为节点 $v in V$。D8 流向给出唯一的下游节点 $d(v)$，据此建立有向边

$ v arrow d(v). $

于是得到有向图 $G=(V,E)$。一个节点最多有一个下游邻居，却可有多个上游邻居。记节点 $v$ 的直接上游集合为

$ cal(U)(v) = {u in V: (u,v) in E}. $

正常流向数据应满足：除出口或洼地节点外，每个节点恰有一个下游节点；图中不存在有向环。若存在环，则无法定义严格的上游—下游计算顺序，运动波汇流也无法按单向依赖推进。

== 拓扑排序

对 $G$ 进行拓扑排序，得到序列

$ pi = (v_1, v_2, dots, v_N), $

并保证对任意边 $(v_i,v_j) in E$ 均有 $i<j$。因此，遍历到节点 $v_j$ 时，其所有直接上游节点已完成本时段更新。该性质是汇流守恒和并行划分的共同基础。

#figure(
  table(
    columns: (1fr, auto, 1fr, auto, 1fr, auto, 1fr),
    align: center,
    inset: 7pt,
    stroke: 0.6pt + rgb("9aadb8"),
    fill: (x, y) => if calc.even(x) { rgb("f2f6f8") },
    [D8 流向栅格], [$arrow.r$], [单元级 DAG], [$arrow.r$], [拓扑排序], [$arrow.r$], [运动波更新],
  ),
  caption: [流向数据到汇流计算的基本链条。],
) <fig:pipeline>

= 运动波汇流

== 控制方程

一维连续方程为

$ (partial A)/(partial t) + (partial Q)/(partial x) = q_l, $ <eq:continuity>

其中，$A$ 为过水断面面积（$"m"^2$），$Q$ 为流量（$"m"^3 slash "s"$），$q_l$ 为单位河长或坡长的侧向入流（$"m"^2 slash "s"$）。

完整圣维南方程还包含局部惯性、平流惯性、压力梯度、底坡和摩阻。运动波近似忽略惯性及压力项，并令摩阻坡降等于河床或坡面坡降：

$ S_f approx S_0. $

在该假设下，动量方程退化为流量与断面面积的单值关系。采用曼宁公式并将湿周 $P$ 视为给定量，可写为

$ Q = 1/n A (A/P)^(2/3) sqrt(S_0), $

等价于

$ A = alpha Q^beta, quad beta = 3/5, $

$ alpha = (n/sqrt(S_0))^(3/5) P^(2/5). $ <eq:rating>

这里 $n$ 为曼宁糙率。Wflow 中 `BETA_KINWAVE = 0.6`，与 $beta=3/5$ 一致。

== 单元离散

对长度为 $Delta x_i$ 的单元 $i$，采用隐式时间离散。其本时段上游入流为

$ Q_("in",i)^(t+Delta t) = sum_(j in cal(U)(i)) Q_j^(t+Delta t). $ <eq:confluence>

将 @eq:rating 代入 @eq:continuity，可得 Wflow 实际求解的代数方程：

$ (Delta t)/(Delta x_i) Q_i^(t+Delta t)
  + alpha_i (Q_i^(t+Delta t))^beta
  = (Delta t)/(Delta x_i) Q_("in",i)^(t+Delta t)
  + alpha_i (Q_i^t)^beta
  + Delta t q_(l,i). $ <eq:discrete>

式 @eq:discrete 同时体现三类水量：上游汇入、单元原有蓄水和侧向补给。由于 $Q_("in",i)^(t+Delta t)$ 使用上游节点的本时段新流量，必须严格按拓扑序计算。

Wflow 令 $u=Q^(1/5)$，把 $beta=3/5$ 的非线性方程改写为

$ f(u) = (Delta t)/(Delta x) u^5 + alpha u^3 - C = 0, $

其中

$ C = (Delta t)/(Delta x) Q_("in") + alpha (Q^t)^(3/5) + Delta t q_l. $

随后采用牛顿迭代

$ u^(r+1) = u^r - f(u^r)/f'(u^r) $

求得 $Q^(t+Delta t)=u^5$ 和 $A^(t+Delta t)=alpha u^3$。

== 汇流点的处理

汇流点不需要单独构造经验汇流公式。对节点 $v$，只需将所有直接上游节点的已更新出流相加：

$ Q_("in",v) = sum_(u in cal(U)(v)) Q_u. $

随后用式 @eq:discrete 更新节点 $v$。因此，汇流过程的正确性取决于两个条件：

1. 上游流量只能在上游节点更新完成后读取；
2. 同一节点的全部上游贡献必须完整求和。

坡面汇流还可按 `flow_fraction_to_river` 将上游出流分为进入河道和继续沿坡面传播的两部分；水库、引调水和河漫滩交换则作为边界项或侧向项进入，但不改变拓扑依赖原则。

= Strahler 河流分级

RiverNetwork 先在单元级 DAG 上计算 Strahler 级别。源头节点级别为 1。对非源头节点 $v$，令

$ m = max_(u in cal(U)(v)) omega(u), $

其中 $omega$ 为河流级别。若最高级别 $m$ 在上游集合中至少出现两次，则

$ omega(v)=m+1; $

否则

$ omega(v)=m. $

因此，只有同级河流汇合才会使级别增加；不同级河流汇合时，下游继承较高级别。算法按拓扑序单次遍历，时间复杂度为 $O(|V|+|E|)$。

Strahler 级别在这里有两项作用：一是表征河网层级；二是用级别阈值控制并行子流域的粒度。

= 运动波子流域划分

== 分界节点

给定最小河流级别 $omega_min$，节点 $v$ 被选为子流域出口，当且仅当满足以下条件之一：

- $omega(v) >= omega_min$，且其级别与下游节点不同；
- $v$ 为整个流域出口或内部终点。

也就是说，算法以高于阈值的*级别跃迁点*和最终出口为分界。阈值以下的细小支流不单独建立任务，而是并入其下游子流域。

== 向上游扩展标签

首先只给分界节点赋予唯一子流域编号，其余节点编号为 0。然后按反向拓扑序，即从下游向上游遍历。若当前节点尚未编号，则继承其唯一直接下游节点的编号。由于下游节点总是先被处理，一次遍历即可把每个单元归入最近的下游分界节点。

该步骤可概括为

$ b(v) = b(d(v)), $

其中 $b(v)$ 为子流域编号。最终，每个节点恰属于一个子流域，各子流域覆盖完整计算域且互不重叠。

== 构建子流域级 DAG

将每个子流域压缩为一个节点。若子流域 $a$ 的出口流向子流域 $b$，则建立边 $a arrow b$，得到子流域级 DAG $G_B$。它保留原河网的依赖关系，却显著减少调度节点数量。

对每个总出口，计算各子流域到出口的图距离。距离较大的子流域位于上游，距离较小者位于下游。RiverNetwork 将相同距离的子流域归为同一并行层，并把无上游邻居的短支流提前到最上游层。最终得到

$ cal(L) = (L_K, L_(K-1), dots, L_0), $

其中 $L_K$ 最靠上游，$L_0$ 含总出口。层内子流域之间不存在直接依赖，可并行计算；层与层之间必须按上游到下游依次执行。

#figure(
  table(
    columns: (auto, 1fr, 1fr, 1fr, 1fr),
    align: center,
    inset: 6pt,
    stroke: 0.6pt + rgb("9aadb8"),
    fill: (x, y) => if y == 0 { rgb("dfeaf0") } else if x == 0 { rgb("f2f6f8") },
    [调度阶段], [线程 1], [线程 2], [线程 3], [同步条件],
    [$L_K$], [上游子流域 A], [上游子流域 B], [源头子流域 C], [层内并行],
    [$L_(K-1)$], [中游子流域 D], [中游子流域 E], [—], [等待 $L_K$],
    [$L_0$], [出口子流域 F], [—], [—], [等待全部上游层],
  ),
  caption: [子流域级并行调度：层内并行、层间串行。],
) <fig:schedule>

== 子流域内部与子流域之间

划分后的计算具有两层循环：

```julia
for upstream_level in order_of_subdomains
    threaded_foreach(upstream_level) do subbasin
        for cell in topological_order[subbasin]
            qin[cell] = sum(q[upstream] for upstream in upstream_nodes[cell])
            q[cell], area[cell] = kinematic_wave(...)
        end
    end
end
```

- *子流域内部*：仍按单元拓扑序串行推进，保证上游流量已知；
- *同一并行层内*：不同子流域相互独立，可由多个线程同时推进；
- *并行层之间*：前一上游层全部结束后，才进入下一下游层。

因此，划分只改变任务调度，不改变单元方程及上下游水量传递。其正确性不依赖线程执行次序，而依赖每层结束时的同步屏障。

= 算法流程

RiverNetwork 的实现可归纳为以下步骤：

1. 根据流向构建单元级 DAG，并检查环路和无效下游方向；
2. 对 DAG 进行上游到下游的拓扑排序；
3. 按拓扑序计算 Strahler 河流级别；
4. 依据 $omega_min$ 识别级别跃迁点和流域出口；
5. 从下游向上游传播编号，形成完整子流域分区；
6. 将单元级 DAG 压缩为子流域级 DAG；
7. 按子流域到总出口的距离构造上游到下游的并行层；
8. 每层内并行、子流域内部按拓扑序串行求解运动波方程。

对于多出口计算域，算法先按出口把单元划为互不相交的流域，再在各流域内独立划分子流域；最后统一子流域编号，并合并相同调度层。单线程运行时，RiverNetwork 直接保留一个覆盖全域的子流域，避免无意义的调度开销。

= 阈值与并行效率

$omega_min$ 决定任务粒度：

- 阈值较低：子流域数量多、任务较细，并行度高，但调度和同步开销增大；
- 阈值较高：子流域数量少、任务较粗，调度开销低，但容易出现线程空闲和负载不均；
- 阈值超过或接近最大河流级别时：整个流域可能退化为单一任务，基本无并行收益。

合理阈值应在任务数量、单个任务工作量和线程数之间折中。它是*性能参数*，原则上不应改变水文方程和汇流路径。实际选择宜以计时结果为准，而非仅依据子流域数量。

= GuanShan 流向数据验证

采用 `GuanShan_flwdir.tif` 验证 RiverNetwork。原始栅格为 ArcGIS D8 编码，尺寸为 $78 times 60$；将编码转换为 RiverNetwork 使用的 PCRaster LDD 后，得到 1837 个有效节点、1836 条边和 1 个出口。河网无环，最大 Strahler 级别为 6，各级节点数如下：

#table(
  columns: (auto,) * 6,
  align: center,
  inset: 5pt,
  stroke: 0.6pt + rgb("9aadb8"),
  fill: (x, y) => if y == 0 { rgb("dfeaf0") },
  [Strahler 级别], [1], [2], [3], [4], [5], [6],
  [节点数], [1274], [316], [122], [74], [33], [18],
)

在 8 线程条件下，不同阈值对应的子流域数和调度层宽度为：

#table(
  columns: (auto, auto, auto, 1fr),
  align: (center, center, center, left),
  inset: 5pt,
  stroke: 0.6pt + rgb("9aadb8"),
  fill: (x, y) => if y == 0 { rgb("dfeaf0") },
  [$omega_min$], [子流域数], [调度层数], [各层子流域数（上游 $arrow$ 下游）],
  [1], [1072], [6], [585, 186, 177, 97, 26, 1],
  [2], [215], [5], [132, 45, 26, 11, 1],
  [3], [45], [4], [30, 10, 4, 1],
  [4], [12], [3], [9, 2, 1],
  [5], [3], [2], [2, 1],
  [6], [1], [1], [1],
)

该结果直观说明：阈值降低会快速增加上游并行任务；但子流域数并非越多越好，最终性能还受任务大小差异、线程数和同步成本控制。

= 适用条件与限制

运动波方法适用于坡度和摩阻主导、回水及惯性效应较弱的单向汇流过程。以下情形应谨慎使用：

- 河网存在显著回水、潮汐、倒流或强惯性效应；
- 流向图含环、分汊或不一致的边界流向；
- 河漫滩横向交换决定主要传播过程；
- 子流域工作量高度不均，静态分层难以充分利用线程。

前两类问题分别属于水动力学假设和河网拓扑约束，不能通过增加子流域数量解决。若回水和双向传播不可忽略，应采用局部惯性或更完整的浅水方程方法。

= 结论

运动波汇流可概括为“拓扑有序的质量守恒”：每个单元先汇总全部上游新流量，再求解本单元的隐式运动波方程。RiverNetwork 的子流域划分则可概括为“依赖保持的图压缩”：以 Strahler 级别跃迁点切分河网，在子流域内部保持串行拓扑序，在同一上游距离层内并行，并按上游到下游设置层间同步。因此，该方法在不改变汇流物理过程的前提下，将单元级长依赖链转化为适合多线程执行的分层任务图。

= 对应实现

本文描述对应以下源文件：

- `Wflow/RiverNetwork.jl/src/network.jl`：流向图构建与拓扑检查；
- `Wflow/RiverNetwork.jl/src/subdomains.jl`：Strahler 分级、子流域标记、子流域图及并行层；
- `Wflow/src/routing/surface/surface_process.jl`：单元运动波非线性方程；
- `Wflow/src/routing/surface/surface_kinwave.jl`：坡面与河道汇流的分层并行执行。

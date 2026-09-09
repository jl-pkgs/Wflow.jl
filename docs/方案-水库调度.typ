#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "", header: "")
#let delta(x) = $Delta #x$
// #v(0.6em)
// #pagebreak()

#h(2em)
// 定义水库控制的集水区域为水库控制单元。
水库控制单元看做是一个特殊的子流域，入库流量 $Q_"in"$ 由集总式水文模型模拟；出库流量 $Q_"out"$ 受调度策略影响。正常情况下，水库在枯水期蓄水以应对干旱，在丰水期降低水位以满足防洪需求。
$chi_"target,max"$、$chi_"target,min"$分别是最高、最低目标蓄水比。*$chi_"target,max"$ 控制防洪*：库容高于 $V_max chi_"target,max"$ 时启动超额下泄，将蓄水压回该目标。汛期降低该上限以腾库迎汛，枯水期抬高以兴利蓄水。*$chi_"target,min"$ 控制保库*：其为下游需水保证率的拐点，库容低于该值时压缩需水，避免水库被放空；库容充足时才按 $"Demand"$ 供水。故前者应对水多，后者应对水少。湖北可按 4–9 月汛期、6 月腾至最低、10 月蓄至最高给定月过程（图@fig:chi-target）。

#figure(
  image("images/chi_target_month.png", height: 7.0cm),
  caption: [湖北水库目标蓄水比月过程示意。蓝色带为 4–9 月汛期；$chi_"target,max"$ 6 月最低、10 月最高。数值仅作形状示意，须按具体水库汛限水位与兴利水位率定。],
) <fig:chi-target>

#h(2em)
Wflow模型中水库调节主要理论如下（对应 Wflow 出流类型 4）。

#h(2em)
*① 首先考虑入库、降水和蒸发对库容 $V$（$m^3$）的影响：*

$ V'_t = max(V_(t-1) + Q_"in" delta(t) + 0.001 (P - "ET"_"water") A_"res", 0) $

其中 $P$ 为时段降水（mm），$"ET"_"water"$ 为时段水面蒸发（mm），$A_"res"$ 为水库面积（$m^2$）。$0.001 (P - "ET"_"water") A_"res"$ 把水深换成体积。$V'_t$ 为扣排泄前的临时库容。

#h(2em)
*② 在当前库容允许的前提下，满足下游用水需求：*

$ f_"guarante" (chi_t) = 1 / (1 + e^(-c (chi_t - a))) $ <eq_scurve>

$ Q_("demand",t) = min(f_"guarante" (chi_t) "Demand"_t, V'_t / delta(t)) $

其中 $"Demand"_t$ 为下游用水需求（$m^3\/s$）。$chi_t = V'_t / V_max$ 是临时库容占总库容之比。$f_"guarante"$ 为用水保证率，值域 $(0,1)$：蓄水不足时接近 0，充足时接近 1。#highlight[Wflow 取 $c=30$]，拐点 $a = chi_"target,min"$（保证率为 50% 时的充满比）。需水上限是当前可放流量 $V'_t \/ delta(t)$，不以死库容截断。

#figure(
  image("images/f_guarante_c.png", height: 7.1cm),
  caption: [用水保证率 $f_"guarante"$ 与充满比 $chi_t$ 的关系（图中 $a=0.3$）。$c$ 越大，在拐点附近转折越陡。],
) <fig:f-guarante>

放完需水后：

$ V_(1,t) = V'_t - Q_("demand",t) delta(t) $

#h(2em)
*③ 在保障目标最高蓄水的前提下，向下游排水：*
超额下泄希望把库容降到 $V_max chi_"target,max"$；超过 $V_max$ 的部分为溢洪。二者均为流量：

$ Q_"want" = max(V_(1,t) - V_max chi_"target,max", 0) / delta(t) $

$ Q_"spill" = max(V_(1,t) - V_max, 0) / delta(t) $ <eq_q_spill>

$Q_"release"^'$ 还受最大泄流能力约束：溢洪道最大排泄 $Q_("release",max)$ 与溢洪 $Q_"spill"$ 之和，并扣掉已放的需水。#highlight[$Q_"release"^'$是想放多少水，以及最快能放多少水的均衡。]

$
  Q_"release"^' = min(Q_"want", Q_"spill" + Q_("release",max) - Q_("demand",t))
$ <eq_q_release>

#h(2em)
最终，总排泄量：

$ Q_"out" = Q_"release" = Q_"release"^' + Q_("demand",t) $

$ V_t = V_(1,t) - Q_"release"^' delta(t) $

得到 $Q_"out"$ 之后，后续演算与其他子流域相同，由汇流模块演算到流域出口。

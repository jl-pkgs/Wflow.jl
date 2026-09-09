#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "", header: "", size: 11pt)

#set document(title: "SBM + 运动波 + 水库：输入数据准备")
#show table: set text(size: 8.8pt)
#show table.cell.where(y: 0): set text(font: "Noto Sans CJK SC", weight: "bold", size: 8.5pt)

#let delta(x) = $Delta #x$
#let req = text(fill: rgb("8b1e1e"))[必填]
#let opt = text(fill: rgb("3d5a40"))[可选]
#let cond = text(fill: rgb("8a5a12"))[条件]

#align(center)[
  #text(font: "Noto Sans CJK SC", size: 18pt, weight: "bold")[
    SBM + 运动波 + 水库的输入数据
  ]
  #v(-0.4em)
  #text(size: 10.5pt)[不含泥沙与地下水（`sbm_gwf`）]
]

#block(
  fill: rgb("f2f6f8"),
  stroke: 0.6pt + rgb("9aadb8"),
  radius: 3pt,
  inset: 10pt,
)[
  *摘要*　本配置只需三类文件：TOML 设置、与强迫同网格的静态 NetCDF、气象强迫 NetCDF。域图给出流向、河道与子流域；SBM 给出土壤与植被；运动波给出坡度、河宽河长与糙率；水库再补位置、面积与出流曲线。积雪、冰川、周期 LAI 与暖启动状态按需增加。不需要泥沙参数、二维地下水参数，也不需要局部惯性法额外高程。
]

== 配置范围 // <!-- omit in toc -->

`sbm` + 运动波（kinematic wave） + 水库：

```toml
[model]
type = "sbm"
land_routing = "kinematic_wave"   # 默认，可省略
river_routing = "kinematic_wave"  # 默认，可省略
reservoir__flag = true
```

侧向壤中流随 `sbm` 一并启用，由土壤厚度、孔隙度和坡度推求。下列内容明确排除：`type = "sbm_gwf"`、`type = "sediment"`、局部惯性法 / 交错网格河道、洪泛区、用水需求。

*文件组成*：路径相对 TOML，或相对 `dir_input`。

#figure(
  table(
    columns: (auto, 1.15fr, 1.6fr),
    align: (left, left, left),
    inset: 6pt,
    stroke: 0.5pt + rgb("9aadb8"),
    fill: (_, y) => if y == 0 { rgb("e6eef2") } else if calc.odd(y) { rgb("f7fafb") },
    [文件], [作用], [格式要点],
    [`settings.toml`], [时段、开关、标准名映射], [唯一入口],
    [`staticmaps.nc`], [域图与空间参数], [与强迫同网格；通量参数*日*为单位],
    [`forcing.nc`], [降水、潜在蒸发（及气温）], [三维 `(x, y, time)`；右标注],
    [`LAI` 等周期场], [年内循环参数], [12 / 365 / 366 层；按当步 DOY 取值],
    [`reservoir_sh_*.csv` / `reservoir_hq_*.csv`], [库容曲线 $S(H)$、出流曲线 $Q(H)$], [仅实测曲线时],
    [`instates.nc`], [暖启动状态], [`cold_start__flag = false` 时],
  ),
  caption: [本配置的输入文件。除水库曲线为 CSV 外，空间数据均为 NetCDF。],
) <tab:files>
#table-note[
  强迫必须右标注：时间戳 `2023-06-15 00:00:00` 的日降水，表示 6 月 14 日 00:00 至 15 日 00:00 的累积量。静态图中的通量参数（如饱和导水率、入渗能力）以日为基准，初始化时换算到 `timestepsecs`。
]

#pagebreak()

= 1 容易获取部分

== 1.1 流域与河网

下列图层不进入 SBM 方程，但决定计算域与汇流拓扑。流向采用 PCRaster LDD（1–9）。

#figure(
  table(
    columns: (auto, 1.7fr, auto, auto, auto),
    align: (left, left, left, center, center),
    inset: 5.5pt,
    stroke: 0.5pt + rgb("9aadb8"),
    fill: (_, y) => if y == 0 { rgb("e6eef2") } else if calc.odd(y) { rgb("f7fafb") },
    [标准名], [含义], [常用变量], [单位], [属性],
    [`basin__local_drain_direction`], [D8 流向], [`wflow_ldd`], [—], [#req],
    [`subbasin_location__count`], [子流域 / 计算域], [`wflow_subcatch`], [—], [#req],
    [`river_location__mask`], [河道掩膜], [`wflow_river`], [—], [#req],
    [`land_surface__elevation`], [地面高程], [`wflow_dem`], [m], [#req],
    [`land_surface__slope`], [坡面坡度], [`Slope`], [m m#super[-1]], [#req],
    [`river_gauge__count`], [站点编号，仅用于点输出], [`wflow_gauges`], [—], [#opt],
  ),
  caption: [域与河网图层。河道掩膜在元数据中可填 0，实际运行必须提供。],
) <tab:domain>

TOML 中这些名称写在 `[input]`，不在 `[input.static]`。

== 1.2 土壤

#let l1 = 0.8cm

土层厚度在 TOML 中给出，默认 `[100, 300, 800]` mm，最后一层由 `soil__thickness` 截断。

#wrap-table(
  table(
    columns: (7.4cm, 2.5cm, auto, 1cm, 2.3cm),
    rows: (1cm, l1, l1, l1, l1, 1.2cm),
    align: (horizon + left, horizon, horizon, horizon, horizon),
    inset: 5.5pt,
    stroke: 0.5pt + rgb("9aadb8"),
    fill: (_, y) => if y == 0 { rgb("e6eef2") } else if calc.odd(y) { rgb("f7fafb") },
    [标准名], [含义], [常用变量], [单位], [属性],
    [`soil__thickness`], [土壤厚度], [`SoilThickness`], [mm], [#req],
    [`soil_water__saturated_volume _fraction`], [饱和含水率 $theta_s$], [`thetaS`], [—], [#req],
    [`soil_water__residual_volume_fraction`], [残余含水率 $theta_r$], [`thetaR`], [—], [#req],
    [`compacted_soil__area_fraction`], [不透水面积], [`PathFrac`], [—], [#req],
    [`soil_layer_water__brooks_ corey_exponent`], [Brooks–Corey 指数], [`c`], [—], [#req],
    [`soil_surface_water__vertical_ saturated_hydraulic_conductivity`],
    [地表垂直饱和导水率],
    [`KsatVer`],
    [mm d#super[-1]],
    [#req],
    
    [`soil_water__vertical_saturated_hydraulic_ conductivity_scale_parameter`],
    [导水率随深度衰减],
    [`f`],
    [mm#super[-1]],
    [#req],
    
    [`compacted_soil_surface_water__ infiltration_capacity`],
    [#highlight[不透水地表入渗能力]],
    [`InfiltCapPath`],
    [mm d#super[-1]],
    [#opt 10],
    
    [`land_water_covered__area_fraction`], [开阔水面比], [`WaterFrac`], [—], [#opt 0],
    [`soil_water_saturated_zone_bottom__ max_leakage_volume_flux`],
    [饱和带最大渗漏],
    [`MaxLeakage`],
    [mm d#super[-1]],
    [#opt 0],
    
    [`soil_water__field_capacity_ volume_fraction`], [田间持水率], [—], [—], [#opt 由 $theta_s,theta_r,c$ 推求],
    [`soil_root__length_density_fraction`], [各层根长密度比], [—], [—], [#opt 由根系深度分配],
  ),
  caption: [土壤参数。未列的 Feddes 临界水头、毛管上升与冻土折减均有默认值。],
) <tab:soil>
#table-note[
  默认导水率剖面为指数型，只需 `KsatVer` 与 `f`。若改为 `exponential_constant`、`layered` 或 `layered_exponential`，再补相应深度或分层 $K$ 图。自然土入渗能力由地表 $K$ 与首层放大系数相乘得到，不必单独给图。
]

== 1.3 气象强迫

#figure(
  table(
    columns: (auto, 1.4fr, auto, auto, auto),
    align: (left, left, left, center, center),
    inset: 5.5pt,
    stroke: 0.5pt + rgb("9aadb8"),
    fill: (_, y) => if y == 0 { rgb("e6eef2") } else if calc.odd(y) { rgb("f7fafb") },
    [标准名], [含义], [常用变量], [单位], [属性],
    [`atmosphere_water__precipitation_volume_flux`], [降水], [`precip`], [mm Δt#super[-1]], [#req],
    [`land_surface_water__potential_evaporation_volume_flux`], [潜在蒸发], [`pet`], [mm Δt#super[-1]], [#req],
    [`atmosphere_air__temperature`], [气温], [`temp`], [°C], [#cond],
  ),
  caption: [强迫变量。单位与模型步长一致，不是静态参数所用的“每日”。],
) <tab:forcing>

*气温仅在打开积雪或冰川时需要。*

= 2 SBM 陆面

SBM 参数写入 `[input.static]`。

== 2.1 植被：两条路径，择一

提供周期 LAI 时，由比叶蓄水、木质部蓄水和消光系数把 LAI 换成冠层间隙率与最大冠层蓄水；否则直接给这两项。

#wrap-table(
  table(
    columns: (6.5cm, 1.35fr, auto, auto, auto),
    align: (left, horizon, horizon, horizon, horizon),
    inset: 5.5pt,
    stroke: 0.5pt + rgb("9aadb8"),
    fill: (_, y) => if y == 0 { rgb("e6eef2") } else if calc.odd(y) { rgb("f7fafb") },
    [标准名], [含义], [常用变量], [单位], [属性],
    [`vegetation_root__depth`], [根系深度], [`RootingDepth`], [mm], [#req],
    [`vegetation__crop_factor`], [作物系数], [—], [—], [#opt 1.0],
    [`vegetation_canopy_water__mean_ evaporation_to_mean_ precipitation_ratio`],
    [Gash 截留比],
    [`EoverR`],
    [—],
    [#opt 0.1],
    table.cell(colspan: 5, fill: rgb("eef4e8"))[路径 A：周期 LAI（推荐）],
    [`vegetation__leaf_area_index`], [叶面积指数], [`LAI`], [m#super[2] m#super[-2]], [#cond],
    [`vegetation__specific_leaf_storage`], [#highlight[比叶蓄水]], [`Sl`], [mm], [#cond],
    [`vegetation_wood_water__ storage_capacity`], [#highlight[木质部蓄水]], [`Swood`], [mm], [#cond],
    [`vegetation_canopy__light_ extinction_coefficient`], [#highlight[消光系数]], [`Kext`], [—], [#cond],
    table.cell(colspan: 5, fill: rgb("eef4e8"))[路径 B：静态冠层],
    [`vegetation_canopy__gap_fraction`], [冠层间隙率], [—], [—], [#cond],
    [`vegetation_water__storage_capacity`], [最大冠层蓄水], [—], [mm], [#cond 默认 1.0],
  ),
  caption: [植被参数。路径 A 将 LAI 写入 `[input.cyclic]`。],
) <tab:veg>

每个时间步，路径 A 仅用下式更新冠层：

$ S_"canopy,max" = S_"leaf" thin "LAI" + S_"wood" $
$ f_"gap" = exp(-k thin "LAI") $

其中 $S_"leaf"$ 为比叶蓄水（`Sl`，mm），$S_"wood"$ 为木质部蓄水（`Swood`，mm），$k$ 为消光系数（`Kext`，无量纲）。三者均为土地利用查表量，不随时间变化；随物候变化的只有 LAI。取值按覆被类型给定（Pitman, 1989；Liu, 1998；van Dijk and Bruijnzeel, 2001）：林地 $S_"leaf"$ 约 0.1–0.3 mm、$S_"wood"$ 约 0.2–2 mm、$k$ 约 0.5–0.7；草地、作物则更小。三者只进入截留，不参与蒸腾或土壤。

= 3 运动波汇流

河道与坡面走运动波；壤中流亦为运动波，参数大多由土壤与坡度派生。

#wrap-table(
  table(
    columns: (5.4cm, 1.4fr, auto, auto, auto),
    align: (left + horizon, horizon, horizon, horizon, horizon),
    inset: 5.5pt,
    stroke: 0.5pt + rgb("9aadb8"),
    fill: (_, y) => if y == 0 { rgb("e6eef2") } else if calc.odd(y) { rgb("f7fafb") },
    [标准名], [含义], [常用变量], [单位], [属性],
    [`river__slope`], [河道坡度], [`RiverSlope`], [m m#super[-1]], [#req],
    [`river__length`], [河段长度], [`wflow_riverlength`], [m], [#req],
    [`river__width`], [河宽], [`wflow_riverwidth`], [m], [#req],
    [`river_water_flow__manning_ n_parameter`], [河道糙率], [`N_River`], [s m#super[-1/3]], [#opt 0.036],
    [`land_surface_water_flow__ manning_n_parameter`], [坡面糙率], [`N`], [s m#super[-1/3]], [#opt 0.072],
    [`subsurface_water__horizontal_ to_vertical_saturated_ hydraulic_conductivity_ratio`],
    [水平 / 垂直导水率比],
    [`KsatHorFrac`],
    [—],
    [#req],
  ),
  caption: [运动波参数。河长、河宽元数据默认为 0，实际必须给有效值。],
) <tab:kw>

#table-note[
  `KsatHorFrac` 把垂直饱和导水率换成水平向，供侧向壤中流使用。满岸水深 `river_bank_water__depth` 仅局部惯性法需要，本配置可省略。
]

= 4 水库

打开 `reservoir__flag` 后，至少需要库区覆盖与出口位置：

```toml
[input]
reservoir_area__count = "wflow_reservoirareas"
reservoir_location__count = "wflow_reservoirlocs"
```

所有水库还要给出水面面积、库容曲线类型、出流曲线类型和初始水位。缺测不允许，即便该类型不用某参数，也需填占位值（如 `-1`）。

#figure(
  table(
    columns: (auto, 1.5fr, auto, auto),
    align: (left, left, left, center),
    inset: 5.5pt,
    stroke: 0.5pt + rgb("9aadb8"),
    fill: (_, y) => if y == 0 { rgb("e6eef2") } else if calc.odd(y) { rgb("f7fafb") },
    [标准名], [含义], [常用变量], [属性],
    [`reservoir_surface__area`], [水面面积 $A$], [`reservoir_area`], [#req],
    [`reservoir_water__storage_curve_type_count`], [库容曲线：\ 1 为 $S=A H$; \ 2 为实测 $S=f(H)$], [`storfunc`], [#req],
    [`reservoir_water__rating_curve_type_count`], [出流曲线类型，见下表], [`outflowfunc`], [#req],
    [`reservoir_water_surface__initial_elevation`], [初始水位], [`waterlevel_reservoir`], [#req],
  ),
  caption: [各类水库共用的静态参数。],
) <tab:res-common>

- *库容曲线* $S=f(H)$：`reservoir_sh_<id>.csv`，两列 `H,S`。第一列 $H$ 为水位（m），第二列 $S$ 为库容（m#super[3]）。
- *出流曲线* $Q=f(H)$：仅出流类型 1 需要 `reservoir_hq_<id>.csv`。第一列必须是水位 $H$（m）；其后 365 列是流量 $Q$（m#super[3] s#super[-1]），依次对应一年中第 1–365 天。代码用当前水位在第一列上插值，再取当天那一列。只有一条 $Q(H)$ 时，把同一列复制 365 遍。CSV 与 TOML 同目录，编号与出口 id 一致。

#figure(
  table(
    columns: (auto, 1.2fr, 1.8fr),
    align: (center, left, left),
    inset: 5.5pt,
    stroke: 0.5pt + rgb("9aadb8"),
    fill: (_, y) => if y == 0 { rgb("e6eef2") } else if calc.odd(y) { rgb("f7fafb") },
    [类型], [出流关系], [额外数据],
    [1], [实测 $Q=f(H)$], [`reservoir_hq_<id>.csv`：第 1 列水位 $H$，后 365 列流量],
    [2], [$Q=b(H-H_0)^e$], [系数 $b$、指数 $e$、阈值水位 $H_0$],
    [3], [Modified Puls，$Q=b(H-H_0)^2$], [系数 $b$、阈值 $H_0$；库容曲线须为 $S=A H$],
    [4], [简单调度规则], [最大库容、目标库容比、下游需水、溢洪道以下最大泄流],
  ),
  caption: [出流曲线类型。库容曲线为类型 2 时另备 `reservoir_sh_<id>.csv`。],
) <tab:res-type>

一般不必用类型 1：类型 2 / 3 给 $b,e,H_0$，类型 4 给调度参数即可。

== 4.1 出流曲线

=== 4.1.1 类型 3（Modified Puls）把水库水量平衡显式解出。

#h(2em)
假定库容$S$与水位$H$成正比，$S = A H$；出流为抛物堰$Q = b (H - H_0)^2$。从而库容$S = A sqrt(Q / b) + A H_0$。其中$A$ 为水面面积（$m^2$），$H_0$ 为出流阈值水位（m），$b$ 为堰流系数，$delta(t)$ 为时间步长（s）。

#h(2em)
代入水量平衡。$"SI"$是本时段来水（流量量纲）：期初库容折成流量，加上入流与库面净雨。$"LF"$ 只含面积与堰流系数，是二次式的系数。

$ "SI" = S(t) / delta(t) + Q_"in" + A (P - E) / delta(t), quad "LF" = A / (delta(t) sqrt(b)) $

其中，$Q$ 为出流（$m^s\/s$），$Q_"in"$ 为入流（$m^s\/s$），$P$、$E$ 为库面降水与蒸发水深（m）。联立解得，出流Q为：
$ Q = cases(
  1/4 (-"LF" + sqrt("LF"^2 + 4 ("SI" - A H_0 / delta(t))))^2 & "SI" > A H_0 / delta(t),
  0 & "SI" <= A H_0 / delta(t),
) $

#highlight[输入只需 $A$、$b$、$H_0$；库容曲线必须是类型 1。与类型 2 的差别是指数固定为 2，一步算出 $Q$。]

实测曲线的 $H$ 多为海拔；库容曲线类型 1（$S=A H$）的 $H$ 是相对库底水深。二者不可混用。

== 4.2 类型 4（简单调度）

#figure(
  table(
    columns: (8cm, 1.5fr, auto),
    align: (left, horizon, horizon),
    inset: 5.5pt,
    stroke: 0.5pt + rgb("9aadb8"),
    fill: (_, y) => if y == 0 { rgb("e6eef2") } else if calc.odd(y) { rgb("f7fafb") },
    [标准名], [含义], [常用变量],
    [`reservoir_water__max_volume`], [最大库容 $S_"max"$], [`ResMaxVolume`],
    [`reservoir_water__target_full_volume_fraction`], [目标充满比], [`ResTargetFullFrac`],
    [`reservoir_water__target_min_volume_fraction`], [目标最小充满比], [`ResTargetMinFrac`],
    [`reservoir_water_demand__required_downstream_ volume_flow_rate`], [下游最小泄流], [`ResDemand`],
    [`reservoir_water_release_below_spillway__max_ volume_flow_rate`], [溢洪道以下最大泄流], [`ResMaxRelease`],
  ),
  caption: [类型 4 参数，可作静态图或周期 / 强迫。],
) <tab:res4>

// = 5 可忽略部分

// == 5.1 积雪、周期场与状态

// 积雪、冰川默认关闭。打开 `snow__flag` 后必须提供气温；度日因子等均有默认。

// 冷启动（`cold_start__flag = true`，默认）不读状态文件。暖启动需 `path_input`，并在 `[state.variables]` 列出冠层蓄水、非饱和层深、饱和带水深、坡面与河道水深流量、壤中流，以及水库水位。

// == 5.2 其他

// - 地下水（`sbm_gwf`）：含水层底板、导水系数二维场、定水头与排水沟。
// - 局部惯性法：满岸高程、下游边界河长、二维地面高程、洪泛区剖面。
// - 用水需求：生活、工业、灌溉与分配区，除非另行打开相应开关。

// 空间图层可用 HydroMT-wflow 从 DEM、河网与土壤植被产品生成，再按标准名写入 TOML。

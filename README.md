# Optimizer.jl 使用指南

Optimizer.jl 是一个小型 Julia 工具包，包含：一维求根（bisection / Brent）、自适应 Gauss–Legendre 积分、数值导数、带 loess 平滑的误差估计（STD），以及拟合用的各类 χ² 函数。

- 版本：`1.1.0-DEV`
- Julia：≥ 1.11
- 依赖：`DataFrames`、`FastGaussQuadrature`、`LinearAlgebra`、`Loess`、`Statistics`
- 许可证：MIT

---

## 1. 安装

```julia
using Pkg
Pkg.add(url="https://github.com/zejianzhuang-uv/Optimizer.jl")
# 或本地开发模式（推荐，边改边用）
Pkg.develop(url="https://github.com/zejianzhuang-uv/Optimizer.jl")
```

```julia
using Optimizer
```

运行测试：

```julia
Pkg.test("Optimizer")
```

## 2. 导出接口一览

| 类别 | 函数 |
|---|---|
| 求根 | `bisectv`, `brentqv` |
| 积分 | `gaussian_quad` |
| 导数 | `derivative` |
| 误差估计 | `STD`, `df_STD`, `loess_smooth` |
| χ² | `χsq`, `weightedχsq`, `red_weightedχsq`, `hyperχsq`, `whyperχsq`, `hyperα` |

> 所有数值参数基本要求 `Float64`（或 `Vector{Float64}`），传 `Int` 会触发 `MethodError`。

---

## 3. 求根

### 3.1 `brentqv(f, init_x, n)`（推荐）

在网格 `init_x` 上扫描变号区间，对每个区间用 Brent 法求根，最多返回 `n` 个根。

```julia
brentqv(f::Function, init_x::AbstractVector{Float64}, n::Int64) -> Vector{Float64}
```

```julia
brentqv(sin, -4π:0.1:4π, 10)
brentqv(x -> x^2 - 1, -2:0.01:2, 10)   # ≈ [-1.0, 1.0]
```

要点：

- `init_x` 可以是 range 或 vector，网格越密越不容易漏根（相邻点之间有两个根时会漏）。
- 返回的根按 `init_x` 顺序排列，数量 ≤ `n`。
- 内部容差：`xtol = rtol = ftol = 1e-8`，最大迭代 1000。
- 收敛后若 `|f(x)| > 0.1`，该根被丢弃（返回 `NaN` 并被过滤）。**这会自动剔除极点处的假变号**，但对函数值本身很陡的根也可能误杀，用的时候留意。
- `|x| < 1e-6` 的根会被取整为 `0.0`。
- 不收敛时会 `println("iteration reaches its maximum.")`。

### 3.2 `bisectv(f, init_x, n)`

接口同上，但用二分法，且 `init_x` 必须是 `Vector{Float64}`（range 要先 `collect`）。

```julia
bisectv(x -> x^2 - 1, collect(-2:0.01:2), 5)
```

更慢但更稳，适合 `brentqv` 表现异常时做交叉验证。

---

## 4. 数值积分：`gaussian_quad`

```julia
gaussian_quad(fcn::Function, a::Float64, b::Float64; npt::Int64=10) -> value
```

自适应 Gauss–Legendre：从 `npt` 点开始，每次加 10 个节点，直到相邻两次结果的相对变化 < `5e-2`；节点数到 100 时停止（实际最高用到 90 点）。

```julia
gaussian_quad(x -> exp(-x^2), -3.0, 3.0)
gaussian_quad(sin, 0.0, 1.0; npt=20)
```

要点：

- 积分限必须是 `Float64`（`0` 会报错，写 `0.0`）。
- 收敛判据 `5e-2` 是写死的，**精度只有 ~5%量级的保证**；光滑被积函数通常远好于此，但奇异/振荡被积函数请自行分段并检查。

---

## 5. 数值导数：`derivative`

```julia
derivative(f::Function, x0::Union{Float64,ComplexF64}; ftol=1e-1, max_iter=50)
```

前向差分，步长从 `eps()*1e12` 起每次缩小 10 倍，直到相邻两次结果的相对变化 < `ftol`。

```julia
derivative(sin, 1.0)             # ≈ cos(1.0)
derivative(z -> z^2, 1.0 + 1.0im)
derivative(sin, 1.0; ftol=1e-4)  # 更严格
```

要点：

- 默认 `ftol=1e-1` 很松，要精确值请调小；调得太小可能被舍入误差卡住并触发 `@warn "Not converged"`。
- `x0` 为 `NaN` 时直接返回 `NaN`。
- 支持复数 `x0`（解析函数）。

---

## 6. 误差估计：`STD` / `df_STD` / `loess_smooth`

### 6.1 `STD(sample; dims, kwargs...)`

沿 `dims` 求标准差，去掉该维度，再对结果做 loess 平滑。结果必须是 1 维或 2 维，否则抛 `ArgumentError`。

```julia
# 例：关联函数样本 (Nt, Nobs, Nconf)，沿组态维求误差
A = rand(48, 5, 100)
err = STD(A; dims=3)          # 48×5，每列沿 t 方向 loess 平滑
err = STD(A; dims=3, corrected=false)   # kwargs 透传给 Statistics.std
```

> 注意：这里是普通样本标准差，**不含 jackknife 的 (N−1) 缩放**。如果传入的是 jackknife 样本，需要自己乘因子。

### 6.2 `df_STD(sample, name; kwargs...)`

同 `STD`，但返回 `DataFrame`，列名由 `name` 给定。`dims` 同样通过 kwargs 传入；结果必须是矩阵（2 维）。

```julia
df = df_STD(A, [:a, :b, :c, :d, :e]; dims=3)
```

### 6.3 `loess_smooth`

```julia
loess_smooth(f::AbstractVector{Float64}; tol=0.0)                          # x = 1:n
loess_smooth(f::AbstractVector{Float64}, x::AbstractVector{Float64}; tol=0.0)
loess_smooth(f::AbstractMatrix{Float64})                                   # 逐列平滑
```

- `|f| ≤ tol` 的点被视为"无效点"：不参与拟合，原值保留（默认 `tol=0.0` 即跳过恰好为 0 的点）。
- 有效点少于 3 个时直接返回原数组的拷贝。

---

## 7. χ² 函数

约定：`yth` 理论值，`yexp` 实验值。

### 7.1 `χsq`

```julia
χsq(yth, yexp, inv_σ::AbstractVector{Float64})   # Σ ((yth-yexp)·inv_σ)²
χsq(yth, yexp, inv_σ::AbstractMatrix{Float64})   # (yth-yexp)ᵀ · inv_σ · (yth-yexp)
```

> **向量版传的是 1/σ（不是 1/σ²）；矩阵版传的是逆协方差矩阵。**

### 7.2 `weightedχsq(yth, yexp, inv_σ, nk)`

`χsq / nk`，`nk` 为该散射过程的观测量数。向量 / 矩阵 `inv_σ` 都支持。

### 7.3 `red_weightedχsq(all_weight_chi2, Ntot, nfit, Nsp)`

多过程联合拟合的约化 χ²：

```
χ²_red = Ntot / (Ntot − nfit) · (1 / Nsp) · all_weight_chi2
```

- `Nsp`：散射过程数（如 K⁻p → {π⁰Λ, …}、π⁻p → {ηn, …}）
- `Ntot`：所有过程的观测量总数
- `nfit`：拟合参数个数
- `all_weight_chi2`：各过程 `weightedχsq` 之和

### 7.4 Hobson–Bridle–Lahav 超参数 χ²（仅向量 `inv_σ`）

对数据集权重 α 以先验 e^(−α) 边缘化后的 −2 ln L：

```julia
hyperχsq(yth, yexp, inv_σ, n)   # (n+2) · log(χ² + 2)
whyperχsq(yth, yexp, inv_σ, n)  # hyperχsq / n（与 weightedχsq 同归一化）
hyperα(yth, yexp, inv_σ, n)     # 有效权重 α = (n+2)/(χ²+2)
```

### 7.5 多数据集拟合示例

```julia
function total_chi2(p)
    c = 0.0
    for (k, d) in enumerate(datasets)
        yth = model(p, d.x)
        c += weightedχsq(yth, d.y, d.inv_σ, length(d.y))
    end
    return c
end

# 约化 χ²
red = red_weightedχsq(total_chi2(p), Ntot, length(p), length(datasets))

# 超参数方案：各数据集自动配权
total_hyper(p) = sum(whyperχsq(model(p, d.x), d.y, d.inv_σ, length(d.y)) for d in datasets)
```

---

## 8. 注意事项 / 已知问题

1. `brentq`、`bisect` 是内部函数，**未导出**，需要时用 `Optimizer.brentq(...)`（签名：`brentq(f, xa, xb, fa, fb; ...)`，要同时传两端函数值）。
2. `brentqv_parallel` / `find_sign_changes_parallel`（多线程版）未导出，且 `brentqv_parallel` 目前调用了不存在的 3 参数 `brentq(f, a, b)`，**暂时不可用**。
3. `bisectv` 要求 `Vector{Float64}`，`brentqv` 接受任意 `AbstractVector{Float64}`。
4. 精度默认较松（`gaussian_quad` 的 5e-2、`derivative` 的 1e-1），做精密计算前请自行验证。
5. README 里的 CI 徽章 URL 写的是 `Optimize.jl`（少了 r），与仓库名 `Optimizer.jl` 不一致，需要修正才能正常显示。

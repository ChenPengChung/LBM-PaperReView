# WENO5 Interpolation 引入 GILBM — 最佳化程式碼初步構思

> **日期**: 2026-03-28
> **目標**: 將現有 7-point Lagrange 插值升級為 WENO5 interpolation
> **對應程式碼**: `gilbm/interpolation_gilbm.h`, `gilbm/evolution_gilbm.h`

---

## 一、現有架構分析

### 1.1 當前插值方式：7-point Lagrange (6th order)

現有 `interpolation_gilbm.h` 使用 dimension-by-dimension 的 7-point Lagrange 插值：

```
位置: gilbm/interpolation_gilbm.h
函式: lagrange_7point_coeffs(t, a[7]) — 計算 Lagrange 權重
函式: lagrange_basis_7(t, idx) — 單基函數 (division-free)
巨集: Intrpl7(...) — 7-point 加權求和
```

在 `evolution_gilbm.h` 的 `gilbm_compute_point()` 中：

```
Step 1: η 方向 — 7×(7×7) 次 Intrpl7 → interpolation1order[7][7]
Step 2: ξ 方向 — 7 次 Intrpl7 → interpolation2order[7]
Step 3: ζ 方向 — 1 次 Intrpl7 → f_streamed
```

每個方向 q (19 directions) 各做一次三維 Lagrange 插值。

### 1.2 問題與限制

1. **Lagrange 無法抑制振盪**: 在分佈函數不連續或劇烈變化處（如壁面附近、hill 分離區），高階 Lagrange 會產生 Gibbs 振盪
2. **無自適應性**: 權重固定，無法根據局部光滑度調整 stencil
3. **可能導致負分佈函數**: f < 0 → 密度/壓力非物理

---

## 二、WENO5 Interpolation 架構設計

### 2.1 核心數學：在任意偏移 σ 處的 WENO interpolation

已知 7 個格點的 point values: {f₀, f₁, ..., f₆}，目標點 x* = x₃ + σ (σ ∈ [0, 1))

**WENO5 使用 3 個候選 3-point stencils (r=0,1,2)**:

```
S₀ = {f₀, f₁, f₂}  → p₀(x*) = 二次多項式
S₁ = {f₁, f₂, f₃}  → p₁(x*) = 二次多項式
S₂ = {f₂, f₃, f₄}  → p₂(x*) = 二次多項式
```

(以 x₃ 為中心，向左移 r=0,1,2)

> 注意：WENO5 reconstruction 用 5 點 → 3 stencils of 3 points = 5th order
> WENO5 interpolation from point values 同樣架構，但 linear weights 隨 σ 變化

### 2.2 Linear Weights d_r(σ) — 核心差異

對於 WENO reconstruction (σ = 1/2)，linear weights 是固定常數：
- d₀ = 1/10, d₁ = 6/10, d₂ = 3/10

對於 WENO interpolation at arbitrary σ，需要解線性系統：
```
P(x*) = d₀(σ)·p₀(x*) + d₁(σ)·p₁(x*) + d₂(σ)·p₂(x*)
```
其中 P(x*) 是大 stencil {f₀,...,f₄} 上的 4th-order 多項式。

**d_r(σ) 是 σ 的有理函數，可能為負！**

→ 需使用 Liu-Shu-Zhang (2009) 正化技巧

### 2.3 Smoothness Indicators β_r

```
β₀ = (13/12)(f₀ - 2f₁ + f₂)² + (1/4)(f₀ - 4f₁ + 3f₂)²
β₁ = (13/12)(f₁ - 2f₂ + f₃)² + (1/4)(f₁ - f₃)²
β₂ = (13/12)(f₂ - 2f₃ + f₄)² + (1/4)(3f₂ - 4f₃ + f₄)²
```

**β_r 不隨 σ 變化** — 只衡量 stencil 上多項式的光滑程度。

### 2.4 Non-linear Weights

```
α_r = d_r(σ) / (ε + β_r)²
ω_r = α_r / Σ α_s
```

ε = 10⁻⁶ (防除零)

---

## 三、程式碼實作策略

### 3.1 新增 `weno5_interpolation_gilbm.h`

```c
// ═══════════════════════════════════════════════════
// WENO5 Interpolation at Arbitrary σ (1D)
// ═══════════════════════════════════════════════════
// Input: f[5] = {f_{i-2}, ..., f_{i+2}} (5 point values)
//        sigma = 目標點相對中心 i 的偏移 (可 ∈ [-0.5, 0.5] 或 [0, 1))
// Output: interpolated value at x_i + sigma

__device__ __forceinline__ double weno5_interp_1d(
    double f0, double f1, double f2, double f3, double f4,
    double sigma
) {
    // Step 1: 候選 stencil 多項式在 σ 處的值
    double p0 = ...; // from {f0, f1, f2}
    double p1 = ...; // from {f1, f2, f3}
    double p2 = ...; // from {f2, f3, f4}

    // Step 2: Linear weights d_r(σ) — 預計算或 runtime 計算
    double d0, d1, d2;
    compute_linear_weights_weno5(sigma, &d0, &d1, &d2);

    // Step 3: Smoothness indicators (σ-independent)
    double beta0 = (13.0/12.0)*SQR(f0-2*f1+f2) + 0.25*SQR(f0-4*f1+3*f2);
    double beta1 = (13.0/12.0)*SQR(f1-2*f2+f3) + 0.25*SQR(f1-f3);
    double beta2 = (13.0/12.0)*SQR(f2-2*f3+f4) + 0.25*SQR(3*f2-4*f3+f4);

    // Step 4: Non-linear weights (Liu-Shu-Zhang positive splitting if needed)
    // ... handle d_r < 0 case ...

    // Step 5: Final interpolation
    return omega0*p0 + omega1*p1 + omega2*p2;
}
```

### 3.2 修改 `evolution_gilbm.h` 的 streaming 插值

目前 dimension-by-dimension 架構保持不變，但每個方向的 Lagrange 替換為 WENO5：

```
原始: lagrange_7point_coeffs(t_eta, coeff) → Intrpl7(f[7], coeff[7])
升級: weno5_interp_1d(f[i-2]...f[i+2], sigma_eta)  (5 個輸入)
```

**注意**: WENO5 只需 5 個 stencil 點 (而非 Lagrange 的 7 個)，因此：
- 記憶體存取量下降: 7→5 per direction
- 但增加了 β_r 計算和非線性加權的運算量

### 3.3 CUDA Kernel 適配

#### Register 壓力分析

現有每個 thread 使用的 local 變數：
```
Lagrangarray_eta[7] + Lagrangarray_xi[7] + Lagrangarray_zeta[7] = 21 doubles
interpolation1order[7][7] = 49 doubles
interpolation2order[7] = 7 doubles
```

WENO5 替換後：
```
每個方向: 5 個 f 值 + 3 個 β + 3 個 d_r + 3 個 ω = 14 doubles
三個方向: ~42 doubles
中間陣列: interpolation1order[5][5] = 25 doubles (5×5, 非 7×7)
```

Register 使用量相近，但中間陣列較小 → **shared memory 壓力降低**

#### Plan C (Fused Kernel) 適配

Plan C 的 fused kernel 使用 `interp1[19*49]` local array。WENO5 可縮減為 `interp1[19*25]`：
```
記憶體節省: 19×49 = 931 → 19×25 = 475 doubles
每 thread 節省: 3648 bytes
```

### 3.4 Stencil 寬度調整

| 項目 | Lagrange-7 (現有) | WENO5 (目標) |
|------|-------------------|-------------|
| 1D stencil 寬度 | 7 點 | 5 點 |
| 3D stencil 體積 | 7³ = 343 | 5³ = 125 |
| Ghost zone 需求 | 3 層 | 2 層 (可降為 2) |
| 精度階 (光滑區) | 6th order | 5th order |
| 非振盪性 | 無 | 有 |
| f_pc 記憶體 | 19×343×GRID_SIZE | 19×125×GRID_SIZE (↓64%) |

**但若保持 ghost=3**: 可混合使用 WENO5 + Lagrange-7 作為 hybrid scheme

---

## 四、關鍵技術挑戰

### 4.1 Linear Weights 的 σ 依賴性

在 GILBM 中，每個方向 q、每個格點 (j,k) 的 departure point 偏移 σ 都不同：
- η 方向: `t_eta = ci - a_local * delta_eta[q]` → σ_η 隨 a_local 變化
- ξ 方向: `t_xi = cj - delta_xi_d[q*NYD6*NZ6 + idx_jk]` → σ_ξ 隨 (j,k,q) 變化
- ζ 方向: `t_zeta = up_k - bk` → σ_ζ 隨 (j,k,q) 變化

每次都需要 runtime 計算 d_r(σ)。

**最佳化策略**:
1. **解析公式**: 直接用 σ 的多項式表達 d_r(σ)，避免線性系統求解
2. **Lookup table**: 預計算 d_r 在 σ ∈ [0,1) 的離散值，GPU texture memory 插值
3. **η 方向特化**: η 方向 sigma 只依賴 a_local (per rank 常數)，可預計算

### 4.2 正化處理 (Liu-Shu-Zhang Splitting)

當 d_r(σ) < 0 時，需要：
```
d_r⁺ = max(d_r, 0) + δ    (正部分)
d_r⁻ = d_r⁺ - d_r         (負部分)
分別做 WENO 加權後相減
```

這會使計算量翻倍。在 GPU 上可用 branch-free 的 `fmax/fmin` 實現：
```c
double dp = fmax(d_r, 0.0);
double dm = dp - d_r;  // = |min(d_r, 0)| ≥ 0
```

### 4.3 Dimension-by-Dimension 的耦合

現有 3D 插值是 dimension-by-dimension: 先 η → 再 ξ → 最後 ζ。
WENO5 同樣可 dimension-by-dimension 執行，但有個微妙問題：

在 η 方向 WENO 後得到的中間值 (interpolation1order) 是否仍是「光滑的 point values」？
如果 η 方向有不連續性，WENO 會在 η 方向抑制它，但 ξ/ζ 方向的 stencil 內看到的是混合了 WENO 非線性效果的值。

**建議**: 首先實作 dimension-by-dimension WENO5 (最簡單)，後續可考慮 multi-dimensional WENO。

### 4.4 壁面附近的特殊處理

在 ζ 方向 (wall-normal)，壁面 k=3 和 k=NZ6-4 附近 stencil 可能不足 5 點。
現有 `bk_precomp_d[k]` 已處理 stencil base 的 wall clamping。

WENO5 需類似處理：
- 壁面附近若只有 3-4 個有效點 → fallback 到低階 WENO 或 Lagrange
- 或使用偏斜 stencil (one-sided WENO)

---

## 五、實作路線圖

### Phase 1: 1D WENO5 核心函式 (1-2 週)

1. 實作 `compute_linear_weights_weno5(sigma, d0, d1, d2)` — σ 依賴的 linear weights
2. 實作 `weno5_interp_1d(f[5], sigma)` — 完整 WENO5 一維插值
3. 驗證: 在均勻網格上比對 Lagrange-7 和 WENO5 的插值精度
4. 測試正化: 確認 Liu-Shu-Zhang splitting 在 d_r < 0 時正確工作

### Phase 2: 整合進 GILBM Streaming (2-3 週)

1. 修改 `gilbm_compute_point()` 中的 interpolation 迴圈
2. 將 7×7×7 stencil 降為 5×5×5 (或保持 7×7×7 但只取 5 點做 WENO)
3. 更新 Plan C fused kernel 的 local array 大小
4. 處理 ghost zone 與 MPI halo exchange 的寬度調整

### Phase 3: 壁面與邊界特殊處理 (1-2 週)

1. ζ 方向壁面附近的 WENO stencil 截斷處理
2. 與 Chapman-Enskog BC 的銜接
3. 周期邊界 (η, ξ 方向) 的 ghost zone WENO 一致性

### Phase 4: 效能最佳化 (2-3 週)

1. **Lookup table**: 將 d_r(σ) 預計算存入 CUDA constant memory 或 texture
2. **Warp divergence**: 正化分支的 branch-free 實現
3. **Register tuning**: 平衡 WENO 額外暫存器與 occupancy
4. **f_pc 記憶體**: 若從 343→125，全場省 ~2.1 GB (以 NX6=39, NYD6=23, NZ6=70 計)
5. Profiling: nvprof/nsight 比較 Lagrange-7 vs WENO5 的 kernel time

### Phase 5: 驗證與基準測試 (2-3 週)

1. Periodic Hill Re=1400 — 對比 Fröhlich DNS benchmark
2. 壁面摩擦係數 Cf 分佈
3. 分離泡再附著點位置
4. Reynolds stress profiles
5. 數值穩定性: 長時間模擬 (100+ FTT) 是否有非物理振盪

---

## 六、效能預估

### 計算量比較 (per direction, per grid point)

| 操作 | Lagrange-7 | WENO5 |
|------|-----------|-------|
| 多項式求值 | 7 乘加 | 3×(3 乘加) = 9 乘加 |
| 權重計算 | 7 個 Lagrange basis | 3 個 β + 3 個 α + 3 個 ω |
| 記憶體讀取 | 7 個 f | 5 個 f |
| Branch | 無 | 正化分支 (可 branch-free) |
| **總 FLOP** | ~14 | ~45 |
| **記憶體** | 7 reads | 5 reads |

WENO5 計算量 ≈ 3× Lagrange-7，但記憶體減少 29%。

在 GPU 上，GILBM streaming 是 **memory-bound** (f_pc 全場讀取佔主導)，因此：
- 記憶體減少 → 可能抵消計算增加
- 預估 kernel 整體慢 10-30%，但物理精度顯著提升

### 記憶體節省

```
f_pc 現有: 19 × 343 × 39 × 23 × 70 × 8 bytes = 3.274 GB (per rank)
f_pc WENO5: 19 × 125 × 39 × 23 × 70 × 8 bytes = 1.194 GB (per rank)
節省: 2.08 GB per rank → 8 GPUs 共省 16.6 GB
```

---

## 七、Hybrid 策略建議

考慮到 Periodic Hill 的流場特性（分離、再附著、回流），建議：

1. **光滑區 (bulk flow)**: WENO5 interpolation
2. **壁面區 (k ≤ 5 or k ≥ NZ6-7)**: 保持 Lagrange-7 或降階 Lagrange-3 (更穩定)
3. **分離點附近**: WENO5 發揮最大優勢 (抑制 Gibbs oscillation)

切換判據可用局部 β 值：若 max(β_r)/min(β_r) > threshold → 啟用 WENO5，否則用 Lagrange。

---

*本文檔為初步構思，將隨實作進展持續更新。*
*最後更新：2026-03-28*

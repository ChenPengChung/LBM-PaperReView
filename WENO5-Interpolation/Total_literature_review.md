# WENO5 Interpolation for Semi-Lagrangian / Off-Lattice LBM — 文獻回顧

> **建立日期**: 2026-03-28  
> **主題**: 將 WENO5 從傳統 FD-LBM 的 reconstruction 框架，擴展到 ISLBM / SLLBM / GILBM 等 interpolation-based LBM 中任意位置（非半格點）的插值問題  
> **關鍵字**: WENO interpolation, semi-Lagrangian LBM, departure point, arbitrary shift, linear weights, smoothness indicators, off-lattice streaming

---

## 目錄

1. [核心概念：Reconstruction vs. Interpolation](#一核心概念reconstruction-vs-interpolation)
2. [理論基礎文獻](#二理論基礎文獻必讀)
3. [Semi-Lagrangian + WENO Interpolation 文獻](#三weno-interpolation-應用於-semi-lagrangian-方法)
4. [直接與 LBM 結合的文獻](#四直接與-lbm-結合的文獻)
5. [WENO Interpolation at Arbitrary Point 專門文獻](#五weno-interpolation-at-arbitrary-point-專門文獻)
6. [數學框架總結](#六數學框架總結如何將-weno-引入任意位置的-interpolation)
7. [建議閱讀順序](#七建議閱讀順序)

---

## 一、核心概念：Reconstruction vs. Interpolation

| 項目 | WENO Reconstruction | WENO Interpolation |
|------|--------------------|--------------------|
| **已知量** | Cell averages $\bar{u}_i$ | Point values $u_i = u(x_i)$ |
| **目標** | 重建半格點 $x_{i+1/2}$ 的 point value | 在 **任意位置** $x^* = x_i + \sigma \Delta x$ 插值 |
| **Linear weights** | 固定常數 (如 $d_0=3/10, d_1=6/10, d_2=1/10$) | **隨 $\sigma$ 變化**，是目標點位置的函數 $d_r(\sigma)$ |
| **典型應用** | FD-LBM 通量重建 | ISLBM/SLLBM/GILBM 的 departure point 插值 |

**核心問題**：在 ISLBM/SLLBM/GILBM 中，streaming step 沿特徵線追蹤到的 departure point 通常不落在格點或半格點上，因此需要一套能在任意偏移量 $\sigma$ 處做高階非振盪插值的 WENO 框架。

---

## 二、理論基礎文獻（必讀）

### 2.1 C.-W. Shu (2020) — WENO 大型 Review
- **標題**: Essentially non-oscillatory and weighted essentially non-oscillatory schemes
- **期刊**: *Acta Numerica*, 2020
- **被引**: 258 次
- **重要性**: 明確區分 WENO reconstruction 和 WENO interpolation 兩套框架。說明如何用 WENO interpolation 在 point values 基礎上逼近任意位置的函數值和導數值，並指出 WENO interpolation 可自然處理非均勻網格。
- **DOI**: 10.1017/S0962492920000057

### 2.2 C.-W. Shu (2009) — 經典 SIAM Review
- **標題**: High order weighted essentially nonoscillatory schemes for convection dominated problems
- **期刊**: *SIAM Review*, 2009
- **被引**: 1243 次
- **重要性**: 提供了 WENO interpolation 的完整框架，說明如何選擇 linear weights 使其保持正性，以及如何將 interpolation 目標點從固定的 $x_{i+1/2}$ 推廣到任意位置。

### 2.3 Y. Liu, C.-W. Shu, M. Zhang (2009) — Linear Weights 正性問題
- **標題**: On the positivity of linear weights in WENO approximations
- **期刊**: *Acta Mathematicae Applicatae Sinica, English Series*, 2009
- **被引**: 42 次
- **重要性**: 專門解決 WENO interpolation 中 linear weights 可能為負的問題。當插值目標點不在半格點時，linear weights 隨位置變化，某些 stencil 組合下可能出現負權重。本文提供了通用的正化處理框架——這是將 WENO 推廣到任意 $\sigma$ 時必須面對的關鍵問題。

### 2.4 Y. Jiang, C.-W. Shu, M. Zhang (2013) — 基於 Interpolation 的替代 WENO 格式
- **標題**: An alternative formulation of finite difference weighted ENO schemes with Lax--Wendroff time discretization for conservation laws
- **期刊**: *SIAM J. Scientific Computing*, 2013
- **被引**: 181 次
- **重要性**: 提出基於 WENO interpolation（而非 reconstruction）的有限差分格式替代形式。詳細描述了如何從 point values 進行 WENO 插值，linear weights $d_r$ 如何根據目標位置調整。

### 2.5 K. Sebastian, C.-W. Shu (2003) — Multidomain WENO with Interpolation
- **標題**: Multidomain WENO finite difference method with interpolation at subdomain interfaces
- **期刊**: *J. Scientific Computing*, 2003
- **被引**: 109 次
- **重要性**: 明確展示了 WENO interpolation 的完整構造過程：從 point values 出發，在任意點做插值。公式可直接推廣到非半格點情形。

---

## 三、WENO Interpolation 應用於 Semi-Lagrangian 方法

> 這些文獻正是將 WENO interpolation 用在 departure point（類似 ISLBM/SLLBM/GILBM 需求）的核心參考。

### 3.1 J.A. Carrillo, F. Vecil (2007) — ⭐ 最關鍵參考
- **標題**: Nonoscillatory interpolation methods applied to Vlasov-based models
- **期刊**: *SIAM J. Scientific Computing*, 2007
- **被引**: 106 次
- **重要性**: 提出基於 WENO interpolation 的 semi-Lagrangian 方法。在 SL 框架下，departure point 不一定落在格點或半格點，因此必須建立能在任意位置做高階非振盪插值的 WENO 方案。**詳細推導了如何根據 departure point 與周圍格點的距離來計算 linear weights 和 smoothness indicators。**

### 3.2 J.-M. Qiu, C.-W. Shu (2011) — Conservative SL-FD-WENO
- **標題**: Conservative semi-Lagrangian finite difference WENO formulations with applications to the Vlasov equation
- **期刊**: *Commun. Comput. Phys.*, 2011
- **被引**: 68 次
- **重要性**: 提供了六階 WENO interpolation 的完整公式，用於 semi-Lagrangian 有限差分格式。核心：departure point 的位置會隨流場變化，WENO interpolation 的 linear weights 和 stencil 選擇都必須相應調整。

### 3.3 D. Yi, H. Kim (2015) — 直接以 SL + WENO Interpolation 為題
- **標題**: A semi-Lagrangian method based on WENO interpolation
- **期刊**: *J. Chungcheong Mathematical Society*, 2015
- **重要性**: 直接討論如何確定 departure point，並用 WENO interpolation 在該點進行值的恢復。包含 smoothness indicator $\beta_r$ 在非半格點位置的定義。

### 3.4 A. Christlieb, M. Link, H. Yang, R. Chang (2023) — 非多項式 WENO
- **標題**: High-order semi-Lagrangian WENO schemes based on non-polynomial space for the Vlasov equation
- **期刊**: *Commun. Appl. Math. Comput.*, 2023
- **被引**: 1 次
- **重要性**: 將 WENO interpolation 推廣到非多項式空間（如指數函數基），用於 semi-Lagrangian Vlasov solver。提供了 generalized WENO interpolation 的新框架。

### 3.5 G.-H. Cottet (2023) — WENO Semi-Lagrangian Particle Methods
- **標題**: WENO semi-Lagrangian particle methods
- **來源**: HAL preprint, 2023
- **重要性**: 設計了 WENO 擴展的 semi-Lagrangian 粒子方法，其中 interpolation 位置取決於速度場和時間推進格式，提供了另一種處理任意位置 WENO interpolation 的視角。

### 3.6 Z. Sun, F. Xiao (2017) — Semi-Lagrangian Multi-Moment + WENO
- **標題**: A semi-Lagrangian multi-moment finite volume method with fourth-order WENO projection
- **期刊**: *Int. J. Numer. Methods Fluids*, 2017
- **被引**: 11 次
- **重要性**: 結合 point value (PV) 和 volume-integrated average (VIA) 資訊在 WENO interpolation 中的應用，departure point 透過 semi-Lagrangian 方法確定。

### 3.7 J.-M. Qiu, G. Russo (2017) — 高階多維 Characteristic Tracing
- **標題**: A high order multi-dimensional characteristic tracing strategy for the Vlasov-Poisson system
- **期刊**: *J. Scientific Computing*, 2017
- **被引**: 27 次
- **重要性**: 提供了六階 WENO interpolation 在多維 characteristic tracing 中的完整公式。

---

## 四、直接與 LBM 結合的文獻

### 4.1 D. Wilde, A. Krämer, D. Reith, H. Foysi (2021) — ⭐ SLLBM 可壓縮湍流
- **標題**: High-order semi-Lagrangian kinetic scheme for compressible turbulence
- **期刊**: *Physical Review E*, 2021
- **被引**: 27 次
- **重要性**: 三維 SLLBM 用於可壓縮湍流（Taylor-Green Vortex）。討論了多種 interpolation 策略的選擇，包括與 WENO/TENO 等方案的比較。展示了 interpolation 引入的數值耗散本身可以起到穩定作用。

### 4.2 D. Wilde, A. Krämer, D. Reith, H. Foysi (2020) — SLLBM 基礎
- **標題**: Semi-Lagrangian lattice Boltzmann method for compressible flows
- **期刊**: *Physical Review E*, 2020
- **被引**: 55 次
- **重要性**: SLLBM 的基礎文獻。使用 finite element interpolation 恢復 off-lattice departure points 的分布函數值，是後續高階 interpolation 方法的起點。

### 4.3 M.H. Saadat, F. Bösch, I.V. Karlin (2020) — 非結構網格 SLLBM
- **標題**: Semi-Lagrangian lattice Boltzmann model for compressible flows on unstructured meshes
- **期刊**: *Physical Review E*, 2020
- **被引**: 33 次
- **重要性**: 在任意非均勻非結構網格上進行 semi-Lagrangian propagation，interpolation 必須處理非均勻間距下的 departure point 問題。

### 4.4 N.G. Kallikounis, B. Dorschner, I.V. Karlin (2021) — Multiscale SLLBM
- **標題**: Multiscale semi-Lagrangian lattice Boltzmann method
- **期刊**: *Physical Review E*, 2021
- **被引**: 18 次
- **重要性**: 多尺度 semi-Lagrangian LBM，使用多個 collocation points 做 interpolation，lifting operation 中涉及非標準位置的插值。

### 4.5 E. Reyhanian (2021) — ETH 博士論文（PonD Framework）
- **標題**: Thermokinetic model for compressible generic fluids
- **來源**: ETH Zurich 博士論文, 2021
- **被引**: 5 次
- **重要性**: 詳細討論 semi-Lagrangian LBM（PonD framework）中的 advection 步驟，包括 WENO 方法和 limiters 如何在 semi-Lagrangian advection 中調整 interpolation weights。

### 4.6 N.G. Kallikounis, B. Dorschner, I.V. Karlin (2022) — Particles on Demand
- **標題**: Particles on demand for flows with strong discontinuities
- **期刊**: *Physical Review E*, 2022
- **被引**: 23 次
- **重要性**: 提出 semi-Lagrangian 和 finite-volume 兩種實現方式，用於處理強不連續性的可壓縮流，討論了 off-lattice propagation 的穩定性。

### 4.7 V.E. Ambruș et al. (2025) — Vielbein LBM on Spherical Surfaces
- **標題**: Vielbein Lattice Boltzmann approach for fluid flows on spherical surfaces
- **期刊**: *Physical Review E*, 2025
- **重要性**: 使用 WENO-5 scheme 對 Boltzmann equation 的 advection 部分做處理，針對球面等非笛卡爾幾何下的 off-lattice streaming。

### 4.8 W. Noh, C. Lee (2025) — Lagrangian Entropic LBM
- **標題**: Lagrangian entropic lattice Boltzmann method for Courant-free supersonic compressible flow simulation
- **來源**: arXiv:2508.06911, 2025
- **重要性**: 提出 Lagrangian entropic LBM (LELBM)，屬 interpolation-based methods，stream off-lattice population，與 GILBM 概念高度相關。

---

## 五、WENO Interpolation at Arbitrary Point 專門文獻

### 5.1 G. Janett et al. (2019) — ⭐ 最推薦的 WENO Interpolation 文獻
- **標題**: A novel fourth-order WENO interpolation technique. A possible new tool designed for radiative transfer
- **期刊**: *Astronomy & Astrophysics*, 2019
- **被引**: 18 次
- **重要性**: **專門討論 WENO interpolation（不是 reconstruction）**，明確處理 off-grid points 的插值問題。文中指出 reconstruction of cell averages 等價於 primitive function 的 point value interpolation，但直接在 point values 上做 WENO interpolation 需要重新推導 linear weights——**這些權重是目標點位置的函數**。提出了一種新的四階 WENO interpolation 技術。

### 5.2 Z. Wang, J. Zhu, L. Tian, N. Zhao (2024) — Interpolation-based vs Reconstruction-based WENO
- **標題**: Assessment of high-order interpolation-based weighted essentially non-oscillatory schemes for compressible Taylor-Green vortex flows
- **期刊**: *Physics of Fluids*, 2024
- **被引**: 2 次
- **重要性**: **直接對比 interpolation-based WENO 和 reconstruction-based WENO** 在可壓縮湍流 (TGV) 中的表現，量化了 interpolation-based 方案的數值耗散特性。

### 5.3 P. Li, W.S. Don, Z. Gao (2020) — WENO Interpolation-Based Schemes
- **標題**: High order well-balanced finite difference WENO interpolation-based schemes for shallow water equations
- **期刊**: *Computers & Fluids*, 2020
- **被引**: 37 次
- **重要性**: 提出基於 WENO interpolation（而非 reconstruction）的有限差分格式，展示了 interpolation-based 方案在保持 well-balanced 性質方面的優勢。

### 5.4 Y. Ha, C.H. Kim, H. Yang, J. Yoon (2024) — 替代 WENO + 指數多項式
- **標題**: A new alternative WENO scheme based on exponential polynomial interpolation with an improved order of accuracy
- **期刊**: *J. Scientific Computing*, 2024
- **重要性**: 提出基於指數多項式插值的新 WENO interpolation 方案，在給定 stencil 上能達到更高的精度階數。

---

## 六、數學框架總結：如何將 WENO 引入任意位置的 Interpolation

### Step 1: 問題定義
已知均勻格點上的 point values $f_i = f(x_i)$，要估計任意位置的值：
$$x^* = x_i + \sigma \Delta x, \quad \sigma \in [0, 1] \text{ (不限於 } \sigma = 1/2 \text{)}$$

### Step 2: 候選 Stencil 上的 Interpolation Polynomial
在每個候選 stencil $S_r = \{x_{i-r}, x_{i-r+1}, \ldots, x_{i-r+k}\}$ 上，建立 Lagrange interpolation polynomial $p_r(x)$，使其通過 stencil 上所有格點的 point values。

### Step 3: Linear Weights 隨 sigma 變化
大 stencil 上的高階多項式 $P(x^*)$ 必須等於各小 stencil 加權和：
$$P(x^*) = \sum_r d_r(\sigma) \cdot p_r(x^*)$$

**Linear weights $d_r(\sigma)$ 不再是常數**，而是目標位置偏移量 $\sigma$ 的函數。

> 例如，對於五階 WENO (WENO5)，三個候選 stencil 的 linear weights 在 $\sigma = 1/2$ 時為經典值 $(1/10, 6/10, 3/10)$，但在其他 $\sigma$ 值時會不同，且可能出現負值。

### Step 4: Smoothness Indicators
Smoothness indicators 的定義與 reconstruction 版本相同：
$$\beta_r = \sum_{l=1}^{k} \int_{x_{i-1/2}}^{x_{i+1/2}} (\Delta x)^{2l-1} \left( \frac{d^l p_r(x)}{dx^l} \right)^2 dx$$

**$\beta_r$ 不隨 $\sigma$ 變化**（它衡量的是 stencil 上多項式的光滑程度，與插值目標點無關）。

### Step 5: Non-linear Weights
$$\omega_r = \frac{\alpha_r}{\sum_s \alpha_s}, \quad \alpha_r = \frac{d_r(\sigma)}{(\epsilon + \beta_r)^p}$$

> **注意**：當 $d_r(\sigma) < 0$ 時，需使用 Liu-Shu-Zhang (2009) 的正化技巧——將負權重拆分為正、負兩組分別做 WENO 加權，最後相減。

### Step 6: 最終插值
$$f(x^*) \approx \sum_r \omega_r \cdot p_r(x^*)$$

---

## 七、建議閱讀順序

| 順序 | 文獻 | 目的 |
|------|------|------|
| 1 | Shu (2020) Acta Numerica §2.5 | 理解 WENO interpolation vs reconstruction 的數學基礎 |
| 2 | Jiang, Shu, Zhang (2013) SIAM | 看 interpolation-based WENO FD 的完整推導 |
| 3 | Liu, Shu, Zhang (2009) | 理解 linear weights 正性問題及處理方法 |
| 4 | Carrillo & Vecil (2007) SIAM | 學習 SL + WENO interpolation 如何處理 departure point |
| 5 | Janett et al. (2019) A&A | 看專門為任意位置設計的 WENO interpolation |
| 6 | Qiu & Shu (2011) CiCP | 六階 WENO interpolation 在 SL-FD 中的完整公式 |
| 7 | Wang et al. (2024) PoF | 對比 interpolation-based 和 reconstruction-based WENO |
| 8 | Wilde et al. (2021) PRE | 看 SLLBM 中高階 interpolation 在可壓縮湍流的實際應用 |

---

## 附錄：相關 LBM 變體縮寫

| 縮寫 | 全稱 | 特點 |
|------|------|------|
| **FD-LBM** | Finite Difference LBM | 用 FD 離散空間導數，WENO reconstruction 用於通量 |
| **ISLBM** | Interpolation-Supplemented LBM | 用插值補充非格點位置的 streaming |
| **SLLBM** | Semi-Lagrangian LBM | 沿特徵線追蹤 departure point，用插值恢復分布函數 |
| **GILBM** | Grid-Interpolation LBM / Generalized Interpolation LBM | 通用 interpolation-based LBM 框架 |
| **PonD** | Particles on Demand | ETH Karlin 組的 off-lattice LBM 框架 |

---

*本文檔由文獻搜尋自動產生，最後更新：2026-03-28*

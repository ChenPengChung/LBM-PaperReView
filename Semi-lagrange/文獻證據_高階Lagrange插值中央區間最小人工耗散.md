# 文獻證據：高階 Lagrange 插值在 central region 使人工耗散最小

> **論點（claim）**：semi-Lagrangian 平流格式中，高階 Lagrange interpolation 的插值 stencil 以落點（departure point）為中心（即落點位於 stencil 的 central interval / central region）時，人工耗散（artificial dissipation / damping / amplitude error）最小。
>
> **驗證方式**：deep-research workflow（101 agents、19 來源、94 條候選引句取前 25 條驗證），每條引句由 3 個獨立驗證者對照**出版社原版 PDF** 投票，25/25 全數 3–0 通過、0 條駁回。標 ✅ 者另經人工逐頁核對頁面影像。
>
> 產生日期：2026-07-13

---

## 結論（直接回答）

1. **沒有任何一篇文獻用 "minimized" 一詞把完整論點寫成單句**。引用時建議用組合引用（見「建議論述」）。
2. 最貼合的單句是 **Bates & McDonald (1982) p. 1834**：限制落點在中央區間後 "complete extinction never occurs and the damping is, in all cases, much less than that given by linear interpolation."
3. centered vs off-centered 的**直接比較句**在 **McDonald (1984) p. 1269**。
4. 「階數越高耗散越低」有多篇直述句，但 **McDonald (1984) 指出此關係非單調**（居中的 quadratic 比 cubic 耗散低）——**居中程度與階數同等重要**，這其實強化本論點。

---

## Tier A：直接命中（centered / central region ⇔ damping）

### A1. Bates & McDonald (1982) — 本資料夾 `4.` ✅

Bates, J. R., and A. McDonald, 1982: Multiply-Upstream, Semi-Lagrangian Advective Schemes: Analysis and Application to a Multi-Level Primitive Equation Model. *Mon. Wea. Rev.*, **110**(12), 1831–1842. DOI: 10.1175/1520-0493(1982)110<1831:MUSLAS>2.0.CO;2

**p. 1833, §2b（centered stencil 定義句）**：
> "The gridpoint nearest the departure point x\* is chosen as the central point of the three interpolation points (I − p + 1), (I − p), (I − p − 1). Thus x\* lies within a half grid interval from (I − p)."

**p. 1834, §2b（最強單句：中央區間 ⇒ 避免完全消滅、耗散大減）**：
> "But our interpolation points are chosen such that −0.5 ≤ α̂ ≤ 0.5. (12) Thus, this scheme is again unconditionally stable. In Fig. 5, |λ|² is plotted as a function of α̂ for various wavelengths L. Since we choose α̂ to lie in the interval defined by (12), we see that complete extinction never occurs and the damping is, in all cases, much less than that given by linear interpolation."

Fig. 5（同頁）：|λ|²–α̂ 曲線顯示 L=2Δx 波在 α̂=±1/√2≈±0.707（中央區間之外）完全消滅（|λ|²=0）；中央點 α̂=0 處 |λ|²=1（零 damping）。

**p. 1836, §2d（2D 版）**：
> "The grid-point nearest the departure point is chosen as the center of the nine points used in the interpolation... As in case (b), complete extinction never occurs and the damping is in all cases much less than that given by bilinear interpolation."

**貼合度**：直接。centered 設計＋耗散結論同段出現；未與同階 off-centered 比較（該比較見 A2）。

### A2. McDonald (1984) — 本資料夾 `5.` ✅

McDonald, A., 1984: Accuracy of Multiply-Upstream, Semi-Lagrangian Advective Schemes. *Mon. Wea. Rev.*, **112**(6), 1267–1275. DOI: 10.1175/1520-0493(1984)112<1267:AOMUSL>2.0.CO;2

**p. 1269, §2c（centered vs off-centered 直接比較）**：
> "Notice that the strictures imposed on α by Eq. (8) mean that the quadratic and quartic cases have amplitudes which are symmetric about α = 0.5. One consequence is a better amplitude representation than if, for instance, i had been chosen such that (i − 1)Δx ≤ x\* ≤ iΔx. Another is that the quadratic amplitude is less damped than the cubic for all wavelengths and all α."

（Eq. (8)（p. 1268）即 (i−½)Δx < x\* ≤ (i+½)Δx 的居中限制；Eqs. (21)–(24) 顯示 amplitude error 為純 damping，故 "better amplitude representation" = 較少 damping，屬原文依據、非推論。）

**p. 1269, §2c（damping 對落點位置的依賴）**：
> "All four schemes are most severely damped when α = 0.5 and undamped when α = 0.0 or 1.0. The damping is heaviest for the shortest resolvable wavelengths and it decreases as the wavelength increases. The linear scheme is the most damped followed by the cubic, quadratic, and quartic."

**p. 1269 總結第 3 條（格點上零誤差）**：
> "3) close to α = 0, the phase error is largest and the amplitude error is least, whereas close to all integer values of α the phase and amplitude errors are zero,"

**p. 1268（居中 ⇒ 無條件穩定）**：
> "It will be shown in Section 2c that these conditions are sufficient to ensure unconditional stability for all four interpolation schemes."

**貼合度**：直接（與一個 off-centered 範例比較，非對所有位移的最優性定理；未用 "minimized" 一詞）。

### A3. Besse & Mehrenberger (2008) — 本資料夾 `2.` ✅

Besse, N., and M. Mehrenberger, 2008: Convergence of classes of high-order semi-Lagrangian schemes for the Vlasov–Poisson system. *Math. Comp.*, **77**(261), 93–123. DOI: 10.1090/S0025-5718-07-01912-6

**pp. 93–94, §1（symmetric = centered Lagrange ⇒ 穩定＋弱耗散）**：
> "The Hermite type interpolation leads to high-order, stable and slightly diffusive schemes, whereas the Lagrange interpolation leads to convergent but too diffusive schemes for orders smaller than two and unstable schemes in any L^p norm for orders greater than two. Nevertheless, on uniform grids, if we use symmetric Lagrange interpolation we can then recover the stability in the L2 discrete norm and keep the high-order and weak diffusion features of our approximation."

**p. 106, §5（階數↑ ⇒ 耗散↓ 的量化句）**：
> "The Lagrange interpolation has the advantage of being local, so that it is well-suited for parallel computing. Nevertheless, it is more diffusive than B-splines. For example, in order to get the same rate of diffusion of the cubic B-splines a 9th order Lagrange interpolation must be used. This has been illustrated numerically in [15]."

（[15] = Filbet, Sonnendrücker & Bertrand, *J. Comput. Phys.* **172** (2001), 166–187。）

**貼合度**：直接（symmetric/centered ⇔ weak diffusion 連結明確；但未用 damping/minimized，對照對象是非結構網格 Lagrange 插值）。

---

## Tier B：高階 ⇒ 低耗散（sub-claim b）

### B1. Staniforth & Côté (1991) 綜述 — 本資料夾 `6.` ✅

Staniforth, A., and J. Côté, 1991: Semi-Lagrangian Integration Schemes for Atmospheric Models—A Review. *Mon. Wea. Rev.*, **119**(9), 2206–2223. DOI: 10.1175/1520-0493(1991)119<2206:SLISFA>2.0.CO;2

**p. 2210, §2d "Interpolation"**：
> "Cubic interpolation gives fourth-order spatial truncation errors with very little damping (it is very scale selective, affecting primarily the smallest scales), whereas linear interpolation (see McDonald 1984 for discussion) has unacceptably large damping (it is also scale selective, but has a much less sharp response)."

**p. 2212, §2e**：
> "The aforementioned stability analyses show that semi-Lagrangian advection schemes have very good phase speeds with little numerical dispersion, but contrary to some Eulerian schemes (e.g., leapfrog-based schemes) there is some damping due to interpolation as discussed in section 2d. This damping is fortunately very scale selective (at least when using high-order interpolators)."

（同頁提到 **McCalpin (1988)** 對 SL 插值耗散與 Laplacian/biharmonic 耗散的定量比較——若需更多耗散專文可追。）

### B2. Durran 教科書 — 本資料夾 `1.` ✅

Durran, D. R., 2010: *Numerical Methods for Fluid Dynamics: With Applications to Geophysics*, 2nd ed., Springer, **p. 362, §7.1.1.3 "Higher-Order Interpolation"**：
> "Upstream differencing generates too much numerical diffusion to be useful in practical computations involving Eulerian problems with smooth solutions. A similar situation holds in the Lagrangian framework, where linearly interpolating the tracer field also generates too much diffusion. Higher-order interpolation is therefore used in most semi-Lagrangian approximations to equations with smooth solutions."

**p. 360, §7.1.1.1（線性插值放大因子精確式）**：|A_k|² = 1 − 2α(1−α)(1−cos kΔx)

**p. 384, §7.5.1（插值 ⇒ damping；non-interpolating 格式的動機）**：
> "Some damping is produced in all the previously described semi-Lagrangian schemes when the prognostic fields are interpolated to the departure point."

### B3. Filbet & Sonnendrücker (2003) — 本資料夾 `8.`

Filbet, F., and E. Sonnendrücker, 2003: Comparison of Eulerian Vlasov solvers. *Comput. Phys. Comm.*, **150**, 247–266. DOI: 10.1016/S0010-4655(02)00694-X

**p. 257, §6.1（放大因子排序）**：
> "Let us first consider the amplification factor for the different methods (see Fig. 1). We observe that methods using a smooth reconstruction (Hermite or spline) are less dissipative than those using only a continuous interpolation. To obtain a similar amplification factor with the Lagrange interpolation as with the spline interpolation, a polynomial of degree nine is required. The dissipation of the conservative method with a quadratic polynomial is identical to the one using cubic Lagrange interpolation. The linear reconstruction used in the (FBM) is the most dissipative."

**pp. 252–253（centered 為設計原則）**：
> "We only choose polynomials of odd degree to have a centered approximation..."
> "...for linear advection with constant coefficients, the use of a centered approximation ensures the conservation of global mass"

### B4. Riishøjgaard et al. (1998) — 本資料夾 `7.`

Riishøjgaard, L. P., S. E. Cohn, Y. Li, and R. Ménard, 1998: The Use of Spline Interpolation in Semi-Lagrangian Transport Models. *Mon. Wea. Rev.*, **126**(7), 2008–2016.

**p. 2012, §3（cubic Lagrange centered 4 點 stencil 的區間內耗散分布）**：
> "Figure 4a shows that the scheme is most diffusive for C = 0.5, where the 2Δx wave is completely eliminated within one time step. The algorithm is moderately diffusive for C close to either 0 or 1. Note that the amplification factor for a given C is equal to that for 1 − C..."

---

## Tier C：數學文獻的「穩定性對應版」（相鄰論述，不用 damping 一詞）

### C1. Falcone & Ferretti (1998) — 本資料夾 `3.` ✅

Falcone, M., and R. Ferretti, 1998: Convergence Analysis for a Class of High-Order Semi-Lagrangian Advection Schemes. *SIAM J. Numer. Anal.*, **35**(3), 909–940.

**p. 925, §5.2**：
> "Therefore, interpolations of order 1 or 2 are L² stable. Interpolations of order 3 and higher are stable only provided z_j ∈ [r/2−1, r/2+1] if r is even, or z_j ∈ [s, s+1] if r is odd."

（即：高階 Lagrange 插值的 von Neumann 條件**只在中央區帶成立**——central region 的穩定性版本。）

### C2. Charles, Després & Mehrenberger (2013) — 本資料夾 `9.`（HAL preprint）

Charles, F., B. Després, and M. Mehrenberger, 2013: Enhanced convergence estimates for semi-Lagrangian schemes: Application to the Vlasov–Poisson equation. *SIAM J. Numer. Anal.*, **51**(2), 840–863. DOI: 10.1137/110851511

**Prop. 7（HAL preprint p. 11, §4.2）**：
> "Proposition 7. The amplification factor satisfies |λ_{ν,k,p}(ψ)| ≤ 1, ∀ψ ∈ ℝ if and only if p ∈ {2k, 2k + 1, 2k + 2}."

（(p+1) 點 Lagrange SL 格式 |放大因子|≤1 ⟺ stencil 居中（p=2k+1）或近居中（p=2k, 2k+2）。全文不使用 damping/dissipation 詞彙。）

### C3. Ferretti & Mehrenberger (2020) — 本資料夾 `10.`（HAL preprint）

Ferretti, R., and M. Mehrenberger, 2020: Stability of semi-Lagrangian schemes of arbitrary odd degree under constant and variable advection speed. *Math. Comp.*, **89**, 1783–1805. DOI: 10.1090/mcom/3494

任意奇數階 **symmetric（centered）Lagrange** 插值 SL 格式的 L² 穩定性證明。全文檢索確認**零次**出現 damping/dissipation——僅以穩定性語言表述。

---

## 建議論述（可直接放入論文，含引用組合）

**英文版**：

> High-order Lagrange interpolation with the stencil centered about the departure point is well documented to keep the interpolation-induced artificial dissipation small. Restricting the departure point to the central interval of the stencil ensures that "complete extinction never occurs and the damping is, in all cases, much less than that given by linear interpolation" (Bates and McDonald 1982, p. 1834), and yields "a better amplitude representation" than an off-centered placement of the same stencil (McDonald 1984, p. 1269); the amplitude error vanishes when the departure point coincides with a grid node and is largest midway between nodes (McDonald 1984). In the rigorous setting, symmetric Lagrange interpolation on uniform grids retains "the high-order and weak diffusion features" of the scheme while guaranteeing discrete L² stability (Besse and Mehrenberger 2008, pp. 93–94), the von Neumann condition for interpolation of order ≥ 3 holding precisely when the departure point lies in the central strip of the stencil (Falcone and Ferretti 1998, p. 925; Charles, Després and Mehrenberger 2013, Prop. 7).

**中文版**：

> 高階 Lagrange 插值以落點為中心取 stencil 可將插值引致的人工耗散壓到很低：限制落點於 stencil 中央區間可保證「完全消滅永不發生，且 damping 在所有情形下都遠小於線性插值」（Bates & McDonald 1982, p. 1834），並比 off-centered 取點得到「更好的振幅表現」（McDonald 1984, p. 1269）；振幅誤差在落點與格點重合時為零、在格點中間時最大（McDonald 1984）。在嚴格數學層面，均勻網格上的 symmetric Lagrange 插值在保證離散 L² 穩定的同時保留「高階與弱耗散特性」（Besse & Mehrenberger 2008, pp. 93–94），且三階以上插值的 von Neumann 條件恰在落點位於 stencil 中央區帶時成立（Falcone & Ferretti 1998, p. 925；Charles, Després & Mehrenberger 2013, Prop. 7）。

---

## 引用時務必注意（誠實限制）

1. **"minimized" 一詞無文獻單句**。McDonald (1984) 的比較是 centered vs「一個」off-centered 範例（"for instance..."），**不是**對所有 stencil 位移的最優性定理；目前文獻中不存在該最優性證明（open question）。
2. **階數–耗散非單調**：McDonald (1984, p. 1269)："The linear scheme is the most damped followed by the cubic, quadratic, and quartic."——居中的 quadratic 比（區間式取點的）cubic 耗散更低。寫「higher order ⇒ less dissipation」須加限定（如 "for centered stencils" 或引成對比較）。
3. **Filbet–Sonnendrücker 引用要分清**：可引的放大因子原句在 **2003 CPC 150** 論文（本資料夾 `8.`）；Besse & Mehrenberger 的 [15] 指 **2001 JCP 172** 論文（數值示例出處）。
4. `9.` 與 `10.` 為 **HAL preprint**（非期刊排版終稿），引用頁碼前請對照期刊版。
5. 未查驗的線索：Falcone & Ferretti 2013 SIAM 專書、McDonald (1987)、Hong & Steinberg (2001)、Sonnendrücker et al. (1999)、Purser & Leslie (1988/1991)、**McCalpin (1988)**（SL 耗散定量分析專文，最值得追）。

---

## BibTeX

```bibtex
@article{BatesMcDonald1982,
  author  = {Bates, J. R. and McDonald, A.},
  title   = {Multiply-Upstream, Semi-Lagrangian Advective Schemes: Analysis and Application to a Multi-Level Primitive Equation Model},
  journal = {Monthly Weather Review},
  volume  = {110},
  number  = {12},
  pages   = {1831--1842},
  year    = {1982},
  doi     = {10.1175/1520-0493(1982)110<1831:MUSLAS>2.0.CO;2}
}

@article{McDonald1984,
  author  = {McDonald, A.},
  title   = {Accuracy of Multiply-Upstream, Semi-Lagrangian Advective Schemes},
  journal = {Monthly Weather Review},
  volume  = {112},
  number  = {6},
  pages   = {1267--1275},
  year    = {1984},
  doi     = {10.1175/1520-0493(1984)112<1267:AOMUSL>2.0.CO;2}
}

@article{StaniforthCote1991,
  author  = {Staniforth, Andrew and C{\^o}t{\'e}, Jean},
  title   = {Semi-Lagrangian Integration Schemes for Atmospheric Models---A Review},
  journal = {Monthly Weather Review},
  volume  = {119},
  number  = {9},
  pages   = {2206--2223},
  year    = {1991},
  doi     = {10.1175/1520-0493(1991)119<2206:SLISFA>2.0.CO;2}
}

@article{Riishojgaard1998,
  author  = {Riish{\o}jgaard, Lars Peter and Cohn, Stephen E. and Li, Yong and M{\'e}nard, Richard},
  title   = {The Use of Spline Interpolation in Semi-Lagrangian Transport Models},
  journal = {Monthly Weather Review},
  volume  = {126},
  number  = {7},
  pages   = {2008--2016},
  year    = {1998},
  doi     = {10.1175/1520-0493(1998)126<2008:TUOSII>2.0.CO;2}
}

@article{BesseMehrenberger2008,
  author  = {Besse, Nicolas and Mehrenberger, Michel},
  title   = {Convergence of classes of high-order semi-Lagrangian schemes for the {V}lasov--{P}oisson system},
  journal = {Mathematics of Computation},
  volume  = {77},
  number  = {261},
  pages   = {93--123},
  year    = {2008},
  doi     = {10.1090/S0025-5718-07-01912-6}
}

@article{FalconeFerretti1998,
  author  = {Falcone, Maurizio and Ferretti, Roberto},
  title   = {Convergence Analysis for a Class of High-Order Semi-Lagrangian Advection Schemes},
  journal = {SIAM Journal on Numerical Analysis},
  volume  = {35},
  number  = {3},
  pages   = {909--940},
  year    = {1998},
  doi     = {10.1137/S0036142994273513}
}

@article{FilbetSonnendrucker2003,
  author  = {Filbet, Francis and Sonnendr{\"u}cker, Eric},
  title   = {Comparison of {E}ulerian {V}lasov solvers},
  journal = {Computer Physics Communications},
  volume  = {150},
  number  = {3},
  pages   = {247--266},
  year    = {2003},
  doi     = {10.1016/S0010-4655(02)00694-X}
}

@article{CharlesDespresMehrenberger2013,
  author  = {Charles, Fr{\'e}d{\'e}rique and Despr{\'e}s, Bruno and Mehrenberger, Michel},
  title   = {Enhanced convergence estimates for semi-Lagrangian schemes: Application to the {V}lasov--{P}oisson equation},
  journal = {SIAM Journal on Numerical Analysis},
  volume  = {51},
  number  = {2},
  pages   = {840--863},
  year    = {2013},
  doi     = {10.1137/110851511}
}

@article{FerrettiMehrenberger2020,
  author  = {Ferretti, Roberto and Mehrenberger, Michel},
  title   = {Stability of semi-Lagrangian schemes of arbitrary odd degree under constant and variable advection speed},
  journal = {Mathematics of Computation},
  volume  = {89},
  number  = {324},
  pages   = {1783--1805},
  year    = {2020},
  doi     = {10.1090/mcom/3494}
}

@book{Durran2010,
  author    = {Durran, Dale R.},
  title     = {Numerical Methods for Fluid Dynamics: With Applications to Geophysics},
  edition   = {2},
  series    = {Texts in Applied Mathematics},
  volume    = {32},
  publisher = {Springer},
  address   = {New York},
  year      = {2010},
  isbn      = {978-1-4419-6411-3}
}
```

---

## 資料夾對照表

| 編號 | 檔案 | 角色 |
|---|---|---|
| 1. | Durran 教科書 | B2：高階插值動機、放大因子公式、non-interpolating 動機 |
| 2. | Besse & Mehrenberger (2008) | A3：symmetric Lagrange ⇒ 穩定＋弱耗散；9 階 ≈ cubic B-spline |
| 3. | Falcone & Ferretti (1998) | C1：中央區帶穩定性條件 |
| 4. | Bates & McDonald (1982) | A1：**最強單句**（中央區間 ⇒ 無完全消滅、耗散大減） |
| 5. | McDonald (1984) | A2：centered vs off-centered 直接比較；耗散排序 |
| 6. | Staniforth & Côté (1991) | B1：cubic "very little damping" 綜述句 |
| 7. | Riishøjgaard et al. (1998) | B4：區間內耗散分布（C=0.5 最耗散） |
| 8. | Filbet & Sonnendrücker (2003) | B3：放大因子排序、degree nine 句 |
| 9. | Charles, Després & Mehrenberger (2013, HAL) | C2：Prop. 7（|λ|≤1 ⟺ 居中） |
| 10. | Ferretti & Mehrenberger (2020, HAL) | C3：任意奇數階 symmetric 插值 L² 穩定 |

# 文獻庫總索引 — WENO5 Interpolation for GILBM

> **建立日期**: 2026-03-28
> **對應 README**: `../README.md`
> **總文獻數**: 24 篇

---

## 資料夾結構

```
LiteratureLibrary/
├── 01_Theoretical_Foundations/          (5 篇) — WENO 理論基礎
│   ├── 2.1_Shu_2020_WENO_Review.md
│   ├── 2.2_Shu_2009_SIAM_Review.md
│   ├── 2.3_Liu_Shu_Zhang_2009_Positivity.md
│   ├── 2.4_Jiang_Shu_Zhang_2013_Alternative_WENO.md
│   └── 2.5_Sebastian_Shu_2003_Multidomain.md
│
├── 02_SemiLagrangian_WENO_Interpolation/ (7 篇) — SL + WENO 插值方法
│   ├── 3.1_Carrillo_Vecil_2007_SL_WENO.md          ⭐ 最關鍵參考
│   ├── 3.2_Qiu_Shu_2011_Conservative_SL.md
│   ├── 3.3_Yi_Kim_2015_SL_WENO_Interp.md
│   ├── 3.4_Christlieb_2023_NonPoly_WENO.md
│   ├── 3.5_Cottet_2023_WENO_SL_Particle.md
│   ├── 3.6_Sun_Xiao_2017_SL_MultiMoment.md
│   └── 3.7_Qiu_Russo_2017_MultiDim_Tracing.md
│
├── 03_LBM_Integration/                  (8 篇) — 與 LBM 直接結合
│   ├── 4.1_Wilde_2021_SLLBM_Turbulence.md          ⭐ SLLBM 可壓縮湍流
│   ├── 4.2_Wilde_2020_SLLBM_Foundation.md
│   ├── 4.3_Saadat_2020_Unstructured_SLLBM.md
│   ├── 4.4_Kallikounis_2021_Multiscale_SLLBM.md
│   ├── 4.5_Reyhanian_2021_ETH_PonD.md
│   ├── 4.6_Kallikounis_2022_PonD_Discontinuities.md
│   ├── 4.7_Ambrus_2025_Vielbein_LBM.md
│   └── 4.8_Noh_Lee_2025_Lagrangian_ELBM.md
│
├── 04_WENO_Arbitrary_Point/             (4 篇) — 任意位置 WENO 插值專論
│   ├── 5.1_Janett_2019_WENO_Interpolation.md        ⭐ 最推薦 WENO 插值文獻
│   ├── 5.2_Wang_2024_Interp_vs_Recon.md
│   ├── 5.3_Li_2020_WENO_Interp_ShallowWater.md
│   └── 5.4_Ha_2024_Alternative_WENO_ExpPoly.md
│
└── 05_Summary/
    ├── INDEX.md                          ← 本檔案
    └── OPTIMIZATION_ANALYSIS.md          ← 程式碼最佳化分析
```

---

## 文獻對照表

| # | 編號 | 作者 | 年份 | 類別 | 關鍵性 | .md 檔案 |
|---|------|------|------|------|--------|----------|
| 1 | 2.1 | Shu | 2020 | 理論 | ★★★ | `01_Theoretical_Foundations/2.1_Shu_2020_WENO_Review.md` |
| 2 | 2.2 | Shu | 2009 | 理論 | ★★★ | `01_Theoretical_Foundations/2.2_Shu_2009_SIAM_Review.md` |
| 3 | 2.3 | Liu, Shu, Zhang | 2009 | 理論 | ★★★ | `01_Theoretical_Foundations/2.3_Liu_Shu_Zhang_2009_Positivity.md` |
| 4 | 2.4 | Jiang, Shu, Zhang | 2013 | 理論 | ★★☆ | `01_Theoretical_Foundations/2.4_Jiang_Shu_Zhang_2013_Alternative_WENO.md` |
| 5 | 2.5 | Sebastian, Shu | 2003 | 理論 | ★★☆ | `01_Theoretical_Foundations/2.5_Sebastian_Shu_2003_Multidomain.md` |
| 6 | 3.1 | Carrillo, Vecil | 2007 | SL+WENO | ⭐⭐⭐ | `02_SemiLagrangian_WENO_Interpolation/3.1_Carrillo_Vecil_2007_SL_WENO.md` |
| 7 | 3.2 | Qiu, Shu | 2011 | SL+WENO | ★★★ | `02_SemiLagrangian_WENO_Interpolation/3.2_Qiu_Shu_2011_Conservative_SL.md` |
| 8 | 3.3 | Yi, Kim | 2015 | SL+WENO | ★★☆ | `02_SemiLagrangian_WENO_Interpolation/3.3_Yi_Kim_2015_SL_WENO_Interp.md` |
| 9 | 3.4 | Christlieb et al. | 2023 | SL+WENO | ★☆☆ | `02_SemiLagrangian_WENO_Interpolation/3.4_Christlieb_2023_NonPoly_WENO.md` |
| 10 | 3.5 | Cottet | 2023 | SL+WENO | ★☆☆ | `02_SemiLagrangian_WENO_Interpolation/3.5_Cottet_2023_WENO_SL_Particle.md` |
| 11 | 3.6 | Sun, Xiao | 2017 | SL+WENO | ★★☆ | `02_SemiLagrangian_WENO_Interpolation/3.6_Sun_Xiao_2017_SL_MultiMoment.md` |
| 12 | 3.7 | Qiu, Russo | 2017 | SL+WENO | ★★☆ | `02_SemiLagrangian_WENO_Interpolation/3.7_Qiu_Russo_2017_MultiDim_Tracing.md` |
| 13 | 4.1 | Wilde et al. | 2021 | LBM | ⭐⭐⭐ | `03_LBM_Integration/4.1_Wilde_2021_SLLBM_Turbulence.md` |
| 14 | 4.2 | Wilde et al. | 2020 | LBM | ★★★ | `03_LBM_Integration/4.2_Wilde_2020_SLLBM_Foundation.md` |
| 15 | 4.3 | Saadat et al. | 2020 | LBM | ★★☆ | `03_LBM_Integration/4.3_Saadat_2020_Unstructured_SLLBM.md` |
| 16 | 4.4 | Kallikounis et al. | 2021 | LBM | ★★☆ | `03_LBM_Integration/4.4_Kallikounis_2021_Multiscale_SLLBM.md` |
| 17 | 4.5 | Reyhanian | 2021 | LBM | ★★☆ | `03_LBM_Integration/4.5_Reyhanian_2021_ETH_PonD.md` |
| 18 | 4.6 | Kallikounis et al. | 2022 | LBM | ★★☆ | `03_LBM_Integration/4.6_Kallikounis_2022_PonD_Discontinuities.md` |
| 19 | 4.7 | Ambruș et al. | 2025 | LBM | ★☆☆ | `03_LBM_Integration/4.7_Ambrus_2025_Vielbein_LBM.md` |
| 20 | 4.8 | Noh, Lee | 2025 | LBM | ★★☆ | `03_LBM_Integration/4.8_Noh_Lee_2025_Lagrangian_ELBM.md` |
| 21 | 5.1 | Janett et al. | 2019 | 任意點 | ⭐⭐⭐ | `04_WENO_Arbitrary_Point/5.1_Janett_2019_WENO_Interpolation.md` |
| 22 | 5.2 | Wang et al. | 2024 | 任意點 | ★★★ | `04_WENO_Arbitrary_Point/5.2_Wang_2024_Interp_vs_Recon.md` |
| 23 | 5.3 | Li et al. | 2020 | 任意點 | ★★☆ | `04_WENO_Arbitrary_Point/5.3_Li_2020_WENO_Interp_ShallowWater.md` |
| 24 | 5.4 | Ha et al. | 2024 | 任意點 | ★☆☆ | `04_WENO_Arbitrary_Point/5.4_Ha_2024_Alternative_WENO_ExpPoly.md` |

---

## 建議閱讀順序（與 README 一致）

1. Shu (2020) §2.5 → 理解 reconstruction vs interpolation
2. Jiang, Shu, Zhang (2013) → interpolation-based WENO FD 完整推導
3. Liu, Shu, Zhang (2009) → linear weights 正性問題
4. Carrillo & Vecil (2007) → SL + WENO interpolation departure point 處理
5. Janett et al. (2019) → 任意位置 WENO interpolation 專門技術
6. Qiu & Shu (2011) → 六階 WENO interpolation 完整公式
7. Wang et al. (2024) → interpolation vs reconstruction 比較
8. Wilde et al. (2021) → SLLBM 高階插值實際應用

---

## 與現有程式碼的對應關係

| 程式碼模組 | 對應文獻類別 | 主要參考 |
|-----------|-------------|---------|
| `interpolation_gilbm.h` (7-point Lagrange) | 待升級 → Cat. 04 | Janett (2019), Wang (2024) |
| `evolution_gilbm.h` (streaming step) | Cat. 03 | Wilde (2020/2021) |
| `metric_terms.h` (curvilinear Jacobian) | Cat. 03 | Saadat (2020), Ambruș (2025) |
| `boundary_conditions.h` (Chapman-Enskog BC) | Cat. 03 | Reyhanian (2021) |
| MRT collision (`gilbm_mrt_collision`) | Cat. 03 | Kallikounis (2021) |
| Local Time Stepping | Cat. 03 | Kallikounis (2022), Reyhanian (2021) |

*本索引自動產生，最後更新：2026-03-28*

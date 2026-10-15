# Hybrid Pixel–Object–Rule-Based Classification Framework with Stability Maps

**Study Area:** Timor-Leste (National Scale)  
**Sensors:** Sentinel-2 L2A (10m & 20m) & Landsat 8/9 CCDC Time-Series  
**Theoretical Baseline:** Malinverni & Sanità (*ISPRS Int. J. Geo-Inf.*, 2026, 15, 241)  
**Schema:** 18-Class Hierarchical LULC Schema (Merged to 17 in final cartography)

---

## 1. Executive Summary & Problem Statement

Traditional land use / land cover (LULC) classification pipelines face two contrasting challenges when applied to heterogeneous tropical island landscapes like Timor-Leste:

1. **Pixel-Based Classifiers (e.g., Random Forest):** Maximize fine spectral and temporal detail from continuous time-series (CCDC harmonics), but suffer from severe **"salt-and-pepper" noise**, terrain shadow artifacts, and high intra-class variance on spectrally identical surfaces (e.g., urban tin/gravel roofs vs. dry braided riverbeds).
2. **Pure Object-Based Classifiers (Standard OBIA):** Enforce spatial regularity through segmentation, but often **over-smooth heterogeneous smallholder mosaics**, erase narrow linear corridors (roads, riparian strips), and propagate classification errors across mixed parcels via rigid majority voting.

To overcome both limitations, this project implements a **Hybrid Pixel–Object–Rule-Based Classification Framework with Stability Maps**. The workflow leverages continuous, soft probabilistic predictions as the bridge between pixel-level spectral nuance and object-level spatial coherence, complemented by an operational **Stability Map** that spatially flags classification uncertainty for GIS-ready decision support.

---

## 2. Theoretical Framework (Malinverni & Sanità, 2026)

Our implementation is directly informed by the hybrid architecture developed by Malinverni & Sanità (2026), structured around five core pillars:

```
┌───────────────────────────┐       ┌───────────────────────────┐
│  Pixel Probability Stack  │       │   SNIC Spatial Segments   │
│ (18 Bands, S2 CCDC + RF)  │       │ (Optical Bands + Indices) │
└─────────────┬─────────────┘       └─────────────┬─────────────┘
              │                                   │
              └─────────────────┬─────────────────┘
                                ▼
              ┌───────────────────────────────────┐
              │     Object Probability Averaging  │
              │     P_obj = Mean(P_pix) in Object │
              └─────────────────┬─────────────────┘
                                ▼
              ┌───────────────────────────────────┐
              │    Intra-Object Stability Map     │
              │   CI = S / W  (Competition Ratio)│
              │   is_stable = (CI <= tau = 0.65)  │
              └─────────────────┬─────────────────┘
                                ▼
              ┌───────────────────────────────────┐
              │     Adaptive Convex Blending      │
              │ P_hyb = w_eff*P_obj + (1-w_eff)*P │
              │   w_eff = w_base * (1 - 0.6*CI)   │
              └─────────────────┬─────────────────┘
                                ▼
              ┌───────────────────────────────────┐
              │      Intrinsic Rule Refinement    │
              │  - False Hydrology Suppression    │
              │  - Post-Merge Riverbed (18 -> 17) │
              └─────────────────┬─────────────────┘
                                ▼
              ┌───────────────────────────────────┐
              │      Final Thematic Deliverables  │
              │   1. LULC Thematic Map (1-17)     │
              │   2. Spatially Explicit Stability │
              │   3. Continuous Confusion Index   │
              └───────────────────────────────────┘
```

---

## 3. Detailed Component Breakdown

### 3.1. Spatial Segmentation via SNIC (`4_1_SNIC_segmentations.js`)
* **Input Features:** Normalized multi-band composite containing `Red (B4)`, `NIR (B8)`, `SWIR1 (B11)`, and `NDVI`.
* **Bilinear Resampling:** Sentinel-2 B11 (natively 20m) is resampled via bilinear interpolation to 10m before computing indices, preventing blocky stair-step artifacts along segment boundaries.
* **Contrast Normalization:** Bands are clamped and normalized to $[0, 1]$ using empirical 2nd/98th percentiles so that Euclidean clustering distances are balanced across bands.
* **SNIC Parameterization:**
  * `size`: Seed spacing along each axis. (Currently 30, with future recommendations of 12–15 for smaller smallholder agricultural parcels).
  * `compactness`: Set to balance regular geometric parcel shapes with flexible meandering boundaries for natural waterways.
  * `connectivity`: 8-neighbor connectivity.
* **Tiling & Global ID Stride:** Timor-Leste is tiled into 5 regional bounding boxes (`central_timor`, `east_box`, `off_dili`, `west_box`, `east_timor`). Segment IDs in each tile are offset by `ID_STRIDE = 300,000,000` to prevent ID collisions when mosaicking across the country.

---

### 3.2. Intra-Object Probability Aggregation & Stability Index
Rather than collapsing pixels into hard discrete classes before OBIA, the workflow aggregates the full **18-class continuous probability vector** within each segment:

$$P_{\text{obj}}(c) = \frac{1}{|R_i|} \sum_{p \in R_i} P_{\text{pix}}(p, c)$$

#### The Confusion Index ($CI$) & Stability Map
To represent classification reliability spatially without relying solely on global confusion matrices, the script extracts the top two class contenders:
* **$W$ (Winner):** Maximum probability in the vector.
* **$S$ (Runner-up):** Second-highest probability in the vector.
* **Confusion Index ($CI$):**
  $$CI = \frac{S}{W} \quad (0.0 \le CI \le 1.0)$$
* **Stability Map (`is_stable`):**
  $$\text{Object is Stable} \iff CI \le \tau \quad (\tau = 0.65)$$
  * **$CI \le 0.65$ (Stable):** Clear class dominance. High confidence, homogeneous parcel.
  * **$CI > 0.65$ (Unstable):** Close competition between classes. Heterogeneous parcel or transitional boundary requiring cautious interpretation.

---

### 3.3. Adaptive Convex Blending
Instead of a rigid Winner-Takes-All (WTA) reclassification, the hybrid model uses convex blending modulated by **class-specific prior smoothness** and **local object stability**:

$$P_{\text{hyb}} = w_{\text{eff}} \cdot P_{\text{obj}} + (1 - w_{\text{eff}}) \cdot P_{\text{pix}}$$

1. **Base Class Weights ($w_{\text{base}}$):** Derived from ecological priors (`smoothnessTL`). Homogeneous, contiguous classes receive high object trust (Primary Forest: $0.90$, Water: $0.80$), while fragmented, fine-grained classes maintain higher pixel sensitivity (Built-up: $0.20$, Dry riverbed: $0.10$).
2. **Dynamic Stability Damping:**
   $$w_{\text{eff}} = w_{\text{base}} \cdot \left(1.0 - 0.60 \times CI_{\text{obj}}\right)$$
   * In **highly stable parcels ($CI \to 0$)**, full object smoothing is enforced to eliminate noise.
   * In **unstable/ambiguous zones ($CI \ge 0.65$)**, object weight automatically drops, allowing the model to fall back on pixel-level spectral measurements to prevent smearing across sharp boundaries.

---

### 3.4. Intrinsic Rule-Based Refinement (`refineRules`)
Applied to resolve systematic spectral confusion identified in Timor-Leste without requiring external DEM or ancillary GIS layers:

1. **False Hydrology Suppression (Urban Roof / Asphalt Shadow Confusion):**
   * *The Problem:* Urban zinc roofs, dark asphalt, and construction clearings have low, flat reflectance that often triggers high probabilities for *Dry riverbed (18)* or *Waterbody (17)*.
   * *The Rule:*
     $$\text{IF } (W \in \{17, 18\}) \text{ AND } (S \in \{15, 16\}) \text{ AND } (CI > 0.65) \implies \text{Reassign to } S$$
     Terrestrial runner-up classes (*Built-up: 16* or *Cleared land: 15*) take precedence over false water/riverbed labels when classification is ambiguous.
2. **Post-Classification Hydrological Merging:**
   * After all stability metrics and refinement rules have executed, *Dry riverbed (18)* is merged into *Waterbody (17)*:
     $$LULC[LULC == 18] = 17$$
   * Keeping them separate until this final step ensures their vastly different seasonal CCDC dynamics do not pollute the pre-classification probability distributions.

---

## 4. Multi-Temporal Architecture: Planned 2020 Mapping

Mapping 2020 requires multi-temporal consistency to prevent spurious change reporting ("classifier churn"). The 2020 workflow is coordinated with 2025 as follows:

```
       2020 S2 CCDC Probabilities               2025 S2 CCDC Probabilities
                  │                                         │
                  └────────────────────┬────────────────────┘
                                       │
                                       ▼
                       Unified Spatial Superpixels (10m)
                     (Prevents Edge Chatter / Slivers)
                                       │
                     ┌─────────────────┴─────────────────┐
                     ▼                                   ▼
          Hybrid Pipeline (2020)              Hybrid Pipeline (2025)
          - Adaptive Blending                 - Adaptive Blending
          - Stability Map 2020                - Stability Map 2025
          - Rule Refinement 2020              - Rule Refinement 2025
                     │                                   │
                     └─────────────────┬─────────────────┘
                                       ▼
                      Vectorized Temporal Consistency Engine
                                (Local Python)
                                       │
               ┌───────────────────────┴───────────────────────┐
               ▼                                               ▼
     Impossible Transitions                          Unlikely Transitions
    (e.g., Built-up -> Forest)                    (Modulated by CCDC Stable Mask)
               │                                               │
               └───────────────────────┬───────────────────────┘
                                       ▼
                   Consistent Multi-Year GIS Deliverables
                    1. Consistent LULC 2020 & 2025 (10m)
                    2. LULC Transition Matrix (Hectares)
                    3. Uncertainty-Weighted Change Confidence
```

### 4.1. Unified Spatial Units (Eliminating Sliver Errors)
* Running an independent SNIC segmentation for 2020 is explicitly avoided.
* Comparing two independent segmentations introduces thousands of artificial "sliver polygons" along parcel boundaries due to slight variations in cloud masking, sun angle, or seed distribution.
* **Strategy:** The 2020 probability stack is processed through the **same spatial segment boundaries**, ensuring that land cover changes represent true attribute transitions within identical parcels.

### 4.2. Local Python Temporal Consistency Engine
To eliminate Google Earth Engine EECU limitations and memory timeouts, temporal consistency rules are evaluated locally in Python using vectorized NumPy arrays:

1. **Ecologically Impossible Transitions (Overridden Everywhere):**
   * Class transitions that cannot occur within a 5-year window (e.g., *Secondary forest $\to$ Primary forest*, *Built-up $\to$ Primary forest*, *Waterbody $\to$ Primary forest*) are automatically corrected:
     $$\text{Transition Pair} = LULC_{2020} \times 100 + LULC_{2025}$$
     $$\text{IF Pair} \in \{201, 701, 1201, 1301, 1401, 1501, 1601, 1701\} \implies LULC_{2020} = LULC_{2025}$$
2. **Unlikely Transitions (Constrained by CCDC Stable Mask):**
   * Questionable transitions (e.g., *Secondary forest $\to$ Shrubland*) are evaluated against the CCDC stable mask (`TL_stable_mask_2020_2025`). If CCDC detected zero breaks and Hansen GFC confirms no loss, the pixel is enforced as stable ($LULC_{2020} = LULC_{2025}$).
3. **Uncertainty-Weighted Change Confidence:**
   * A pixel marked as "changed" is only reliable if both dates were classified with high stability:
     $$\text{Change Confidence} = (1.0 - CI_{2020}) \times (1.0 - CI_{2025})$$
   * Prevents classifier confusion along transitional boundaries from being falsely counted as deforestation or land degradation in national statistics.

---

## 5. Repository File Map & Execution Workflow

| File | Environment | Role | Key Inputs & Outputs |
| :--- | :--- | :--- | :--- |
| [`Sentinel-2/4_1_SNIC_segmentations.js`](file:///c:/Users/AFahrezi/Documents/GitHub/Temporal_segment_gee/Sentinel-2/4_1_SNIC_segmentations.js) | GEE (JS) | Generates and exports the 5 regional 10m SNIC segment tiles. | **In:** S2 2025 Median Composite<br>**Out:** `projects/gee-v1-510606/assets/snic/seg_*` |
| [`Sentinel-2/6_hybrid_pixel_object_assembly.js`](file:///c:/Users/AFahrezi/Documents/GitHub/Temporal_segment_gee/Sentinel-2/6_hybrid_pixel_object_assembly.js) | GEE (JS) | Interactive testing script for small AOIs with visualization layers. | **In:** 2025 Probabilities + SNIC assets<br>**Out:** Interactive Map layers & console diagnostic metrics |
| [`Sentinel-2/7_export_probabilities_and_segments.js`](file:///c:/Users/AFahrezi/Documents/GitHub/Temporal_segment_gee/Sentinel-2/7_export_probabilities_and_segments.js) | GEE (JS) | Optimized exporter to Google Drive for 2020 & 2025 probabilities and segments. | **In:** Probability assets + SNIC assets<br>**Out:** Multi-tile Cloud-Optimized GeoTIFFs in Drive |
| [`python/hybrid_assembly.py`](file:///c:/Users/AFahrezi/Documents/GitHub/Temporal_segment_gee/python/hybrid_assembly.py) | Local Python | High-performance batch processor for blending, stability, and rules. | **In:** Downloaded GeoTIFF tiles<br>**Out:** Final 10m LULC, Stability Map, and Confusion Index rasters |

---

## 6. Target 18-Class Schema Reference

| ID | Class Name | Base Smoothness ($w$) | Final Post-Merge ID | Ecological Group |
| :-: | :--- | :-: | :-: | :--- |
| **1** | Primary_forest | 4.5 ($w=0.90$) | 1 | Natural Forest |
| **2** | Secondary_forest | 3.0 ($w=0.60$) | 2 | Natural Forest / Regeneration |
| **3** | Mangrove_forest | 3.0 ($w=0.60$) | 3 | Coastal / Tidal Wetland |
| **4** | Coastal_forest | 3.5 ($w=0.70$) | 4 | Coastal Littoral Forest |
| **5** | Teak_plantation | 1.2 ($w=0.24$) | 5 | Tree Plantation (Deciduous) |
| **6** | Eucalyptus_plantation | 1.2 ($w=0.24$) | 6 | Tree Plantation (Evergreen) |
| **7** | Coffe_agroforestry | 1.5 ($w=0.30$) | 7 | Montane Agroforestry |
| **8** | Cocoa_agroforestry | 1.5 ($w=0.30$) | 8 | Lowland Agroforestry |
| **9** | Coconut_agroforestry | 1.5 ($w=0.30$) | 9 | Coastal Agroforestry |
| **10** | Mixed_garden | 1.2 ($w=0.24$) | 10 | Rural Settlement Mosaic |
| **11** | Paddy_field | 1.0 ($w=0.20$) | 11 | Irrigated / Lowland Agriculture |
| **12** | Other_cropland | 1.0 ($w=0.20$) | 12 | Rainfed / Upland Cropland |
| **13** | Grassland | 1.5 ($w=0.30$) | 13 | Savanna / Herbaceous |
| **14** | Shrubland | 1.5 ($w=0.30$) | 14 | Savanna / Woody Shrub |
| **15** | Cleared_land | 1.5 ($w=0.30$) | 15 | Bare Ground / Clearing |
| **16** | Built-up | 1.0 ($w=0.20$) | 16 | Impervious / Settlements |
| **17** | Waterbody | 4.0 ($w=0.80$) | 17 | Hydrological (Water & Riverbed) |
| **18** | Dry_riverbed | 0.1 ($w=0.20$) | **Merged into 17** | Hydrological (Seasonal Braided Riverbed) |

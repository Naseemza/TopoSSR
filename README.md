# TopoSSR: Topological Data Analysis for Protein Mutations

---

## 📐 Mathematical Foundation

### Persistent Homology
TopoSSR uses **persistent homology** to extract topological features from protein structures. For each protein block, we compute a **Vietoris-Rips complex** in 9-dimensional torsion angle space:

**Filtration:** $C_r = \{(p_i, p_j) : d(p_i, p_j) \leq r\}$

where $p_i$ = torsion angle vector for residue $i$, and $d$ = Euclidean distance.

### Maximum Persistent Betti-1 (MaxPB-1)
The **Betti-1 score** $\hat{\beta}_1$ counts 1-dimensional topological features (loops/cycles):

$$\hat{\beta}_1 = \max_{r \in [0, r_{\max}]} \beta_1(r)$$

**Block Deformation Detection:**
$$\Delta\beta_1 = |\hat{\beta}_1(\text{PDB1}) - \hat{\beta}_1(\text{PDB2})|$$

- **Match (Conserved):** $\Delta\beta_1 \leq 2$
- **No Match (Deformed):** $\Delta\beta_1 > 2$

### Torsion Angle Features
9-dimensional feature vector for each residue:

$$\mathbf{x}_i = [\phi_i, \psi_i, \chi_1^{(i)}, \chi_2^{(i)}, \chi_3^{(i)}, \chi_4^{(i)}, \chi_5^{(i)}, \omega_i, \alpha_i]$$

where each angle is **sigmoid-normalized:**

$$\theta' = \frac{1}{1 + e^{-\theta/60}}$$

---

## 🔬 Methodology

### Step 1: Extract Torsion Angles
- Download PDB structure from RCSB
- Use R/bio3d to compute dihedral angles
- Apply sigmoid normalization to map angles to [0,1]

### Step 2: Segment by Secondary Structure
- Identify α-helices (AH), β-sheets (BS), loops (TL) from PDB header
- Group residues into contiguous blocks of same SSR type

### Step 3: Gaussian Process Upsampling
For each block with SSR type $s$, upsample by factor $k$ using GP regression:

$$\text{GP kernel} = \begin{cases} 
\text{ExpSineSquared}(3.6, 1.0) & \text{if } s = \text{AH} \\
\text{Matern}(1.0, \nu=2.5) & \text{if } s = \text{BS} \\
\text{RBF}(0.5) & \text{if } s = \text{TL}
\end{cases}$$

### Step 4: Compute Persistent Homology
- **Method 1:** Vietoris-Rips complex on 9-D Euclidean distances
- **Method 2:** Alpha complex on 2-D PCA projection

Calculate persistent intervals $[b_i, d_i)$ for 1-cycles.

### Step 5: Compute MaxPB-1
From persistence diagram, extract:

$$\hat{\beta}_1 = \max_r \{i : b_i \leq r < d_i\}$$

### Step 6: Compare Blocks
Compute $\Delta\beta_1$ between wild-type and mutant for each block.

---

## 📋 Prerequisites

### 1. Install Python 3.8+
```bash
python3 --version  # Should show 3.8 or higher
```

### 2. Install R
**Linux:**
```bash
sudo apt-get install r-base r-base-dev
```

**macOS:**
```bash
brew install r
```

### 3. Install bio3d R package
```bash
Rscript -e 'install.packages("bio3d")'
```

### 4. Install Python packages
```bash
pip install pandas numpy scikit-learn matplotlib ripser gudhi rpy2
```

---

## ▶️ How to Run

### Option 1: Simple Copy-Paste (Recommended)

Save this as `run_analysis.py`:

```python
from tda_analysis_v3 import run_full_analysis

# Analyze T4 Lysozyme: Wild-type (2LZM) vs L99A Mutant (1L63)
torsion_df, coord_df, raw_df1, raw_df2 = run_full_analysis(
    pdb1="2LZM",
    pdb2="1L63",
    out_dir="results"
)

print("\n=== TORSION ANGLE ANALYSIS ===")
print(torsion_df)

print("\n=== COORDINATE ANALYSIS ===")
print(coord_df)

print("\n✅ Analysis complete! Check 'results/' folder for plots.")
```

**Run it:**
```bash
python3 run_analysis.py
```

---

### Option 2: Interactive Python

```python
from tda_analysis_v3 import run_full_analysis

# Analyze SARS-CoV-2 Spike: Wuhan (6VSB) vs Omicron (7TGW)
torsion_df, coord_df, raw_df1, raw_df2 = run_full_analysis(
    pdb1="6VSB",
    pdb2="7TGW",
    c1="A",
    c2="A",
    multiplier=10,
    out_dir="spike_analysis"
)
```

---

## 📊 Example: T4 Lysozyme

**Compare:** Wild-type (2LZM) vs L99A cavity mutant (1L63)

```python
from tda_analysis_v3 import run_full_analysis

torsion_df, coord_df, raw_df1, raw_df2 = run_full_analysis(
    pdb1="2LZM",
    pdb2="1L63",
    out_dir="lysozyme_results"
)

# Expected result:
# Block 18 (α-helix, residues 93-106) shows Δβ₁ = 5 → DEFORMED
# This matches the engineered L99A cavity location
```

---

## 📊 Example: SARS-CoV-2 Spike Protein

**Compare:** Wuhan reference (6VSB) vs Omicron variant (7TGW)

```python
from tda_analysis_v3 import run_full_analysis

torsion_df, coord_df, raw_df1, raw_df2 = run_full_analysis(
    pdb1="6VSB",
    pdb2="7TGW",
    c1="A",
    c2="A",
    multiplier=10,
    out_dir="spike_results"
)

# Expected result:
# ~30 blocks show Δβ₁ > 2 (22% deformation rate)
# Major changes in NTD and RBD domains
```

---

## 📤 Output Files

After running, check the output directory for:

1. **`residue_vs_betti_torsion.png`** — Comparative β₁ profiles
2. **`residue_vs_betti_coords.png`** — Coordinate-based analysis
3. **`tda_torsion_block_*.png`** — Detailed block analysis (persistence diagrams + β₁ curves)
4. **Summary tables** — Console output with deformation scores

---

## 🔑 Key Parameters

| Parameter | Default | What It Does |
|-----------|---------|--------------|
| `pdb1` | Required | First PDB ID (e.g., "2LZM") |
| `pdb2` | Required | Second PDB ID (e.g., "1L63") |
| `c1`, `c2` | `None` | Chain ID (optional; uses first chain if not specified) |
| `multiplier` | 10 | GP upsampling factor (increase for smoother curves) |
| `out_dir` | "." | Output directory for plots |

---

## 📋 Interpretation

### Δβ₁ Scores

| Score | Meaning | Example |
|-------|---------|---------|
| 0–2 | ✅ Topologically conserved | Minor structural adjustments |
| 3–5 | ⚠️ Moderate deformation | Local loop rearrangement |
| 6–10 | ❌ Significant disruption | Helix collapse/reorganization |
| >10 | ❌ Severe deformation | Complete topological change |

---



**Contact:** achakraborty@cus.ac.in

---

## 📚 References

1. Edelsbrunner, H., & Harer, J. (2010). *Computational Topology*. AMS.
2. Carlsson, G. (2009). Topology and data. *Bulletin of AMS*, 46(2), 255–308.
3. Zomorodian, A. (2005). Fast algorithm for computing persistence homology. *Algorithms*, 2005.
4. Bauer, U., Kerber, M., Reininghaus, J., & Wagner, H. (2014). Phat. *Software*, github.com/compute-homology/phat.
5. Grant, B. J., et al. (2006). bio3d: An R package for the comparative analysis of protein structures. *Bioinformatics*, 22(21), 2695–2696.

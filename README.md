# CFTR Variant Analysis with ClinVar

A beginner-friendly bioinformatics project that explores variants in the **CFTR** gene, the gene mutated in cystic fibrosis, using the public [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) database.

## Background

Cystic fibrosis (CF) is an inherited disorder caused by mutations in the *CFTR* gene, which encodes a chloride channel in epithelial cells. Faulty channels lead to thick, sticky mucus, mainly affecting the lungs and digestive system. Chronic lung damage can also strain the heart (pulmonary hypertension and right-heart strain), which links this gene to cardiopulmonary research.

Hundreds of CFTR variants have been reported. This project asks:

1. How many CFTR variants are classified as pathogenic, benign or of uncertain significance?
2. What molecular types of variants are most common?
3. Where along the gene do pathogenic and benign variants cluster?

## Data

- **Source:** ClinVar `variant_summary.txt.gz` (NCBI, public, updated regularly)
- **Filter:** records for gene *CFTR*, assembly GRCh38
- The raw file is large and **not stored in this repo**; the script downloads it automatically.

## Methods

- Python (pandas) filters the ClinVar table to CFTR.
- ClinVar's many significance labels are grouped into five categories.
- Matplotlib produces the figures.

## How to run

```bash
git clone https://github.com/<your-username>/cftr-variant-analysis.git
cd cftr-variant-analysis
conda env create -f environment.yml
conda activate cftr-variants
python scripts/analyze_cftr.py
```

Outputs:

- `data/cftr_variants.csv`: filtered CFTR table
- `results/summary_by_category.csv`: counts per category
- `results/fig1_clinical_significance.png`, `fig2_variant_types.png`, `fig3_position_along_gene.png`

## Results

_Add your figures here after running the script, then write 2-3 sentences describing what you see._

```markdown
![Clinical significance](results/fig1_clinical_significance.png)
```

## Limitations

- ClinVar contains submitter-reported classifications, which can change or conflict.
- Counts reflect what has been submitted, not how common a variant is in the population.
- This is an educational analysis, not clinical advice.

## Possible extensions

- Highlight the well-known F508del variant.
- Map variants onto CFTR protein domains using UniProt.
- Compare CFTR with another gene.
- Turn the script into a Jupyter notebook.

## Data and licence

ClinVar data are public domain from NCBI. Code released under the MIT licence (add a LICENSE file via GitHub).

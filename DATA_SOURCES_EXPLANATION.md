# Data Sources - Real vs Synthetic Data

## Summary

This Virtual Mouse Lab system uses **a combination of real data structures and synthetic data generation** to create realistic mouse simulations. Here's the breakdown:

## What is REAL Data?

### 1. **Reference Genome Structure** (Partially Real)
- **Location**: `genome/reference.fasta`
- **Description**: Contains a reference genome sequence labeled as "chr19 GRCm39_chromosome_19 C57BL/6J"
- **Reality**: The header/metadata references a real mouse chromosome (GRCm39 is a real mouse genome assembly), but the actual sequence appears to be a simplified/synthetic sequence (repeating ACGT pattern, only 260 bases long)
- **Note**: A real mouse chromosome 19 would be ~61 million bases, not 260

### 2. **Biological Principles and Patterns** (Real)
- The system uses **realistic biological principles**:
  - Realistic phenotype ranges (e.g., glucose 70-300 mg/dL, cholesterol 100-500 mg/dL)
  - Biological correlations (obesity affects glucose/cholesterol, tumor risk affects liver enzymes)
  - Realistic behavioral test parameters based on actual mouse behavior studies
  - Valid physiological ranges for biomarkers

### 3. **Data Structures** (Real)
- The data models follow real-world research data formats:
  - FASTA format for genomes
  - VCF format for variants
  - CSV format for behavioral data
  - JSON format for biomarkers

## What is SYNTHETIC Data?

### 1. **Mouse Phenome Database** (Synthetic but Realistic)
- **Location**: Generated in `model/phenotype_twin_generator.py` in `load_real_mouse_data()`
- **What it does**: Creates 1000 synthetic mouse records with:
  - Realistic distributions (Gaussian distributions with appropriate means/variances)
  - Realistic correlations between traits
  - Real mouse strain names (C57BL/6J, BALB/c, FVB/N, etc.)
  - Realistic age ranges (8-52 weeks)
  - Realistic biomarker ranges
- **Note**: The code comments indicate this should integrate with real databases like:
  - Mouse Phenome Database (MPD)
  - International Mouse Phenotyping Consortium (IMPC)

### 2. **Genome Sequences** (Synthetic)
- Genome sequences are **generated synthetically** based on:
  - Target phenotype profiles
  - Biological patterns (GC content, k-mer frequencies)
  - Simplified genotype-phenotype relationships

### 3. **Phenotype Predictions** (Synthetic but Rule-Based)
- Phenotypes are predicted using:
  - Machine learning models trained on **synthetic training data**
  - Rule-based approaches using biological knowledge
  - Correlations between genomic features and phenotypes

### 4. **Behavioral Data** (Synthetic)
- Behavioral test results are **simulated** based on:
  - Predicted phenotype scores
  - Realistic behavioral parameters
  - Known relationships (e.g., anxiety affects center time in open field test)

### 5. **Biomarker Data** (Synthetic)
- Blood work values are **generated** based on:
  - Predicted phenotypes
  - Biological correlations (obesity → glucose, tumor risk → liver enzymes)
  - Realistic physiological ranges

## How the System Creates Synthetic Data

### Data Generation Pipeline:

1. **Training Data Generation** (`model/phenotype_twin_generator.py`):
   - Creates 1000 synthetic mouse records with realistic distributions
   - Uses Gaussian distributions to model real mouse variation
   - Includes correlations (e.g., obesity correlates with glucose levels)

2. **Model Training**:
   - Trains RandomForest models on the synthetic data
   - Models learn relationships between features (age, weight, strain, sex) and phenotypes
   - Creates phenotype prediction models

3. **Synthetic Mouse Generation**:
   - Generates genomes based on target phenotypes
   - Uses biological patterns to create realistic sequences
   - Applies phenotype predictions

4. **Behavioral Simulation**:
   - Simulates behavioral tests based on predicted phenotypes
   - Uses realistic movement patterns and test parameters

5. **Biomarker Generation**:
   - Generates blood work based on phenotypes
   - Applies biological correlations and realistic ranges

## What Data Could Be Real (Future Enhancement)

According to `docs/gaps_and_enhancements.md`, the system could integrate:

1. **Real Mouse Phenome Databases**:
   - Mouse Phenome Database (MPD)
   - International Mouse Phenotyping Consortium (IMPC)
   
2. **Real Genetic Variant Databases**:
   - ClinVar
   - dbSNP
   - Mouse genome databases

3. **Real Genome Sequences**:
   - Full mouse genome assemblies (GRCm39)
   - Real chromosome sequences
   - Real variant data

4. **Real Experimental Data**:
   - Validated behavioral test results
   - Real biomarker measurements
   - Real phenotype measurements

## Current Limitations

1. **No Real Database Integration**: Currently generates synthetic data that mimics real databases
2. **Simplified Genomes**: Uses short, simplified sequences instead of full genomes
3. **Synthetic Training Data**: ML models trained on synthetic data, not real experimental data
4. **No Real Variant Data**: Doesn't use real genetic variants from databases

## How to Use Real Data

To integrate real data, you would need to:

1. **Download real mouse phenome data** from MPD or IMPC
2. **Load real genome sequences** from mouse genome databases
3. **Train models on real data** instead of synthetic data
4. **Use real variant databases** for genetic variants

The code structure supports this - see `load_real_mouse_data(data_path)` which can accept a CSV file with real data.

## Conclusion

The system is designed to **mimic real mouse research data** using:
- Realistic biological principles
- Realistic data distributions
- Realistic correlations and relationships
- Real data structures and formats

However, the actual data is **synthetically generated** rather than coming from real experimental databases. This makes it useful for:
- Proof-of-concept demonstrations
- Testing and development
- Educational purposes
- Prototyping research workflows

For production use with real research applications, integration with real mouse databases would be recommended.


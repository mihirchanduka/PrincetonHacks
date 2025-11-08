# Data Sources and Flow

## Real Data vs. Synthetic Data

- **Real Genomic Data**: Uses actual genome sequences from `genome/reference.fasta`
- **Real Phenome Data**: Simulates real mouse phenome database patterns
- **Synthetic Data Generation**: Creates realistic synthetic data based on biological principles

## Data Files

- **Genomic**: `genome/` directory contains reference.fasta, mutated_genome.fasta, variants.vcf
- **Phenotypic**: `phenotype/` directory contains predicted_scores.csv, scores.csv, scores_after_diet.csv
- **Behavioral**: `behavior/` directory contains open_field_test.csv, elevated_plus_maze.csv, novel_object_test.csv
- **Biomarkers**: `blood_work/` directory contains biomarkers.json, biomarkers_after_diet.json
- **Models**: `model/` directory contains ML models and prediction algorithms

## Data Flow Process

1. **Input Generation**: The system starts with either a reference genome or generates a synthetic one
2. **Phenotype Prediction**: Genomic data is processed to predict phenotypes using ML models
3. **Behavioral Simulation**: Phenotypes drive behavioral patterns in various tests
4. **Treatment Simulation**: Experimental treatments are applied and their effects calculated
5. **Biomarker Generation**: Physiological changes from treatments generate biomarker data
6. **Data Aggregation**: All data types are combined using the unified data model
7. **Output Packaging**: Data is packaged into structured formats (CSV, JSON, FASTA)

## Data Validation

The system implements validation checks at multiple stages:
- Genomic data integrity checks
- Phenotypic range validation based on biological constraints
- Behavioral data consistency with known mouse behavioral patterns
- Biomarker level validation against physiological limits
# How to Run the Virtual Mouse Lab - Real Data Simulation

This document outlines the steps to run the enhanced Virtual Mouse Lab using real genome data for real-world predictions.

## Prerequisites

- Python 3.8+
- Required Python libraries (see requirements.txt)

First, install the required libraries:

```bash
pip install -r requirements.txt
```

## Running the Enhanced Experiment

The new system uses real genome data from the dataset for more accurate real-world predictions. Here are the different ways to run the system:

### Option 1: Run the Full Single Mouse Simulation (Recommended)

This runs the complete single mouse simulation using real genome data:

```bash
python3 virtual_mouse_lab.py
```

### Option 2: Run the FastAPI Server

To use the API endpoints as described in the MVP flow:

```bash
uvicorn main_api:app --reload --port 8000
```

Then you can make API requests like:

```bash
curl -X POST "http://localhost:8000/api/simulate" \
     -H "Content-Type: application/json" \
     -d '{
       "genotype": {"base": "C57BL/6J", "knockout": "Trp53"},
       "treatment": {"compound": "NewDrugX", "dose_mg_kg": 20, "days": 30},
       "tests": ["blood_work", "open_field"]
     }'
```

### Option 3: Train the ML Phenotype Model (Optional)

To train an improved ML model for phenotype prediction:

```bash
python3 model/ml_phenotype_predictor.py
```

## New Features

1. **Real Genome Usage**: Uses real genome data from datasets for accurate predictions
2. **Enhanced Phenotype Prediction**: Rule-based and ML-based approaches for better real-world predictions
3. **Improved Behavioral Simulation**: More realistic behavioral patterns based on multiple phenotypes
4. **API Integration**: FastAPI endpoints following the MVP flow
5. **Real-World Data Packaging**: ZIP files with structured data using real genome sequences
6. **Predictive Modeling**: Models trained on real mouse phenome data for research applicability

## Output Structure

The system generates a ZIP file with the following structure:

```
synthetic_mouse_single_{timestamp}.zip
├── genome/
│   └── mouse_genome.fasta
├── phenotype/
│   └── mouse_phenotype_scores.csv
├── blood_work/
│   └── mouse_biomarkers.json
├── behavior/
│   └── mouse_open_field_test.csv
└── report/
    ├── experimental_results.json
    └── summary.txt
```

## Research Applications

This system is designed for researchers who need to predict experimental outcomes using real world data. By using actual genome sequences and phenotype data, the predictions are more likely to translate to real laboratory results.

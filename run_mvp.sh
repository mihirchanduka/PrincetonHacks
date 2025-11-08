#!/bin/bash

# This script runs the full MVP workflow for the synthetic mouse experiment.

# Ensure the script is run from the 'synthetic_mouse_001' directory
if [ ! -f "generate_genome.py" ]; then
    echo "Error: This script must be run from the 'synthetic_mouse_001' directory."
    exit 1
fi

echo "--- Step 1: Generating Mutated Genome ---"
python3 generate_genome.py

echo ""
echo "--- Step 2: Predicting Phenotype ---"
python3 predict_phenotype.py

echo ""
echo "--- Step 3: Running Behavioral Test ---"
python3 run_behavioral_test.py

echo ""
echo "--- Step 4: Simulating High-Fat Diet ---"
# Copy predicted scores to the location expected by the diet simulation script
cp phenotype/predicted_scores.csv phenotype/scores.csv
python3 simulate_experiment.py

echo ""
echo "--- MVP Workflow Complete ---"


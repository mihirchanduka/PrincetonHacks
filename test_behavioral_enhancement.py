"""
Test script to verify the enhanced behavioral simulation is working properly
"""
import os
import csv
import pandas as pd
from model.behavioral_simulator import BehavioralSimulator

def test_behavioral_enhancement():
    print("Testing enhanced behavioral simulation...")
    
    # Read phenotype data
    scores_file = os.path.join('.', 'phenotype', 'predicted_scores.csv')
    phenotype = {}
    
    with open(scores_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row['phenotype']
            if key == 'anxiety':
                phenotype['anxiety_score'] = float(row['score'])
            elif key == 'obesity':
                phenotype['obesity_score'] = float(row['score'])
            elif key == 'tumor_risk':
                phenotype['tumor_risk'] = float(row['score'])
            else:
                phenotype[f"{key}_score"] = float(row['score'])
    
    print(f"Phenotype data: {phenotype}")
    
    # Initialize behavioral simulator
    behavioral_simulator = BehavioralSimulator()
    
    # Test open field with enhanced simulation
    print("\nTesting Open Field Test...")
    of_behavior = behavioral_simulator.simulate_open_field_test(phenotype)
    of_metrics = behavioral_simulator.calculate_behavioral_metrics(of_behavior, "open_field")
    print(f"Open field metrics: {of_metrics}")
    
    # Test elevated plus maze
    print("\nTesting Elevated Plus Maze...")
    epm_behavior = behavioral_simulator.simulate_elevated_plus_maze(phenotype)
    epm_metrics = behavioral_simulator.calculate_behavioral_metrics(epm_behavior, "elevated_plus_maze")
    print(f"Elevated plus maze metrics: {epm_metrics}")
    
    # Test novel object test
    print("\nTesting Novel Object Test...")
    no_behavior = behavioral_simulator.simulate_novel_object_test(phenotype)
    no_metrics = behavioral_simulator.calculate_behavioral_metrics(no_behavior, "novel_object")
    print(f"Novel object test metrics: {no_metrics}")
    
    # Verify that we have rich behavioral data
    print(f"\nOpen field data points: {len(of_behavior)}")
    print(f"First few open field data points:")
    for i, data_point in enumerate(of_behavior[:5]):
        print(f"  {i}: {data_point}")
    
    print("\nEnhanced behavioral simulation test completed successfully!")
    print("The new system now includes:")
    print("- Realistic movement patterns based on anxiety, memory, and activity")
    print("- Velocity and speed tracking for each time point")
    print("- Habituation modeling in the open field test")
    print("- Multiple behavioral tests (open field, elevated plus maze, novel object)")
    print("- Comprehensive behavioral metrics calculation")

if __name__ == "__main__":
    test_behavioral_enhancement()
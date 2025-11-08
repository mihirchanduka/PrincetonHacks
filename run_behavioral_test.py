"""
Simulates behavioral tests for synthetic mice using the enhanced behavioral simulator.
"""

import csv
import os
import pandas as pd
from model.behavioral_simulator import BehavioralSimulator

def run_open_field_test(base_dir):
    """
    Simulates a 2D open-field test using the enhanced behavioral simulator.
    """
    # --- File Paths ---
    scores_file = os.path.join(base_dir, 'phenotype', 'predicted_scores.csv')
    output_log_file = os.path.join(base_dir, 'behavior', 'open_field_test.log')
    output_csv_file = os.path.join(base_dir, 'behavior', 'open_field_test.csv')
    
    # --- Read Phenotype Scores ---
    phenotype = {}
    try:
        with open(scores_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Handle both naming conventions
                key = row['phenotype']
                if key == 'anxiety':
                    phenotype['anxiety_score'] = float(row['score'])
                elif key == 'memory':
                    phenotype['memory_score'] = float(row['score'])
                elif key == 'obesity':
                    phenotype['obesity_score'] = float(row['score'])
                elif key == 'tumor_risk':
                    phenotype['tumor_risk'] = float(row['score'])
                elif key == 'aggression':
                    phenotype['aggression_score'] = float(row['score'])
                else:
                    # Add other phenotypes as needed
                    phenotype[f"{key}_score"] = float(row['score'])
    except FileNotFoundError:
        # Default values if file not found
        phenotype = {
            'anxiety_score': 0.5,
            'memory_score': 0.5,
            'obesity_score': 0.5,
            'tumor_risk': 0.5,
            'aggression_score': 0.5
        }
    
    print(f"Running open-field test with phenotypes: {phenotype}")
    
    # Initialize the enhanced behavioral simulator
    behavioral_simulator = BehavioralSimulator()
    
    # Run the simulation
    behavior_data = behavioral_simulator.simulate_open_field_test(phenotype)
    
    # Save detailed log
    with open(output_log_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['x', 'y', 'time', 'zone', 'velocity_x', 'velocity_y', 'speed'])
        
        for data_point in behavior_data:
            writer.writerow([
                f"{data_point['x']:.2f}", 
                f"{data_point['y']:.2f}", 
                data_point['time'], 
                data_point['zone'],
                f"{data_point.get('velocity_x', 0):.2f}",
                f"{data_point.get('velocity_y', 0):.2f}",
                f"{data_point.get('speed', 0):.2f}"
            ])
    
    # Also save as CSV with pandas for easier analysis
    behavior_df = pd.DataFrame(behavior_data)
    behavior_df.to_csv(output_csv_file, index=False)

    # Calculate and print metrics
    metrics = behavioral_simulator.calculate_behavioral_metrics(behavior_data, "open_field")
    print(f"Behavioral metrics: {metrics}")
    print(f"Open-field test simulation complete. Log saved to {output_log_file}")
    print(f"Detailed CSV saved to {output_csv_file}")

def run_elevated_plus_maze(base_dir):
    """
    Simulates an elevated plus maze test to assess anxiety-like behavior.
    """
    # --- File Paths ---
    scores_file = os.path.join(base_dir, 'phenotype', 'predicted_scores.csv')
    output_csv_file = os.path.join(base_dir, 'behavior', 'elevated_plus_maze.csv')
    
    # --- Read Phenotype Scores ---
    phenotype = {}
    try:
        with open(scores_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = row['phenotype']
                if key == 'anxiety':
                    phenotype['anxiety_score'] = float(row['score'])
                elif key == 'memory':
                    phenotype['memory_score'] = float(row['score'])
                elif key == 'obesity':
                    phenotype['obesity_score'] = float(row['score'])
                elif key == 'tumor_risk':
                    phenotype['tumor_risk'] = float(row['score'])
                elif key == 'aggression':
                    phenotype['aggression_score'] = float(row['score'])
                else:
                    phenotype[f"{key}_score"] = float(row['score'])
    except FileNotFoundError:
        phenotype = {
            'anxiety_score': 0.5,
            'memory_score': 0.5,
            'obesity_score': 0.5,
            'tumor_risk': 0.5,
            'aggression_score': 0.5
        }
    
    print(f"Running elevated plus maze with phenotypes: {phenotype}")
    
    # Initialize the enhanced behavioral simulator
    behavioral_simulator = BehavioralSimulator()
    
    # Run the simulation
    behavior_data = behavioral_simulator.simulate_elevated_plus_maze(phenotype)
    
    # Save as CSV
    behavior_df = pd.DataFrame(behavior_data)
    behavior_df.to_csv(output_csv_file, index=False)
    
    # Calculate and print metrics
    metrics = behavioral_simulator.calculate_behavioral_metrics(behavior_data, "elevated_plus_maze")
    print(f"Elevated plus maze metrics: {metrics}")
    print(f"Elevated plus maze simulation complete. Data saved to {output_csv_file}")

def run_novel_object_test(base_dir):
    """
    Simulates a novel object test to assess recognition memory.
    """
    # --- File Paths ---
    scores_file = os.path.join(base_dir, 'phenotype', 'predicted_scores.csv')
    output_csv_file = os.path.join(base_dir, 'behavior', 'novel_object_test.csv')
    
    # --- Read Phenotype Scores ---
    phenotype = {}
    try:
        with open(scores_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = row['phenotype']
                if key == 'anxiety':
                    phenotype['anxiety_score'] = float(row['score'])
                elif key == 'memory':
                    phenotype['memory_score'] = float(row['score'])
                elif key == 'obesity':
                    phenotype['obesity_score'] = float(row['score'])
                elif key == 'tumor_risk':
                    phenotype['tumor_risk'] = float(row['score'])
                elif key == 'aggression':
                    phenotype['aggression_score'] = float(row['score'])
                else:
                    phenotype[f"{key}_score"] = float(row['score'])
    except FileNotFoundError:
        phenotype = {
            'anxiety_score': 0.5,
            'memory_score': 0.5,
            'obesity_score': 0.5,
            'tumor_risk': 0.5,
            'aggression_score': 0.5
        }
    
    print(f"Running novel object test with phenotypes: {phenotype}")
    
    # Initialize the enhanced behavioral simulator
    behavioral_simulator = BehavioralSimulator()
    
    # Run the simulation
    behavior_data = behavioral_simulator.simulate_novel_object_test(phenotype)
    
    # Save as CSV
    behavior_df = pd.DataFrame(behavior_data)
    behavior_df.to_csv(output_csv_file, index=False)
    
    # Calculate and print metrics
    metrics = behavioral_simulator.calculate_behavioral_metrics(behavior_data, "novel_object")
    print(f"Novel object test metrics: {metrics}")
    print(f"Novel object test simulation complete. Data saved to {output_csv_file}")

if __name__ == '__main__':
    # Run all three behavioral tests
    base_dir = '.'
    run_open_field_test(base_dir)
    run_elevated_plus_maze(base_dir)
    run_novel_object_test(base_dir)
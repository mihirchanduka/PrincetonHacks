
import csv
import json
import os

def simulate_high_fat_diet(base_dir):
    """
    Simulates the effect of a high-fat diet on phenotype scores and biomarkers.
    """
    # Define file paths
    scores_file = os.path.join(base_dir, 'phenotype', 'scores.csv')
    biomarkers_file = os.path.join(base_dir, 'blood_work', 'biomarkers.json')
    
    output_scores_file = os.path.join(base_dir, 'phenotype', 'scores_after_diet.csv')
    output_biomarkers_file = os.path.join(base_dir, 'blood_work', 'biomarkers_after_diet.json')

    # --- Process Phenotype Scores ---
    phenotypes = {}
    with open(scores_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            phenotypes[row[0]] = float(row[1])

    # Apply rule: Increase obesity score by 20%
    if 'obesity' in phenotypes:
        phenotypes['obesity'] *= 1.2

    # Write updated scores
    with open(output_scores_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for phenotype, score in phenotypes.items():
            writer.writerow([phenotype, f"{score:.2f}"])

    print(f"Updated phenotype scores saved to {output_scores_file}")

    # --- Process Biomarkers ---
    with open(biomarkers_file, 'r') as f:
        biomarkers = json.load(f)

    # Apply rule: Increase cholesterol by 50%
    if 'cholesterol' in biomarkers:
        biomarkers['cholesterol']['value'] *= 1.5

    # Write updated biomarkers
    with open(output_biomarkers_file, 'w') as f:
        json.dump(biomarkers, f, indent=2)
        
    print(f"Updated biomarkers saved to {output_biomarkers_file}")

if __name__ == '__main__':
    # Assuming the script is run from the parent directory of 'synthetic_mouse_001'
    # or we can just specify the path directly.
    # For this execution, we'll assume the script is in 'synthetic_mouse_001'
    # and the base directory is the current directory.
    simulate_high_fat_diet('.')


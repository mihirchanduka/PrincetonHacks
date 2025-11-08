
import csv
from Bio import SeqIO

def predict_phenotype(base_dir):
    """
    Predicts phenotype scores based on a mutated genome using a simple rule-based model.
    """
    # Define file paths
    mutated_genome_file = f"{base_dir}/genome/mutated_genome.fasta"
    output_scores_file = f"{base_dir}/phenotype/predicted_scores.csv"

    # --- Baseline Phenotype Scores ---
    # In a real model, these would be predicted, but we'll start with a baseline.
    phenotypes = {
        "obesity": 0.5,
        "tumor_risk": 0.5,
        "anxiety": 0.5
    }

    # --- Read Mutated Genome ---
    mutated_genome = SeqIO.read(mutated_genome_file, "fasta")
    sequence = mutated_genome.seq

    # --- Apply Rules Based on Genotype ---
    # Rule 1: If mutation at position 10 is 'G', increase anxiety.
    # VCF was 1-based, sequence is 0-based.
    if sequence[9] == 'G':
        phenotypes["anxiety"] += 0.3
        print("Rule applied: Mutation at position 10 (C->G) detected. Increasing anxiety score.")

    # Rule 2: If mutation at position 20 is 'A', decrease tumor risk.
    if sequence[19] == 'A':
        phenotypes["tumor_risk"] -= 0.15
        print("Rule applied: Mutation at position 20 (T->A) detected. Decreasing tumor_risk score.")

    # --- Write Predicted Scores ---
    with open(output_scores_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['phenotype', 'score'])
        for phenotype, score in phenotypes.items():
            writer.writerow([phenotype, f"{score:.2f}"])

    print(f"\\nPredicted phenotype scores saved to {output_scores_file}")

if __name__ == '__main__':
    predict_phenotype('.')

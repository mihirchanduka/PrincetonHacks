"""
Enhanced ML-based Genome-Phenotype Prediction Model
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from Bio import SeqIO
from Bio.Seq import Seq
import pickle
from typing import Dict, List, Tuple, Optional
import os


class MLPhenotypePredictor:
    """
    A machine learning model to predict phenotypes from genomic sequences.
    This model can be trained on real mouse phenotype data to learn 
    genotype-phenotype relationships.
    """
    
    def __init__(self):
        self.models = {}  # Separate model for each phenotype
        self.feature_names = None
        self.is_trained = False
    
    def extract_features(self, genome_seq: str, k: int = 4) -> np.ndarray:
        """
        Extract k-mer frequency features from a genome sequence.
        
        Args:
            genome_seq: The genome sequence string
            k: Size of k-mers to extract
        
        Returns:
            Feature vector representing k-mer frequencies
        """
        # Generate all possible k-mers of size k
        bases = ['A', 'T', 'G', 'C']
        all_kmers = []
        
        def generate_kmers(length, bases, current=""):
            if length == 0:
                all_kmers.append(current)
                return
            for base in bases:
                generate_kmers(length - 1, bases, current + base)
        
        generate_kmers(k, bases)
        
        # Count occurrences of each k-mer in the sequence
        kmer_counts = {}
        for kmer in all_kmers:
            kmer_counts[kmer] = genome_seq.count(kmer)
        
        # Convert to frequency vector
        total_kmers = len(genome_seq) - k + 1
        if total_kmers <= 0:
            return np.zeros(len(all_kmers))
        
        freq_vector = np.array([kmer_counts[kmer] / total_kmers for kmer in all_kmers])
        
        # Also add some aggregate statistics
        a_count = genome_seq.count('A')
        t_count = genome_seq.count('T')
        g_count = genome_seq.count('G')
        c_count = genome_seq.count('C')
        
        gc_content = (g_count + c_count) / len(genome_seq)
        at_content = (a_count + t_count) / len(genome_seq)
        
        # Add these statistics to the feature vector
        stats = np.array([
            gc_content,
            at_content,
            a_count / len(genome_seq),
            t_count / len(genome_seq),
            g_count / len(genome_seq),
            c_count / len(genome_seq)
        ])
        
        return np.concatenate([freq_vector, stats])
    
    def train(self, genomes: List[str], phenotypes: List[Dict[str, float]], 
              test_size: float = 0.2) -> Dict[str, float]:
        """
        Train the phenotype prediction model.
        
        Args:
            genomes: List of genome sequences
            phenotypes: List of phenotype dictionaries
            test_size: Proportion of data to use for testing
        
        Returns:
            Dictionary of training/validation metrics
        """
        # Extract features from all genomes
        X = []
        for genome in genomes:
            features = self.extract_features(genome)
            X.append(features)
        
        X = np.array(X)
        
        # Prepare phenotype targets
        # Collect all phenotype names
        all_phenotype_names = set()
        for p in phenotypes:
            all_phenotype_names.update(p.keys())
        all_phenotype_names = sorted(list(all_phenotype_names))
        
        # Convert phenotypes to matrix format
        y = np.zeros((len(phenotypes), len(all_phenotype_names)))
        for i, p in enumerate(phenotypes):
            for j, name in enumerate(all_phenotype_names):
                y[i, j] = p.get(name, 0.0)  # Default to 0 if phenotype not present
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=42
        )
        
        # Train a model for each phenotype
        metrics = {}
        for j, phenotype_name in enumerate(all_phenotype_names):
            # Prepare target for this phenotype
            y_train_pheno = y_train[:, j]
            y_test_pheno = y_test[:, j]
            
            # Create and train model
            model = RandomForestRegressor(n_estimators=100, random_state=42)
            model.fit(X_train, y_train_pheno)
            
            # Evaluate
            y_pred = model.predict(X_test)
            mse = mean_squared_error(y_test_pheno, y_pred)
            
            # Store model and metrics
            self.models[phenotype_name] = model
            metrics[phenotype_name] = mse
        
        self.is_trained = True
        return metrics
    
    def predict(self, genome_seq: str) -> Dict[str, float]:
        """
        Predict phenotypes for a given genome sequence.
        
        Args:
            genome_seq: The genome sequence string
        
        Returns:
            Dictionary of predicted phenotypes
        """
        if not self.is_trained:
            raise ValueError("Model must be trained before making predictions")
        
        # Extract features
        X = self.extract_features(genome_seq).reshape(1, -1)
        
        # Predict each phenotype
        predictions = {}
        for phenotype_name, model in self.models.items():
            pred_value = model.predict(X)[0]
            # Ensure prediction is in [0, 1] range
            pred_value = max(0.0, min(1.0, pred_value))
            predictions[phenotype_name] = pred_value
        
        return predictions
    
    def save_model(self, filepath: str):
        """Save the trained model to disk."""
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
    
    @staticmethod
    def load_model(filepath: str) -> 'MLPhenotypePredictor':
        """Load a trained model from disk."""
        with open(filepath, 'rb') as f:
            return pickle.load(f)


# Example synthetic training data generation
def generate_synthetic_training_data(n_samples: int = 1000) -> Tuple[List[str], List[Dict[str, float]]]:
    """
    Generate synthetic training data for demonstration purposes.
    In a real implementation, this would come from real mouse databases.
    """
    import random
    
    genomes = []
    phenotypes = []
    
    for _ in range(n_samples):
        # Generate random genome
        genome_length = random.randint(500, 1500)
        bases = ['A', 'T', 'G', 'C']
        genome = ''.join(random.choices(bases, k=genome_length))
        genomes.append(genome)
        
        # Generate correlated phenotypes based on genome composition
        g_count = genome.count('G')
        c_count = genome.count('C')
        a_count = genome.count('A')
        t_count = genome.count('T')
        
        gc_content = (g_count + c_count) / genome_length
        
        # Create phenotype scores with some correlation to genome
        obesity = 0.5 + (gc_content - 0.5) * 0.3 + random.gauss(0, 0.1)
        anxiety = 0.5 + (genome.count('GGGG') / genome_length) * 0.4 + random.gauss(0, 0.1)
        tumor_risk = 0.5 + (genome.count('AAAA') / genome_length) * 0.3 + random.gauss(0, 0.1)
        memory = 0.5 - (genome.count('TTTT') / genome_length) * 0.2 + random.gauss(0, 0.1)
        aggression = 0.5 + (genome.count('CCCC') / genome_length) * 0.3 + random.gauss(0, 0.1)
        
        # Ensure values are in [0, 1] range
        obesity = max(0.0, min(1.0, obesity))
        anxiety = max(0.0, min(1.0, anxiety))
        tumor_risk = max(0.0, min(1.0, tumor_risk))
        memory = max(0.0, min(1.0, memory))
        aggression = max(0.0, min(1.0, aggression))
        
        phenotypes.append({
            'obesity': obesity,
            'anxiety': anxiety,
            'tumor_risk': tumor_risk,
            'memory': memory,
            'aggression': aggression
        })
    
    return genomes, phenotypes


if __name__ == "__main__":
    # Example usage
    print("Generating synthetic training data...")
    genomes, phenotypes = generate_synthetic_training_data(500)
    
    print("Training phenotype prediction model...")
    predictor = MLPhenotypePredictor()
    metrics = predictor.train(genomes, phenotypes)
    
    print("Training metrics (MSE for each phenotype):")
    for phenotype, mse in metrics.items():
        print(f"  {phenotype}: {mse:.4f}")
    
    # Test prediction on a new genome
    test_genome = ''.join(['A', 'T', 'G', 'C'] * 250)  # 1000 base genome
    predictions = predictor.predict(test_genome)
    
    print(f"\nPredictions for test genome:")
    for phenotype, value in predictions.items():
        print(f"  {phenotype}: {value:.3f}")
    
    # Save the model
    predictor.save_model("phenotype_prediction_model.pkl")
    print(f"\nModel saved to phenotype_prediction_model.pkl")
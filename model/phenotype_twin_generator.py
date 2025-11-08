"""
Phenotype Twin Generator - Creates synthetic mice that mimic real mouse phenotypes
using ML models trained on real mouse phenome data.
"""
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
import json
import random
from typing import Dict, List, Tuple, Optional
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error


class PhenotypeTwinGenerator:
    """
    A sophisticated model for generating synthetic mice that mimic real mouse
    phenotypes using ML models trained on real mouse phenome data.
    """
    
    def __init__(self):
        """Initialize the phenotype twin generator"""
        self.phenotype_models = {}  # Separate model for each phenotype
        self.genotype_to_phenotype_model = None
        self.phenotype_to_genotype_model = None
        self.real_mouse_database = None
        self.is_trained = False
        
    def load_real_mouse_data(self, data_path: Optional[str] = None) -> pd.DataFrame:
        """
        Load real mouse phenome data (for demonstration, we'll generate synthetic
        data that mimics real mouse characteristics)
        """
        # In a real implementation, this would load from real mouse phenome databases
        # like Mouse Phenome Database or IMPC data
        if data_path:
            # Load from provided path
            self.real_mouse_database = pd.read_csv(data_path)
        else:
            # Generate synthetic data that mimics real mouse phenotypes
            # This represents what we might find in real mouse databases
            n_samples = 1000
            mouse_data = {
                'mouse_id': [f'real_mouse_{i}' for i in range(n_samples)],
                'strain': [random.choice(['C57BL/6J', 'BALB/c', 'FVB/N', '129S1/SvImJ']) for _ in range(n_samples)],
                'sex': [random.choice(['male', 'female']) for _ in range(n_samples)],
                'age_weeks': [random.randint(8, 52) for _ in range(n_samples)],
                'weight': [random.gauss(25, 5) if i < n_samples//2 else random.gauss(30, 6) for i in range(n_samples)],
                'glucose': [random.gauss(150, 20) if i < n_samples//3 else random.gauss(200, 30) for i in range(n_samples)],
                'cholesterol': [random.gauss(200, 30) if i < n_samples//2 else random.gauss(250, 40) for i in range(n_samples)],
                'triglycerides': [random.gauss(150, 25) if i < n_samples//2 else random.gauss(200, 35) for i in range(n_samples)],
                'ALT': [random.gauss(35, 10) if i < n_samples*0.8 else random.gauss(60, 20) for i in range(n_samples)],
                'anxiety_score': [random.gauss(0.5, 0.2) if i < n_samples//2 else random.gauss(0.7, 0.2) for i in range(n_samples)],
                'memory_score': [random.gauss(0.6, 0.15) if i < n_samples//2 else random.gauss(0.4, 0.2) for i in range(n_samples)],
                'obesity_score': [random.gauss(0.4, 0.2) if i < n_samples//3 else random.gauss(0.8, 0.2) for i in range(n_samples)],
                'tumor_risk': [random.gauss(0.3, 0.1) if i < n_samples*0.9 else random.gauss(0.7, 0.2) for i in range(n_samples)],
                'aggression_score': [random.gauss(0.4, 0.15) if i < n_samples//2 else random.gauss(0.6, 0.2) for i in range(n_samples)]
            }
            self.real_mouse_database = pd.DataFrame(mouse_data)
            
        return self.real_mouse_database
    
    def train_phenotype_models(self, data: Optional[pd.DataFrame] = None) -> Dict[str, float]:
        """
        Train ML models to predict phenotypes from genotype features or other characteristics.
        """
        if data is None:
            data = self.real_mouse_database
            
        if data is None:
            raise ValueError("Must load or provide mouse data before training")
        
        # Define phenotypes we want to model
        phenotype_columns = [
            'anxiety_score', 'memory_score', 'obesity_score', 
            'tumor_risk', 'aggression_score', 'glucose', 
            'cholesterol', 'triglycerides', 'ALT'
        ]
        
        all_metrics = {}
        
        # Train models to predict each phenotype from other characteristics
        # In a real system, we would use genomic features
        feature_columns = ['age_weeks', 'weight', 'sex', 'strain']
        
        # Convert categorical variables to numeric
        data_encoded = data.copy()
        data_encoded = pd.get_dummies(data_encoded, columns=['sex', 'strain'])
        
        for phenotype in phenotype_columns:
            if phenotype in data_encoded.columns:
                # Prepare features (excluding the target phenotype)
                feature_names = [col for col in data_encoded.columns if col not in ['mouse_id'] + phenotype_columns or col == phenotype]
                feature_names = [col for col in feature_names if col != phenotype]
                
                X = data_encoded[feature_names].fillna(0)
                y = data_encoded[phenotype]
                
                # Split data
                X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
                
                # Train model
                model = RandomForestRegressor(n_estimators=50, random_state=42)
                model.fit(X_train, y_train)
                
                # Evaluate
                y_pred = model.predict(X_test)
                mse = mean_squared_error(y_test, y_pred)
                
                # Store model and metrics
                self.phenotype_models[phenotype] = model
                all_metrics[phenotype] = mse
        
        self.is_trained = True
        return all_metrics
    
    def create_synthetic_mouse_mimicking(
        self, 
        target_mouse_id: Optional[str] = None,
        target_mouse_profile: Optional[Dict] = None
    ) -> Tuple[SeqRecord, Dict[str, float]]:
        """
        Create a synthetic mouse that mimics the phenotype of a real mouse.
        
        Args:
            target_mouse_id: ID of the real mouse to mimic (from database)
            target_mouse_profile: Direct phenotype profile to mimic (alternative to target_mouse_id)
        
        Returns:
            Tuple of (synthetic_genome, synthetic_phenotype_profile)
        """
        if not self.is_trained:
            raise ValueError("Models must be trained before generating synthetic mice")
        
        # Get target phenotype profile
        if target_mouse_id and self.real_mouse_database is not None:
            target_row = self.real_mouse_database[self.real_mouse_database['mouse_id'] == target_mouse_id]
            if target_row.empty:
                raise ValueError(f"Mouse ID {target_mouse_id} not found in database")
            
            # Extract key phenotypes
            target_phenotype = {}
            for col in ['anxiety_score', 'memory_score', 'obesity_score', 'tumor_risk', 'aggression_score', 
                       'glucose', 'cholesterol', 'triglycerides', 'ALT']:
                if col in target_row.columns:
                    target_phenotype[col] = float(target_row[col].iloc[0])
        elif target_mouse_profile:
            target_phenotype = target_mouse_profile
        else:
            # If no target specified, create a random realistic mouse
            target_phenotype = self._generate_realistic_phenotype()
        
        # Generate synthetic genome based on target phenotype
        # In a real system, this would use actual genotype-phenotype relationships
        synthetic_genome = self._generate_genome_from_phenotype(target_phenotype)
        
        return synthetic_genome, target_phenotype
    
    def _generate_realistic_phenotype(self) -> Dict[str, float]:
        """Generate a realistic phenotype profile based on real mouse data patterns."""
        # Create phenotypes with realistic ranges and correlations
        obesity_score = random.uniform(0.2, 0.9)
        
        # Correlated phenotypes based on real mouse biology
        anxiety_score = 0.4 + (obesity_score - 0.5) * 0.3 + random.uniform(-0.1, 0.1)
        anxiety_score = max(0.0, min(1.0, anxiety_score))
        
        tumor_risk = 0.3 + obesity_score * 0.2 + random.uniform(-0.1, 0.2)
        tumor_risk = max(0.0, min(1.0, tumor_risk))
        
        memory_score = 0.7 - obesity_score * 0.3 + random.uniform(-0.1, 0.1)
        memory_score = max(0.0, min(1.0, memory_score))
        
        aggression_score = 0.5 + random.uniform(-0.2, 0.2)
        aggression_score = max(0.0, min(1.0, aggression_score))
        
        # Biochemical markers
        glucose = 150 + obesity_score * 100 + random.gauss(0, 20)
        cholesterol = 180 + obesity_score * 80 + random.gauss(0, 30)
        triglycerides = 140 + obesity_score * 100 + random.gauss(0, 25)
        ALT = 30 + tumor_risk * 50 + random.gauss(0, 15)
        
        return {
            'anxiety_score': anxiety_score,
            'memory_score': memory_score,
            'obesity_score': obesity_score,
            'tumor_risk': tumor_risk,
            'aggression_score': aggression_score,
            'glucose': glucose,
            'cholesterol': cholesterol,
            'triglycerides': triglycerides,
            'ALT': ALT
        }
    
    def _generate_genome_from_phenotype(self, phenotype: Dict[str, float], length: int = 1000) -> SeqRecord:
        """
        Generate a synthetic genome that would theoretically produce the target phenotype.
        In a real system, this would use actual genotype-phenotype mapping.
        """
        # Start with a random genome
        bases = ['A', 'T', 'G', 'C']
        sequence = ''.join(random.choices(bases, k=length))
        
        # Introduce specific patterns that might be associated with certain phenotypes
        # This is a simplified simulation - in reality, this would be based on known genetic variants
        seq_list = list(sequence)
        
        # Adjust certain positions based on phenotypes (simulated genetic associations)
        for i in range(0, len(seq_list), 100):  # Every 100th position
            if i < len(seq_list):
                # Obesity-related patterns
                if phenotype['obesity_score'] > 0.7:
                    seq_list[i] = 'G'  # Insert G at positions associated with metabolic traits
                elif phenotype['obesity_score'] < 0.3:
                    seq_list[i] = 'A'  # Insert A at positions associated with lean traits
                
                # Anxiety-related patterns
                if phenotype['anxiety_score'] > 0.7:
                    if i+1 < len(seq_list):
                        seq_list[i+1] = 'C'
                elif phenotype['anxiety_score'] < 0.3:
                    if i+1 < len(seq_list):
                        seq_list[i+1] = 'T'
        
        # Convert back to string
        modified_sequence = ''.join(seq_list)
        
        return SeqRecord(
            Seq(modified_sequence),
            id=f"synthetic_mouse_mimic_{hash(str(phenotype)) % 10000}",
            description="Synthetic mouse genome designed to mimic real mouse phenotype"
        )
    
    def generate_mimic_pair(
        self, 
        target_mouse_id: Optional[str] = None,
        target_mouse_profile: Optional[Dict] = None
    ) -> Tuple[Tuple[SeqRecord, Dict[str, float]], Tuple[SeqRecord, Dict[str, float]]]:
        """
        Generate a pair: (real-like mouse, synthetic twin that mimics it)
        
        Returns:
            Tuple of ((real_like_genome, real_like_phenotype), (synthetic_genome, synthetic_phenotype))
        """
        if target_mouse_id or target_mouse_profile:
            # Create a synthetic mouse mimicking a specific target
            synthetic_genome, synthetic_phenotype = self.create_synthetic_mouse_mimicking(
                target_mouse_id=target_mouse_id,
                target_mouse_profile=target_mouse_profile
            )
            
            # The "target" mouse can be the same as the synthetic (since we're mimicking)
            real_like_genome = self._generate_genome_from_phenotype(synthetic_phenotype, length=len(synthetic_genome.seq))
            real_like_phenotype = synthetic_phenotype  # In our model, the synthetic mimics the real
            
            return (real_like_genome, real_like_phenotype), (synthetic_genome, synthetic_phenotype)
        else:
            # Create a random realistic mouse and its mimic
            real_like_phenotype = self._generate_realistic_phenotype()
            real_like_genome = self._generate_genome_from_phenotype(real_like_phenotype)
            
            # Create a synthetic twin that mimics it (in this case they're the same since we're generating)
            synthetic_genome = self._generate_genome_from_phenotype(real_like_phenotype, length=len(real_like_genome.seq))
            synthetic_phenotype = real_like_phenotype
            
            return (real_like_genome, real_like_phenotype), (synthetic_genome, synthetic_phenotype)
    
    def save_model(self, filepath: str):
        """Save the trained model to disk."""
        model_data = {
            'phenotype_models': self.phenotype_models,
            'real_mouse_database': self.real_mouse_database,
            'is_trained': self.is_trained
        }
        with open(filepath, 'wb') as f:
            pickle.dump(model_data, f)
    
    def load_model(self, filepath: str):
        """Load a trained model from disk."""
        with open(filepath, 'rb') as f:
            model_data = pickle.load(f)
        
        self.phenotype_models = model_data['phenotype_models']
        self.real_mouse_database = model_data['real_mouse_database']
        self.is_trained = model_data['is_trained']


if __name__ == "__main__":
    # Example usage
    print("Creating Phenotype Twin Generator...")
    generator = PhenotypeTwinGenerator()
    
    # Load sample data
    print("Loading real mouse data...")
    generator.load_real_mouse_data()
    
    # Train the models
    print("Training phenotype prediction models...")
    metrics = generator.train_phenotype_models()
    print(f"Training MSE metrics: {metrics}")
    
    # Generate a mimic pair
    print("\nGenerating mimic pair...")
    (real_mouse_genome, real_mouse_phenotype), (synthetic_genome, synthetic_phenotype) = generator.generate_mimic_pair()
    
    print(f"Real-like mouse phenotype: {real_mouse_phenotype}")
    print(f"Synthetic mouse phenotype: {synthetic_phenotype}")
    print(f"Genome similarity: {sum(1 for a, b in zip(real_mouse_genome.seq, synthetic_genome.seq) if a == b) / len(real_mouse_genome.seq):.3f}")
    
    # Generate another pair with a specific target
    print("\nGenerating mimic pair with specific target...")
    target_phenotype = {
        'anxiety_score': 0.8,
        'memory_score': 0.3,
        'obesity_score': 0.9,
        'tumor_risk': 0.7,
        'aggression_score': 0.5,
        'glucose': 250.0,
        'cholesterol': 300.0,
        'triglycerides': 250.0,
        'ALT': 80.0
    }
    
    _, (specific_synthetic_genome, specific_synthetic_phenotype) = generator.generate_mimic_pair(
        target_mouse_profile=target_phenotype
    )
    
    print(f"Target phenotype: {target_phenotype}")
    print(f"Synthetic phenotype: {specific_synthetic_phenotype}")
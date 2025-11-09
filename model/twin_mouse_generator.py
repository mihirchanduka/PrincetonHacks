"""
Twin Mouse Generator - Core model for creating synthetic twin mice
"""
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
# import vcf  # Not currently used - commented out to avoid dependency issues
import json
import random
from typing import Dict, List, Tuple, Optional
import pickle
from .behavioral_simulator import BehavioralSimulator


class TwinMouseGenerator:
    """
    A sophisticated model for generating synthetic twin mice with realistic 
    genotype-phenotype relationships and behavioral patterns.
    """
    
    def __init__(self):
        """Initialize the twin mouse generator"""
        self.phenotype_model = None
        self.twin_similarity_matrix = {}
        self.known_mutations = {}
        self.behavioral_simulator = BehavioralSimulator()
        
    def generate_reference_genome(self, length: int = 1000, base_name: str = "chr19") -> SeqRecord:
        """Generate a random reference genome sequence"""
        bases = ['A', 'T', 'G', 'C']
        sequence = ''.join(random.choices(bases, k=length))
        return SeqRecord(
            Seq(sequence),
            id=base_name,
            description="Synthetic reference genome"
        )
    
    def create_twin_genomes(
        self, 
        reference_genome: SeqRecord, 
        mutation_rate: float = 0.001,
        twin_similarity: float = 0.95
    ) -> Tuple[SeqRecord, SeqRecord]:
        """
        Create twin genomes with controlled genetic similarity
        
        Args:
            reference_genome: The base genome to modify
            mutation_rate: Overall mutation rate
            twin_similarity: How similar the twins should be (0-1)
            
        Returns:
            A tuple of two SeqRecord objects representing the twin genomes
        """
        # Convert to mutable sequences
        seq1 = MutableSeq(reference_genome.seq)
        seq2 = MutableSeq(reference_genome.seq)
        
        genome_length = len(seq1)
        
        # Generate shared mutations (for both twins)
        shared_mutation_positions = set()
        for i in range(genome_length):
            if random.random() < mutation_rate * twin_similarity:
                shared_mutation_positions.add(i)
        
        # Apply shared mutations to both
        for pos in shared_mutation_positions:
            original_base = seq1[pos]
            new_base = random.choice([b for b in ['A', 'T', 'G', 'C'] if b != original_base])
            seq1[pos] = new_base
            seq2[pos] = new_base  # Same mutation for both twins
        
        # Generate individual mutations (for each twin separately)
        individual_mutation_rate = mutation_rate * (1 - twin_similarity)
        
        # Twin 1 individual mutations
        for i in range(genome_length):
            if i not in shared_mutation_positions and random.random() < individual_mutation_rate:
                original_base = seq1[i]
                new_base = random.choice([b for b in ['A', 'T', 'G', 'C'] if b != original_base])
                seq1[i] = new_base
        
        # Twin 2 individual mutations
        for i in range(genome_length):
            if i not in shared_mutation_positions and random.random() < individual_mutation_rate:
                original_base = seq2[i]
                new_base = random.choice([b for b in ['A', 'T', 'G', 'C'] if b != original_base])
                seq2[i] = new_base
        
        # Create SeqRecord objects for each twin
        twin1 = SeqRecord(Seq(seq1), id="twin1_genome", description="Twin 1 genome")
        twin2 = SeqRecord(Seq(seq2), id="twin2_genome", description="Twin 2 genome")
        
        return twin1, twin2
    
    def predict_phenotype_from_genotype(
        self, 
        genome: SeqRecord,
        trained_model: Optional = None
    ) -> Dict[str, float]:
        """
        Predict phenotype from genotype using a trained model or rule-based approach.
        
        Args:
            genome: The genome sequence
            trained_model: Optional trained ML model for prediction
            
        Returns:
            A dictionary of phenotypes with their scores
        """
        if trained_model is not None:
            return self._predict_with_model(genome, trained_model)
        else:
            return self._predict_with_rules(genome)
    
    def _predict_with_rules(self, genome: SeqRecord) -> Dict[str, float]:
        """
        Rule-based phenotype prediction for demonstration purposes.
        In a real implementation, this would use machine learning.
        """
        sequence = str(genome.seq)
        genome_length = len(sequence)
        
        # Initialize baseline phenotypes
        phenotypes = {
            "obesity": 0.5,
            "tumor_risk": 0.5,
            "anxiety": 0.5,
            "aggression": 0.5,
            "memory": 0.5
        }
        
        # Count specific motifs in the genome that might affect phenotypes
        # This is a simplified example - in reality, you'd train on real data
        g_content = sequence.count('G') / genome_length
        c_content = sequence.count('C') / genome_length
        gc_content = (g_content + c_content) / 2.0
        
        # Rule: Higher GC content might affect anxiety
        phenotypes["anxiety"] += (gc_content - 0.5) * 0.3
        phenotypes["anxiety"] = max(0.0, min(1.0, phenotypes["anxiety"]))
        
        # Look for specific k-mers that might affect phenotypes
        anxiety_kmers = ['GGGG', 'CCCC', 'GCGC', 'ATAT']
        anxiety_score_modifier = 0
        
        for kmer in anxiety_kmers:
            anxiety_score_modifier += sequence.count(kmer) * 0.01
        
        phenotypes["anxiety"] += anxiety_score_modifier
        phenotypes["anxiety"] = max(0.0, min(1.0, phenotypes["anxiety"]))
        
        # Obesity-related patterns
        obesity_kmers = ['AAAA', 'TTTT', 'AATT', 'CCAA']
        obesity_score_modifier = 0
        
        for kmer in obesity_kmers:
            obesity_score_modifier += sequence.count(kmer) * 0.01
        
        phenotypes["obesity"] += obesity_score_modifier
        phenotypes["obesity"] = max(0.0, min(1.0, phenotypes["obesity"]))
        
        # Memory-related patterns
        memory_kmers = ['AGCT', 'TCGA', 'GCTA', 'ATCG']
        memory_score_modifier = 0
        
        for kmer in memory_kmers:
            memory_score_modifier += sequence.count(kmer) * 0.01
        
        phenotypes["memory"] -= memory_score_modifier  # Higher score = better memory
        phenotypes["memory"] = max(0.0, min(1.0, phenotypes["memory"]))
        
        return phenotypes
    
    def _predict_with_model(self, genome: SeqRecord, trained_model) -> Dict[str, float]:
        """
        Use a trained ML model to predict phenotypes.
        This is a placeholder for when we have a real trained model.
        """
        # In a real implementation, this would:
        # 1. Extract features from the genome sequence
        # 2. Apply the trained model to predict phenotypes
        # 3. Return the predicted phenotypes
        return self._predict_with_rules(genome)
    
    def simulate_behavior(
        self, 
        phenotype: Dict[str, float], 
        duration: int = 300,
        test_type: str = "open_field"
    ) -> List[Dict[str, float]]:
        """
        Simulate behavioral data based on phenotype scores using the advanced behavioral simulator.
        
        Args:
            phenotype: Dictionary of phenotype scores
            duration: Duration of simulation in seconds
            test_type: Type of behavioral test to run ("open_field", "elevated_plus_maze", "novel_object")
            
        Returns:
            List of behavioral data points with x, y, time, and zone information
        """
        # Use the new BehavioralSimulator for more realistic behavior
        return self.behavioral_simulator.simulate_behavior(phenotype, duration, test_type)
    
    def simulate_experiment(
        self, 
        phenotype: Dict[str, float], 
        experiment_type: str = "high_fat_diet",
        treatment_params: Dict = None
    ) -> Dict[str, float]:
        """
        Simulate an experiment (e.g., drug treatment, diet change) and modify phenotypes.
        
        Args:
            phenotype: Original phenotype scores
            experiment_type: Type of experiment to simulate
            treatment_params: Parameters for the treatment
        
        Returns:
            Modified phenotype scores after experiment
        """
        if treatment_params is None:
            treatment_params = {}
        
        # Create a copy of the original phenotypes
        modified_phenotype = phenotype.copy()
        
        # Handle both naming conventions: 'obesity_score' and 'obesity'
        obesity_key = 'obesity' if 'obesity' in modified_phenotype else 'obesity_score' if 'obesity_score' in modified_phenotype else None
        anxiety_key = 'anxiety' if 'anxiety' in modified_phenotype else 'anxiety_score' if 'anxiety_score' in modified_phenotype else None
        memory_key = 'memory' if 'memory' in modified_phenotype else 'memory_score' if 'memory_score' in modified_phenotype else None
        tumor_key = 'tumor_risk' if 'tumor_risk' in modified_phenotype else None
        aggression_key = 'aggression' if 'aggression' in modified_phenotype else 'aggression_score' if 'aggression_score' in modified_phenotype else None
        
        if experiment_type == "high_fat_diet":
            # Increase obesity score
            obesity_increase = treatment_params.get('obesity_increase', 0.2)
            if obesity_key and obesity_key in modified_phenotype:
                modified_phenotype[obesity_key] = min(1.0, modified_phenotype[obesity_key] + obesity_increase)
            
        elif experiment_type == "drug_treatment":
            drug_name = treatment_params.get('drug_name', 'unknown')
            dose = treatment_params.get('dose_mg_kg', 10.0)
            
            if drug_name == "NewDrugX":
                # Example: NewDrugX affects tumor risk and anxiety
                if tumor_key and tumor_key in modified_phenotype:
                    modified_phenotype[tumor_key] = max(0.0, modified_phenotype[tumor_key] - (dose * 0.01))
                if anxiety_key and anxiety_key in modified_phenotype:
                    modified_phenotype[anxiety_key] = min(1.0, modified_phenotype[anxiety_key] + (dose * 0.005))
                
        elif experiment_type == "environmental":
            # Example: Environmental enrichment affects memory
            enrichment_factor = treatment_params.get('enrichment_factor', 0.1)
            if memory_key and memory_key in modified_phenotype:
                modified_phenotype[memory_key] = min(1.0, modified_phenotype[memory_key] + enrichment_factor)
        
        # Ensure all values are within [0, 1] range
        for key in modified_phenotype:
            if isinstance(modified_phenotype[key], (int, float)):
                modified_phenotype[key] = max(0.0, min(1.0, modified_phenotype[key]))
        
        return modified_phenotype
    
    def generate_blood_work(
        self, 
        phenotype: Dict[str, float]
    ) -> Dict[str, Dict[str, float]]:
        """
        Generate realistic blood work based on phenotype.
        
        Args:
            phenotype: Phenotype scores
            
        Returns:
            Blood work data in appropriate format
        """
        # Base values
        blood_work = {
            "glucose": {"value": 120.0, "units": "mg/dL"},
            "cholesterol": {"value": 200.0, "units": "mg/dL"},
            "triglycerides": {"value": 150.0, "units": "mg/dL"},
            "ALT": {"value": 35.0, "units": "U/L"},
            "AST": {"value": 40.0, "units": "U/L"}
        }
        
        # Modify based on phenotypes (using original naming for backward compatibility)
        obesity = phenotype.get('obesity', phenotype.get('obesity_score', 0.5))
        tumor_risk = phenotype.get('tumor_risk', 0.5)
        
        # Obesity affects glucose and cholesterol
        blood_work['glucose']['value'] += (obesity - 0.5) * 40.0
        blood_work['cholesterol']['value'] += (obesity - 0.5) * 60.0
        blood_work['triglycerides']['value'] += (obesity - 0.5) * 50.0
        
        # Tumor risk affects liver enzymes
        blood_work['ALT']['value'] += tumor_risk * 30.0
        blood_work['AST']['value'] += tumor_risk * 25.0
        
        # Ensure realistic ranges
        for key, value_info in blood_work.items():
            # Simple bounds checking
            if key == "glucose":
                value_info["value"] = max(70.0, min(300.0, value_info["value"]))
            elif key == "cholesterol":
                value_info["value"] = max(100.0, min(500.0, value_info["value"]))
            elif key == "triglycerides":
                value_info["value"] = max(50.0, min(1000.0, value_info["value"]))
            elif key in ["ALT", "AST"]:
                value_info["value"] = max(10.0, min(200.0, value_info["value"]))
        
        return blood_work


# Example usage
if __name__ == "__main__":
    # Create a twin mouse generator
    generator = TwinMouseGenerator()
    
    # Generate a reference genome
    ref_genome = generator.generate_reference_genome(length=1000)
    
    # Create twin genomes with high similarity
    twin1_genome, twin2_genome = generator.create_twin_genomes(
        ref_genome, 
        mutation_rate=0.01, 
        twin_similarity=0.9
    )
    
    # Predict phenotypes for both twins
    twin1_phenotype = generator.predict_phenotype_from_genotype(twin1_genome)
    twin2_phenotype = generator.predict_phenotype_from_genotype(twin2_genome)
    
    print("Twin 1 Phenotype:", twin1_phenotype)
    print("Twin 2 Phenotype:", twin2_phenotype)
    
    # Simulate behavior for both twins
    twin1_behavior = generator.simulate_behavior(twin1_phenotype, duration=50)
    twin2_behavior = generator.simulate_behavior(twin2_phenotype, duration=50)
    
    print(f"Twin 1 Behavior data points: {len(twin1_behavior)}")
    print(f"Twin 2 Behavior data points: {len(twin2_behavior)}")
    
    # Generate blood work for both twins
    twin1_blood = generator.generate_blood_work(twin1_phenotype)
    twin2_blood = generator.generate_blood_work(twin2_phenotype)
    
    print("Twin 1 Blood work:", twin1_blood)
    print("Twin 2 Blood work:", twin2_blood)
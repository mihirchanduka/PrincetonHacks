"""
Twin Comparison and Analysis Module
"""
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt
import seaborn as sns
from .twin_mouse_generator import TwinMouseGenerator


class TwinComparator:
    """
    Module for comparing twin mice and analyzing their similarities and differences.
    """
    
    def __init__(self):
        pass
    
    def calculate_genetic_similarity(
        self, 
        genome1: str, 
        genome2: str
    ) -> float:
        """
        Calculate genetic similarity between two genomes.
        
        Args:
            genome1: First genome sequence
            genome2: Second genome sequence (must be same length)
        
        Returns:
            Similarity score between 0 and 1
        """
        if len(genome1) != len(genome2):
            raise ValueError("Genomes must be the same length for comparison")
        
        if len(genome1) == 0:
            return 1.0  # Empty genomes are considered identical
        
        matches = sum(1 for a, b in zip(genome1, genome2) if a == b)
        return matches / len(genome1)
    
    def calculate_phenotype_similarity(
        self, 
        phenotype1: Dict[str, float], 
        phenotype2: Dict[str, float]
    ) -> float:
        """
        Calculate similarity between two phenotype profiles.
        
        Args:
            phenotype1: First phenotype dictionary
            phenotype2: Second phenotype dictionary
        
        Returns:
            Similarity score between 0 and 1
        """
        all_keys = set(phenotype1.keys()) | set(phenotype2.keys())
        
        if not all_keys:
            return 1.0  # If both are empty, they are identical
        
        # Calculate the difference for each phenotype and average
        differences = []
        for key in all_keys:
            val1 = phenotype1.get(key, 0.0)
            val2 = phenotype2.get(key, 0.0)
            diff = abs(val1 - val2)
            differences.append(diff)
        
        avg_diff = sum(differences) / len(differences)
        
        # Convert difference to similarity (0 diff = 1 similarity, 1 diff = 0 similarity)
        return 1.0 - avg_diff
    
    def compare_behavior(
        self, 
        behavior1: List[Dict], 
        behavior2: List[Dict]
    ) -> Dict[str, float]:
        """
        Compare behavioral patterns between two mice.
        
        Args:
            behavior1: First mouse's behavior data
            behavior2: Second mouse's behavior data
        
        Returns:
            Dictionary of behavioral similarity metrics
        """
        if not behavior1 or not behavior2:
            return {
                'center_time_similarity': 1.0,
                'distance_similarity': 1.0,
                'speed_similarity': 1.0,
                'zone_transitions_similarity': 1.0
            }
        
        if len(behavior1) != len(behavior2):
            # If different lengths, compare up to shorter length
            min_len = min(len(behavior1), len(behavior2))
            behavior1 = behavior1[:min_len]
            behavior2 = behavior2[:min_len]
        
        # Calculate similarity metrics
        center_time_similarity = self._compare_center_time(behavior1, behavior2)
        distance_similarity = self._compare_distance_traveled(behavior1, behavior2)
        speed_similarity = self._compare_movement_speed(behavior1, behavior2)
        
        # Calculate zone transitions similarity (exploratory behavior)
        zone_transitions_similarity = self._compare_zone_transitions(behavior1, behavior2)
        
        return {
            'center_time_similarity': center_time_similarity,
            'distance_similarity': distance_similarity,
            'speed_similarity': speed_similarity,
            'zone_transitions_similarity': zone_transitions_similarity
        }
    
    def _compare_center_time(self, behavior1: List[Dict], behavior2: List[Dict]) -> float:
        """Compare time spent in center between two behaviors."""
        center_time1 = sum(1 for b in behavior1 if b['zone'] == 'center')
        center_time2 = sum(1 for b in behavior2 if b['zone'] == 'center')
        
        max_time = max(len(behavior1), len(behavior2))
        if max_time == 0:
            return 1.0
        
        prop1 = center_time1 / len(behavior1)
        prop2 = center_time2 / len(behavior2)
        
        return 1.0 - abs(prop1 - prop2)
    
    def _compare_distance_traveled(self, behavior1: List[Dict], behavior2: List[Dict]) -> float:
        """Compare total distance traveled between two behaviors."""
        def calculate_distance(behavior_list):
            if len(behavior_list) < 2:
                return 0.0
            
            total_dist = 0.0
            for i in range(1, len(behavior_list)):
                x1, y1 = behavior_list[i-1]['x'], behavior_list[i-1]['y']
                x2, y2 = behavior_list[i]['x'], behavior_list[i]['y']
                dist = ((x2-x1)**2 + (y2-y1)**2)**0.5
                total_dist += dist
            return total_dist
        
        dist1 = calculate_distance(behavior1)
        dist2 = calculate_distance(behavior2)
        
        # Normalize by the maximum distance
        max_dist = max(dist1, dist2)
        if max_dist == 0:
            return 1.0
        
        # Similarity based on relative difference
        return 1.0 - abs(dist1 - dist2) / max_dist
    
    def _compare_movement_speed(self, behavior1: List[Dict], behavior2: List[Dict]) -> float:
        """Compare average movement speed between two behaviors."""
        def calculate_avg_speed(behavior_list):
            if len(behavior_list) < 2:
                return 0.0
            
            total_speed = 0.0
            for i in range(1, len(behavior_list)):
                x1, y1 = behavior_list[i-1]['x'], behavior_list[i-1]['y']
                x2, y2 = behavior_list[i]['x'], behavior_list[i]['y']
                dist = ((x2-x1)**2 + (y2-y1)**2)**0.5
                total_speed += dist  # Assuming 1 time unit between points
            return total_speed / (len(behavior_list) - 1)
        
        speed1 = calculate_avg_speed(behavior1)
        speed2 = calculate_avg_speed(behavior2)
        
        # Normalize by the maximum speed
        max_speed = max(speed1, speed2)
        if max_speed == 0:
            return 1.0
        
        # Similarity based on relative difference
        return 1.0 - abs(speed1 - speed2) / max_speed
    
    def _compare_zone_transitions(self, behavior1: List[Dict], behavior2: List[Dict]) -> float:
        """Compare zone transition patterns between two behaviors (exploratory behavior)."""
        def count_zone_transitions(behavior_list):
            if len(behavior_list) < 2:
                return 0
            
            transitions = 0
            for i in range(1, len(behavior_list)):
                if behavior_list[i].get('zone', '') != behavior_list[i-1].get('zone', ''):
                    transitions += 1
            return transitions
        
        transitions1 = count_zone_transitions(behavior1)
        transitions2 = count_zone_transitions(behavior2)
        
        # Normalize by maximum possible transitions
        max_transitions = max(len(behavior1), len(behavior2))
        if max_transitions == 0:
            return 1.0
        
        # Calculate similarity - if both have similar transition rates, similarity should be high
        max_possible_transitions = max_transitions - 1  # Max possible transitions in a sequence
        if max_possible_transitions == 0:
            return 1.0
        
        # Normalize transition counts
        norm_transitions1 = transitions1 / max_possible_transitions
        norm_transitions2 = transitions2 / max_possible_transitions
        
        # Similarity based on normalized transition rates
        return 1.0 - abs(norm_transitions1 - norm_transitions2)
    
    def generate_comparison_report(
        self, 
        twin1_data: Dict, 
        twin2_data: Dict
    ) -> Dict[str, any]:
        """
        Generate a comprehensive comparison report between twin mice.
        
        Args:
            twin1_data: Complete data for first twin (genome, phenotype, behavior, etc.)
            twin2_data: Complete data for second twin (genome, phenotype, behavior, etc.)
        
        Returns:
            Dictionary containing comprehensive comparison results
        """
        # Extract components
        genome1 = str(twin1_data['genome'].seq) if 'genome' in twin1_data else ""
        genome2 = str(twin2_data['genome'].seq) if 'genome' in twin2_data else ""
        
        phenotype1 = twin1_data.get('phenotype', {})
        phenotype2 = twin2_data.get('phenotype', {})
        
        behavior1 = twin1_data.get('behavior', [])
        behavior2 = twin2_data.get('behavior', [])
        
        # Calculate similarities
        genetic_similarity = self.calculate_genetic_similarity(genome1, genome2) if genome1 and genome2 else 1.0
        phenotype_similarity = self.calculate_phenotype_similarity(phenotype1, phenotype2)
        
        behavior_similarity = self.compare_behavior(behavior1, behavior2) if behavior1 and behavior2 else {}
        
        # Create comprehensive report
        report = {
            'genetic_similarity': genetic_similarity,
            'phenotype_similarity': phenotype_similarity,
            'behavior_similarity': behavior_similarity,
            'overall_similarity': (genetic_similarity + phenotype_similarity) / 2,  # Simple average
            'twin1_summary': {
                'genome_length': len(genome1) if genome1 else 0,
                'phenotypes': phenotype1,
                'behavior_duration': len(behavior1)
            },
            'twin2_summary': {
                'genome_length': len(genome2) if genome2 else 0,
                'phenotypes': phenotype2,
                'behavior_duration': len(behavior2)
            }
        }
        
        return report
    
    def visualize_comparison(self, report: Dict[str, any], save_path: str = None):
        """
        Create visualizations comparing the twins.
        
        Args:
            report: The comparison report generated by generate_comparison_report
            save_path: Optional path to save the visualization
        """
        # Set up the plot
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Twin Mouse Comparison Report', fontsize=16)
        
        # 1. Similarity metrics
        metrics = ['Genetic', 'Phenotype', 'Overall']
        values = [
            report['genetic_similarity'],
            report['phenotype_similarity'],
            report['overall_similarity']
        ]
        
        axes[0, 0].bar(metrics, values)
        axes[0, 0].set_title('Similarity Metrics')
        axes[0, 0].set_ylabel('Similarity (0-1)')
        axes[0, 0].set_ylim([0, 1])
        
        # Add value labels on bars
        for i, v in enumerate(values):
            axes[0, 0].text(i, v + 0.01, f'{v:.2f}', ha='center')
        
        # 2. Phenotype comparison
        pheno1 = report['twin1_summary']['phenotypes']
        pheno2 = report['twin2_summary']['phenotypes']
        
        phenotypes = list(set(pheno1.keys()) | set(pheno2.keys()))
        twin1_values = [pheno1.get(p, 0) for p in phenotypes]
        twin2_values = [pheno2.get(p, 0) for p in phenotypes]
        
        x = np.arange(len(phenotypes))
        width = 0.35
        
        axes[0, 1].bar(x - width/2, twin1_values, width, label='Twin 1', alpha=0.8)
        axes[0, 1].bar(x + width/2, twin2_values, width, label='Twin 2', alpha=0.8)
        axes[0, 1].set_xlabel('Phenotypes')
        axes[0, 1].set_ylabel('Score (0-1)')
        axes[0, 1].set_title('Phenotype Comparison')
        axes[0, 1].set_xticks(x)
        axes[0, 1].set_xticklabels(phenotypes, rotation=45)
        axes[0, 1].legend()
        
        # 3. Behavior comparison (Center time)
        behavior_sim = report.get('behavior_similarity', {})
        if behavior_sim:
            behavior_metrics = list(behavior_sim.keys())
            behavior_values = list(behavior_sim.values())
            
            axes[1, 0].bar(behavior_metrics, behavior_values)
            axes[1, 0].set_title('Behavior Similarity Metrics')
            axes[1, 0].set_ylabel('Similarity (0-1)')
            axes[1, 0].set_ylim([0, 1])
            axes[1, 0].tick_params(axis='x', rotation=45)
            
            # Add value labels on bars
            for i, v in enumerate(behavior_values):
                axes[1, 0].text(i, v + 0.01, f'{v:.2f}', ha='center')
        
        # 4. Combined similarity radar chart
        # Prepare data for radar chart
        categories = ['Genetic', 'Phenotype', 'Behavior']
        values_radar = [
            report['genetic_similarity'],
            report['phenotype_similarity'],
            sum(behavior_sim.values()) / len(behavior_sim) if behavior_sim else 0.5  # Avg behavior
        ]
        
        # Add the first value to close the circular graph:
        values_radar += values_radar[:1]
        angles = [n / float(len(categories)) * 2 * np.pi for n in range(len(categories))]
        angles += angles[:1]
        
        ax_radar = plt.subplot(2, 2, 4, projection='polar')
        ax_radar.plot(angles, values_radar, 'o-', linewidth=2)
        ax_radar.fill(angles, values_radar, alpha=0.25)
        ax_radar.set_xticks(angles[:-1])
        ax_radar.set_xticklabels(categories)
        ax_radar.set_ylim(0, 1)
        ax_radar.set_title('Overall Similarity Radar')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()


if __name__ == "__main__":
    # Example usage
    print("Creating twin comparison example...")
    
    # Create twins using the generator
    generator = TwinMouseGenerator()
    ref_genome = generator.generate_reference_genome(length=1000)
    twin1_genome, twin2_genome = generator.create_twin_genomes(
        ref_genome, 
        mutation_rate=0.005, 
        twin_similarity=0.85
    )
    
    # Predict phenotypes
    twin1_phenotype = generator.predict_phenotype_from_genome(twin1_genome)
    twin2_phenotype = generator.predict_phenotype_from_genome(twin2_genome)
    
    # Generate behavior
    twin1_behavior = generator.simulate_behavior(twin1_phenotype, duration=100)
    twin2_behavior = generator.simulate_behavior(twin2_phenotype, duration=100)
    
    # Create data dictionaries
    twin1_data = {
        'genome': twin1_genome,
        'phenotype': twin1_phenotype,
        'behavior': twin1_behavior
    }
    
    twin2_data = {
        'genome': twin2_genome,
        'phenotype': twin2_phenotype,
        'behavior': twin2_behavior
    }
    
    # Compare twins
    comparator = TwinComparator()
    report = comparator.generate_comparison_report(twin1_data, twin2_data)
    
    print("Comparison Report:")
    for key, value in report.items():
        print(f"  {key}: {value}")
    
    print("\nBehavior Similarity Details:")
    for key, value in report['behavior_similarity'].items():
        print(f"  {key}: {value:.3f}")
"""
Virtual Mouse Lab - Main Simulation Orchestrator
This module orchestrates the entire phenotypic twin mouse simulation process.
"""
import os
import json
import zipfile
import pandas as pd
from datetime import datetime
from typing import Dict, List, Optional
from Bio import SeqIO

from model.phenotype_twin_generator import PhenotypeTwinGenerator
from model.twin_mouse_generator import TwinMouseGenerator
from model.ml_phenotype_predictor import MLPhenotypePredictor
from model.twin_comparator import TwinComparator
from data.unified_data_model import MouseSubject, GenomicData, PhenotypicData, BehavioralData, BiomarkerData, Treatment


class VirtualMouseLab:
    """
    Main orchestrator for the virtual mouse lab that handles:
    - Phenotypic twin generation (synthetic mouse mimicking real mouse)
    - Phenotype prediction
    - Behavioral simulation
    - Experiment simulation
    - Data packaging
    """
    
    def __init__(self):
        self.phenotype_generator = PhenotypeTwinGenerator()
        self.behavior_generator = TwinMouseGenerator()  # Keep for behavior simulation
        self.predictor = MLPhenotypePredictor()
        self.comparator = TwinComparator()
        self.output_dir = "synthetic_mimic_mice"
    
    def run_single_mouse_simulation(
        self,
        target_mouse_id: Optional[str] = None,
        target_mouse_profile: Optional[Dict] = None,
        experiments: List[Dict] = None,
        include_behavior: bool = True,
        use_real_genome: bool = True
    ) -> str:
        """
        Run a complete single mouse simulation using real genome data.
        This allows researchers to predict experiment results using real world data.
        
        Args:
            target_mouse_id: ID of real mouse to mimic (from database)
            target_mouse_profile: Direct phenotype profile to mimic
            experiments: List of experiments to run
            include_behavior: Whether to include behavioral testing
            use_real_genome: Whether to use real genome data from dataset
        
        Returns:
            Path to the generated data package
        """
        if experiments is None:
            experiments = []
        
        # Load real mouse data and train models
        print("Loading real mouse data and training models...")
        self.phenotype_generator.load_real_mouse_data()
        self.phenotype_generator.train_phenotype_models()
        
        if use_real_genome:
            # Use real genome data from the dataset
            from Bio import SeqIO
            import os
            
            # Load real genome from dataset
            real_genome_path = "genome/reference.fasta"
            if os.path.exists(real_genome_path):
                real_genome = next(SeqIO.parse(real_genome_path, "fasta"))
                print(f"Using real genome data from {real_genome_path}")
            else:
                # Fallback: generate a genome from the target phenotype
                print("Real genome file not found, generating from phenotype...")
                real_genome, target_phenotype = self.phenotype_generator.create_synthetic_mouse_mimicking(
                    target_mouse_id=target_mouse_id,
                    target_mouse_profile=target_mouse_profile
                )
        else:
            # Generate synthetic genome as before
            real_genome, target_phenotype = self.phenotype_generator.create_synthetic_mouse_mimicking(
                target_mouse_id=target_mouse_id,
                target_mouse_profile=target_mouse_profile
            )
        
        # If we don't have a target phenotype yet (e.g. when using real genome file), generate it
        if 'target_phenotype' not in locals() or target_phenotype is None:
            if target_mouse_id or target_mouse_profile:
                # Use the phenotype from the target
                if target_mouse_id:
                    target_row = self.phenotype_generator.real_mouse_database[
                        self.phenotype_generator.real_mouse_database['mouse_id'] == target_mouse_id
                    ]
                    if not target_row.empty:
                        target_phenotype = {}
                        for col in ['anxiety_score', 'memory_score', 'obesity_score', 'tumor_risk', 'aggression_score', 
                                   'glucose', 'cholesterol', 'triglycerides', 'ALT']:
                            if col in target_row.columns:
                                target_phenotype[col] = float(target_row[col].iloc[0])
                elif target_mouse_profile:
                    target_phenotype = target_mouse_profile
            else:
                # Generate a realistic phenotype if no target is specified
                target_phenotype = self.phenotype_generator._generate_realistic_phenotype()
        
        # Create unified mouse subject
        mouse_subject = MouseSubject()
        
        # Populate genomic data
        genomic_data = GenomicData()
        genomic_data.reference_genome = real_genome
        mouse_subject.genomic_data = genomic_data
        
        # Populate phenotypic data
        phenotypic_data = PhenotypicData()
        phenotypic_data.anxiety_score = target_phenotype.get('anxiety_score')
        phenotypic_data.memory_score = target_phenotype.get('memory_score')
        phenotypic_data.obesity_score = target_phenotype.get('obesity_score')
        phenotypic_data.tumor_risk = target_phenotype.get('tumor_risk')
        phenotypic_data.aggression_score = target_phenotype.get('aggression_score')
        mouse_subject.phenotypic_data = phenotypic_data
        
        # Run experiments on the mouse
        experimental_results = {}
        mouse_phenotype_after_exp = target_phenotype.copy()
        
        for i, experiment in enumerate(experiments):
            exp_type = experiment.get("type", "control")
            exp_params = experiment.get("params", {})
            
            # Apply experiment to the mouse
            mouse_phenotype_after_exp = self.behavior_generator.simulate_experiment(
                mouse_phenotype_after_exp, exp_type, exp_params
            )
            
            # Add treatment to mouse subject
            treatment = Treatment(
                treatment_type=exp_type,
                compound=exp_params.get('drug_name', exp_type),
                dose=exp_params.get('dose_mg_kg', 0.0),
                dose_unit='mg/kg',
                duration=exp_params.get('days', 30)
            )
            mouse_subject.add_treatment(treatment)
            
            experimental_results[f"experiment_{i}"] = {
                "type": exp_type,
                "params": exp_params
            }
        
        # Update phenotypic data after experiments
        mouse_subject.phenotypic_data.anxiety_score = mouse_phenotype_after_exp.get('anxiety_score')
        mouse_subject.phenotypic_data.memory_score = mouse_phenotype_after_exp.get('memory_score')
        mouse_subject.phenotypic_data.obesity_score = mouse_phenotype_after_exp.get('obesity_score')
        mouse_subject.phenotypic_data.tumor_risk = mouse_phenotype_after_exp.get('tumor_risk')
        mouse_subject.phenotypic_data.aggression_score = mouse_phenotype_after_exp.get('aggression_score')
        
        # Generate behavioral data based on phenotype
        if include_behavior:
            # Add behavioral test results
            behavioral_data = BehavioralData()
            
            # Simulate open field test
            open_field_test = self.behavior_generator.simulate_behavior(mouse_phenotype_after_exp)
            behavioral_data.open_field_test = open_field_test
            behavioral_data.calculate_behavioral_metrics()
            
            mouse_subject.behavioral_data = behavioral_data
        
        # Generate blood work
        mouse_blood = self.behavior_generator.generate_blood_work(mouse_phenotype_after_exp)
        
        # Populate biomarker data
        biomarker_data = BiomarkerData()
        for key, value_info in mouse_blood.items():
            if isinstance(value_info, dict) and 'value' in value_info:
                biomarker_data.add_blood_work(key, value_info['value'], value_info.get('units', ''))
        
        mouse_subject.biomarker_data = biomarker_data
        
        # Link all data types together
        mouse_subject.link_data_types()
        
        # Package all results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = f"{self.output_dir}_single_mouse_{timestamp}"
        package_path = self._package_results_single_unified(
            mouse_subject, experimental_results, output_dir
        )
        
        return package_path
    
    def _package_results(
        self, 
        real_like_mouse_data: Dict, 
        synthetic_mouse_data: Dict, 
        comparison_report: Dict, 
        experimental_results: Dict,
        output_dir: str
    ) -> str:
        """
        Package all simulation results into a downloadable format (for twin simulation).
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Create subdirectories
        for subdir in ['genome', 'phenotype', 'behavior', 'blood_work', 'report']:
            os.makedirs(os.path.join(output_dir, subdir), exist_ok=True)
        
        # Save real-like mouse data
        SeqIO.write(real_like_mouse_data['genome'], 
                   os.path.join(output_dir, 'genome', 'real_like_mouse_genome.fasta'), 
                   'fasta')
        
        # Save real-like mouse phenotype
        phenotype_df = pd.DataFrame(
            list(real_like_mouse_data['phenotype'].items()), 
            columns=['phenotype', 'score']
        )
        phenotype_df.to_csv(
            os.path.join(output_dir, 'phenotype', 'real_like_mouse_scores.csv'), 
            index=False
        )
        
        # Save real-like mouse behavior
        if real_like_mouse_data['behavior']:
            behavior_df = pd.DataFrame(real_like_mouse_data['behavior'])
            behavior_df.to_csv(
                os.path.join(output_dir, 'behavior', 'real_like_mouse_open_field_test.csv'), 
                index=False
            )
        
        # Save real-like mouse blood work
        with open(os.path.join(output_dir, 'blood_work', 'real_like_mouse_biomarkers.json'), 'w') as f:
            json.dump(real_like_mouse_data['blood_work'], f, indent=2)
        
        # Save synthetic mouse data
        SeqIO.write(synthetic_mouse_data['genome'], 
                   os.path.join(output_dir, 'genome', 'synthetic_mouse_genome.fasta'), 
                   'fasta')
        
        # Save synthetic mouse phenotype
        phenotype_df = pd.DataFrame(
            list(synthetic_mouse_data['phenotype'].items()), 
            columns=['phenotype', 'score']
        )
        phenotype_df.to_csv(
            os.path.join(output_dir, 'phenotype', 'synthetic_mouse_scores.csv'), 
            index=False
        )
        
        # Save synthetic mouse behavior
        if synthetic_mouse_data['behavior']:
            behavior_df = pd.DataFrame(synthetic_mouse_data['behavior'])
            behavior_df.to_csv(
                os.path.join(output_dir, 'behavior', 'synthetic_mouse_open_field_test.csv'), 
                index=False
            )
        
        # Save synthetic mouse blood work
        with open(os.path.join(output_dir, 'blood_work', 'synthetic_mouse_biomarkers.json'), 'w') as f:
            json.dump(synthetic_mouse_data['blood_work'], f, indent=2)
        
        # Save comparison report
        with open(os.path.join(output_dir, 'report', 'phenotype_comparison_report.json'), 'w') as f:
            json.dump(comparison_report, f, indent=2, default=str)
        
        # Save experimental results
        with open(os.path.join(output_dir, 'report', 'experimental_results.json'), 'w') as f:
            json.dump(experimental_results, f, indent=2)
        
        # Create summary report
        summary = self._create_summary_report(
            real_like_mouse_data, synthetic_mouse_data, comparison_report, experimental_results
        )
        with open(os.path.join(output_dir, 'report', 'summary.txt'), 'w') as f:
            f.write(summary)
        
        # Create ZIP package
        zip_path = f"{output_dir}.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    arc_path = os.path.relpath(file_path, os.path.dirname(output_dir))
                    zipf.write(file_path, arc_path)
        
        # Clean up directory after zipping
        import shutil
        shutil.rmtree(output_dir)
        
        return zip_path

    def _package_results_single_unified(
        self, 
        mouse_subject: MouseSubject, 
        experimental_results: Dict,
        output_dir: str
    ) -> str:
        """
        Package single mouse simulation results using the unified data model.
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Create subdirectories
        for subdir in ['genome', 'phenotype', 'behavior', 'blood_work', 'report', 'unified_data']:
            os.makedirs(os.path.join(output_dir, subdir), exist_ok=True)
        
        # Save mouse genome using unified data model
        if mouse_subject.genomic_data.reference_genome:
            SeqIO.write(mouse_subject.genomic_data.reference_genome, 
                       os.path.join(output_dir, 'genome', 'mouse_genome.fasta'), 
                       'fasta')
        
        # Save mouse phenotype using unified data model
        phenotype_dict = mouse_subject.phenotypic_data.to_dict()
        phenotype_df = pd.DataFrame(
            [(k, v) for k, v in phenotype_dict.items() if v is not None], 
            columns=['phenotype', 'value']
        )
        phenotype_df.to_csv(
            os.path.join(output_dir, 'phenotype', 'mouse_phenotype_scores.csv'), 
            index=False
        )
        
        # Save mouse behavior using unified data model
        behavioral_data = mouse_subject.behavioral_data
        if behavioral_data.open_field_test:
            behavior_df = pd.DataFrame(behavioral_data.open_field_test)
            behavior_df.to_csv(
                os.path.join(output_dir, 'behavior', 'mouse_open_field_test.csv'), 
                index=False
            )
        
        # Also save other behavioral tests if available
        if behavioral_data.elevated_plus_maze:
            eplus_df = pd.DataFrame(behavioral_data.elevated_plus_maze)
            eplus_df.to_csv(
                os.path.join(output_dir, 'behavior', 'mouse_elevated_plus_maze.csv'),
                index=False
            )
            
        if behavioral_data.novel_object_test:
            novel_df = pd.DataFrame(behavioral_data.novel_object_test)
            novel_df.to_csv(
                os.path.join(output_dir, 'behavior', 'mouse_novel_object_test.csv'),
                index=False
            )
        
        # Save mouse biomarkers using unified data model
        with open(os.path.join(output_dir, 'blood_work', 'mouse_biomarkers.json'), 'w') as f:
            json.dump(mouse_subject.biomarker_data.to_dict(), f, indent=2)
        
        # Save the full unified data model
        with open(os.path.join(output_dir, 'unified_data', 'mouse_subject.json'), 'w') as f:
            json.dump(mouse_subject.to_dict(), f, indent=2)
        
        # Save experimental results
        with open(os.path.join(output_dir, 'report', 'experimental_results.json'), 'w') as f:
            json.dump(experimental_results, f, indent=2)
        
        # Create summary report using unified data
        summary = self._create_summary_report_single_unified(
            mouse_subject, experimental_results
        )
        with open(os.path.join(output_dir, 'report', 'summary.txt'), 'w') as f:
            f.write(summary)
        
        # Create ZIP package
        zip_path = f"{output_dir}.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    arc_path = os.path.relpath(file_path, os.path.dirname(output_dir))
                    zipf.write(file_path, arc_path)
        
        # Clean up directory after zipping
        import shutil
        shutil.rmtree(output_dir)
        
        return zip_path
    
    def _package_results_single(
        self, 
        mouse_data: Dict, 
        experimental_results: Dict,
        output_dir: str
    ) -> str:
        """
        Package single mouse simulation results into a downloadable format (legacy format).
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Create subdirectories
        for subdir in ['genome', 'phenotype', 'behavior', 'blood_work', 'report']:
            os.makedirs(os.path.join(output_dir, subdir), exist_ok=True)
        
        # Save mouse genome
        SeqIO.write(mouse_data['genome'], 
                   os.path.join(output_dir, 'genome', 'mouse_genome.fasta'), 
                   'fasta')
        
        # Save mouse phenotype
        phenotype_df = pd.DataFrame(
            list(mouse_data['phenotype'].items()), 
            columns=['phenotype', 'score']
        )
        phenotype_df.to_csv(
            os.path.join(output_dir, 'phenotype', 'mouse_phenotype_scores.csv'), 
            index=False
        )
        
        # Save mouse behavior
        if mouse_data['behavior']:
            behavior_df = pd.DataFrame(mouse_data['behavior'])
            behavior_df.to_csv(
                os.path.join(output_dir, 'behavior', 'mouse_open_field_test.csv'), 
                index=False
            )
        
        # Save mouse blood work
        with open(os.path.join(output_dir, 'blood_work', 'mouse_biomarkers.json'), 'w') as f:
            json.dump(mouse_data['blood_work'], f, indent=2)
        
        # Save experimental results
        with open(os.path.join(output_dir, 'report', 'experimental_results.json'), 'w') as f:
            json.dump(experimental_results, f, indent=2)
        
        # Create summary report
        summary = self._create_summary_report_single(
            mouse_data, experimental_results
        )
        with open(os.path.join(output_dir, 'report', 'summary.txt'), 'w') as f:
            f.write(summary)
        
        # Create ZIP package
        zip_path = f"{output_dir}.zip"
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    arc_path = os.path.relpath(file_path, os.path.dirname(output_dir))
                    zipf.write(file_path, arc_path)
        
        # Clean up directory after zipping
        import shutil
        shutil.rmtree(output_dir)
        
        return zip_path
    
    def _create_summary_report(
        self, 
        real_like_mouse_data: Dict, 
        synthetic_mouse_data: Dict, 
        comparison_report: Dict, 
        experimental_results: Dict
    ) -> str:
        """Create a human-readable summary report."""
        summary = f"""
Virtual Mouse Lab - Phenotype Twin Simulation Summary
=====================================================

Simulation Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

PHENOTYPIC TWIN COMPARISON SUMMARY:
-----------------------------------
Genetic Similarity: {comparison_report['genetic_similarity']:.3f}
Phenotype Similarity: {comparison_report['phenotype_similarity']:.3f}
Overall Similarity: {comparison_report['overall_similarity']:.3f}

REAL-LIKE MOUSE PHENOTYPES:
---------------------------
"""
        for phenotype, score in real_like_mouse_data['phenotype'].items():
            summary += f"  {phenotype}: {score:.3f}\n"
        
        summary += f"""
SYNTHETIC MOUSE PHENOTYPES (MIMIC):
-----------------------------------
"""
        for phenotype, score in synthetic_mouse_data['phenotype'].items():
            summary += f"  {phenotype}: {score:.3f}\n"
        
        if experimental_results:
            summary += f"""
EXPERIMENTAL RESULTS:
---------------------
"""
            for exp_name, exp_data in experimental_results.items():
                summary += f"  {exp_name}: {exp_data['type']} with {exp_data['params']}\n"
        
        summary += f"""
BEHAVIOR ANALYSIS:
------------------
"""
        behavior_sim = comparison_report.get('behavior_similarity', {})
        for metric, value in behavior_sim.items():
            summary += f"  {metric}: {value:.3f}\n"
        
        return summary

    def _create_summary_report_single_unified(
        self, 
        mouse_subject: MouseSubject, 
        experimental_results: Dict
    ) -> str:
        """Create a human-readable summary report for a single mouse using unified data model."""
        summary = f"""
Virtual Mouse Lab - Single Mouse Simulation Summary (Unified Data Model)
========================================================================

Simulation Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Mouse ID: {mouse_subject.mouse_id}
Strain: {mouse_subject.strain.value}
Sex: {mouse_subject.sex.value}
Age: {mouse_subject.age_weeks} weeks

UNIFIED DATA MODEL COMPONENTS:
------------------------------
- Genomic Data: Available
- Phenotypic Data: Available  
- Behavioral Data: {len(mouse_subject.behavioral_data.open_field_test) > 0}
- Biomarker Data: {len(mouse_subject.biomarker_data.blood_work) > 0}

MOUSE PHENOTYPES:
-----------------
"""
        phenotypic_dict = mouse_subject.phenotypic_data.to_dict()
        for key, value in phenotypic_dict.items():
            if value is not None and key not in ['phenotype_timeline', 'measurement_date', 'measurement_method']:
                summary += f"  {key}: {value:.3f}\n"
        
        if mouse_subject.treatments:
            summary += f"""
TREATMENTS APPLIED:
-------------------
"""
            for i, treatment in enumerate(mouse_subject.treatments):
                summary += f"  Treatment {i+1}: {treatment.compound} ({treatment.dose} {treatment.dose_unit})\n"
        
        if experimental_results:
            summary += f"""
EXPERIMENTAL RESULTS:
---------------------
"""
            for exp_name, exp_data in experimental_results.items():
                summary += f"  {exp_name}: {exp_data['type']} with {exp_data['params']}\n"
        
        behavioral_metrics = mouse_subject.behavioral_data.behavioral_metrics
        if behavioral_metrics:
            summary += f"""
BEHAVIORAL METRICS:
-------------------
"""
            for metric, value in behavioral_metrics.items():
                summary += f"  {metric}: {value:.3f}\n"
        
        if mouse_subject.biomarker_data.blood_work:
            summary += f"""
BIOMARKER DATA:
---------------
"""
            for biomarker, data in mouse_subject.biomarker_data.blood_work.items():
                summary += f"  {biomarker}: {data['value']} {data['units']}\n"
        
        summary += f"""
DATA INTEGRATION NOTES:
-----------------------
{len(mouse_subject.notes)} notes recorded during data integration.
"""
        
        return summary

    def _create_summary_report_single(
        self, 
        mouse_data: Dict, 
        experimental_results: Dict
    ) -> str:
        """Create a human-readable summary report for a single mouse."""
        summary = f"""
Virtual Mouse Lab - Single Mouse Simulation Summary
===================================================

Simulation Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Genome: {mouse_data['genome'].id}
Description: {mouse_data['genome'].description}

MOUSE PHENOTYPES:
-----------------
"""
        for phenotype, score in mouse_data['phenotype'].items():
            summary += f"  {phenotype}: {score:.3f}\n"
        
        if experimental_results:
            summary += f"""
EXPERIMENTAL RESULTS:
---------------------
"""
            for exp_name, exp_data in experimental_results.items():
                summary += f"  {exp_name}: {exp_data['type']} with {exp_data['params']}\n"
        
        if mouse_data.get('behavior'):
            summary += f"""
BEHAVIOR ANALYSIS:
------------------
Total behavior data points: {len(mouse_data['behavior'])}
Time span: {max([b['time'] for b in mouse_data['behavior']]) if mouse_data['behavior'] else 0} seconds
Center zone visits: {len([b for b in mouse_data['behavior'] if b['zone'] == 'center'])}
Periphery zone visits: {len([b for b in mouse_data['behavior'] if b['zone'] == 'periphery'])}
"""
        
        return summary


# Example usage
if __name__ == "__main__":
    print("Initializing Virtual Mouse Lab...")
    lab = VirtualMouseLab()
    
    # Define experiments to run
    experiments = [
        {
            "type": "high_fat_diet",
            "params": {"obesity_increase": 0.15}
        },
        {
            "type": "drug_treatment", 
            "params": {"drug_name": "NewDrugX", "dose_mg_kg": 20.0}
        }
    ]
    
    print("Running single mouse simulation using real genome data...")
    package_path = lab.run_single_mouse_simulation(
        target_mouse_profile={
            'anxiety_score': 0.6,
            'memory_score': 0.5,
            'obesity_score': 0.4,
            'tumor_risk': 0.3,
            'aggression_score': 0.4,
            'glucose': 180.0,
            'cholesterol': 220.0,
            'triglycerides': 160.0,
            'ALT': 40.0
        },
        experiments=experiments,
        include_behavior=True,
        use_real_genome=True  # Use real genome from dataset
    )
    
    print(f"Simulation complete! Results packaged in: {package_path}")
    
    # Unpack and show the zipped contents for verification
    import zipfile
    with zipfile.ZipFile(package_path, 'r') as zipf:
        print("\nContents of the package:")
        for file in zipf.namelist():
            print(f"  {file}")
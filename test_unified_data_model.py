"""
Test script for the unified data model implementation.
This script tests the integration of genomic, phenotypic, behavioral, and biomarker data.
"""
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from data.unified_data_model import MouseSubject, GenomicData, PhenotypicData, BehavioralData, BiomarkerData, Treatment, MouseStrain, Sex
from virtual_mouse_lab import VirtualMouseLab
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json

def test_unified_data_model():
    """Test the unified data model with sample data."""
    print("Testing unified data model...")
    
    # Create a mouse subject
    mouse = MouseSubject(
        mouse_id="test_mouse_001",
        strain=MouseStrain.C57BL_6J,
        sex=Sex.MALE,
        age_weeks=12
    )
    
    # Add genomic data
    from data.unified_data_model import GenomicVariantType
    
    genomic_seq = SeqRecord(
        Seq("ATCGATCGATCGATCGATCGATCGATCGATCGATCG"),
        id="test_genome_001",
        description="Test genome sequence"
    )
    mouse.genomic_data.reference_genome = genomic_seq
    mouse.genomic_data.add_variant(
        position=100,
        reference_allele="A",
        alternative_allele="G",
        variant_type=GenomicVariantType.SNP,
        quality_score=0.98
    )
    
    # Add phenotypic data
    mouse.phenotypic_data.anxiety_score = 0.6
    mouse.phenotypic_data.memory_score = 0.7
    mouse.phenotypic_data.obesity_score = 0.3
    mouse.phenotypic_data.tumor_risk = 0.2
    mouse.phenotypic_data.aggression_score = 0.4
    
    # Add behavioral data
    for i in range(5):
        mouse.behavioral_data.open_field_test.append({
            'x': 50.0 + i,
            'y': 50.0 + i,
            'time': i,
            'zone': 'center' if i > 2 else 'periphery',
            'velocity_x': 1.0,
            'velocity_y': 1.0,
            'speed': 1.41
        })
    
    # Calculate behavioral metrics
    mouse.behavioral_data.calculate_behavioral_metrics()
    
    # Add biomarker data
    mouse.biomarker_data.add_blood_work("glucose", 145.0, "mg/dL")
    mouse.biomarker_data.add_blood_work("cholesterol", 195.0, "mg/dL")
    mouse.biomarker_data.add_blood_work("triglycerides", 140.0, "mg/dL")
    mouse.biomarker_data.add_blood_work("ALT", 38.0, "U/L")
    mouse.biomarker_data.add_protein_level("BDNF", 1.1, "ng/mL")
    
    # Add a treatment
    treatment = Treatment(
        treatment_type="drug_treatment",
        compound="TestDrugX",
        dose=10.0,
        dose_unit="mg/kg",
        duration=14
    )
    mouse.add_treatment(treatment)
    
    # Link data types
    mouse.link_data_types()
    
    # Print the unified data
    print("Unified mouse data created successfully!")
    print(f"Mouse ID: {mouse.mouse_id}")
    print(f"Genomic variants: {len(mouse.genomic_data.variants)}")
    print(f"Phenotypic scores: {mouse.phenotypic_data.anxiety_score}, {mouse.phenotypic_data.memory_score}")
    print(f"Behavioral data points: {len(mouse.behavioral_data.open_field_test)}")
    print(f"Biomarkers: {list(mouse.biomarker_data.blood_work.keys())}")
    print(f"Treatments: {len(mouse.treatments)}")
    print(f"Notes: {len(mouse.notes)}")
    
    # Test serialization
    json_data = mouse.to_json()
    print(f"JSON serialization successful, length: {len(json_data)}")
    
    # Test deserialization
    reconstructed_mouse = MouseSubject.from_dict(json.loads(json_data))
    print(f"Reconstructed mouse ID: {reconstructed_mouse.mouse_id}")
    print(f"Reconstructed anxiety score: {reconstructed_mouse.phenotypic_data.anxiety_score}")
    print(f"Reconstructed treatments: {len(reconstructed_mouse.treatments)}")
    
    return True

def test_virtual_mouse_lab_integration():
    """Test integration with the virtual mouse lab."""
    print("\nTesting Virtual Mouse Lab integration...")
    
    try:
        # Initialize the lab
        lab = VirtualMouseLab()
        
        # Run a simulation with the unified data model
        package_path = lab.run_single_mouse_simulation(
            target_mouse_profile={
                'anxiety_score': 0.5,
                'memory_score': 0.6,
                'obesity_score': 0.4,
                'tumor_risk': 0.3,
                'aggression_score': 0.4,
                'glucose': 150.0,
                'cholesterol': 200.0,
                'triglycerides': 150.0,
                'ALT': 35.0
            },
            experiments=[{
                "type": "drug_treatment",
                "params": {
                    "drug_name": "TestCompound",
                    "dose_mg_kg": 15.0
                }
            }],
            include_behavior=True,
            use_real_genome=True
        )
        
        print(f"Simulation completed successfully! Package: {package_path}")
        return True
        
    except Exception as e:
        print(f"Error in Virtual Mouse Lab integration: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("Starting tests for unified data model implementation...")
    
    success1 = test_unified_data_model()
    success2 = test_virtual_mouse_lab_integration()
    
    if success1 and success2:
        print("\nAll tests passed! Unified data model implementation is working correctly.")
    else:
        print("\nSome tests failed. Please check the implementation.")
    
    print("\nTest completed.")
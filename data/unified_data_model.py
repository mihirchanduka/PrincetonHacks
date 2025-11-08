"""
Unified Data Model for Mouse Research Data
This module provides a comprehensive data model that connects genomic, phenotypic, 
behavioral, and biomarker data in a coherent way.
"""
from datetime import datetime
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from enum import Enum
import json
import uuid
from Bio.SeqRecord import SeqRecord


class MouseStrain(Enum):
    """Common mouse strains used in research"""
    C57BL_6J = "C57BL/6J"
    BALB_C = "BALB/c"
    FVB_N = "FVB/N"
    C57BL_6NJ = "C57BL/6NJ"
    C3H_HEJ = "C3H/HeJ"
    DBA_2J = "DBA/2J"
    OTHER = "Other"


class Sex(Enum):
    """Biological sex of the mouse"""
    MALE = "male"
    FEMALE = "female"
    UNKNOWN = "unknown"


class GenomicVariantType(Enum):
    """Types of genomic variants"""
    SNP = "SNP"
    INDEL = "INDEL"
    CNV = "CNV"
    SV = "SV"
    OTHER = "Other"


@dataclass
class GenomicData:
    """
    Represents genomic data for a mouse subject.
    """
    # Primary sequence data
    reference_genome: Optional[SeqRecord] = None
    alternate_genomes: List[SeqRecord] = field(default_factory=list)
    
    # Variants and mutations
    variants: List[Dict[str, Any]] = field(default_factory=list)
    specific_mutations: List[Dict[str, Any]] = field(default_factory=list)
    
    # Quality metrics
    coverage_metrics: Dict[str, float] = field(default_factory=dict)
    assembly_quality: str = ""
    
    # Metadata
    sequencing_date: Optional[datetime] = None
    sequencing_technology: str = ""
    
    def add_variant(self, position: int, reference_allele: str, 
                   alternative_allele: str, variant_type: GenomicVariantType, 
                   quality_score: float = 0.0):
        """Add a genomic variant to the mouse."""
        self.variants.append({
            "position": position,
            "reference_allele": reference_allele,
            "alternative_allele": alternative_allele,
            "type": variant_type.value,
            "quality_score": quality_score,
            "id": str(uuid.uuid4())
        })
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format for serialization."""
        result = {
            "variants": self.variants,
            "specific_mutations": self.specific_mutations,
            "coverage_metrics": self.coverage_metrics,
            "assembly_quality": self.assembly_quality,
            "sequencing_date": self.sequencing_date.isoformat() if self.sequencing_date else None,
            "sequencing_technology": self.sequencing_technology
        }
        
        if self.reference_genome:
            result["reference_genome"] = {
                "id": self.reference_genome.id,
                "description": self.reference_genome.description,
                "sequence": str(self.reference_genome.seq),
                "length": len(self.reference_genome.seq)
            }
        
        result["alternate_genomes"] = [
            {
                "id": genome.id,
                "description": genome.description,
                "sequence": str(genome.seq),
                "length": len(genome.seq)
            } for genome in self.alternate_genomes
        ]
        
        return result


@dataclass
class PhenotypicData:
    """
    Represents phenotypic data for a mouse subject.
    Includes both behavioral and physical phenotypes.
    """
    # Behavioral phenotypes
    anxiety_score: Optional[float] = None
    memory_score: Optional[float] = None
    aggression_score: Optional[float] = None
    depression_score: Optional[float] = None
    locomotor_activity: Optional[float] = None
    
    # Physical phenotypes
    obesity_score: Optional[float] = None
    tumor_risk: Optional[float] = None
    weight: Optional[float] = None
    length: Optional[float] = None
    coat_color: Optional[str] = None
    
    # Time-dependent phenotypes (for tracking changes over time)
    phenotype_timeline: List[Dict[str, Any]] = field(default_factory=list)
    
    # Metadata
    measurement_date: Optional[datetime] = None
    measurement_method: str = ""
    
    def add_phenotype_measurement(self, phenotype_name: str, value: float, 
                                date: Optional[datetime] = None, 
                                method: str = "", notes: str = ""):
        """Add a phenotypic measurement to the timeline."""
        if date is None:
            date = datetime.now()
        
        self.phenotype_timeline.append({
            "phenotype": phenotype_name,
            "value": value,
            "date": date.isoformat(),
            "method": method,
            "notes": notes
        })
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format for serialization."""
        return {
            "anxiety_score": self.anxiety_score,
            "memory_score": self.memory_score,
            "aggression_score": self.aggression_score,
            "depression_score": self.depression_score,
            "locomotor_activity": self.locomotor_activity,
            "obesity_score": self.obesity_score,
            "tumor_risk": self.tumor_risk,
            "weight": self.weight,
            "length": self.length,
            "coat_color": self.coat_color,
            "phenotype_timeline": self.phenotype_timeline,
            "measurement_date": self.measurement_date.isoformat() if self.measurement_date else None,
            "measurement_method": self.measurement_method
        }


@dataclass
class BehavioralData:
    """
    Represents behavioral test data for a mouse subject.
    """
    # Test results
    open_field_test: List[Dict[str, Any]] = field(default_factory=list)
    elevated_plus_maze: List[Dict[str, Any]] = field(default_factory=list)
    novel_object_test: List[Dict[str, Any]] = field(default_factory=list)
    morris_water_maze: List[Dict[str, Any]] = field(default_factory=list)
    
    # Behavioral metrics derived from tests
    behavioral_metrics: Dict[str, float] = field(default_factory=dict)
    
    # Metadata
    test_date: Optional[datetime] = None
    testing_conditions: Dict[str, str] = field(default_factory=dict)
    
    def add_behavioral_test_result(self, test_type: str, data: List[Dict[str, Any]]):
        """Add results for a behavioral test."""
        if test_type == "open_field":
            self.open_field_test.extend(data)
        elif test_type == "elevated_plus_maze":
            self.elevated_plus_maze.extend(data)
        elif test_type == "novel_object":
            self.novel_object_test.extend(data)
        elif test_type == "morris_water_maze":
            self.morris_water_maze.extend(data)
    
    def calculate_behavioral_metrics(self) -> Dict[str, float]:
        """Calculate various behavioral metrics from test data."""
        metrics = {}
        
        # Open field test metrics
        if self.open_field_test:
            center_time = sum(1 for d in self.open_field_test if d.get('zone') == 'center')
            total_time = len(self.open_field_test)
            if total_time > 0:
                metrics['center_time_ratio'] = center_time / total_time
                # Calculate total distance traveled
                total_distance = 0.0
                if len(self.open_field_test) > 1:
                    for i in range(1, len(self.open_field_test)):
                        dx = self.open_field_test[i]['x'] - self.open_field_test[i-1]['x']
                        dy = self.open_field_test[i]['y'] - self.open_field_test[i-1]['y']
                        total_distance += (dx**2 + dy**2)**0.5
                metrics['total_distance'] = total_distance
        
        # Elevated plus maze metrics
        if any('summary' in d for d in self.elevated_plus_maze):
            summary_record = next(d for d in self.elevated_plus_maze if d.get('summary', False))
            metrics['time_in_open'] = summary_record.get('time_in_open', 0)
            metrics['entries_to_open'] = summary_record.get('entries_to_open', 0)
            metrics['open_time_ratio'] = summary_record.get('open_time_ratio', 0)
        
        # Novel object test metrics
        if self.novel_object_test:
            final_record = self.novel_object_test[-1] if self.novel_object_test else {}
            if 'object_0_interactions' in final_record:
                obj0_interactions = final_record.get('object_0_interactions', 0)
                obj1_interactions = final_record.get('object_1_interactions', 0)
                total_interactions = obj0_interactions + obj1_interactions
                metrics['object_0_interactions'] = obj0_interactions
                metrics['object_1_interactions'] = obj1_interactions
                metrics['total_interactions'] = total_interactions
                if total_interactions > 0:
                    metrics['discrimination_ratio'] = (obj1_interactions - obj0_interactions) / total_interactions
        
        self.behavioral_metrics = metrics
        return metrics
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format for serialization."""
        return {
            "open_field_test": self.open_field_test,
            "elevated_plus_maze": self.elevated_plus_maze,
            "novel_object_test": self.novel_object_test,
            "morris_water_maze": self.morris_water_maze,
            "behavioral_metrics": self.behavioral_metrics,
            "test_date": self.test_date.isoformat() if self.test_date else None,
            "testing_conditions": self.testing_conditions
        }


@dataclass
class BiomarkerData:
    """
    Represents biomarker data for a mouse subject.
    Includes blood work, protein levels, metabolites, etc.
    """
    # Blood chemistry
    blood_work: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    
    # Protein levels
    protein_levels: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    
    # Metabolites
    metabolites: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    
    # Inflammatory markers
    inflammatory_markers: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    
    # Hormones
    hormone_levels: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    
    # Metadata
    collection_date: Optional[datetime] = None
    collection_method: str = ""
    
    def add_blood_work(self, test_name: str, value: float, units: str = "", 
                      reference_range: Optional[List[float]] = None):
        """Add a blood work result."""
        self.blood_work[test_name] = {
            "value": value,
            "units": units,
            "reference_range": reference_range,
            "timestamp": datetime.now().isoformat()
        }
    
    def add_protein_level(self, protein_name: str, value: float, 
                         units: str = "", method: str = ""):
        """Add a protein level measurement."""
        self.protein_levels[protein_name] = {
            "value": value,
            "units": units,
            "method": method,
            "timestamp": datetime.now().isoformat()
        }
    
    def add_metabolite(self, metabolite_name: str, value: float, 
                      units: str = "", method: str = ""):
        """Add a metabolite measurement."""
        self.metabolites[metabolite_name] = {
            "value": value,
            "units": units,
            "method": method,
            "timestamp": datetime.now().isoformat()
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format for serialization."""
        return {
            "blood_work": self.blood_work,
            "protein_levels": self.protein_levels,
            "metabolites": self.metabolites,
            "inflammatory_markers": self.inflammatory_markers,
            "hormone_levels": self.hormone_levels,
            "collection_date": self.collection_date.isoformat() if self.collection_date else None,
            "collection_method": self.collection_method
        }


@dataclass
class Treatment:
    """
    Represents a treatment or experimental intervention.
    """
    treatment_type: str
    compound: str = ""
    dose: float = 0.0
    dose_unit: str = ""
    duration: int = 0  # in days
    start_date: Optional[datetime] = None
    end_date: Optional[datetime] = None
    route: str = "oral"  # oral, i.p., s.c., etc.
    frequency: str = "once"  # once, daily, twice_daily, etc.
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format for serialization."""
        return {
            "treatment_type": self.treatment_type,
            "compound": self.compound,
            "dose": self.dose,
            "dose_unit": self.dose_unit,
            "duration": self.duration,
            "start_date": self.start_date.isoformat() if self.start_date else None,
            "end_date": self.end_date.isoformat() if self.end_date else None,
            "route": self.route,
            "frequency": self.frequency
        }


@dataclass
class MouseSubject:
    """
    The main unified data model that connects all mouse research data.
    This class combines genomic, phenotypic, behavioral, and biomarker data
    with metadata to provide a complete view of the subject.
    """
    # Unique identifiers
    mouse_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    cohort_id: Optional[str] = None
    
    # Basic subject information
    strain: MouseStrain = MouseStrain.C57BL_6J
    sex: Sex = Sex.UNKNOWN
    age_weeks: int = 0
    generation: str = "F0"  # F0, F1, etc.
    
    # Primary data components
    genomic_data: GenomicData = field(default_factory=GenomicData)
    phenotypic_data: PhenotypicData = field(default_factory=PhenotypicData)
    behavioral_data: BehavioralData = field(default_factory=BehavioralData)
    biomarker_data: BiomarkerData = field(default_factory=BiomarkerData)
    
    # Experimental information
    treatments: List[Treatment] = field(default_factory=list)
    experimental_conditions: Dict[str, str] = field(default_factory=dict)
    
    # Metadata
    created_date: datetime = field(default_factory=datetime.now)
    last_modified: datetime = field(default_factory=datetime.now)
    notes: List[str] = field(default_factory=list)
    
    def add_treatment(self, treatment: Treatment):
        """Add a treatment to the mouse subject."""
        self.treatments.append(treatment)
        self.last_modified = datetime.now()
    
    def add_note(self, note: str):
        """Add a note about the mouse subject."""
        self.notes.append(f"[{datetime.now().isoformat()}] {note}")
        self.last_modified = datetime.now()
    
    def link_data_types(self):
        """
        Establish connections between different data types based on 
        biological relationships and temporal associations.
        """
        # Example: Link anxiety phenotype to behavioral tests
        if self.phenotypic_data.anxiety_score is not None:
            # This score might be derived from elevated plus maze performance
            # or open field test center time
            if self.behavioral_data.behavioral_metrics:
                anxiety_corr = self.phenotypic_data.anxiety_score
                center_time = self.behavioral_data.behavioral_metrics.get('center_time_ratio', 0)
                self.add_note(f"Anxiety score ({anxiety_corr}) correlates with center time ({center_time})")
        
        # Example: Link obesity to biomarkers
        if self.phenotypic_data.obesity_score is not None:
            obesity = self.phenotypic_data.obesity_score
            if 'glucose' in self.biomarker_data.blood_work:
                glucose_val = self.biomarker_data.blood_work['glucose']['value']
                self.add_note(f"Obesity score ({obesity}) correlates with glucose level ({glucose_val})")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the entire mouse subject to a dictionary format."""
        return {
            "mouse_id": self.mouse_id,
            "cohort_id": self.cohort_id,
            "strain": self.strain.value,
            "sex": self.sex.value,
            "age_weeks": self.age_weeks,
            "generation": self.generation,
            "genomic_data": self.genomic_data.to_dict(),
            "phenotypic_data": self.phenotypic_data.to_dict(),
            "behavioral_data": self.behavioral_data.to_dict(),
            "biomarker_data": self.biomarker_data.to_dict(),
            "treatments": [t.to_dict() for t in self.treatments],
            "experimental_conditions": self.experimental_conditions,
            "created_date": self.created_date.isoformat(),
            "last_modified": self.last_modified.isoformat(),
            "notes": self.notes
        }
    
    def to_json(self) -> str:
        """Convert the mouse subject to JSON string."""
        return json.dumps(self.to_dict(), indent=2, default=str)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'MouseSubject':
        """Create a MouseSubject from a dictionary."""
        # Create the basic object
        mouse = cls(
            mouse_id=data.get('mouse_id', str(uuid.uuid4())),
            cohort_id=data.get('cohort_id'),
            strain=MouseStrain(data.get('strain', 'C57BL/6J')),
            sex=Sex(data.get('sex', 'unknown')),
            age_weeks=data.get('age_weeks', 0),
            generation=data.get('generation', 'F0'),
            created_date=datetime.fromisoformat(data['created_date']) if data.get('created_date') else datetime.now(),
            last_modified=datetime.fromisoformat(data['last_modified']) if data.get('last_modified') else datetime.now(),
            experimental_conditions=data.get('experimental_conditions', {}),
            notes=data.get('notes', [])
        )
        
        # Populate genomic data
        genomic_data = data.get('genomic_data', {})
        mouse.genomic_data = GenomicData()
        mouse.genomic_data.variants = genomic_data.get('variants', [])
        mouse.genomic_data.specific_mutations = genomic_data.get('specific_mutations', [])
        mouse.genomic_data.coverage_metrics = genomic_data.get('coverage_metrics', {})
        mouse.genomic_data.assembly_quality = genomic_data.get('assembly_quality', '')
        mouse.genomic_data.sequencing_technology = genomic_data.get('sequencing_technology', '')
        if genomic_data.get('sequencing_date'):
            mouse.genomic_data.sequencing_date = datetime.fromisoformat(genomic_data['sequencing_date'])
        
        # Populate phenotypic data
        phenotypic_data = data.get('phenotypic_data', {})
        mouse.phenotypic_data = PhenotypicData(
            anxiety_score=phenotypic_data.get('anxiety_score'),
            memory_score=phenotypic_data.get('memory_score'),
            aggression_score=phenotypic_data.get('aggression_score'),
            depression_score=phenotypic_data.get('depression_score'),
            locomotor_activity=phenotypic_data.get('locomotor_activity'),
            obesity_score=phenotypic_data.get('obesity_score'),
            tumor_risk=phenotypic_data.get('tumor_risk'),
            weight=phenotypic_data.get('weight'),
            length=phenotypic_data.get('length'),
            coat_color=phenotypic_data.get('coat_color'),
            phenotype_timeline=phenotypic_data.get('phenotype_timeline', []),
            measurement_method=phenotypic_data.get('measurement_method', '')
        )
        if phenotypic_data.get('measurement_date'):
            mouse.phenotypic_data.measurement_date = datetime.fromisoformat(phenotypic_data['measurement_date'])
        
        # Populate behavioral data
        behavioral_data = data.get('behavioral_data', {})
        mouse.behavioral_data = BehavioralData(
            open_field_test=behavioral_data.get('open_field_test', []),
            elevated_plus_maze=behavioral_data.get('elevated_plus_maze', []),
            novel_object_test=behavioral_data.get('novel_object_test', []),
            morris_water_maze=behavioral_data.get('morris_water_maze', []),
            behavioral_metrics=behavioral_data.get('behavioral_metrics', {}),
            testing_conditions=behavioral_data.get('testing_conditions', {})
        )
        if behavioral_data.get('test_date'):
            mouse.behavioral_data.test_date = datetime.fromisoformat(behavioral_data['test_date'])
        
        # Populate biomarker data
        biomarker_data = data.get('biomarker_data', {})
        mouse.biomarker_data = BiomarkerData(
            blood_work=biomarker_data.get('blood_work', {}),
            protein_levels=biomarker_data.get('protein_levels', {}),
            metabolites=biomarker_data.get('metabolites', {}),
            inflammatory_markers=biomarker_data.get('inflammatory_markers', {}),
            hormone_levels=biomarker_data.get('hormone_levels', {}),
            collection_method=biomarker_data.get('collection_method', '')
        )
        if biomarker_data.get('collection_date'):
            mouse.biomarker_data.collection_date = datetime.fromisoformat(biomarker_data['collection_date'])
        
        # Populate treatments
        for treatment_data in data.get('treatments', []):
            treatment = Treatment(
                treatment_type=treatment_data['treatment_type'],
                compound=treatment_data.get('compound', ''),
                dose=treatment_data.get('dose', 0.0),
                dose_unit=treatment_data.get('dose_unit', ''),
                duration=treatment_data.get('duration', 0),
                route=treatment_data.get('route', 'oral'),
                frequency=treatment_data.get('frequency', 'once')
            )
            if treatment_data.get('start_date'):
                treatment.start_date = datetime.fromisoformat(treatment_data['start_date'])
            if treatment_data.get('end_date'):
                treatment.end_date = datetime.fromisoformat(treatment_data['end_date'])
            mouse.treatments.append(treatment)
        
        return mouse


# Example usage and testing
if __name__ == "__main__":
    # Create a sample mouse subject with unified data
    mouse = MouseSubject(
        mouse_id="mouse_001",
        strain=MouseStrain.C57BL_6J,
        sex=Sex.MALE,
        age_weeks=12
    )
    
    # Add genomic data
    mouse.genomic_data.add_variant(
        position=1000,
        reference_allele="A",
        alternative_allele="G",
        variant_type=GenomicVariantType.SNP,
        quality_score=0.95
    )
    
    # Add phenotypic data
    mouse.phenotypic_data.anxiety_score = 0.7
    mouse.phenotypic_data.obesity_score = 0.4
    mouse.phenotypic_data.memory_score = 0.6
    
    # Add behavioral test results
    # Simulate some open field test data
    for i in range(10):
        mouse.behavioral_data.open_field_test.append({
            'x': 50.0 + i,
            'y': 50.0 + i,
            'time': i,
            'zone': 'center' if i > 3 else 'periphery',
            'velocity_x': 1.0,
            'velocity_y': 1.0,
            'speed': 1.41
        })
    
    # Add biomarker data
    mouse.biomarker_data.add_blood_work("glucose", 150.0, "mg/dL")
    mouse.biomarker_data.add_blood_work("cholesterol", 200.0, "mg/dL")
    mouse.biomarker_data.add_protein_level("BDNF", 1.2, "ng/mL")
    
    # Add treatment
    treatment = Treatment(
        treatment_type="high_fat_diet",
        compound="high fat diet",
        duration=30,
        start_date=datetime.now()
    )
    mouse.add_treatment(treatment)
    
    # Link the data types together
    mouse.link_data_types()
    
    # Print the unified data
    print("Unified Mouse Data Model Example:")
    print(mouse.to_json())
    
    # Test serialization/deserialization
    json_str = mouse.to_json()
    reconstructed_mouse = MouseSubject.from_dict(json.loads(json_str))
    
    print("\nReconstructed mouse anxiety score:", reconstructed_mouse.phenotypic_data.anxiety_score)
    print("Number of treatments:", len(reconstructed_mouse.treatments))
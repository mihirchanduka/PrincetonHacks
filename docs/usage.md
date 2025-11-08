# Usage Guide

## Installation and Setup

1. Install Python 3.8+
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Running the System

### Option 1: Full Simulation (Recommended)
```bash
python3 virtual_mouse_lab.py
```

### Option 2: API Server
```bash
uvicorn main_api:app --reload --port 8000
```

Then send API requests:
```bash
curl -X POST "http://localhost:8000/api/simulate" \
     -H "Content-Type: application/json" \
     -d '{
       "genotype": {"base": "C57BL/6J", "knockout": "Trp53"},
       "treatment": {"compound": "NewDrugX", "dose_mg_kg": 20, "days": 30},
       "tests": ["blood_work", "open_field"]
     }'
```

### Option 3: Legacy MVP Workflow
```bash
./run_mvp.sh
```
This runs the four-step workflow: generate genome → predict phenotype → run behavioral test → simulate experiment.

## Main Components Usage

### Virtual Mouse Lab
```python
from virtual_mouse_lab import VirtualMouseLab

lab = VirtualMouseLab()
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
    experiments=[{
        "type": "high_fat_diet",
        "params": {"obesity_increase": 0.15}
    }],
    include_behavior=True,
    use_real_genome=True
)
```

### Data Model Usage
```python
from data.unified_data_model import MouseSubject, Treatment

mouse = MouseSubject(
    mouse_id="mouse_001",
    strain=MouseStrain.C57BL_6J,
    sex=Sex.MALE,
    age_weeks=12
)

# Add phenotypic data
mouse.phenotypic_data.anxiety_score = 0.6
mouse.phenotypic_data.memory_score = 0.7

# Add behavioral data
mouse.behavioral_data.open_field_test.append({
    'x': 50.0,
    'y': 50.0,
    'time': 0,
    'zone': 'center',
    'velocity_x': 1.0,
    'velocity_y': 1.0,
    'speed': 1.41
})

# Add biomarker data
mouse.biomarker_data.add_blood_work("glucose", 150.0, "mg/dL")

# Add treatment
treatment = Treatment(
    treatment_type="drug_treatment",
    compound="TestCompound",
    dose=10.0,
    dose_unit="mg/kg",
    duration=14
)
mouse.add_treatment(treatment)
```
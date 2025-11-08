# Output Structure and Format

## Package Structure

The system generates comprehensive ZIP packages with this structure:
```
synthetic_mouse_single_{timestamp}.zip
├── genome/
│   └── mouse_genome.fasta
├── phenotype/
│   └── mouse_phenotype_scores.csv
├── blood_work/
│   └── mouse_biomarkers.json
├── behavior/
│   ├── mouse_open_field_test.csv
│   ├── mouse_elevated_plus_maze.csv
│   └── mouse_novel_object_test.csv
├── report/
│   ├── experimental_results.json
│   └── summary.txt
└── unified_data/
    └── mouse_subject.json
```

## File Format Specifications

### Mouse Genome (FASTA)
```
>mouse_genome_sequence
ATCGATCGATCG...
```

### Phenotype Scores (CSV)
```
phenotype_type,score,confidence,description
anxiety,0.6,0.85,Mice show moderate anxiety-like behavior
memory,0.7,0.78,Recognition memory performance
obesity,0.45,0.81,Body weight relative to controls
tumor_risk,0.3,0.92,Predisposition to tumor formation
aggression,0.2,0.75,Aggressive behavior score
```

### Biomarkers (JSON)
```json
{
  "mouse_id": "mouse_001",
  "timestamp": "2025-11-08T14:30:00Z",
  "measurements": [
    {
      "name": "glucose",
      "value": 150.0,
      "unit": "mg/dL",
      "reference_range": {"min": 120.0, "max": 180.0},
      "status": "normal"
    },
    {
      "name": "cholesterol",
      "value": 220.0,
      "unit": "mg/dL",
      "reference_range": {"min": 100.0, "max": 200.0},
      "status": "high"
    }
  ]
}
```

### Behavioral Data (CSV)

#### Open Field Test
```
timestamp,x,y,zone,velocity_x,velocity_y,speed,acceleration,rotation
0,50.0,50.0,center,1.0,0.0,1.0,0.0,0
1,51.0,50.0,center,1.1,0.1,1.1,0.1,5
...
```

#### Elevated Plus Maze
```
timestamp,zone,entry_time,exit_time,duration,velocity,speed
0,open_arm,0,15,15,1.2,1.2
15,closed_arm,15,30,15,0.8,0.8
...
```

#### Novel Object Test
```
timestamp,trial,object_type,duration,frequency,distance,velocity
0,habituation,novel,300,12,45.5,1.2
300,test,familiar,150,8,32.1,0.9
450,test,novel,150,15,52.3,1.1
```

### Experimental Results (JSON)
```json
{
  "experiment_id": "exp_001",
  "mouse_id": "mouse_001",
  "treatment": {
    "type": "high_fat_diet",
    "start_time": "2025-11-08T10:00:00Z",
    "duration_days": 30
  },
  "pre_treatment": {
    "weight": 25.0,
    "anxiety_score": 0.6,
    "memory_score": 0.7
  },
  "post_treatment": {
    "weight": 32.5,
    "anxiety_score": 0.65,
    "memory_score": 0.65
  },
  "changes": {
    "weight_change": 7.5,
    "weight_percent_change": 30.0,
    "anxiety_change": 0.05,
    "memory_change": -0.05
  }
}
```

### Summary Report (TXT)
```
Virtual Mouse Laboratory - Experimental Summary

Mouse ID: mouse_001
Strain: C57BL/6J
Sex: Male
Age: 12 weeks

Phenotype Profile:
- Anxiety Score: 0.6 (Moderate anxiety-like behavior)
- Memory Score: 0.7 (Above average recognition memory)
- Obesity Score: 0.45 (Normal body weight)
- Tumor Risk: 0.3 (Low risk)
- Aggression Score: 0.2 (Low aggressive behavior)

Treatment Applied: High fat diet for 30 days
Weight Change: +30% (25.0g to 32.5g)
Phenotypic Changes: Slight increase in anxiety, minor memory decline
Biomarker Changes: Cholesterol increased to 220 mg/dL (above normal)

Conclusion: Mouse showed expected metabolic responses to high-fat diet.
```

### Unified Mouse Subject (JSON)
```json
{
  "mouse_subject": {
    "id": "mouse_001",
    "strain": "C57BL/6J",
    "sex": "MALE",
    "age_weeks": 12,
    "genotype": {
      "reference_genome": "GRCm39",
      "variants": [
        {
          "chromosome": "chr1",
          "position": 123456,
          "reference": "A",
          "alternative": "T",
          "type": "SNP",
          "effect": "synonymous"
        }
      ]
    },
    "phenotypic_data": {
      "anxiety_score": 0.6,
      "memory_score": 0.7,
      "obesity_score": 0.45,
      "tumor_risk": 0.3,
      "aggression_score": 0.2
    },
    "behavioral_data": {
      "open_field_test": [...],
      "elevated_plus_maze": [...],
      "novel_object_test": [...]
    },
    "biomarker_data": {...},
    "experimental_data": [...]
  }
}
```

## Data Validation and Quality Assurance

All output files include:
- Timestamp information
- Quality metrics where applicable
- Reference ranges for biomarkers
- Confidence scores for predictions
- Validation flags for data integrity
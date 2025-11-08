# API Documentation

## Base URL
```
http://localhost:8000
```

## Endpoints

### 1. Simulate Mouse Experiment
- **Endpoint**: `POST /api/simulate`
- **Description**: Run a complete mouse simulation with specified parameters
- **Request Body**:
```json
{
  "genotype": {
    "base": "C57BL/6J",
    "knockout": "Trp53",
    "transgene": null
  },
  "phenotype_profile": {
    "anxiety_score": 0.6,
    "memory_score": 0.7,
    "obesity_score": 0.4,
    "tumor_risk": 0.3,
    "aggression_score": 0.4
  },
  "treatment": {
    "type": "drug_treatment",
    "compound": "NewDrugX",
    "dose_mg_kg": 20,
    "days": 30
  },
  "tests": ["blood_work", "open_field", "elevated_plus_maze"],
  "include_behavior": true
}
```
- **Response**: ZIP file containing complete simulation results

### 2. Generate Synthetic Mouse
- **Endpoint**: `POST /api/mouse/generate`
- **Description**: Create a synthetic mouse with specified characteristics
- **Request Body**:
```json
{
  "strain": "C57BL/6J",
  "sex": "MALE",
  "age_weeks": 12,
  "phenotype_profile": {
    "anxiety_score": 0.6,
    "memory_score": 0.7
  }
}
```
- **Response**:
```json
{
  "mouse_id": "mouse_001",
  "genome_path": "/genome/mouse_001.fasta",
  "phenotype_path": "/phenotype/mouse_001.csv",
  "status": "created"
}
```

### 3. Run Behavioral Tests
- **Endpoint**: `POST /api/behavior/test`
- **Description**: Run behavioral tests on a synthetic mouse
- **Request Body**:
```json
{
  "mouse_id": "mouse_001",
  "tests": ["open_field", "elevated_plus_maze", "novel_object"],
  "duration_minutes": 30
}
```
- **Response**:
```json
{
  "test_results": {
    "open_field": "/behavior/mouse_001_open_field.csv",
    "elevated_plus_maze": "/behavior/mouse_001_elevated_plus_maze.csv",
    "novel_object": "/behavior/mouse_001_novel_object.csv"
  },
  "analysis": {
    "anxiety_score": 0.62,
    "locomotion_activity": "high",
    "exploration_behavior": "balanced"
  }
}
```

### 4. Apply Treatment
- **Endpoint**: `POST /api/treatment/apply`
- **Description**: Apply a treatment to a synthetic mouse and observe changes
- **Request Body**:
```json
{
  "mouse_id": "mouse_001",
  "treatment": {
    "type": "high_fat_diet",
    "params": {
      "duration_days": 30,
      "fat_content_percent": 60
    }
  }
}
```
- **Response**:
```json
{
  "treatment_id": "trt_001",
  "pre_treatment_phenotype": {
    "weight": 25.0,
    "obesity_score": 0.4
  },
  "post_treatment_phenotype": {
    "weight": 32.0,
    "obesity_score": 0.7
  },
  "biomarker_changes": {
    "glucose": {"before": 150, "after": 180},
    "cholesterol": {"before": 200, "after": 280}
  }
}
```

### 5. Get Phenotype Prediction
- **Endpoint**: `POST /api/phenotype/predict`
- **Description**: Predict phenotypes from genomic data
- **Request Body**:
```json
{
  "genome_path": "/genome/mouse_001.fasta",
  "phenotypes": ["anxiety", "memory", "obesity", "tumor_risk"]
}
```
- **Response**:
```json
{
  "predictions": {
    "anxiety_score": 0.6,
    "memory_score": 0.7,
    "obesity_score": 0.4,
    "tumor_risk": 0.3,
    "confidence_intervals": {
      "anxiety": [0.55, 0.65],
      "memory": [0.65, 0.75]
    }
  }
}
```

### 6. Run Complete Experiment
- **Endpoint**: `POST /api/experiment/run`
- **Description**: Run a complete experiment with multiple mice and treatments
- **Request Body**:
```json
{
  "experiment_name": "DrugX_Efficacy_Study",
  "mice": [
    {
      "mouse_id": "mouse_001",
      "treatment_group": "control"
    },
    {
      "mouse_id": "mouse_002", 
      "treatment_group": "treatment",
      "treatment_params": {
        "compound": "DrugX",
        "dose_mg_kg": 20
      }
    }
  ],
  "timeline": {
    "baseline_measurements": 0,
    "treatment_start": 1,
    "midpoint_measurements": 15,
    "endpoint_measurements": 30
  },
  "endpoints": ["behavior", "biomarkers", "phenotype"]
}
```
- **Response**:
```json
{
  "experiment_id": "exp_001",
  "status": "completed",
  "results_path": "/report/exp_001_results.zip",
  "summary": {
    "mice_count": 2,
    "treatment_effect": "significant",
    "p_value": 0.023,
    "effect_size": 0.75
  }
}
```

## Error Responses

All endpoints return standardized error responses:
```json
{
  "error": "error_type",
  "message": "Human-readable error message",
  "details": "Additional error details if available"
}
```

## Status Codes

- `200`: Success
- `201`: Created
- `400`: Bad request (invalid parameters)
- `404`: Resource not found
- `500`: Internal server error

## Authentication

For production deployments, endpoints may require authentication headers:
```
Authorization: Bearer <token>
```

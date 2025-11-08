# Architecture and How the Code Works

## Core Architecture

The system is built around a **unified data model** that connects genomic, phenotypic, behavioral, and biomarker data into a coherent framework. The main components include:

1. **VirtualMouseLab Class** (`virtual_mouse_lab.py`): The main orchestrator that coordinates the entire simulation process
2. **Unified Data Model** (`data/unified_data_model.py`): Connects all data types using class structures
3. **Genomic Processing** (`model/twin_mouse_generator.py`): Handles genome creation and manipulation
4. **Phenotype Prediction** (`model/phenotype_twin_generator.py`): Predicts phenotypes from genomic data
5. **Behavioral Simulation** (`model/behavioral_simulator.py`): Simulates realistic mouse behaviors
6. **API Layer** (`main_api.py`): Provides REST endpoints for the system

## Key Functionalities

### 1. Twin Generation
- Creates synthetic mice that mimic real mouse phenotypes
- Uses ML models trained on real mouse phenome data
- Generates highly similar genomic sequences with controlled variations

### 2. Phenotype Prediction
- Rule-based and ML-based approaches
- Predicts multiple phenotypes (anxiety, memory, obesity, tumor risk, aggression)
- Correlates genomic features with phenotypic outcomes

### 3. Behavioral Simulation
- Implements three major behavioral tests:
  - Open Field Test (anxiety, locomotion)
  - Elevated Plus Maze (anxiety-like behavior)
  - Novel Object Test (recognition memory)
- Simulates realistic movement patterns based on phenotype

### 4. Experimental Simulation
- Simulates various treatments (diet, drugs, environmental factors)
- Tracks phenotype changes over time
- Generates biomarker data (blood work, protein levels)

## Data Flow

1. **Genomic Input**: Starts with reference genome or creates synthetic genome
2. **Phenotype Prediction**: ML models predict phenotypes from genomic data
3. **Behavioral Simulation**: Phenotypes drive realistic behavioral patterns
4. **Experimental Treatment**: Simulates experiments and updates phenotypes
5. **Biomarker Generation**: Creates realistic biomarker data based on phenotypes
6. **Data Packaging**: Packages all data in structured format

## Architecture and Design Patterns

### Unified Data Model
The core innovation is the unified data model that connects different data types:
- **GenomicData**: Stores sequence data and variants
- **PhenotypicData**: Stores behavioral and physical phenotypes
- **BehavioralData**: Stores behavioral test results and metrics
- **BiomarkerData**: Stores blood work, protein levels, metabolites
- **MouseSubject**: Integrates all data types with metadata

### Modular Design
The system follows a modular architecture:
- **Model Layer**: Contains business logic and algorithms
- **Data Layer**: Handles data structures and models
- **API Layer**: Provides REST endpoints
- **Utility Layer**: Contains helper functions and scripts
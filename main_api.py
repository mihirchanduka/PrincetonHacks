"""
Virtual Mouse Lab - FastAPI Server
Implements the API for the virtual mouse lab as per the MVP flow.
"""
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import Dict, List, Optional
import uuid
import os
from datetime import datetime

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from virtual_mouse_lab import VirtualMouseLab
from fastapi.middleware.cors import CORSMiddleware


app = FastAPI(
    title="Virtual Mouse Lab API",
    description="An in-silico platform for creating synthetic mice, running experiments, and generating research data",
    version="1.0.0"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://127.0.0.1:3000",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class MouseGenotype(BaseModel):
    base: str = "C57BL/6J"
    knockout: Optional[str] = None
    custom_mutations: Optional[List[Dict]] = None


class Treatment(BaseModel):
    compound: str = "control"
    dose_mg_kg: Optional[float] = 10.0
    days: Optional[int] = 30


class ExperimentRequest(BaseModel):
    genotype: MouseGenotype
    treatment: Optional[Treatment] = None
    experiments: Optional[List[Dict]] = []
    tests: List[str] = ["blood_work", "open_field"]
    twin_params: Optional[Dict] = {
        "mutation_rate": 0.005,
        "twin_similarity": 0.90,
        "genome_length": 1000
    }


class SimulationResponse(BaseModel):
    simulation_id: str
    download_url: str
    summary: Dict
    created_at: str


# Store for simulation results (in production, use a database)
simulations_store = {}


@app.post("/api/simulate", response_model=SimulationResponse)
async def run_simulation(request: ExperimentRequest):
    """
    Run a complete mouse simulation with the specified parameters.
    
    This endpoint creates synthetic mice, runs experiments, and returns
    a downloadable package with research-ready data.
    """
    try:
        # Generate a unique simulation ID
        simulation_id = str(uuid.uuid4())
        
        # Initialize the virtual mouse lab
        lab = VirtualMouseLab()
        
        # Prepare experiments from request
        experiments = []
        if request.treatment:
            experiments.append({
                "type": "drug_treatment",
                "params": {
                    "drug_name": request.treatment.compound,
                    "dose_mg_kg": request.treatment.dose_mg_kg
                }
            })
        
        # Add any additional experiments from the request
        if request.experiments:
            experiments.extend(request.experiments)
        
        # Run the single mouse simulation using real genome data
        include_behavior = "open_field" in request.tests
        
        # Use the single mouse simulation method with real genome data
        package_path = lab.run_single_mouse_simulation(
            target_mouse_profile={
                'anxiety_score': 0.5,
                'memory_score': 0.5,
                'obesity_score': 0.5,
                'tumor_risk': 0.3,
                'aggression_score': 0.4,
                'glucose': 150.0,
                'cholesterol': 200.0,
                'triglycerides': 150.0,
                'ALT': 35.0
            },
            experiments=experiments,
            include_behavior=include_behavior,
            use_real_genome=True  # Use real genome from dataset for real world predictions
        )
        
        # Create download URL (in a real implementation, this would be a temporary URL)
        download_url = f"/api/download/{simulation_id}"
        
        # Store simulation result (in a real implementation, use a database)
        simulation_result = {
            "package_path": package_path,
            "created_at": datetime.now().isoformat(),
            "genotype": request.genotype,
            "treatment": request.treatment,
            "tests": request.tests
        }
        simulations_store[simulation_id] = simulation_result
        
        # For now, return basic summary
        response = SimulationResponse(
            simulation_id=simulation_id,
            download_url=download_url,
            summary={
                "status": "completed",
                "phenotype_mimicry": True,
                "experiments_run": len(experiments),
                "tests_conducted": request.tests
            },
            created_at=simulation_result["created_at"]
        )
        
        return response
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Simulation failed: {str(e)}")


@app.get("/api/download/{simulation_id}")
async def download_results(simulation_id: str):
    """
    Download the results package for a completed simulation.
    """
    if simulation_id not in simulations_store:
        raise HTTPException(status_code=404, detail="Simulation not found")
    
    simulation = simulations_store[simulation_id]
    package_path = simulation["package_path"]
    
    if not os.path.exists(package_path):
        raise HTTPException(status_code=404, detail="Results package not found")
    
    return FileResponse(
        path=package_path,
        filename=os.path.basename(package_path),
        media_type='application/zip'
    )


@app.get("/")
def read_root():
    """
    Root endpoint providing API information.
    """
    return {
        "service": "Virtual Mouse Lab",
        "version": "1.0.0",
        "description": "An in-silico platform for creating synthetic mice, running experiments, and generating research data",
        "endpoints": {
            "POST /api/simulate": "Run a mouse simulation",
            "GET /api/download/{id}": "Download simulation results"
        }
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
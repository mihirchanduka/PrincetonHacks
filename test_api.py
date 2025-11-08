"""
Test script for the Virtual Mouse Lab API
"""
import requests
import json

def test_api():
    """Test the API endpoint"""
    # Define the payload for the simulation
    payload = {
        "genotype": {"base": "C57BL/6J", "knockout": "Trp53"},
        "treatment": {"compound": "NewDrugX", "dose_mg_kg": 20, "days": 30},
        "tests": ["blood_work", "open_field"],
        "twin_params": {
            "mutation_rate": 0.005,
            "twin_similarity": 0.90,
            "genome_length": 1000
        }
    }
    
    print("Testing the Virtual Mouse Lab API...")
    print(f"Payload: {json.dumps(payload, indent=2)}")
    
    # In a real scenario, you would send this to the running server:
    # response = requests.post("http://localhost:8000/api/simulate", json=payload)
    # print(f"Response: {response.status_code}, {response.json()}")
    
    print("\nAPI endpoint structure ready!")
    print("To test the full API:")
    print("1. Start the server: uvicorn main_api:app --reload --port 8000")
    print("2. Send a POST request to http://localhost:8000/api/simulate with the above payload")
    print("3. Download results using the returned download URL")

if __name__ == "__main__":
    test_api()
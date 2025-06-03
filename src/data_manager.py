import json
import os
import uuid
from datetime import datetime
from typing import Dict, Optional
import logging

logger = logging.getLogger(__name__)

class DataManager:
    """Service for managing analysis data and results"""
    
    def __init__(self):
        self.data_dir = "data"
        self.results_dir = os.path.join(self.data_dir, "results")
        self._ensure_directories()
    
    def _ensure_directories(self):
        """Ensure data directories exist"""
        os.makedirs(self.results_dir, exist_ok=True)
    
    def save_analysis_result(self, result_data: Dict) -> str:
        """Save analysis result and return unique ID"""
        try:
            result_id = str(uuid.uuid4())
            
            # Add metadata
            result_data["id"] = result_id
            result_data["created_at"] = datetime.now().isoformat()
            result_data["version"] = "1.0"
            
            # Save to file
            file_path = os.path.join(self.results_dir, f"{result_id}.json")
            with open(file_path, 'w') as f:
                json.dump(result_data, f, indent=2)
            
            logger.info(f"Saved analysis result with ID: {result_id}")
            return result_id
            
        except Exception as e:
            logger.error(f"Error saving analysis result: {str(e)}")
            raise
    
    def get_analysis_result(self, result_id: str) -> Optional[Dict]:
        """Get analysis result by ID"""
        try:
            file_path = os.path.join(self.results_dir, f"{result_id}.json")
            
            if not os.path.exists(file_path):
                return None
            
            with open(file_path, 'r') as f:
                result_data = json.load(f)
            
            return result_data
            
        except Exception as e:
            logger.error(f"Error getting analysis result: {str(e)}")
            return None
    
    def list_analysis_results(self) -> list:
        """List all saved analysis results"""
        try:
            results = []
            
            for filename in os.listdir(self.results_dir):
                if filename.endswith('.json'):
                    result_id = filename[:-5]  # Remove .json extension
                    file_path = os.path.join(self.results_dir, filename)
                    
                    try:
                        with open(file_path, 'r') as f:
                            data = json.load(f)
                        
                        results.append({
                            "id": result_id,
                            "created_at": data.get("created_at"),
                            "species_count": data.get("species_count", 0),
                            "sequences_processed": data.get("sequences_processed", 0)
                        })
                    except Exception as e:
                        logger.warning(f"Error reading result file {filename}: {str(e)}")
                        continue
            
            # Sort by creation date (newest first)
            results.sort(key=lambda x: x.get("created_at", ""), reverse=True)
            
            return results
            
        except Exception as e:
            logger.error(f"Error listing analysis results: {str(e)}")
            return []
    
    def delete_analysis_result(self, result_id: str) -> bool:
        """Delete analysis result by ID"""
        try:
            file_path = os.path.join(self.results_dir, f"{result_id}.json")
            
            if os.path.exists(file_path):
                os.remove(file_path)
                logger.info(f"Deleted analysis result: {result_id}")
                return True
            else:
                return False
                
        except Exception as e:
            logger.error(f"Error deleting analysis result: {str(e)}")
            return False

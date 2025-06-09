import json
import os
import glob
from typing import List, Dict, Optional
import logging

logger = logging.getLogger(__name__)

class IUCNService:
    """Simplified service for reading IUCN data from multiple JSON files"""
    
    def __init__(self, data_folder_path: str = "iucn/red_list"):
        self.data_folder_path = data_folder_path
        self.species_data = {}
        self.species_index = {}
        self._load_all_iucn_data()
        
    def _load_all_iucn_data(self) -> None:
        """Load all IUCN data from JSON files in the folder"""
        try:
            if not os.path.exists(self.data_folder_path):
                logger.warning(f"IUCN data folder not found: {self.data_folder_path}")
                return
            
            json_files = glob.glob(os.path.join(self.data_folder_path, "*.json"))
            
            if not json_files:
                logger.warning(f"No JSON files found in: {self.data_folder_path}")
                return
            
            total_species = 0
            
            for json_file in json_files:
                logger.info(f"Loading IUCN data from: {json_file}")
                
                try:
                    with open(json_file, 'r', encoding='utf-8') as f:
                        data = json.load(f)
                    
                    if isinstance(data, list):
                        file_species_count = self._process_species_array(data)
                        total_species += file_species_count
                        logger.info(f"Loaded {file_species_count} species from {os.path.basename(json_file)}")
                    else:
                        logger.warning(f"Expected array in {json_file}, got {type(data)}")
                        
                except json.JSONDecodeError as e:
                    logger.error(f"Invalid JSON in {json_file}: {str(e)}")
                except Exception as e:
                    logger.error(f"Error loading {json_file}: {str(e)}")
            
            self._build_search_index()
            logger.info(f"Successfully loaded {total_species} species from {len(json_files)} files")
            
        except Exception as e:
            logger.error(f"Error loading IUCN data: {str(e)}")
    
    def _process_species_array(self, species_array: List[Dict]) -> int:
        """Process array of species data - only extract existing attributes"""
        processed_count = 0
        
        for species_entry in species_array:
            try:
                scientific_name = species_entry.get("taxon_scientific_name")
                
                # Skip if no scientific name
                if not scientific_name:
                    continue
                
                # Only extract attributes that exist in the JSON
                species_data = {}
                
                # Copy all existing attributes from the JSON
                for key, value in species_entry.items():
                    species_data[key] = value
                
                # Store with scientific name as key
                self.species_data[scientific_name] = species_data
                processed_count += 1
                
            except Exception as e:
                logger.warning(f"Error processing species entry: {str(e)}")
                continue
        
        return processed_count
    
    def _build_search_index(self):
        """Build search index for faster searching"""
        self.species_index = {
            'by_category': {},
            'by_extinction_risk': {},
            'by_name_parts': {}
        }
        
        for scientific_name, species_data in self.species_data.items():
            # Index by red list category
            category = species_data.get('red_list_category_code') or species_data.get('code')
            if category:
                if category not in self.species_index['by_category']:
                    self.species_index['by_category'][category] = []
                self.species_index['by_category'][category].append(scientific_name)
            
            # Index by extinction risk
            extinction_risk = 'high_risk' if (
                species_data.get('possibly_extinct', False) or 
                species_data.get('possibly_extinct_in_the_wild', False) or
                category in ['CR', 'EN', 'EX', 'EW']
            ) else 'lower_risk'
            
            if extinction_risk not in self.species_index['by_extinction_risk']:
                self.species_index['by_extinction_risk'][extinction_risk] = []
            self.species_index['by_extinction_risk'][extinction_risk].append(scientific_name)
            
            # Index by name parts for partial matching
            if scientific_name:
                name_parts = scientific_name.lower().split()
                for part in name_parts:
                    if part not in self.species_index['by_name_parts']:
                        self.species_index['by_name_parts'][part] = []
                    self.species_index['by_name_parts'][part].append(scientific_name)
    
    def search_species(self, query: str, limit: int = 30) -> List[Dict]:
        """Search for species by name or partial match"""
        try:
            query_lower = query.lower().strip()
            results = []
            seen_species = set()
            
            # 1. Exact scientific name match
            for scientific_name in self.species_data.keys():
                if scientific_name.lower() == query_lower:
                    if scientific_name not in seen_species:
                        results.append(self.species_data[scientific_name])
                        seen_species.add(scientific_name)
            
            # 2. Starts with query
            for scientific_name in self.species_data.keys():
                if scientific_name.lower().startswith(query_lower):
                    if scientific_name not in seen_species:
                        results.append(self.species_data[scientific_name])
                        seen_species.add(scientific_name)
            
            # 3. Contains query
            for scientific_name in self.species_data.keys():
                if query_lower in scientific_name.lower():
                    if scientific_name not in seen_species:
                        results.append(self.species_data[scientific_name])
                        seen_species.add(scientific_name)
            
            # 4. Search by name parts using index
            if query_lower in self.species_index['by_name_parts']:
                for scientific_name in self.species_index['by_name_parts'][query_lower]:
                    if scientific_name not in seen_species:
                        results.append(self.species_data[scientific_name])
                        seen_species.add(scientific_name)
            
            # Sort by relevance (exact matches first, then alphabetical)
            results.sort(key=lambda x: (
                not x.get('taxon_scientific_name', '').lower().startswith(query_lower),
                x.get('taxon_scientific_name', '')
            ))
            
            return results[:limit]
            
        except Exception as e:
            logger.error(f"Error searching species: {str(e)}")
            return []
    
    def get_species_details(self, species_name: str) -> Optional[Dict]:
        """Get detailed information about a specific species"""
        try:
            # Try exact match first
            if species_name in self.species_data:
                return self.species_data[species_name]
            
            # Try case-insensitive search
            species_name_lower = species_name.lower()
            for key, species_info in self.species_data.items():
                if key.lower() == species_name_lower:
                    return species_info
            
            logger.warning(f"Species not found in IUCN data: {species_name}")
            return None
            
        except Exception as e:
            logger.error(f"Error getting species details: {str(e)}")
            return None
    
    def get_species_by_category(self, category: str) -> List[str]:
        """Get all species by conservation category"""
        return self.species_index['by_category'].get(category, [])
    
    def get_high_risk_species(self) -> List[str]:
        """Get species with high extinction risk"""
        return self.species_index['by_extinction_risk'].get('high_risk', [])
    
    def get_possibly_extinct_species(self) -> List[Dict]:
        """Get species that are possibly extinct"""
        extinct_species = []
        
        for scientific_name, species_data in self.species_data.items():
            if (species_data.get('possibly_extinct', False) or 
                species_data.get('possibly_extinct_in_the_wild', False)):
                extinct_species.append(species_data)
        
        return extinct_species
    
    def get_latest_assessments(self) -> List[Dict]:
        """Get species with latest assessments"""
        latest_species = []
        
        for species_data in self.species_data.values():
            if species_data.get('latest', False):
                latest_species.append(species_data)
        
        return latest_species
    
    def get_statistics(self) -> Dict:
        """Get statistics about loaded IUCN data"""
        stats = {
            'total_species': len(self.species_data),
            'by_category': {},
            'possibly_extinct': 0,
            'possibly_extinct_in_wild': 0,
            'latest_assessments': 0,
            'by_year': {}
        }
        
        for species_data in self.species_data.values():
            # Count by category
            category = species_data.get('red_list_category_code') or species_data.get('code')
            if category:
                stats['by_category'][category] = stats['by_category'].get(category, 0) + 1
            
            # Count extinction risks
            if species_data.get('possibly_extinct', False):
                stats['possibly_extinct'] += 1
            if species_data.get('possibly_extinct_in_the_wild', False):
                stats['possibly_extinct_in_wild'] += 1
            if species_data.get('latest', False):
                stats['latest_assessments'] += 1
            
            # Count by year
            year = species_data.get('year_published')
            if year:
                stats['by_year'][year] = stats['by_year'].get(year, 0) + 1
        
        return stats
    
    def get_all_species(self) -> Dict:
        """Get all species data"""
        return self.species_data
    
    def reload_data(self):
        """Reload all IUCN data from files"""
        logger.info("Reloading IUCN data...")
        self.species_data = {}
        self.species_index = {}
        self._load_all_iucn_data()
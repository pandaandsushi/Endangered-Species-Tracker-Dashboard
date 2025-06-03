import json
import os
import glob
from typing import List, Dict, Optional
import logging

logger = logging.getLogger(__name__)

class IUCNService:
    """Service for reading IUCN data from multiple JSON files in a folder"""
    
    def __init__(self, data_folder_path: str = "iucn/red_list"):
        self.data_folder_path = data_folder_path
        self.species_data = {}
        self.species_index = {}  # For faster searching
        self._load_all_iucn_data()
        
    def _load_all_iucn_data(self) -> None:
        """Load all IUCN data from JSON files in the folder"""
        try:
            if not os.path.exists(self.data_folder_path):
                logger.warning(f"IUCN data folder not found: {self.data_folder_path}")
                return
            
            # Get all JSON files in the folder
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
                    
                    # Process array of species data
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
            
            # Build search index
            self._build_search_index()
            
            logger.info(f"Successfully loaded {total_species} species from {len(json_files)} files")
            
        except Exception as e:
            logger.error(f"Error loading IUCN data: {str(e)}")
    
    def _process_species_array(self, species_array: List[Dict]) -> int:
        """Process array of species data and extract required fields"""
        processed_count = 0
        
        for species_entry in species_array:
            try:
                # Extract required fields
                scientific_name = species_entry.get("taxon_scientific_name")
                red_list_code = species_entry.get("code")  # or "red_list_category_code"
                possibly_extinct = species_entry.get("possibly_extinct", False)
                possibly_extinct_wild = species_entry.get("possibly_extinct_in_the_wild", False)
                
                # Additional useful fields
                year_published = species_entry.get("year_published")
                latest = species_entry.get("latest", False)
                assessment_id = species_entry.get("assessment_id")
                url = species_entry.get("url")
                sis_taxon_id = species_entry.get("sis_taxon_id")
                
                # Skip if no scientific name
                if not scientific_name:
                    continue
                
                # Convert red list code to standard format
                conservation_status = self._convert_red_list_code(red_list_code)
                
                # Create species data structure
                species_data = {
                    'scientific_name': scientific_name,
                    'category': conservation_status,
                    'red_list_code': red_list_code,
                    'possibly_extinct': possibly_extinct,
                    'possibly_extinct_in_the_wild': possibly_extinct_wild,
                    'year_published': year_published,
                    'latest': latest,
                    'assessment_id': assessment_id,
                    'url': url,
                    'sis_taxon_id': sis_taxon_id,
                    
                    # Extract taxonomic info from scientific name
                    'genus': self._extract_genus(scientific_name),
                    'species': self._extract_species(scientific_name),
                    
                    # Default values for missing fields
                    'kingdom': 'ANIMALIA',
                    'phylum': 'UNKNOWN',
                    'class': 'UNKNOWN',
                    'order': 'UNKNOWN',
                    'family': 'UNKNOWN',
                    'population_trend': 'Unknown',
                    'marine_system': False,
                    'freshwater_system': False,
                    'terrestrial_system': True
                }
                
                # Store with scientific name as key
                self.species_data[scientific_name] = species_data
                processed_count += 1
                
            except Exception as e:
                logger.warning(f"Error processing species entry: {str(e)}")
                continue
        
        return processed_count
    
    def _convert_red_list_code(self, code: str) -> str:
        """Convert IUCN red list code to standard conservation status"""
        if not code:
            return 'NE'
        
        # Mapping dari kode IUCN ke status standar
        code_mapping = {
            'EX': 'EX',    # Extinct
            'EW': 'EW',    # Extinct in the Wild
            'CR': 'CR',    # Critically Endangered
            'EN': 'EN',    # Endangered
            'VU': 'VU',    # Vulnerable
            'NT': 'NT',    # Near Threatened
            'LC': 'LC',    # Least Concern
            'DD': 'DD',    # Data Deficient
            'NE': 'NE',    # Not Evaluated
            'I': 'DD',     # Indeterminate (old category) -> Data Deficient
            'K': 'DD',     # Insufficiently Known (old category) -> Data Deficient
            'R': 'NT',     # Rare (old category) -> Near Threatened
            'V': 'VU',     # Vulnerable (old notation)
            'E': 'EN',     # Endangered (old notation)
        }
        
        return code_mapping.get(code.upper(), 'NE')
    
    def _extract_genus(self, scientific_name: str) -> str:
        """Extract genus from scientific name"""
        if not scientific_name:
            return 'Unknown'
        
        parts = scientific_name.split()
        return parts[0] if parts else 'Unknown'
    
    def _extract_species(self, scientific_name: str) -> str:
        """Extract species from scientific name"""
        if not scientific_name:
            return 'Unknown'
        
        parts = scientific_name.split()
        return parts[1] if len(parts) > 1 else 'sp.'
    
    def _build_search_index(self):
        """Build search index for faster searching"""
        self.species_index = {
            'by_genus': {},
            'by_category': {},
            'by_extinction_risk': {},
            'by_name_parts': {}
        }
        
        for scientific_name, species_data in self.species_data.items():
            genus = species_data['genus']
            category = species_data['category']
            
            # Index by genus
            if genus not in self.species_index['by_genus']:
                self.species_index['by_genus'][genus] = []
            self.species_index['by_genus'][genus].append(scientific_name)
            
            # Index by category
            if category not in self.species_index['by_category']:
                self.species_index['by_category'][category] = []
            self.species_index['by_category'][category].append(scientific_name)
            
            # Index by extinction risk
            extinction_risk = 'high_risk' if (
                species_data['possibly_extinct'] or 
                species_data['possibly_extinct_in_the_wild'] or
                category in ['CR', 'EN', 'EX', 'EW']
            ) else 'lower_risk'
            
            if extinction_risk not in self.species_index['by_extinction_risk']:
                self.species_index['by_extinction_risk'][extinction_risk] = []
            self.species_index['by_extinction_risk'][extinction_risk].append(scientific_name)
            
            # Index by name parts for partial matching
            name_parts = scientific_name.lower().split()
            for part in name_parts:
                if part not in self.species_index['by_name_parts']:
                    self.species_index['by_name_parts'][part] = []
                self.species_index['by_name_parts'][part].append(scientific_name)
    
    def search_species(self, query: str, limit: int = 20) -> List[Dict]:
        """Search for species by name, genus, or partial match"""
        try:
            query_lower = query.lower().strip()
            results = []
            seen_species = set()
            
            # 1. Exact scientific name match
            for scientific_name in self.species_data.keys():
                if scientific_name.lower() == query_lower:
                    if scientific_name not in seen_species:
                        results.append(self._format_search_result(scientific_name))
                        seen_species.add(scientific_name)
            
            # 2. Starts with query
            for scientific_name in self.species_data.keys():
                if scientific_name.lower().startswith(query_lower):
                    if scientific_name not in seen_species:
                        results.append(self._format_search_result(scientific_name))
                        seen_species.add(scientific_name)
            
            # 3. Contains query
            for scientific_name in self.species_data.keys():
                if query_lower in scientific_name.lower():
                    if scientific_name not in seen_species:
                        results.append(self._format_search_result(scientific_name))
                        seen_species.add(scientific_name)
            
            # 4. Search by genus using index
            if query_lower.title() in self.species_index['by_genus']:
                for scientific_name in self.species_index['by_genus'][query_lower.title()]:
                    if scientific_name not in seen_species:
                        results.append(self._format_search_result(scientific_name))
                        seen_species.add(scientific_name)
            
            # 5. Search by name parts using index
            if query_lower in self.species_index['by_name_parts']:
                for scientific_name in self.species_index['by_name_parts'][query_lower]:
                    if scientific_name not in seen_species:
                        results.append(self._format_search_result(scientific_name))
                        seen_species.add(scientific_name)
            
            # Sort by relevance (exact matches first, then alphabetical)
            results.sort(key=lambda x: (
                not x['scientific_name'].lower().startswith(query_lower),
                x['scientific_name']
            ))
            
            return results[:limit]
            
        except Exception as e:
            logger.error(f"Error searching species: {str(e)}")
            return []
    
    def _format_search_result(self, scientific_name: str) -> Dict:
        """Format species data for search results"""
        species_data = self.species_data[scientific_name]
        
        return {
            "taxonid": species_data.get('sis_taxon_id'),
            "scientific_name": scientific_name,
            "kingdom": species_data.get('kingdom'),
            "phylum": species_data.get('phylum'),
            "class": species_data.get('class'),
            "order": species_data.get('order'),
            "family": species_data.get('family'),
            "genus": species_data.get('genus'),
            "category": species_data.get('category'),
            "red_list_code": species_data.get('red_list_code'),
            "possibly_extinct": species_data.get('possibly_extinct'),
            "possibly_extinct_in_the_wild": species_data.get('possibly_extinct_in_the_wild'),
            "year_published": species_data.get('year_published'),
            "url": species_data.get('url')
        }
    
    def get_species_details(self, species_name: str) -> Dict:
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
            
            # If not found, return default structure
            logger.warning(f"Species not found in IUCN data: {species_name}")
            return self._get_default_species_details(species_name)
            
        except Exception as e:
            logger.error(f"Error getting species details: {str(e)}")
            return self._get_default_species_details(species_name)
    
    def get_species_by_category(self, category: str) -> List[str]:
        """Get all species by conservation category"""
        return self.species_index['by_category'].get(category, [])
    
    def get_species_by_genus(self, genus: str) -> List[str]:
        """Get all species by genus"""
        return self.species_index['by_genus'].get(genus, [])
    
    def get_high_risk_species(self) -> List[str]:
        """Get species with high extinction risk"""
        return self.species_index['by_extinction_risk'].get('high_risk', [])
    
    def get_possibly_extinct_species(self) -> List[Dict]:
        """Get species that are possibly extinct"""
        extinct_species = []
        
        for scientific_name, species_data in self.species_data.items():
            if (species_data.get('possibly_extinct') or 
                species_data.get('possibly_extinct_in_the_wild')):
                extinct_species.append({
                    'scientific_name': scientific_name,
                    'category': species_data.get('category'),
                    'possibly_extinct': species_data.get('possibly_extinct'),
                    'possibly_extinct_in_the_wild': species_data.get('possibly_extinct_in_the_wild'),
                    'year_published': species_data.get('year_published')
                })
        
        return extinct_species
    
    def get_statistics(self) -> Dict:
        """Get statistics about loaded IUCN data"""
        stats = {
            'total_species': len(self.species_data),
            'by_category': {},
            'possibly_extinct': 0,
            'possibly_extinct_in_wild': 0,
            'by_genus_count': len(self.species_index['by_genus']),
            'latest_assessments': 0
        }
        
        for species_data in self.species_data.values():
            # Count by category
            category = species_data.get('category', 'NE')
            stats['by_category'][category] = stats['by_category'].get(category, 0) + 1
            
            # Count extinction risks
            if species_data.get('possibly_extinct'):
                stats['possibly_extinct'] += 1
            if species_data.get('possibly_extinct_in_the_wild'):
                stats['possibly_extinct_in_wild'] += 1
            if species_data.get('latest'):
                stats['latest_assessments'] += 1
        
        return stats
    
    def _get_default_species_details(self, species_name: str) -> Dict:
        """Return default species details structure"""
        return {
            'scientific_name': species_name,
            'category': 'NE',
            'red_list_code': 'NE',
            'possibly_extinct': False,
            'possibly_extinct_in_the_wild': False,
            'genus': self._extract_genus(species_name),
            'species': self._extract_species(species_name),
            'kingdom': 'ANIMALIA',
            'phylum': 'UNKNOWN',
            'class': 'UNKNOWN',
            'order': 'UNKNOWN',
            'family': 'UNKNOWN',
            'population_trend': 'Unknown',
            'marine_system': False,
            'freshwater_system': False,
            'terrestrial_system': True,
            'year_published': None,
            'latest': False,
            'assessment_id': None,
            'url': None,
            'sis_taxon_id': None
        }
    
    def get_all_species(self) -> Dict:
        """Get all species data"""
        return self.species_data
    
    def reload_data(self):
        """Reload all IUCN data from files"""
        logger.info("Reloading IUCN data...")
        self.species_data = {}
        self.species_index = {}
        self._load_all_iucn_data()

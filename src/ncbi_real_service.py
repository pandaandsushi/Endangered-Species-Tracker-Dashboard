"""
SYNCHRONIZED NCBI Real Service - properly integrated with IUCN service
"""

import os
import time
import logging
import traceback
from contextlib import contextmanager
from typing import Dict, List, Any, Optional, Union
import json
import re

# Biopython imports
from Bio import Entrez, SeqIO, AlignIO, Align
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class NCBIRealServiceFixed:
    """
    SYNCHRONIZED NCBI service with proper IUCN integration
    """
    
    def __init__(self, email: str = "raffaelsiahaan@gmail.com"):
        """Initialize NCBI service with IUCN integration"""
        self.email = email
        Entrez.email = email
        
        # Initialize IUCN service
        try:
            from iucn import IUCNService
            self.iucn_service = IUCNService()
            logger.info("✅ IUCN Service integrated successfully")
        except Exception as e:
            logger.warning(f"⚠️ IUCN Service initialization failed: {str(e)}")
            self.iucn_service = None

        # Set timeouts
        self.search_timeout = 300  # seconds
        self.blast_timeout = 300  # seconds
        self.fetch_timeout = 300   # seconds
        
        # Rate limiting
        self.last_request_time = 0
        self.min_request_interval = 0.5  # seconds
        
        logger.info(f"✅ SYNCHRONIZED NCBI Service initialized with email: {email}")
        logger.info(f"   Timeouts: search={self.search_timeout}s, blast={self.blast_timeout}s, fetch={self.fetch_timeout}s")
    
    @contextmanager
    def _timeout_context(self, operation: str, timeout: int):
        """Context manager for timeout handling"""
        start_time = time.time()
        try:
            # Rate limiting
            elapsed = time.time() - self.last_request_time
            if elapsed < self.min_request_interval:
                time.sleep(self.min_request_interval - elapsed)
            
            self.last_request_time = time.time()
            yield
            
        except Exception as e:
            elapsed = time.time() - start_time
            if elapsed >= timeout:
                logger.error(f"⏱️ {operation} timed out after {elapsed:.1f}s")
                raise TimeoutError(f"{operation} timed out after {elapsed:.1f} seconds")
            else:
                logger.error(f"❌ {operation} failed: {str(e)}")
                raise
    
    def _clean_species_name(self, species_name: str) -> str:
        """Clean species name for IUCN lookup"""
        if not species_name:
            return ""
        
        # Remove common prefixes that interfere with IUCN lookup
        cleaned = species_name.strip()
        
        # Remove "mitochondrion" prefix
        if cleaned.startswith("mitochondrion "):
            cleaned = cleaned.replace("mitochondrion ", "")
        
        # Remove subspecies (keep only genus + species)
        parts = cleaned.split()
        if len(parts) >= 2:
            # Keep only first two parts (genus + species)
            cleaned = f"{parts[0]} {parts[1]}"
        
        # Remove common suffixes in parentheses
        if "(" in cleaned:
            cleaned = cleaned.split("(")[0].strip()
        
        return cleaned
    
    def get_conservation_status_from_iucn(self, species_name: str) -> Dict[str, Any]:
        """Get conservation status from IUCN service with better name matching"""
        try:
            if not self.iucn_service:
                logger.warning("IUCN service not available, using default status")
                return {
                    "status": "DD",
                    "source": "default",
                    "year_published": None,
                    "possibly_extinct": False,
                    "possibly_extinct_in_the_wild": False
                }
        
            # Clean the species name for better matching
            original_name = species_name
            cleaned_name = self._clean_species_name(species_name)
        
            logger.info(f"Looking up IUCN data: '{original_name}' -> '{cleaned_name}'")
        
            # Try with cleaned name first
            species_details = self.iucn_service.get_species_details(cleaned_name)
        
            if species_details and species_details.get('taxon_scientific_name'):
                conservation_info = {
                    "status": species_details.get('red_list_category_code', 'DD'),
                    "source": "IUCN Red List",
                    "year_published": species_details.get('year_published'),
                    "possibly_extinct": species_details.get('possibly_extinct', False),
                    "possibly_extinct_in_the_wild": species_details.get('possibly_extinct_in_the_wild', False),
                    "assessment_id": species_details.get('assessment_id'),
                    "url": species_details.get('url'),
                    "sis_taxon_id": species_details.get('sis_taxon_id'),
                    "latest": species_details.get('latest', False),
                    "matched_name": species_details.get('taxon_scientific_name')
                }
            
                logger.info(f"✅ Found IUCN data for {cleaned_name}: {conservation_info['status']}")
                return conservation_info
            else:
                # Try searching for similar species names
                search_results = self.iucn_service.search_species(cleaned_name, limit=3)
            
                if search_results:
                    best_match = search_results[0]
                    conservation_info = {
                        "status": best_match.get('red_list_category_code', 'DD'),
                        "source": "IUCN Red List (search match)",
                        "year_published": best_match.get('year_published'),
                        "possibly_extinct": best_match.get('possibly_extinct', False),
                        "possibly_extinct_in_the_wild": best_match.get('possibly_extinct_in_the_wild', False),
                        "matched_name": best_match.get('taxon_scientific_name'),
                        "assessment_id": best_match.get('assessment_id'),
                        "url": best_match.get('url'),
                        "original_query": original_name,
                        "cleaned_query": cleaned_name
                    }
                
                    logger.info(f"✅ Found IUCN match for {cleaned_name} -> {best_match.get('taxon_scientific_name')}: {conservation_info['status']}")
                    return conservation_info
                else:
                    logger.warning(f"⚠️ No IUCN data found for {cleaned_name}, defaulting to DD")
                    return {
                        "status": "DD",
                        "source": "IUCN Red List (not found)",
                        "year_published": None,
                        "possibly_extinct": False,
                        "possibly_extinct_in_the_wild": False,
                        "original_query": original_name,
                        "cleaned_query": cleaned_name
                    }
                
        except Exception as e:
            logger.error(f"❌ Error getting IUCN status for {species_name}: {str(e)}")
            return {
                "status": "DD",
                "source": "error",
                "year_published": None,
                "possibly_extinct": False,
                "possibly_extinct_in_the_wild": False,
                "error": str(e),
                "original_query": species_name
            }
    
    def check_connection_safe(self) -> bool:
        """Check if NCBI connection is working"""
        try:
            with self._timeout_context("NCBI connection check", 10):
                # Simple search to test connection
                handle = Entrez.esearch(db="nucleotide", term="Homo sapiens[Organism] AND COI[Gene]", retmax=1)
                record = Entrez.read(handle)
                handle.close()
                
                return int(record["Count"]) > 0
        except Exception as e:
            logger.warning(f"NCBI connection check failed: {str(e)}")
            return False
    
    def get_sequence_info_real_safe(self, species_name: str, gene: str = "COI") -> Dict[str, Any]:
        """
        Get sequence information with proper error handling and timeout
        """
        try:
            logger.info(f"🔍 Getting sequence info for {species_name} ({gene})")
            
            # Step 1: Search for the sequence
            search_term = f"{species_name}[Organism] AND {gene}[Gene]"
            
            with self._timeout_context("NCBI search", self.search_timeout):
                search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1)
                search_results = Entrez.read(search_handle)
                search_handle.close()
            
            if not search_results["IdList"]:
                logger.warning(f"No {gene} sequences found for {species_name}")
                
                # Try broader search
                broader_term = f"{species_name}[Organism]"
                with self._timeout_context("NCBI broader search", self.search_timeout):
                    search_handle = Entrez.esearch(db="nucleotide", term=broader_term, retmax=1)
                    search_results = Entrez.read(search_handle)
                    search_handle.close()
                
                if not search_results["IdList"]:
                    logger.warning(f"No sequences found for {species_name}")
                    return {}
            
            # Step 2: Fetch the sequence
            seq_id = search_results["IdList"][0]
            
            with self._timeout_context("NCBI fetch", self.fetch_timeout):
                fetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
                record = SeqIO.read(fetch_handle, "genbank")
                fetch_handle.close()
            
            # Step 3: Extract and format information (JSON safe)
            sequence_info = {
                "accession": str(record.id),
                "length": int(len(record.seq)),
                "description": str(record.description),
                "sequence": str(record.seq)[:100] + "..." if len(record.seq) > 100 else str(record.seq),
                "organism": str(record.annotations.get("organism", species_name)),
                "gene": str(gene),
                "authors": str(", ".join(record.annotations.get("references", [{}])[0].authors.split(",")[:3])) 
                          if record.annotations.get("references") else "Unknown",
                "journal": str(record.annotations.get("references", [{}])[0].journal) 
                          if record.annotations.get("references") else "Unknown",
                "source": "NCBI Nucleotide Database"
            }
            
            logger.info(f"✅ Retrieved sequence info for {species_name}: {record.id}")
            return sequence_info
            
        except Exception as e:
            logger.error(f"❌ Error getting sequence info: {str(e)}")
            logger.error(traceback.format_exc())
            
            # Return minimal info
            return {
                "accession": "Error",
                "length": 0,
                "description": f"Error retrieving sequence: {str(e)}",
                "organism": species_name,
                "gene": gene,
                "error": str(e)
            }
    
    def search_similar_species_real(self, 
                                   species_name: str, 
                                   gene: str = "COI", 
                                   max_results: int = 20,
                                   min_similarity: float = 0.7) -> Dict[str, Any]:
        """
        Search for similar species using REAL NCBI BLAST with IUCN integration
        """
        try:
            logger.info(f"🔍 Searching for similar species to {species_name} ({gene})")
            
            # Step 1: Get query sequence
            query_info = self.get_sequence_info_real_safe(species_name, gene)
            
            accession = query_info['accession']
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = handle.read()
            import re
            match = re.search(r"SOURCE\s+(.*)", record)
            
            # PERBAIKAN: Simpan nama organisme dari GenBank sebagai variabel terpisah
            genbank_organism = match.group(1).strip() if match else species_name
            
            # Gunakan species_name asli untuk konsistensi
            logger.info(f"📋 Input species: {species_name}")
            logger.info(f"📋 GenBank organism: {genbank_organism}")

            if not query_info or "error" in query_info:
                logger.warning(f"❌ Could not find {gene} sequence for {organism}")
                return {
                    "query_species": species_name,
                    "genbank_organism": genbank_organism,
                    "query_gene": gene,
                    "query_sequence": "",
                    "query_accession": "Not found",
                    "similar_species": [],
                    "total_found": 0,
                    "search_method": "NCBI BLAST (failed)",
                    "error": f"Could not find {gene} sequence for {species_name}"
                }
            
            query_sequence = query_info.get("sequence", "")
            
            # Remove ellipsis if present
            if "..." in query_sequence:
                # Get full sequence
                with self._timeout_context("NCBI fetch full sequence", self.fetch_timeout):
                    fetch_handle = Entrez.efetch(
                        db="nucleotide", 
                        id=query_info["accession"], 
                        rettype="fasta", 
                        retmode="text"
                    )
                    record = SeqIO.read(fetch_handle, "fasta")
                    fetch_handle.close()
                    query_sequence = str(record.seq)
            
            if not query_sequence:
                logger.error(f"❌ Empty sequence for {organism}")
                return {
                    "query_species": organism,
                    "query_gene": gene,
                    "query_sequence": "",
                    "query_accession": query_info.get("accession", "Not found"),
                    "similar_species": [],
                    "total_found": 0,
                    "search_method": "NCBI BLAST (failed)",
                    "error": "Empty sequence"
                }
            
            try:
                with self._timeout_context("NCBI BLAST", self.blast_timeout):
                    blast_handle = NCBIWWW.qblast(
                        program="blastn",
                        database="nt",
                        sequence=query_sequence,
                        entrez_query=f"{gene}[Gene] NOT {species_name}[Organism]",
                        expect=0.001,
                        hitlist_size=max_results * 2  # Get more for filtering
                    )
                    blast_results = NCBIXML.read(blast_handle)
                    blast_handle.close()
            except Exception as e:
                logger.error(f"❌ BLAST search failed: {str(e)}")
                return {
                    "query_species": species_name,
                    "query_gene": gene,
                    "query_sequence": query_sequence[:50] + "..." if len(query_sequence) > 50 else query_sequence,
                    "query_accession": query_info.get("accession", "Not found"),
                    "similar_species": [],
                    "total_found": 0,
                    "search_method": "NCBI BLAST (failed)",
                    "error": f"BLAST search failed: {str(e)}"
                }
            
            # Step 3: Process BLAST results
            similar_species = []
            organism_hits = {}  # Track hits per organism
            
            for alignment in blast_results.alignments:
                for hsp in alignment.hsps:
                    # Calculate similarity
                    similarity = hsp.identities / hsp.align_length
                    
                    if similarity >= min_similarity:
                        # Extract organism name from accession
                        try:
                            fetch_handle = Entrez.efetch(
                                db="nucleotide", 
                                id=alignment.accession, 
                                rettype="gb", 
                                retmode="text"
                            )
                            record = fetch_handle.read()
                            fetch_handle.close()
                            
                            # Extract organism from SOURCE field
                            import re
                            match = re.search(r"SOURCE\s+(.*)", record)
                            if match:
                                organism = match.group(1).strip()
                            else:
                                # Fallback: extract from title
                                match = re.search(r'\[(.*?)\]', alignment.title)
                                organism = match.group(1) if match else alignment.title[:50]
                            
                        except Exception as e:
                            logger.warning(f"Could not extract organism for {alignment.accession}: {e}")
                            # Fallback to title parsing
                            match = re.search(r'\[(.*?)\]', alignment.title)
                            organism = match.group(1) if match else alignment.title[:50]
                        
                        # Skip if same as query species
                        if (organism.lower() == species_name.lower() or 
                            organism.lower() == genbank_organism.lower()):
                            continue
                        
                        # Get sequence
                        try:
                            fetch_handle = Entrez.efetch(
                                db="nucleotide", 
                                id=alignment.accession, 
                                rettype="fasta", 
                                retmode="text"
                            )
                            seq_record = SeqIO.read(fetch_handle, "fasta")
                            fetch_handle.close()
                            hit_sequence = str(seq_record.seq)
                        except Exception as e:
                            logger.warning(f"Could not fetch sequence for {alignment.accession}: {e}")
                            hit_sequence = hsp.sbjct
                        
                        # Create species data
                        species_data = {
                            "organism": organism,
                            "accession": str(alignment.accession),
                            "similarity": float(similarity),
                            "alignment_score": float(hsp.score),
                            "e_value": float(hsp.expect),
                            "identity": float(hsp.identities / hsp.align_length * 100),
                            "sequence": str(hit_sequence),
                            "status": "Unknown"  
                        }
                        
                        # Check if organism already exists
                        if organism not in organism_hits:
                            organism_hits[organism] = species_data
                        else:
                            # Compare with existing hit for same organism
                            current_best = organism_hits[organism]
                            
                            # Use better hit based on multiple criteria
                            if (species_data["similarity"] > current_best["similarity"] or
                                (species_data["similarity"] == current_best["similarity"] and
                                species_data["alignment_score"] > current_best["alignment_score"]) or
                                (species_data["similarity"] == current_best["similarity"] and
                                species_data["alignment_score"] == current_best["alignment_score"] and
                                species_data["e_value"] < current_best["e_value"])):
                                
                                logger.info(f"Replacing {organism} with better hit: "
                                        f"sim {current_best['similarity']:.3f}→{species_data['similarity']:.3f}")
                                organism_hits[organism] = species_data
            
            similar_species = list(organism_hits.values())
            similar_species.sort(key=lambda x: (
                -x["similarity"],
                -x["alignment_score"], 
                x["e_value"]
            ))
            similar_species = similar_species[:max_results]

            # Step 4: Add REAL conservation status from IUCN
            logger.info("🔍 Getting conservation status from IUCN for all species...")
            
            for species in similar_species:
                conservation_info = self.get_conservation_status_from_iucn(species["organism"])
                species.update({
                    "status": conservation_info["status"],
                    "conservation_source": conservation_info["source"],
                    "iucn_year_published": conservation_info.get("year_published"),
                    "possibly_extinct": conservation_info.get("possibly_extinct", False),
                    "possibly_extinct_in_the_wild": conservation_info.get("possibly_extinct_in_the_wild", False),
                    "iucn_url": conservation_info.get("url"),
                    "assessment_id": conservation_info.get("assessment_id"),
                    "sis_taxon_id": conservation_info.get("sis_taxon_id")
                })
            
            # Step 5: Return results (JSON safe)
            results = {
                "query_species": str(species_name),
                "genbank_organism": str(genbank_organism),
                "query_gene": str(gene),
                "query_sequence": str(query_sequence),
                "query_accession": str(query_info.get("accession", "Not found")),
                "similar_species": similar_species,
                "total_found": int(len(similar_species)),
                "search_method": "NCBI BLAST + IUCN Red List (real data)"
            }
            
            logger.info(f"✅ Found {len(similar_species)} similar species to {organism} with IUCN conservation data")
            return results
            
        except Exception as e:
            logger.error(f"❌ Error in similarity search: {str(e)}")
            logger.error(traceback.format_exc())
            
            # Return error results
            return {
                "query_species": str(species_name),
                "query_gene": str(gene),
                "query_sequence": "",
                "query_accession": "Error",
                "similar_species": [],
                "total_found": 0,
                "search_method": "NCBI BLAST + IUCN (error)",
                "error": str(e)
            }

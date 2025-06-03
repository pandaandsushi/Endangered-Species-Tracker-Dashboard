"""
FIXED NCBI Real Service - dengan proper error handling dan JSON serialization
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
    FIXED NCBI service for real data retrieval with proper error handling
    """
    
    def __init__(self, email: str = "raffaelsiahaan@gmail.com"):
        """Initialize NCBI service with proper timeouts"""
        self.email = email
        Entrez.email = email
        
        # Set timeouts
        self.search_timeout = 240  # seconds
        self.blast_timeout = 120  # seconds
        self.fetch_timeout = 240   # seconds
        
        # Rate limiting
        self.last_request_time = 0
        self.min_request_interval = 0.5  # seconds
        
        logger.info(f"✅ FIXED NCBI Service initialized with email: {email}")
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
                                   max_results: int = 8,
                                   min_similarity: float = 0.7) -> Dict[str, Any]:
        """
        Search for similar species using REAL NCBI BLAST with proper error handling
        """
        try:
            logger.info(f"🔍 Searching for similar species to {species_name} ({gene})")
            
            # Step 1: Get query sequence
            query_info = self.get_sequence_info_real_safe(species_name, gene)
            
            if not query_info or "error" in query_info:
                logger.warning(f"❌ Could not find {gene} sequence for {species_name}")
                return {
                    "query_species": species_name,
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
                logger.error(f"❌ Empty sequence for {species_name}")
                return {
                    "query_species": species_name,
                    "query_gene": gene,
                    "query_sequence": "",
                    "query_accession": query_info.get("accession", "Not found"),
                    "similar_species": [],
                    "total_found": 0,
                    "search_method": "NCBI BLAST (failed)",
                    "error": "Empty sequence"
                }
            
            # Step 2: Perform BLAST search
            logger.info(f"🔍 Running BLAST search for {species_name} ({len(query_sequence)} bp)")
            
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
            
            for alignment in blast_results.alignments:
                for hsp in alignment.hsps:
                    # Calculate similarity
                    similarity = hsp.identities / hsp.align_length
                    
                    if similarity >= min_similarity:
                        # Extract species name from hit title
                        hit_title = alignment.title
                        match = re.search(r'\[(.*?)\]', hit_title)
                        if match:
                            hit_species = match.group(1)
                        else:
                            # Try to extract from title
                            parts = hit_title.split()
                            if len(parts) >= 2:
                                hit_species = f"{parts[0]} {parts[1]}"
                            else:
                                hit_species = hit_title[:30]
                        
                        # Get sequence
                        try:
                            with self._timeout_context("NCBI fetch hit", self.fetch_timeout):
                                fetch_handle = Entrez.efetch(
                                    db="nucleotide", 
                                    id=alignment.accession, 
                                    rettype="fasta", 
                                    retmode="text"
                                )
                                record = SeqIO.read(fetch_handle, "fasta")
                                fetch_handle.close()
                                hit_sequence = str(record.seq)
                        except Exception as e:
                            logger.warning(f"Could not fetch sequence for {alignment.accession}: {e}")
                            hit_sequence = hsp.sbjct
                        
                        # Add to results (JSON safe)
                        similar_species.append({
                            "name": str(hit_species),
                            "accession": str(alignment.accession),
                            "similarity": float(similarity),
                            "alignment_score": float(hsp.score),
                            "e_value": float(hsp.expect),
                            "identity": float(hsp.identities / hsp.align_length * 100),
                            "sequence": str(hit_sequence),
                            "status": "Unknown"  # Will be updated with IUCN data
                        })
            
            # Sort by similarity
            similar_species.sort(key=lambda x: x["similarity"], reverse=True)
            
            # Limit results
            similar_species = similar_species[:max_results]
            
            # Step 4: Add conservation status (mock for now)
            conservation_statuses = {
                "Panthera leo": "VU",
                "Panthera tigris": "EN",
                "Panthera pardus": "VU",
                "Panthera onca": "NT",
                "Acinonyx jubatus": "VU",
                "Puma concolor": "LC"
            }
            
            for species in similar_species:
                species_name = species["name"]
                if species_name in conservation_statuses:
                    species["status"] = conservation_statuses[species_name]
                else:
                    species["status"] = "DD"  # Data Deficient
            
            # Step 5: Return results (JSON safe)
            results = {
                "query_species": str(species_name),
                "query_gene": str(gene),
                "query_sequence": str(query_sequence),
                "query_accession": str(query_info.get("accession", "Not found")),
                "similar_species": similar_species,
                "total_found": int(len(similar_species)),
                "search_method": "NCBI BLAST (real data)"
            }
            
            logger.info(f"✅ Found {len(similar_species)} similar species to {species_name}")
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
                "search_method": "NCBI BLAST (error)",
                "error": str(e)
            }

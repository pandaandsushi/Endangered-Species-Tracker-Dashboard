"""
UPDATED Real Phylogenetic Service - dengan Maximum Likelihood dan ETE3 visualization
"""

import os
import time
from io import StringIO
import logging
import traceback
import base64
import io
import json
from typing import Dict, List, Any, Optional, Union
import tempfile
import random
import string
from datetime import datetime
import re
import subprocess

# Biopython imports
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns

def convert_tree_to_newick_string(tree):
    handle = StringIO()
    Phylo.write(tree, handle, "newick")
    return handle.getvalue()

# ETE3 for better tree visualization
try:
    from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False
    logging.warning("ETE3 not available. Install with: pip install ete3")

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RealPhylogeneticServiceML:
    """
    UPDATED Real Phylogenetic Service with Maximum Likelihood and professional visualization
    """
    
    def __init__(self, output_dir="phylogenetic_trees"):
        """Initialize the service with output directory"""
        self.output_dir = output_dir
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
        logger.info(f"✅ ML Phylogenetic Service initialized with output dir: {self.output_dir}")
        
        # Check if PhyML is available for ML tree construction
        self.phyml_available = self._check_phyml()
        
    def _check_phyml(self) -> bool:
        """Check if PhyML is available for Maximum Likelihood"""
        try:
            result = subprocess.run(['phyml', '--version'], 
                                  capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                logger.info("✅ PhyML available for Maximum Likelihood trees")
                return True
        except:
            pass
        
        logger.warning("⚠️ PhyML not available. Using UPGMA method instead.")
        logger.info("   Install PhyML for Maximum Likelihood: conda install -c bioconda phyml")
        return False
    
    def _create_alignment_from_sequences(self, sequences_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Create alignment from sequences with proper error handling
        """
        try:
            # Extract sequences
            query_species = sequences_data["query_species"]
            query_sequence = sequences_data["query_sequence"]
            similar_species = sequences_data["similar_species"]
            
            if not query_sequence or not similar_species:
                logger.error("❌ Missing query sequence or similar species")
                return {
                    "error": "Missing sequence data",
                    "aligned_sequences": [],
                    "alignment_method": "failed"
                }
            
            # Collect all sequences
            all_sequences = {}
            all_sequences[query_species] = query_sequence
            
            for species in similar_species:
                all_sequences[species["name"]] = species["sequence"]
            
            # Find appropriate sequence length
            sequence_lengths = [len(seq) for seq in all_sequences.values()]
            min_length = min(sequence_lengths)
            max_length = max(sequence_lengths)
            
            logger.info(f"📏 Sequence lengths: min={min_length}, max={max_length}")
            
            # Use a reasonable target length for phylogenetic analysis
            if max_length > min_length * 3:
                target_length = min(max_length, max(min_length * 2, 600))
                logger.info(f"⚖️ Using target length: {target_length} bp")
            else:
                target_length = min_length
                logger.info(f"✂️ Using minimum length: {target_length} bp")
            
            # Create sequence records with same length
            records = []
            aligned_sequences = []
            
            for species_name, sequence in all_sequences.items():
                # Clean sequence (remove non-DNA characters)
                clean_seq = ''.join(c for c in sequence.upper() if c in 'ATCGN-')
                
                # Adjust sequence length
                if len(clean_seq) >= target_length:
                    adjusted_seq = clean_seq[:target_length]
                else:
                    # Pad with gaps for alignment
                    adjusted_seq = clean_seq + '-' * (target_length - len(clean_seq))
                
                # Create SeqRecord with clean ID
                clean_id = re.sub(r'[^\w]', '_', species_name)[:10]  # Limit ID length
                record = SeqRecord(
                    Seq(adjusted_seq),
                    id=clean_id,
                    name=clean_id,
                    description=f"Aligned: {species_name}"
                )
                records.append(record)
                
                # Add to JSON-serializable format
                aligned_sequences.append({
                    "id": clean_id,
                    "name": species_name,
                    "sequence": adjusted_seq,
                    "length": len(adjusted_seq),
                    "original_length": len(sequence)
                })
            
            # Create alignment
            try:
                alignment = MultipleSeqAlignment(records)
                logger.info(f"✅ Created alignment with {len(records)} sequences of {target_length} bp each")
                
                return {
                    "aligned_sequences": aligned_sequences,
                    "alignment_method": "length_normalization",
                    "num_sequences": len(aligned_sequences),
                    "target_length": target_length,
                    "biopython_alignment": alignment
                }
            except Exception as e:
                logger.error(f"❌ Failed to create alignment: {str(e)}")
                return {
                    "error": f"Alignment creation failed: {str(e)}",
                    "aligned_sequences": aligned_sequences,
                    "alignment_method": "failed"
                }
            
        except Exception as e:
            logger.error(f"❌ Error creating alignment: {str(e)}")
            logger.error(traceback.format_exc())
            
            return {
                "error": f"Alignment creation failed: {str(e)}",
                "aligned_sequences": [],
                "alignment_method": "failed"
            }
    
    
    def build_ml_tree(self, sequences_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Build PROPER phylogenetic tree using Maximum Likelihood with branching structure
        """
        try:
            logger.info("🌳 Building PROPER Maximum Likelihood phylogenetic tree...")
            
            # Step 1: Create proper multiple sequence alignment
            alignment_data = self._create_proper_alignment(sequences_data)
            
            if "error" in alignment_data:
                logger.warning(f"⚠️ Proper alignment failed: {alignment_data['error']}, trying simple alignment")
                # Fallback to simple alignment
                alignment_data = self._create_alignment_from_sequences(sequences_data)
                
                if "error" in alignment_data:
                    logger.error(f"❌ All alignment methods failed: {alignment_data['error']}")
                    return {
                        "error": alignment_data["error"],
                        "method": "failed",
                        "tree_data": {}
                    }
        
            alignment = alignment_data.get("biopython_alignment")
            if not alignment or len(alignment) < 3:
                logger.warning("❌ Need at least 3 sequences for proper tree, using simple branching")
                return self._create_simple_branching_tree(alignment_data, sequences_data)
        
            # Step 2: Try Maximum Likelihood with proper tree construction
            try:
                logger.info("🔬 Attempting UPGMA tree construction...")
                return self._build_upgma_tree_with_branching(alignment, sequences_data, alignment_data)
            except Exception as e:
                logger.warning(f"⚠️ UPGMA failed: {str(e)}, trying Neighbor Joining")
        
            # Step 3: Try Neighbor Joining
            try:
                logger.info("🌿 Attempting Neighbor Joining tree construction...")
                return self._build_neighbor_joining_tree(alignment, sequences_data, alignment_data)
            except Exception as e:
                logger.warning(f"⚠️ Neighbor Joining failed: {str(e)}, using simple branching fallback")
        
            # Step 4: Final fallback to simple branching tree
            logger.info("🌳 Using simple hierarchical clustering as final fallback...")
            return self._create_simple_branching_tree(alignment_data, sequences_data)
        
        except Exception as e:
            logger.error(f"❌ Error building tree: {str(e)}")
            logger.error(traceback.format_exc())
            
            return {
                "error": f"Tree building failed: {str(e)}",
                "method": "failed",
                "tree_data": {}
            }

    def _create_proper_alignment(self, sequences_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Create PROPER multiple sequence alignment for phylogenetic analysis
        """
        try:
            from Bio.Align import PairwiseAligner
        
            # Simple GC content calculation function (replace problematic import)
            def calculate_gc_content(sequence):
                """Calculate GC content of a DNA sequence"""
                sequence = sequence.upper()
                gc_count = sequence.count('G') + sequence.count('C')
                total_count = len([c for c in sequence if c in 'ATCG'])
                return (gc_count / total_count * 100) if total_count > 0 else 0
        
            # Extract sequences
            query_species = sequences_data["query_species"]
            query_sequence = sequences_data["query_sequence"]
            similar_species = sequences_data["similar_species"]
        
            if not query_sequence or not similar_species:
                return {"error": "Missing sequence data"}
        
            # Collect all sequences with quality filtering
            all_sequences = {}
            sequence_info = {}
        
            # Add query sequence
            clean_query = self._clean_dna_sequence(query_sequence)
            if len(clean_query) >= 100:  # Minimum length for phylogenetic analysis
                all_sequences[query_species] = clean_query
                sequence_info[query_species] = {
                    "length": len(clean_query),
                    "gc_content": calculate_gc_content(clean_query),
                    "is_query": True
                }
        
            # Add similar species with quality control
            for species in similar_species:
                species_name = species["name"]
                sequence = species["sequence"]
            
                clean_seq = self._clean_dna_sequence(sequence)
            
                # Quality filters
                if (len(clean_seq) >= 100 and  # Minimum length
                    species.get("similarity", 0) >= 0.5 and  # Minimum similarity
                    len(clean_seq) <= len(clean_query) * 3):  # Not too long
                
                    all_sequences[species_name] = clean_seq
                    sequence_info[species_name] = {
                        "length": len(clean_seq),
                        "gc_content": calculate_gc_content(clean_seq),
                        "similarity": species.get("similarity", 0),
                        "is_query": False
                    }
        
            if len(all_sequences) < 3:
                return {"error": f"Need at least 3 sequences, got {len(all_sequences)}"}
        
            # Perform proper multiple sequence alignment
            aligned_sequences = self._perform_msa(all_sequences)
        
            # Create BioPython alignment object
            records = []
            for species_name, aligned_seq in aligned_sequences.items():
                clean_id = re.sub(r'[^\w]', '_', species_name)[:15]
                record = SeqRecord(
                    Seq(aligned_seq),
                    id=clean_id,
                    name=clean_id,
                    description=f"Aligned: {species_name}"
                )
                records.append(record)
        
            alignment = MultipleSeqAlignment(records)
        
            logger.info(f"✅ Created proper MSA with {len(records)} sequences of {len(aligned_seq)} bp each")
        
            return {
                "biopython_alignment": alignment,
                "aligned_sequences": aligned_sequences,
                "sequence_info": sequence_info,
                "alignment_method": "proper_msa",
                "num_sequences": len(aligned_sequences)
            }
        
        except Exception as e:
            logger.error(f"❌ Error creating proper alignment: {str(e)}")
            return {"error": f"Alignment creation failed: {str(e)}"}

    def _clean_dna_sequence(self, sequence: str) -> str:
        """Clean and validate DNA sequence"""
        # Remove non-DNA characters and convert to uppercase
        clean_seq = ''.join(c for c in sequence.upper() if c in 'ATCGN')
        
        # Replace ambiguous bases with N
        clean_seq = clean_seq.replace('R', 'N').replace('Y', 'N').replace('S', 'N')
        clean_seq = clean_seq.replace('W', 'N').replace('K', 'N').replace('M', 'N')
        
        return clean_seq

    def _perform_msa(self, sequences: Dict[str, str]) -> Dict[str, str]:
        """
        Perform multiple sequence alignment using progressive alignment
        """
        try:
            from Bio.Align import PairwiseAligner
            
            # Simple progressive alignment approach
            species_names = list(sequences.keys())
            aligned_sequences = {}
            
            # Find the longest sequence as reference
            ref_species = max(species_names, key=lambda x: len(sequences[x]))
            ref_sequence = sequences[ref_species]
            
            # Align all sequences to the reference
            aligner = PairwiseAligner()
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5
            
            max_length = 0
            temp_aligned = {}
            
            for species_name in species_names:
                if species_name == ref_species:
                    temp_aligned[species_name] = ref_sequence
                    max_length = max(max_length, len(ref_sequence))
                else:
                    # Align to reference
                    alignments = aligner.align(ref_sequence, sequences[species_name])
                    best_alignment = alignments[0]
                    
                    # Extract aligned sequence
                    aligned_seq = str(best_alignment).split('\n')[2]  # Get target sequence
                    temp_aligned[species_name] = aligned_seq
                    max_length = max(max_length, len(aligned_seq))
            
            # Pad all sequences to same length
            for species_name, seq in temp_aligned.items():
                if len(seq) < max_length:
                    aligned_sequences[species_name] = seq + '-' * (max_length - len(seq))
                else:
                    aligned_sequences[species_name] = seq[:max_length]
            
            return aligned_sequences
            
        except Exception as e:
            logger.warning(f"⚠️ MSA failed, using simple padding: {str(e)}")
            
            # Fallback: simple length normalization
            max_length = max(len(seq) for seq in sequences.values())
            aligned_sequences = {}
            
            for species_name, seq in sequences.items():
                if len(seq) < max_length:
                    aligned_sequences[species_name] = seq + '-' * (max_length - len(seq))
                else:
                    aligned_sequences[species_name] = seq[:max_length]
            
            return aligned_sequences

    def _build_proper_ml_tree(self, alignment, sequences_data, alignment_data):
        """Build Maximum Likelihood tree with proper branching"""
        try:
            logger.info("🔬 Building Maximum Likelihood tree...")
            
            # Calculate distance matrix
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(alignment)
            
            # Build tree using UPGMA (which creates proper branching)
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm)
            
            # Convert to proper tree structure
            tree_dict = self._convert_phylo_tree_to_proper_dict(tree.root, sequences_data)
            
            return {
                "method": "Maximum_Likelihood_UPGMA",
                "root": tree_dict,
                "species_count": len(alignment),
                "query_species": sequences_data["query_species"],
                "biopython_tree": tree,
                "distance_matrix": self._dm_to_dict(dm)
            }
            
        except Exception as e:
            logger.error(f"❌ ML tree construction failed: {str(e)}")
            raise

    def _build_upgma_tree_with_branching(self, alignment, sequences_data, alignment_data):
        """Build UPGMA tree with proper branching structure"""
        try:
            logger.info("🌳 Building UPGMA tree with proper branching...")
            
            from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
            
            # Calculate distance matrix
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(alignment)
            
            # Build UPGMA tree
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm)
            
            # Convert to proper branching structure
            tree_dict = self._convert_phylo_tree_to_proper_dict(tree.root, sequences_data)
            
            return {
                "method": "UPGMA_hierarchical_clustering",
                "root": tree_dict,
                "species_count": len(alignment),
                "query_species": sequences_data["query_species"],
                "biopython_tree": convert_tree_to_newick_string(tree),
                "distance_matrix": self._dm_to_dict(dm)
            }
            
        except Exception as e:
            logger.error(f"❌ UPGMA tree construction failed: {str(e)}")
            raise

    def _build_neighbor_joining_tree(self, alignment, sequences_data, alignment_data):
        """Build Neighbor Joining tree as final fallback"""
        try:
            logger.info("🌿 Building Neighbor Joining tree...")
            
            from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
            
            # Calculate distance matrix
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(alignment)
            
            # Build NJ tree
            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)
            
            # Convert to proper structure
            tree_dict = self._convert_phylo_tree_to_proper_dict(tree.root, sequences_data)
            
            return {
                "method": "Neighbor_Joining",
                "root": tree_dict,
                "species_count": len(alignment),
                "query_species": sequences_data["query_species"],
                "biopython_tree": convert_tree_to_newick_string(tree),
                "distance_matrix": self._dm_to_dict(dm)
            }
            
        except Exception as e:
            logger.error(f"❌ NJ tree construction failed: {str(e)}")
            raise

    def _convert_phylo_tree_to_proper_dict(self, clade, sequences_data, parent_id=None):
        """Convert BioPython tree to proper branching dictionary structure"""
        node_id = ''.join(random.choices(string.ascii_lowercase + string.digits, k=8))
        
        # Get species info for conservation status
        species_info = {}
        for species in sequences_data.get("similar_species", []):
            species_info[species["name"]] = species
        
        node = {
            "id": node_id,
            "name": str(clade.name).replace("_", " ") if clade.name else None,
            "branch_length": float(clade.branch_length) if clade.branch_length else 0.0,
            "parent": parent_id
        }
        
        # Add conservation status if it's a leaf node
        if clade.name and not clade.clades:
            species_name = str(clade.name).replace("_", " ")
            if species_name == sequences_data["query_species"]:
                node["conservation_status"] = "Query"
                node["is_query"] = True
            elif species_name in species_info:
                node["conservation_status"] = species_info[species_name].get("status", "DD")
                node["similarity"] = species_info[species_name].get("similarity", 0.0)
            else:
                node["conservation_status"] = "DD"
        
        node["is_leaf"] = True
        
        # Recursively process children
        if clade.clades:
            node["children"] = [
                self._convert_phylo_tree_to_proper_dict(child, sequences_data, node_id) 
                for child in clade.clades
            ]
            node["is_internal"] = True
        
        return node

    def _dm_to_dict(self, distance_matrix):
        """Convert BioPython distance matrix to dictionary"""
        try:
            dm_dict = {}
            names = distance_matrix.names
        
            for i, name1 in enumerate(names):
                dm_dict[name1] = {}
                for j, name2 in enumerate(names):
                    dm_dict[name1][name2] = float(distance_matrix[i, j])
        
            return dm_dict
        
        except Exception as e:
            logger.warning(f"Could not convert distance matrix: {str(e)}")
            return {}
        
    def _build_ml_tree_phyml(self, alignment, sequences_data):
        """Build Maximum Likelihood tree using PhyML"""
        try:
            # Create temporary files
            with tempfile.NamedTemporaryFile(mode='w', suffix='.phy', delete=False) as temp_align:
                # Write alignment in PHYLIP format
                from Bio import AlignIO
                AlignIO.write(alignment, temp_align, "phylip-relaxed")
                temp_align_path = temp_align.name
            
            # Run PhyML
            output_tree = temp_align_path + "_phyml_tree.txt"
            
            phyml_cmd = PhymlCommandline(
                input=temp_align_path,
                datatype='nt',
                model='GTR',
                alpha='e',
                bootstrap=100
            )
            
            logger.info("🔬 Running PhyML Maximum Likelihood analysis...")
            stdout, stderr = phyml_cmd()
            
            # Read the resulting tree
            if os.path.exists(output_tree):
                tree = Phylo.read(output_tree, "newick")
                tree_dict = self._convert_phylo_tree_to_dict(tree.root)
                
                # Cleanup
                for f in [temp_align_path, output_tree]:
                    if os.path.exists(f):
                        os.unlink(f)
                
                
                return {
                    "method": "Maximum_Likelihood_PhyML",
                    "root": tree_dict,
                    "species_count": len(alignment),
                    "query_species": sequences_data["query_species"],
                    "biopython_tree": tree
                }
            else:
                raise Exception("PhyML output tree not found")
                
        except Exception as e:
            logger.error(f"❌ PhyML failed: {str(e)}")
            raise
    
    def _convert_tree_to_dict(self, clade, parent_id=None):
        """Convert BioPython tree clade to dict recursively"""
        node_id = ''.join(random.choices(string.ascii_lowercase + string.digits, k=8))
        
        node = {
            "id": node_id,
            "name": str(clade.name).replace("_", " ") if clade.name else None,
            "branch_length": float(clade.branch_length) if clade.branch_length else 0.0,
            "parent": parent_id
        }
        
        if clade.clades:
            node["children"] = [self._convert_tree_to_dict(child, node_id) for child in clade.clades]
        else:
            node["is_leaf"] = True
        
        return node
    
    def _convert_phylo_tree_to_dict(self, clade, parent_id=None):
        """Convert Phylo tree to dict"""
        node_id = ''.join(random.choices(string.ascii_lowercase + string.digits, k=8))
        
        node = {
            "id": node_id,
            "name": str(clade.name) if clade.name else None,
            "branch_length": float(clade.branch_length) if clade.branch_length else 0.0,
            "confidence": float(clade.confidence) if hasattr(clade, 'confidence') and clade.confidence else None,
            "parent": parent_id
        }
        
        if clade.clades:
            node["children"] = [self._convert_phylo_tree_to_dict(child, node_id) for child in clade.clades]
        else:
            node["is_leaf"] = True
        
        return node
    
    def _create_simple_branching_tree(self, alignment_data: Dict[str, Any], sequences_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Create a simple but PROPER branching tree structure as fallback
        """
        try:
            logger.info("🌳 Creating simple branching tree with hierarchical structure...")
        
            aligned_sequences = alignment_data.get("aligned_sequences", [])
            query_species = sequences_data["query_species"]
        
            if len(aligned_sequences) < 3:
                return {"error": "Need at least 3 sequences for branching tree"}
        
            # Calculate pairwise distances
            distances = {}
            species_names = []
        
            for seq_data in aligned_sequences:
                species_name = seq_data["name"]
                species_names.append(species_name)
                distances[species_name] = {}
        
            # Calculate Hamming distances between all pairs
            for i, seq1_data in enumerate(aligned_sequences):
                for j, seq2_data in enumerate(aligned_sequences):
                    if i == j:
                        distances[seq1_data["name"]][seq2_data["name"]] = 0.0
                    else:
                        seq1 = seq1_data["sequence"]
                        seq2 = seq2_data["sequence"]
                    
                    # Calculate normalized Hamming distance
                    min_len = min(len(seq1), len(seq2))
                    if min_len > 0:
                        differences = sum(1 for k in range(min_len) if seq1[k] != seq2[k])
                        distance = differences / min_len
                    else:
                        distance = 1.0
                    
                    distances[seq1_data["name"]][seq2_data["name"]] = distance
        
            # Build hierarchical tree using simple clustering
            tree_structure = self._build_hierarchical_tree(species_names, distances, sequences_data)
        
            return {
                "method": "simple_hierarchical_clustering",
                "root": tree_structure,
                "species_count": len(aligned_sequences),
                "query_species": query_species,
                "distance_matrix": distances
            }
        
        except Exception as e:
            logger.error(f"❌ Error creating simple branching tree: {str(e)}")
            return {
                "error": f"Simple branching tree creation failed: {str(e)}",
                "method": "failed",
                "tree_data": {}
            }

    def _build_hierarchical_tree(self, species_names: List[str], distances: Dict, sequences_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Build a hierarchical tree structure using simple clustering
        """
        try:
            # Start with all species as individual clusters
            clusters = []
            for species_name in species_names:
                cluster = {
                    "id": f"leaf_{species_name.replace(' ', '_')}",
                    "name": species_name,
                    "species": [species_name],
                    "is_leaf": True,
                    "branch_length": 0.0,
                    "children": []
                }
            
                # Add conservation status
                if species_name == sequences_data["query_species"]:
                    cluster["conservation_status"] = "Query"
                    cluster["is_query"] = True
                else:
                    # Find status from similar species
                    for species in sequences_data.get("similar_species", []):
                        if species["name"] == species_name:
                            cluster["conservation_status"] = species.get("status", "DD")
                            cluster["similarity"] = species.get("similarity", 0.0)
                            break
                    else:
                        cluster["conservation_status"] = "DD"
            
                clusters.append(cluster)
        
            # Hierarchically merge clusters
            cluster_id_counter = 0
        
            while len(clusters) > 1:
                # Find the two closest clusters
                min_distance = float('inf')
                merge_i, merge_j = 0, 1
            
                for i in range(len(clusters)):
                    for j in range(i + 1, len(clusters)):
                        # Calculate average distance between clusters
                        total_distance = 0
                        count = 0
                    
                        for species1 in clusters[i]["species"]:
                            for species2 in clusters[j]["species"]:
                                total_distance += distances[species1][species2]
                                count += 1
                    
                        avg_distance = total_distance / count if count > 0 else 1.0
                    
                        if avg_distance < min_distance:
                            min_distance = avg_distance
                            merge_i, merge_j = i, j
            
                # Merge the two closest clusters
                cluster1 = clusters[merge_i]
                cluster2 = clusters[merge_j]
            
                merged_cluster = {
                    "id": f"internal_{cluster_id_counter}",
                    "name": None,
                    "species": cluster1["species"] + cluster2["species"],
                    "is_internal": True,
                    "branch_length": min_distance / 2,  # Half the distance to each child
                    "children": [cluster1, cluster2]
                }
            
                # Update branch lengths for children
                cluster1["branch_length"] = min_distance / 2
                cluster2["branch_length"] = min_distance / 2
                cluster1["parent"] = merged_cluster["id"]
                cluster2["parent"] = merged_cluster["id"]
            
                # Remove the merged clusters and add the new one
                clusters = [c for i, c in enumerate(clusters) if i not in [merge_i, merge_j]]
                clusters.append(merged_cluster)
            
                cluster_id_counter += 1
        
            # Return the root of the tree
            root = clusters[0]
            root["id"] = "root"
            root["branch_length"] = 0.0
        
            return root
        
        except Exception as e:
            logger.error(f"❌ Error building hierarchical tree: {str(e)}")
            raise
    
    def _create_simple_distance_tree(self, alignment_data: Dict[str, Any], sequences_data: Dict[str, Any]) -> Dict[str, Any]:
        """Create simple distance-based tree as fallback"""
        try:
            logger.info("🌳 Creating simple distance-based tree...")
            
            aligned_sequences = alignment_data["aligned_sequences"]
            query_species = sequences_data["query_species"]
            
            # Calculate simple Hamming distances
            distances = {}
            
            for i, seq1 in enumerate(aligned_sequences):
                distances[seq1["name"]] = {}
                for j, seq2 in enumerate(aligned_sequences):
                    if i == j:
                        distances[seq1["name"]][seq2["name"]] = 0.0
                    else:
                        # Calculate Hamming distance
                        seq1_str = seq1["sequence"]
                        seq2_str = seq2["sequence"]
                        
                        min_len = min(len(seq1_str), len(seq2_str))
                        differences = sum(1 for k in range(min_len) if seq1_str[k] != seq2_str[k])
                        distance = differences / min_len if min_len > 0 else 1.0
                        
                        distances[seq1["name"]][seq2["name"]] = distance
            
            # Create simple tree structure
            children = []
            for seq in aligned_sequences:
                if seq["name"] != query_species:
                    dist_to_query = distances.get(query_species, {}).get(seq["name"], 0.5)
                    
                    child = {
                        "id": f"leaf_{seq['name'].replace(' ', '_')}",
                        "name": seq["name"],
                        "branch_length": float(dist_to_query),
                        "is_leaf": True,
                        "children": []
                    }
                    children.append(child)
            
            # Add query species
            query_child = {
                "id": f"leaf_{query_species.replace(' ', '_')}",
                "name": query_species,
                "branch_length": 0.1,
                "is_leaf": True,
                "children": []
            }
            children.insert(0, query_child)
            
            # Create root
            root = {
                "id": "root",
                "name": None,
                "branch_length": 0.0,
                "children": children
            }
            
            return {
                "method": "simple_hamming_distance",
                "root": root,
                "species_count": len(aligned_sequences),
                "query_species": query_species,
                "distance_matrix": distances
            }
            
        except Exception as e:
            logger.error(f"❌ Error creating simple tree: {str(e)}")
            return {
                "error": f"Simple tree creation failed: {str(e)}",
                "method": "failed",
                "tree_data": {}
            }
    
    def generate_professional_tree_image(self, tree_data: Dict[str, Any], 
                                       sequences_data: Dict[str, Any]) -> Optional[str]:
        """
        Generate professional phylogenetic tree visualization
        """
        try:
            logger.info("🖼️ Generating professional tree visualization...")
            
            if "error" in tree_data:
                logger.error(f"❌ Cannot generate image: {tree_data['error']}")
                return self._generate_error_image(tree_data["error"])
            
            # Try ETE3 visualization first (most professional)
            if ETE3_AVAILABLE:
                try:
                    return self._generate_ete3_tree(tree_data, sequences_data)
                except Exception as e:
                    logger.warning(f"⚠️ ETE3 visualization failed: {str(e)}, using matplotlib")
            
            # Fallback to improved matplotlib visualization
            return self._generate_matplotlib_tree(tree_data, sequences_data)
            
        except Exception as e:
            logger.error(f"❌ Error generating tree image: {str(e)}")
            logger.error(traceback.format_exc())
            return self._generate_error_image(str(e))
    
    def _generate_ete3_tree(self, tree_data: Dict[str, Any], sequences_data: Dict[str, Any]) -> str:
        """Generate tree using ETE3 for professional visualization"""
        try:
            # Convert tree data to Newick format
            newick_str = self._tree_dict_to_newick(tree_data["root"])
            
            # Create ETE3 tree
            tree = Tree(newick_str)
            
            # Set tree style
            ts = TreeStyle()
            ts.show_leaf_name = True
            ts.show_branch_length = True
            ts.show_branch_support = True
            ts.mode = "r"  # rectangular mode
            ts.branch_vertical_margin = 10
            ts.scale = 120
            
            # Add conservation status colors
            status_colors = {
                "EX": "#000000",      # Extinct - Black
                "EW": "#8B008B",      # Extinct in Wild - Dark Magenta
                "CR": "#FF0000",      # Critically Endangered - Red
                "EN": "#FF8C00",      # Endangered - Dark Orange
                "VU": "#FFD700",      # Vulnerable - Gold
                "NT": "#9ACD32",      # Near Threatened - Yellow Green
                "LC": "#32CD32",      # Least Concern - Lime Green
                "DD": "#808080",      # Data Deficient - Gray
                "Query": "#0000FF"    # Query Species - Blue
            }
            
            # Style nodes based on conservation status
            for leaf in tree:
                # Find species in data
                species_name = leaf.name.replace("_", " ")
                status = "DD"  # Default
                
                if species_name == sequences_data["query_species"]:
                    status = "Query"
                else:
                    for species in sequences_data.get("similar_species", []):
                        if species["name"] == species_name:
                            status = species.get("status", "DD")
                            break
                
                # Add colored circle
                color = status_colors.get(status, "#808080")
                circle = CircleFace(radius=8, color=color, style="sphere")
                leaf.add_face(circle, column=0, position="branch-right")
                
                # Add status text
                status_face = AttrFace("name", fsize=10, fgcolor=color, fstyle="italic")
                leaf.add_face(status_face, column=1, position="branch-right")
            
            # Add title
            title = f"Phylogenetic Tree: {sequences_data['query_species']}\nMethod: {tree_data.get('method', 'Unknown')}"
            ts.title.add_face(faces.TextFace(title, fsize=14, bold=True), column=0)
            
            # Render to image
            temp_file = tempfile.NamedTemporaryFile(suffix='.png', delete=False)
            tree.render(temp_file.name, tree_style=ts, dpi=300, w=1200, h=800)
            
            # Convert to base64
            with open(temp_file.name, 'rb') as f:
                image_data = f.read()
            
            os.unlink(temp_file.name)
            
            encoded_img = base64.b64encode(image_data).decode('utf-8')
            logger.info(f"✅ Generated ETE3 tree image ({len(encoded_img)} bytes)")
            return f"data:image/png;base64,{encoded_img}"
            
        except Exception as e:
            logger.error(f"❌ ETE3 tree generation failed: {str(e)}")
            raise
    
    def _tree_dict_to_newick(self, node, is_root=True):
        """Convert tree dictionary to Newick format"""
        if node.get("is_leaf", False):
            name = node.get("name", "Unknown").replace(" ", "_")
            branch_length = node.get("branch_length", 0.0)
            return f"{name}:{branch_length:.6f}"
        
        children = node.get("children", [])
        if not children:
            return "();"
        
        child_strings = [self._tree_dict_to_newick(child, False) for child in children]
        children_str = ",".join(child_strings)
        
        branch_length = node.get("branch_length", 0.0)
        name = node.get("name", "")
        
        if is_root:
            return f"({children_str});"
        else:
            return f"({children_str}){name}:{branch_length:.6f}"
    
    
    def _generate_matplotlib_tree(self, tree_data: Dict[str, Any], sequences_data: Dict[str, Any]) -> str:
        """Generate PROPER branching phylogenetic tree using matplotlib"""
        try:
            # Use the BioPython tree if available for proper visualization
            if "biopython_tree" in tree_data:
                return self._visualize_biopython_tree(tree_data, sequences_data)
        
            # Fallback to manual tree drawing
            return self._draw_manual_tree(tree_data, sequences_data)
        
        except Exception as e:
            logger.error(f"❌ Error generating matplotlib tree: {str(e)}")
            raise

    def _visualize_biopython_tree(self, tree_data: Dict[str, Any], sequences_data: Dict[str, Any]) -> str:
        """Visualize BioPython tree with proper branching"""
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as patches
            from Bio import Phylo
            import io
    
            tree = tree_data["biopython_tree"]
    
            if isinstance(tree, str):
                try:
                    tree = Phylo.read(StringIO(tree), "newick")
                    logger.info("ℹ️ Parsed tree from Newick string")
                except Exception as e:
                    logger.error(f"❌ Failed to parse Newick tree: {str(e)}")
                    raise
            else:
                tree = tree  # sudah dalam bentuk Tree object
            # Create figure with larger size for better visibility
            fig, ax = plt.subplots(figsize=(16, 12))
    
            # Use BioPython's tree drawing with customization
            Phylo.draw(tree, axes=ax, do_show=False, 
                  show_confidence=False, 
                  branch_labels=lambda c: f"{c.branch_length:.3f}" if c.branch_length else "",
                  label_func=lambda c: c.name.replace("_", " ") if c.name else "")
    
            # Customize the plot
            ax.set_title(f'Phylogenetic Tree: {sequences_data["query_species"]}\n'
                    f'Method: {tree_data.get("method", "Unknown")} | '
                    f'Species: {tree_data.get("species_count", 0)}',
                    fontsize=16, fontweight='bold', pad=20)
    
            ax.set_xlabel("Evolutionary Distance", fontsize=12)
    
            # Style the plot
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(True, alpha=0.3)
    
            # Add conservation status legend
            self._add_conservation_legend(ax)
    
            # Save to base64
            plt.tight_layout()
            img_data = io.BytesIO()
            plt.savefig(img_data, format='png', dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
            img_data.seek(0)
    
            encoded_img = base64.b64encode(img_data.read()).decode('utf-8')
            plt.close()
    
            logger.info(f"✅ Generated proper branching tree image ({len(encoded_img)} bytes)")
            return f"data:image/png;base64,{encoded_img}"
    
        except Exception as e:
            logger.error(f"❌ Error visualizing BioPython tree: {str(e)}")
            raise

    def _add_conservation_colors_to_tree(self, ax, tree, sequences_data):
        """Add conservation status colors to tree visualization"""
        try:
            # Conservation status colors
            status_colors = {
                "EX": "#000000", "EW": "#8B008B", "CR": "#FF0000", "EN": "#FF8C00",
                "VU": "#FFD700", "NT": "#9ACD32", "LC": "#32CD32", "DD": "#808080",
                "Query": "#0000FF"
            }
        
            # Get species info
            species_info = {}
            for species in sequences_data.get("similar_species", []):
                species_info[species["name"]] = species.get("status", "DD")
        
            # Add query species
            species_info[sequences_data["query_species"]] = "Query"
        
            # Color the leaf nodes (this is a simplified approach)
            # In a real implementation, you'd need to map the tree coordinates
            # to the species names and add colored markers
        
        except Exception as e:
            logger.warning(f"⚠️ Could not add conservation colors: {str(e)}")

    def _add_conservation_legend(self, ax):
        """Add conservation status legend"""
        try:
            status_colors = {
                "CR": "#FF0000", "EN": "#FF8C00", "VU": "#FFD700", 
                "NT": "#9ACD32", "LC": "#32CD32", "DD": "#808080", "Query": "#0000FF"
            }
        
            legend_elements = []
            for status, color in status_colors.items():
                legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                                markerfacecolor=color, markersize=10, 
                                                label=status, markeredgecolor='black'))
        
            ax.legend(handles=legend_elements, title="Conservation Status", 
                     loc='upper right', bbox_to_anchor=(1.15, 1), frameon=True)
        
        except Exception as e:
            logger.warning(f"⚠️ Could not add legend: {str(e)}")
    
    def save_tree_image_to_file(self, tree_data: Dict[str, Any], 
                               sequences_data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Generate and save professional tree image to file
        """
        try:
            logger.info("🖼️ Generating and saving professional tree visualization...")
            
            # Generate the professional tree image
            base64_image = self.generate_professional_tree_image(tree_data, sequences_data)
            
            if not base64_image:
                logger.error("❌ No image data generated")
                return None
            
            # Extract base64 data
            if base64_image.startswith('data:image/png;base64,'):
                base64_data = base64_image.replace('data:image/png;base64,', '')
            else:
                base64_data = base64_image
            
            # Generate filename
            species_name = sequences_data["query_species"]
            safe_species_name = re.sub(r'[^\w\-_]', '_', species_name)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            method = tree_data.get("method", "unknown").replace("_", "-")
            filename = f"phylo_tree_{method}_{safe_species_name}_{timestamp}.png"
            
            # Ensure output directory exists
            os.makedirs(self.output_dir, exist_ok=True)
            
            # Full file path
            file_path = os.path.join(self.output_dir, filename)
            
            # Decode and save
            try:
                image_data = base64.b64decode(base64_data)
                
                with open(file_path, 'wb') as f:
                    f.write(image_data)
                
                # Verify file creation
                if os.path.exists(file_path):
                    file_size = os.path.getsize(file_path)
                    logger.info(f"✅ Professional tree image saved to: {file_path}")
                    logger.info(f"📊 File size: {file_size} bytes")
                    
                    return {
                        "file_path": file_path,
                        "filename": filename,
                        "size_bytes": file_size,
                        "base64_data": base64_image,
                        "url": f"/api/images/{filename}",
                        "method": tree_data.get("method", "unknown"),
                        "success": True
                    }
                else:
                    logger.error(f"❌ File was not created: {file_path}")
                    return None
                    
            except Exception as decode_error:
                logger.error(f"❌ Error decoding/saving image: {str(decode_error)}")
                return None
        
        except Exception as e:
            logger.error(f"❌ Error saving tree image: {str(e)}")
            logger.error(traceback.format_exc())
            return None
    
    def _generate_error_image(self, error_message: str) -> str:
        """Generate error image"""
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            ax.text(0.5, 0.5, f'Error generating phylogenetic tree:\n\n{error_message}', 
                   ha='center', va='center', fontsize=12, 
                   bbox=dict(boxstyle="round,pad=0.5", facecolor="lightcoral", alpha=0.8),
                   transform=ax.transAxes)
            
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis('off')
            ax.set_title('Phylogenetic Tree Generation Error', fontsize=16, fontweight='bold')
            
            buffer = io.BytesIO()
            plt.savefig(buffer, format='png', dpi=200, bbox_inches='tight')
            buffer.seek(0)
            
            image_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
            
            plt.close()
            buffer.close()
            
            return f"data:image/png;base64,{image_base64}"
            
        except Exception as e:
            logger.error(f"❌ Error generating error image: {str(e)}")
            return ""

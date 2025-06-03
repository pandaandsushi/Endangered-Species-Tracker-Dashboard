"""
FIXED Real Phylogenetic Service - dengan proper error handling dan JSON serialization
"""

import os
import time
import logging
import traceback
import base64
import io
import json
from typing import Dict, List, Any, Optional, Union
import tempfile
import random
import string

# Biopython imports
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RealPhylogeneticServiceFixed:
    """
    FIXED Real Phylogenetic Service with proper error handling and JSON serialization
    """
    
    def __init__(self):
        """Initialize the service"""
        logger.info("✅ FIXED Phylogenetic Service initialized")
    
    def _create_alignment_from_sequences(self, sequences_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Create alignment from sequences with proper error handling
        Returns a JSON-serializable dict instead of Bio.Align.MultipleSeqAlignment
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
            
            # FIXED: Make all sequences the same length
            # Find the minimum length among all sequences
            sequence_lengths = [len(seq) for seq in all_sequences.values()]
            min_length = min(sequence_lengths)
            max_length = max(sequence_lengths)
            
            logger.info(f"📏 Sequence lengths: min={min_length}, max={max_length}")
            
            # If sequences are too different in length, use a reasonable common length
            if max_length > min_length * 3:  # If max is 3x longer than min
                # Use a reasonable middle ground
                target_length = min(max_length, max(min_length * 2, 500))
                logger.info(f"⚖️ Using target length: {target_length} bp")
            else:
                target_length = min_length
                logger.info(f"✂️ Truncating all sequences to: {target_length} bp")
            
            # Create sequence records with same length
            records = []
            aligned_sequences = []
            
            for species_name, sequence in all_sequences.items():
                # Clean sequence (remove non-DNA characters)
                clean_seq = ''.join(c for c in sequence.upper() if c in 'ATCGN')
                
                # Adjust sequence length
                if len(clean_seq) >= target_length:
                    # Truncate if too long
                    adjusted_seq = clean_seq[:target_length]
                else:
                    # Pad with N's if too short
                    adjusted_seq = clean_seq + 'N' * (target_length - len(clean_seq))
                
                # Create SeqRecord
                record = SeqRecord(
                    Seq(adjusted_seq),
                    id=species_name.replace(" ", "_"),  # Remove spaces for compatibility
                    name=species_name.replace(" ", "_"),
                    description=f"Aligned: {species_name}"
                )
                records.append(record)
                
                # Add to JSON-serializable format
                aligned_sequences.append({
                    "id": species_name,
                    "name": species_name,
                    "sequence": adjusted_seq,
                    "length": len(adjusted_seq),
                    "original_length": len(sequence)
                })
            
            # Now create alignment with same-length sequences
            try:
                alignment = MultipleSeqAlignment(records)
                logger.info(f"✅ Created alignment with {len(records)} sequences of {target_length} bp each")
            except Exception as e:
                logger.error(f"❌ Still failed to create alignment: {str(e)}")
                # Return without BioPython alignment object
                return {
                    "aligned_sequences": aligned_sequences,
                    "alignment_method": "manual_truncation",
                    "num_sequences": len(aligned_sequences),
                    "target_length": target_length,
                    "note": "Created without BioPython alignment due to compatibility issues"
                }
            
            return {
                "aligned_sequences": aligned_sequences,
                "alignment_method": "truncation_and_padding",
                "num_sequences": len(aligned_sequences),
                "target_length": target_length,
                "biopython_alignment": alignment  # Keep for distance calculation
            }
            
        except Exception as e:
            logger.error(f"❌ Error creating alignment: {str(e)}")
            logger.error(traceback.format_exc())
            
            # Return error result
            return {
                "error": f"Alignment creation failed: {str(e)}",
                "aligned_sequences": [],
                "alignment_method": "failed"
            }
    
    def build_real_tree_fixed(self, sequences_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Build phylogenetic tree from real sequence data with proper error handling
        Returns a JSON-serializable dict
        """
        try:
            logger.info("🌳 Building real phylogenetic tree...")
            
            # Step 1: Create alignment
            alignment_data = self._create_alignment_from_sequences(sequences_data)
            
            if "error" in alignment_data:
                logger.error(f"❌ Alignment failed: {alignment_data['error']}")
                return {
                    "error": alignment_data["error"],
                    "method": "failed",
                    "tree_data": {}
                }
            
            # Step 2: Calculate distance matrix and build tree
            try:
                # If we have BioPython alignment, use it for distance calculation
                if "biopython_alignment" in alignment_data:
                    alignment = alignment_data["biopython_alignment"]
                    
                    # Calculate distance matrix
                    calculator = DistanceCalculator('identity')
                    dm = calculator.get_distance(alignment)
                    
                    # Build tree
                    constructor = DistanceTreeConstructor()
                    tree = constructor.upgma(dm)
                    
                    # Convert tree to JSON-serializable format
                    tree_dict = self._convert_tree_to_dict(tree.root)
                    
                    tree_data = {
                        "method": "UPGMA_with_biopython",
                        "root": tree_dict,
                        "species_count": len(alignment_data["aligned_sequences"]),
                        "query_species": sequences_data["query_species"]
                    }
                    
                else:
                    # Fallback: create simple distance-based tree
                    tree_data = self._create_simple_distance_tree(alignment_data, sequences_data)
                
            except Exception as e:
                logger.error(f"❌ Error building tree with BioPython: {str(e)}")
                # Fallback to simple tree
                tree_data = self._create_simple_distance_tree(alignment_data, sequences_data)
            
            logger.info(f"✅ Built phylogenetic tree with {tree_data.get('species_count', 0)} species")
            return tree_data
            
        except Exception as e:
            logger.error(f"❌ Error building tree: {str(e)}")
            logger.error(traceback.format_exc())
            
            # Return error result
            return {
                "error": f"Tree building failed: {str(e)}",
                "method": "failed",
                "tree_data": {}
            }
    
    def _convert_tree_to_dict(self, clade, parent_id=None):
        """Convert tree clade to dict recursively"""
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
                    # Get distance to query
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
    
    def generate_real_tree_image_fixed(self, tree_data: Dict[str, Any], 
                                      sequences_data: Dict[str, Any]) -> Optional[str]:
        """
        Generate tree image with proper error handling
        Returns base64 encoded PNG image
        """
        try:
            logger.info("🖼️ Generating tree visualization...")
            
            if "error" in tree_data:
                logger.error(f"❌ Cannot generate image: {tree_data['error']}")
                return self._generate_error_image(tree_data["error"])
            
            # Create a simple visualization based on tree_data
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Extract species names and conservation status
            species_names = [sequences_data["query_species"]]
            statuses = ["Query"]
            similarities = [1.0]  # Query has 100% similarity to itself
            
            for species in sequences_data.get("similar_species", []):
                species_names.append(species["name"])
                statuses.append(species.get("status", "DD"))
                similarities.append(species.get("similarity", 0.0))
            
            # Map conservation status to colors
            status_colors = {
                "EX": "black",      # Extinct
                "EW": "purple",     # Extinct in the Wild
                "CR": "red",        # Critically Endangered
                "EN": "orange",     # Endangered
                "VU": "yellow",     # Vulnerable
                "NT": "yellowgreen", # Near Threatened
                "LC": "green",      # Least Concern
                "DD": "gray",       # Data Deficient
                "Query": "blue"     # Query Species
            }
            
            # Create colors list
            colors = [status_colors.get(status, "gray") for status in statuses]
            
            # Create a simple tree visualization
            y_positions = list(range(len(species_names)))
            
            # Draw branches based on similarity (closer = more similar)
            for i, (species, y, similarity) in enumerate(zip(species_names, y_positions, similarities)):
                # Branch length inversely proportional to similarity
                branch_length = 1.0 - similarity
                
                # Draw horizontal line
                ax.plot([0, branch_length], [y, y], 'k-', linewidth=2)
                
                # Draw vertical connector to common ancestor
                if i > 0:
                    ax.plot([branch_length, 1.0], [y, y], 'k-', linewidth=1)
            
            # Draw vertical line connecting all species
            if len(y_positions) > 1:
                ax.plot([1.0, 1.0], [min(y_positions), max(y_positions)], 'k-', linewidth=2)
            
            # Draw species names and status indicators
            for i, (species, y, color, similarity) in enumerate(zip(species_names, y_positions, colors, similarities)):
                # Species name
                ax.text(1.1, y, f"{species} ({similarity:.1%})", fontsize=10, va='center')
                
                # Status indicator
                ax.scatter(1.0 - similarity, y, c=color, s=150, zorder=10, edgecolors='black')
            
            # Add legend for conservation status
            legend_elements = []
            unique_statuses = list(set(statuses))
            for status in unique_statuses:
                color = status_colors.get(status, "gray")
                legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                                markerfacecolor=color, markersize=10, label=status))
            
            ax.legend(handles=legend_elements, title="Conservation Status", 
                     loc='upper right', bbox_to_anchor=(1.3, 1))
            
            # Set labels and title
            ax.set_title(f'Phylogenetic Tree: {sequences_data["query_species"]}\n'
                        f'Method: {tree_data.get("method", "Unknown")} | '
                        f'Species: {tree_data.get("species_count", 0)}',
                        fontsize=14, fontweight='bold')
            ax.set_xlabel("Evolutionary Distance (based on sequence similarity)")
            ax.set_ylabel("Species")
            
            # Set axis limits
            ax.set_xlim(-0.1, 2.0)
            ax.set_ylim(-0.5, len(species_names) - 0.5)
            
            # Style axes
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_yticks([])
            
            # Save to base64
            plt.tight_layout()
            img_data = io.BytesIO()
            plt.savefig(img_data, format='png', dpi=300, bbox_inches='tight')
            img_data.seek(0)
            
            # Encode to base64
            encoded_img = base64.b64encode(img_data.read()).decode('utf-8')
            plt.close()
            
            logger.info(f"✅ Generated tree image ({len(encoded_img)} bytes)")
            return f"data:image/png;base64,{encoded_img}"
            
        except Exception as e:
            logger.error(f"❌ Error generating tree image: {str(e)}")
            logger.error(traceback.format_exc())
            return self._generate_error_image(str(e))
    
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

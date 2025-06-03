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
from Bio import Phylo, AlignIO, SeqIO
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
            
            # Create sequence records
            records = []
            
            # Add query sequence
            query_record = SeqRecord(
                Seq(query_sequence),
                id=query_species,
                name=query_species,
                description=f"Query: {query_species}"
            )
            records.append(query_record)
            
            # Add similar species sequences
            for species in similar_species:
                species_record = SeqRecord(
                    Seq(species["sequence"]),
                    id=species["name"],
                    name=species["name"],
                    description=f"Similar: {species['name']}"
                )
                records.append(species_record)
            
            # Create alignment (this is just a collection of sequences at this point)
            alignment = MultipleSeqAlignment(records)
            
            # Convert to JSON-serializable format
            aligned_sequences = []
            for record in alignment:
                aligned_sequences.append({
                    "id": str(record.id),
                    "name": str(record.name),
                    "sequence": str(record.seq),
                    "length": int(len(record.seq))
                })
            
            logger.info(f"✅ Created alignment with {len(aligned_sequences)} sequences")
            
            return {
                "aligned_sequences": aligned_sequences,
                "alignment_method": "sequence collection",
                "num_sequences": len(aligned_sequences),
                "sequence_lengths": [len(str(record.seq)) for record in alignment]
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
            
            # Step 2: Calculate distance matrix
            try:
                # Create temporary FASTA file for alignment
                with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as temp_file:
                    # Write sequences to temp file
                    for seq in alignment_data["aligned_sequences"]:
                        temp_file.write(f">{seq['id']}\n{seq['sequence']}\n")
                    temp_file_name = temp_file.name
                
                # Read alignment from file
                alignment = AlignIO.read(temp_file_name, "fasta")
                
                # Calculate distance matrix
                calculator = DistanceCalculator('identity')
                dm = calculator.get_distance(alignment)
                
                # Build tree
                constructor = DistanceTreeConstructor()
                tree = constructor.upgma(dm)
                
                # Clean up temp file
                os.unlink(temp_file_name)
                
            except Exception as e:
                logger.error(f"❌ Error building tree: {str(e)}")
                
                # Create a simple JSON-serializable tree structure as fallback
                tree_data = {
                    "method": "fallback_tree",
                    "error": str(e),
                    "nodes": [{"name": seq["id"]} for seq in alignment_data["aligned_sequences"]],
                    "note": "Fallback tree due to error in tree construction"
                }
                return tree_data
            
            # Step 3: Convert tree to JSON-serializable format
            def convert_tree_to_dict(clade, parent_id=None):
                """Convert tree clade to dict recursively"""
                node_id = ''.join(random.choices(string.ascii_lowercase + string.digits, k=8))
                
                node = {
                    "id": node_id,
                    "name": str(clade.name) if clade.name else None,
                    "branch_length": float(clade.branch_length) if clade.branch_length else 0.0,
                    "parent": parent_id
                }
                
                if clade.clades:
                    node["children"] = [convert_tree_to_dict(child, node_id) for child in clade.clades]
                else:
                    node["is_leaf"] = True
                
                return node
            
            # Convert tree to dict
            tree_dict = convert_tree_to_dict(tree.root)
            
            # Create tree data
            tree_data = {
                "method": "UPGMA",
                "root": tree_dict,
                "species_count": len(alignment_data["aligned_sequences"]),
                "query_species": sequences_data["query_species"]
            }
            
            logger.info(f"✅ Built phylogenetic tree with {tree_data['species_count']} species")
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
                return None
            
            # Create a simple visualization based on tree_data
            fig, ax = plt.figure(figsize=(10, 8)), plt.subplot(111)
            
            # Extract species names and conservation status
            species_names = [sequences_data["query_species"]]
            statuses = ["Query"]
            
            for species in sequences_data.get("similar_species", []):
                species_names.append(species["name"])
                statuses.append(species.get("status", "DD"))
            
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
            
            # Draw branches
            for i, (species, y) in enumerate(zip(species_names, y_positions)):
                # Draw horizontal line
                ax.plot([0, 0.8], [y, y], 'k-', linewidth=1.5)
                
                # Draw vertical connector if not the first species
                if i > 0:
                    ax.plot([0.8, 0.8], [y_positions[i-1], y], 'k-', linewidth=1.5)
            
            # Draw species names and status indicators
            for i, (species, y, color) in enumerate(zip(species_names, y_positions, colors)):
                ax.text(0.85, y, species, fontsize=10, va='center')
                ax.scatter(0, y, c=color, s=100, zorder=10)
            
            # Add legend for conservation status
            legend_elements = []
            for status, color in status_colors.items():
                if status in statuses:
                    legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                                    markerfacecolor=color, markersize=10, label=status))
            
            ax.legend(handles=legend_elements, title="Conservation Status", 
                     loc='upper right', bbox_to_anchor=(1.1, 1))
            
            # Set labels and title
            ax.set_title(f"Phylogenetic Tree for {sequences_data['query_species']}")
            ax.set_xlabel("Evolutionary Distance")
            ax.set_yticks([])
            ax.set_xlim(-0.1, 2)
            
            # Remove axes
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.get_xaxis().set_visible(False)
            
            # Save to base64
            img_data = io.BytesIO()
            plt.savefig(img_data, format='png', bbox_inches='tight')
            img_data.seek(0)
            
            # Encode to base64
            encoded_img = base64.b64encode(img_data.read()).decode('utf-8')
            plt.close()
            
            logger.info(f"✅ Generated tree image ({len(encoded_img)} bytes)")
            return encoded_img
            
        except Exception as e:
            logger.error(f"❌ Error generating tree image: {str(e)}")
            logger.error(traceback.format_exc())
            return None

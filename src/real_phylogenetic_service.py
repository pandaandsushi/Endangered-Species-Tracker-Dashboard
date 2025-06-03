
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
        Build phylogenetic tree using Maximum Likelihood or UPGMA
        """
        try:
            logger.info("🌳 Building Maximum Likelihood phylogenetic tree...")
            
            # Step 1: Create alignment
            alignment_data = self._create_alignment_from_sequences(sequences_data)
            
            if "error" in alignment_data:
                logger.error(f"❌ Alignment failed: {alignment_data['error']}")
                return {
                    "error": alignment_data["error"],
                    "method": "failed",
                    "tree_data": {}
                }
            
            alignment = alignment_data.get("biopython_alignment")
            if not alignment:
                logger.error("❌ No alignment available")
                return {
                    "error": "No alignment available",
                    "method": "failed",
                    "tree_data": {}
                }
            
            # Step 2: Try Maximum Likelihood if PhyML is available
            if self.phyml_available:
                try:
                    tree_data = self._build_ml_tree_phyml(alignment, sequences_data)
                    if tree_data and "error" not in tree_data:
                        return tree_data
                except Exception as e:
                    logger.warning(f"⚠️ PhyML failed: {str(e)}, falling back to UPGMA")
            
            # Step 3: Fallback to UPGMA
            try:
                logger.info("🌳 Using UPGMA method...")
                calculator = DistanceCalculator('identity')
                dm = calculator.get_distance(alignment)
                
                constructor = DistanceTreeConstructor()
                tree = constructor.upgma(dm)
                
                # Convert tree to JSON-serializable format
                tree_dict = self._convert_tree_to_dict(tree.root)
                
                return {
                    "method": "UPGMA_distance_based",
                    "root": tree_dict,
                    "species_count": len(alignment_data["aligned_sequences"]),
                    "query_species": sequences_data["query_species"],
                    "biopython_tree": tree  # Keep for visualization
                }
                
            except Exception as e:
                logger.error(f"❌ UPGMA also failed: {str(e)}")
                return self._create_simple_distance_tree(alignment_data, sequences_data)
            
        except Exception as e:
            logger.error(f"❌ Error building tree: {str(e)}")
            logger.error(traceback.format_exc())
            
            return {
                "error": f"Tree building failed: {str(e)}",
                "method": "failed",
                "tree_data": {}
            }
    
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
        """Generate tree using improved matplotlib visualization"""
        try:
            # Create figure with better styling
            plt.style.use('seaborn-v0_8')
            fig, ax = plt.subplots(figsize=(14, 10))
            
            # Extract species data
            species_names = [sequences_data["query_species"]]
            statuses = ["Query"]
            similarities = [1.0]
            
            for species in sequences_data.get("similar_species", []):
                species_names.append(species["name"])
                statuses.append(species.get("status", "DD"))
                similarities.append(species.get("similarity", 0.0))
            
            # Conservation status colors
            status_colors = {
                "EX": "#000000", "EW": "#8B008B", "CR": "#FF0000", "EN": "#FF8C00",
                "VU": "#FFD700", "NT": "#9ACD32", "LC": "#32CD32", "DD": "#808080",
                "Query": "#0000FF"
            }
            
            colors = [status_colors.get(status, "#808080") for status in statuses]
            
            # Create dendrogram-style tree
            y_positions = list(range(len(species_names)))
            
            # Draw tree branches based on evolutionary distance
            for i, (species, y, similarity) in enumerate(zip(species_names, y_positions, similarities)):
                # Branch length based on genetic distance
                branch_length = 1.0 - similarity
                
                # Main horizontal branch
                ax.plot([0, branch_length], [y, y], 'k-', linewidth=2, alpha=0.8)
                
                # Vertical connector
                if i > 0:
                    ax.plot([branch_length, 1.0], [y, y], 'k-', linewidth=1, alpha=0.6)
            
            # Main vertical trunk
            if len(y_positions) > 1:
                ax.plot([1.0, 1.0], [min(y_positions), max(y_positions)], 'k-', linewidth=3, alpha=0.8)
            
            # Add species labels and status indicators
            for i, (species, y, color, similarity) in enumerate(zip(species_names, y_positions, colors, similarities)):
                # Species name with similarity
                label = f"{species}\n({similarity:.1%} similarity)"
                ax.text(1.1, y, label, fontsize=11, va='center', ha='left', 
                       style='italic' if i > 0 else 'normal',
                       weight='bold' if i == 0 else 'normal')
                
                # Conservation status circle
                ax.scatter(1.0 - similarity, y, c=color, s=200, zorder=10, 
                          edgecolors='black', linewidth=2, alpha=0.9)
                
                # Branch length annotation
                if branch_length > 0:
                    ax.text(branch_length/2, y+0.1, f"{branch_length:.3f}", 
                           fontsize=8, ha='center', alpha=0.7)
            
            # Add legend
            legend_elements = []
            unique_statuses = list(set(statuses))
            for status in unique_statuses:
                color = status_colors.get(status, "#808080")
                legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                                markerfacecolor=color, markersize=12, 
                                                label=status, markeredgecolor='black'))
            
            ax.legend(handles=legend_elements, title="Conservation Status", 
                     loc='upper right', bbox_to_anchor=(1.4, 1), frameon=True)
            
            # Styling
            ax.set_title(f'Phylogenetic Tree: {sequences_data["query_species"]}\n'
                        f'Method: {tree_data.get("method", "Unknown")} | '
                        f'Species: {tree_data.get("species_count", 0)}',
                        fontsize=16, fontweight='bold', pad=20)
            
            ax.set_xlabel("Evolutionary Distance (genetic divergence)", fontsize=12)
            ax.set_ylabel("Species", fontsize=12)
            
            # Set limits and remove y-axis ticks
            ax.set_xlim(-0.1, 2.2)
            ax.set_ylim(-0.5, len(species_names) - 0.5)
            ax.set_yticks([])
            
            # Style the plot
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.grid(True, alpha=0.3, axis='x')
            
            # Save to base64
            plt.tight_layout()
            img_data = io.BytesIO()
            plt.savefig(img_data, format='png', dpi=300, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            img_data.seek(0)
            
            encoded_img = base64.b64encode(img_data.read()).decode('utf-8')
            plt.close()
            
            logger.info(f"✅ Generated matplotlib tree image ({len(encoded_img)} bytes)")
            return f"data:image/png;base64,{encoded_img}"
            
        except Exception as e:
            logger.error(f"❌ Error generating matplotlib tree: {str(e)}")
            raise
    
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

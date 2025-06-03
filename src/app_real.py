"""
UPDATED Flask Application - dengan Maximum Likelihood phylogenetic trees
"""

from flask import Flask, jsonify, request, send_from_directory
from flask_cors import CORS
import os
import sys
import logging
import traceback

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from ncbi_real_service import NCBIRealServiceFixed
from real_phylogenetic_service import RealPhylogeneticServiceML
from data_manager import DataManager

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)

IMAGE_OUTPUT_DIR = "phylogenetic_trees"
os.makedirs(IMAGE_OUTPUT_DIR, exist_ok=True)

def initialize_services():
    """Initialize all services"""
    services = {}
    
    ncbi_email = os.getenv('NCBI_EMAIL', 'raffaelsiahaan@gmail.com')
    
    try:
        services['ncbi'] = NCBIRealServiceFixed(email=ncbi_email)
        logger.info("✅ NCBI Service initialized")
    except Exception as e:
        logger.error(f"❌ NCBI Service failed: {e}")
        services['ncbi'] = None
    
    try:
        services['phylogenetic'] = RealPhylogeneticServiceML(output_dir=IMAGE_OUTPUT_DIR)
        logger.info("✅ ML Phylogenetic Service initialized")
    except Exception as e:
        logger.error(f"❌ Phylogenetic Service failed: {e}")
        services['phylogenetic'] = None
    
    try:
        services['data_manager'] = DataManager()
        logger.info("✅ Data Manager initialized")
    except Exception as e:
        logger.error(f"❌ Data Manager failed: {e}")
        services['data_manager'] = None
    
    return services

# Initialize services
services = initialize_services()
ncbi_service = services['ncbi']
phylogenetic_service = services['phylogenetic']
data_manager = services['data_manager']

@app.route('/api/images/<filename>')
def serve_image(filename):
    """Serve phylogenetic tree images"""
    try:
        return send_from_directory(IMAGE_OUTPUT_DIR, filename)
    except Exception as e:
        logger.error(f"Error serving image {filename}: {str(e)}")
        return jsonify({"error": "Image not found"}), 404

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check"""
    services_status = {
        "ncbi_service": ncbi_service is not None,
        "phylogenetic_ml_service": phylogenetic_service is not None,
        "data_manager": data_manager is not None
    }
    
    return jsonify({
        "status": "healthy" if all(services_status.values()) else "partial",
        "message": "IUCN Red List Phylogenetic Tree API with Maximum Likelihood",
        "services": services_status,
        "features": {
            "maximum_likelihood": phylogenetic_service.phyml_available if phylogenetic_service else False,
            "ete3_visualization": True,  # Will be checked at runtime
            "professional_trees": True
        }
    })

@app.route('/api/search/phylogenetic', methods=['POST'])
def phylogenetic_search():
    """
    Main endpoint for Maximum Likelihood phylogenetic search
    """
    if not all([ncbi_service, phylogenetic_service, data_manager]):
        return jsonify({
            "success": False,
            "error": "Required services not available",
            "services_available": {
                "ncbi": ncbi_service is not None,
                "phylogenetic": phylogenetic_service is not None,
                "data_manager": data_manager is not None
            }
        }), 503
    
    try:
        data = request.get_json()
        
        if not data or not data.get('species_name'):
            return jsonify({
                "success": False,
                "error": "species_name is required"
            }), 400
        
        species_name = data['species_name'].strip()
        gene = data.get('gene', 'COI')
        max_results = min(data.get('max_results', 20), 20)
        min_similarity = data.get('min_similarity', 0.7)
        
        logger.info(f"🔍 ML PHYLOGENETIC SEARCH: {species_name} ({gene})")
        
        # Step 1: NCBI similarity search
        logger.info("Step 1: Performing NCBI BLAST search...")
        ncbi_results = ncbi_service.search_similar_species_real(
            species_name, gene, max_results, min_similarity
        )
        
        if ncbi_results["total_found"] == 0:
            return jsonify({
                "success": False,
                "error": f"No similar species found for {species_name}",
                "suggestion": "Try lowering the min_similarity parameter"
            })
        
        # Step 2: Build Maximum Likelihood phylogenetic tree
        logger.info("Step 2: Building Maximum Likelihood phylogenetic tree...")
        tree_data = phylogenetic_service.build_ml_tree(ncbi_results)
        
        if "error" in tree_data:
            return jsonify({
                "success": False,
                "error": tree_data["error"],
                "ncbi_results": ncbi_results
            })
        
        # Step 3: Generate professional tree visualization
        logger.info("Step 3: Generating professional tree visualization...")
        image_result = phylogenetic_service.save_tree_image_to_file(tree_data, ncbi_results)
        
        # Step 4: Format results
        search_results = []
        for species_data in ncbi_results["similar_species"]:
            search_results.append({
                "name": str(species_data["name"]),
                "status": str(species_data.get("status", "DD")),
                "alignment": f"{float(species_data['similarity']):.1%}",
                "metadata": {
                    "accession": str(species_data["accession"]),
                    "bit_score": float(species_data["alignment_score"]),
                    "e_value": f"{float(species_data['e_value']):.2e}",
                    "identity": f"{float(species_data['identity']):.1f}%",
                    "sequence_length": int(len(species_data["sequence"]))
                }
            })
        
        # Step 5: Create response
        response_data = {
            "success": True,
            "query_species": {
                "name": str(species_name),
                "status": "Query Species",
                "accession": str(ncbi_results.get("query_accession", "N/A"))
            },
            "search_results": search_results,
            "phylogenetic_tree": {
                "tree_data": tree_data,
                "tree_image": image_result["base64_data"] if image_result else None,
                "image_info": image_result,
            },
            "metadata": {
                "total_found": int(ncbi_results["total_found"]),
                "search_method": str(ncbi_results["search_method"]),
                "tree_method": str(tree_data.get("method", "unknown")),
                "gene": str(gene),
                "results_summary": {
                    "total_found": int(ncbi_results["total_found"]),
                    "search_method": str(ncbi_results["search_method"]),
                    "tree_method": str(tree_data.get("method", "unknown"))
                }
            }
        }
        
        # Step 6: Save results
        try:
            result_id = data_manager.save_analysis_result(response_data)
            response_data["result_id"] = result_id
        except Exception as e:
            logger.error(f"Error saving results: {str(e)}")
            response_data["result_id"] = "save_failed"
        
        logger.info(f"✅ ML Phylogenetic search completed successfully!")
        return jsonify(response_data)
        
    except Exception as e:
        logger.error(f"❌ ML Phylogenetic search failed: {str(e)}")
        logger.error(traceback.format_exc())
        return jsonify({
            "success": False,
            "error": f"Search failed: {str(e)}"
        }), 500

# Compatibility endpoints
@app.route('/api/search/fixed', methods=['POST'])
@app.route('/api/search/real', methods=['POST'])
def compatibility_search():
    """Compatibility endpoint that forwards to ML search"""
    return phylogenetic_search()

if __name__ == '__main__':
    print("🧬 MAXIMUM LIKELIHOOD PHYLOGENETIC TREE API 🧬")
    print("=" * 60)
    print("✅ Maximum Likelihood tree construction")
    print("✅ Professional ETE3 visualization")
    print("✅ Conservation status integration")
    print("✅ Real NCBI BLAST data")
    print(f"📁 Image output: {IMAGE_OUTPUT_DIR}")
    print("🔍 Main endpoint: POST /api/search/phylogenetic")
    
    app.run(debug=True, host='0.0.0.0', port=5000)

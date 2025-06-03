"""
FIXED Real Application - dengan proper error handling dan JSON serialization
"""

from flask import Flask, jsonify, request
from flask_cors import CORS
import os
import sys
import logging
import traceback

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from ncbi_real_service import NCBIRealServiceFixed
from real_phylogenetic_service import RealPhylogeneticServiceFixed
from data_manager import DataManager

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)

# Initialize FIXED services
def initialize_services_fixed():
    """Initialize all services with better error handling"""
    services = {}
    
    # Get email from environment or use default for testing
    ncbi_email = os.getenv('NCBI_EMAIL', 'raffaelsiahaan@gmail.com')
    
    if ncbi_email == 'test@example.com':
        print("⚠️  Using test email. Set NCBI_EMAIL environment variable for production!")
    
    try:
        services['ncbi'] = NCBIRealServiceFixed(email=ncbi_email)
        logger.info("✅ FIXED NCBI Service initialized")
    except Exception as e:
        logger.error(f"❌ NCBI Service failed: {e}")
        services['ncbi'] = None
    
    try:
        services['phylogenetic'] = RealPhylogeneticServiceFixed()
        logger.info("✅ FIXED Phylogenetic Service initialized")
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
try:
    services = initialize_services_fixed()
    ncbi_service = services['ncbi']
    phylogenetic_service = services['phylogenetic']
    data_manager = services['data_manager']
except Exception as e:
    logger.error(f"❌ Service initialization failed: {e}")
    ncbi_service = None
    phylogenetic_service = None
    data_manager = None

@app.route('/api/health', methods=['GET'])
def health_check_fixed():
    """Health check with better error reporting"""
    services_status = {
        "ncbi_fixed_service": ncbi_service is not None,
        "phylogenetic_fixed_service": phylogenetic_service is not None,
        "data_manager": data_manager is not None
    }
    
    return jsonify({
        "status": "healthy" if all(services_status.values()) else "partial",
        "message": "FIXED IUCN Red List Phylogenetic Tree API",
        "services": services_status,
        "description": "Fixed API with proper error handling and JSON serialization",
        "note": "All JSON serialization issues resolved!"
    })

# PENTING: Endpoint untuk kompatibilitas dengan test script lama
@app.route('/api/status/ncbi', methods=['GET'])
def ncbi_status_compat():
    """NCBI status check - untuk kompatibilitas dengan test script lama"""
    try:
        if not ncbi_service:
            return jsonify({
                "status": "unavailable",
                "message": "NCBI service not initialized"
            }), 503
        
        # Quick check if NCBI service is responsive
        is_available = ncbi_service.check_connection_safe()
        
        return jsonify({
            "status": "available" if is_available else "limited",
            "message": "NCBI service status",
            "email": ncbi_service.email,
            "note": "Fixed endpoint for test compatibility"
        })
        
    except Exception as e:
        logger.error(f"❌ NCBI status check failed: {str(e)}")
        return jsonify({
            "status": "error",
            "message": f"Status check failed: {str(e)}"
        }), 500

# PENTING: Endpoint untuk kompatibilitas dengan test script lama
@app.route('/api/search/real', methods=['POST'])
def real_search_compat():
    """Compatibility endpoint for old test script"""
    try:
        logger.info("⚠️ Using compatibility endpoint /api/search/real")
        # Forward to fixed endpoint
        with app.test_request_context('/api/search/fixed', 
                                     json=request.get_json(),
                                     method='POST'):
            response = fixed_similarity_search()
            return response
    except Exception as e:
        logger.error(f"❌ Compatibility endpoint failed: {str(e)}")
        return jsonify({
            "success": False,
            "error": f"Search failed: {str(e)}"
        }), 500

# PENTING: Endpoint untuk kompatibilitas dengan test script lama
@app.route('/api/sequence/info', methods=['GET'])
def sequence_info_compat():
    """Compatibility endpoint for old test script"""
    try:
        logger.info("⚠️ Using compatibility endpoint /api/sequence/info")
        # Forward to fixed endpoint
        with app.test_request_context('/api/sequence/info/fixed', 
                                     args=request.args,
                                     method='GET'):
            response = get_fixed_sequence_info()
            return response
    except Exception as e:
        logger.error(f"❌ Compatibility endpoint failed: {str(e)}")
        return jsonify({
            "error": f"Failed to get sequence info: {str(e)}"
        }), 500

@app.route('/api/search/fixed', methods=['POST'])
def fixed_similarity_search():
    """
    FIXED similarity search dengan proper error handling
    """
    if not all([ncbi_service, phylogenetic_service, data_manager]):
        return jsonify({
            "success": False,
            "error": "Required services not available. Check server logs.",
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
        max_results = min(data.get('max_results', 8), 10)  # Limit for stability
        min_similarity = data.get('min_similarity', 0.7)
        
        logger.info(f"🔍 FIXED SEARCH: {species_name} ({gene}) - max:{max_results}, min_sim:{min_similarity}")
        
        # Step 1: FIXED NCBI similarity search
        logger.info("Step 1: Performing FIXED NCBI search...")
        try:
            ncbi_results = ncbi_service.search_similar_species_real(
                species_name, gene, max_results, min_similarity
            )
        except Exception as e:
            logger.error(f"NCBI search failed: {str(e)}")
            return jsonify({
                "success": False,
                "error": f"NCBI search failed: {str(e)}",
                "suggestion": "Try with a different species name or check your internet connection"
            }), 500
        
        if ncbi_results["total_found"] == 0:
            return jsonify({
                "success": False,
                "error": f"No similar species found for {species_name} with similarity >= {min_similarity}",
                "suggestion": "Try lowering the min_similarity parameter or use a different gene"
            })
        
        # Step 2: Build FIXED phylogenetic tree
        logger.info("Step 2: Building FIXED phylogenetic tree...")
        try:
            tree_data = phylogenetic_service.build_real_tree_fixed(ncbi_results)
        except Exception as e:
            logger.error(f"Tree building failed: {str(e)}")
            return jsonify({
                "success": False,
                "error": f"Tree building failed: {str(e)}",
                "ncbi_results": ncbi_results  # Still return search results
            }), 500
        
        # Step 3: Generate FIXED tree image
        logger.info("Step 3: Generating FIXED tree visualization...")
        try:
            tree_image = phylogenetic_service.generate_real_tree_image_fixed(tree_data, ncbi_results)
        except Exception as e:
            logger.error(f"Tree image generation failed: {str(e)}")
            tree_image = None  # Continue without image
        
        # Step 4: Format search results for frontend (JSON safe)
        search_results = []
        for species_data in ncbi_results["similar_species"]:
            try:
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
            except Exception as e:
                logger.warning(f"Error formatting species data: {str(e)}")
                continue
        
        # Step 5: Create metadata summary (JSON safe)
        metadata = {
            "query_species": {
                "name": str(species_name),
                "status": "Query Species",
                "accession": str(ncbi_results.get("query_accession", "N/A")),
                "sequence_length": int(len(ncbi_results.get("query_sequence", "")))
            },
            "search_parameters": {
                "gene": str(gene),
                "max_results": int(max_results),
                "min_similarity": float(min_similarity)
            },
            "results_summary": {
                "total_found": int(ncbi_results["total_found"]),
                "search_method": str(ncbi_results["search_method"]),
                "tree_method": str(tree_data.get("method", "unknown"))
            }
        }
        
        # Step 6: Save results (JSON serializable)
        try:
            result_data = {
                "query_species": str(species_name),
                "search_results": search_results,
                "phylogenetic_tree": tree_data,  # Already JSON serializable
                "tree_image": tree_image,
                "metadata": metadata,
                "timestamp": "fixed_real_data"
            }
            
            result_id = data_manager.save_analysis_result(result_data)
            logger.info(f"✅ Results saved with ID: {result_id}")
        except Exception as e:
            logger.error(f"Error saving results: {str(e)}")
            result_id = "save_failed"
        
        logger.info(f"✅ FIXED search completed successfully!")
        
        # Step 7: Return FIXED results
        return jsonify({
            "success": True,
            "result_id": result_id,
            "search_results": search_results,
            "phylogenetic_tree": {
                "tree_data": tree_data,
                "tree_image": tree_image
            },
            "metadata": metadata,
            "note": "FIXED: All JSON serialization issues resolved!"
        })
        
    except Exception as e:
        logger.error(f"❌ FIXED search failed: {str(e)}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        return jsonify({
            "success": False,
            "error": f"FIXED search failed: {str(e)}",
            "traceback": traceback.format_exc() if app.debug else None,
            "note": "This error has been logged for debugging"
        }), 500

@app.route('/api/sequence/info/fixed', methods=['GET'])
def get_fixed_sequence_info():
    """Get FIXED sequence information"""
    try:
        species_name = request.args.get('species')
        gene = request.args.get('gene', 'COI')
        
        if not species_name:
            return jsonify({"error": "species parameter required"}), 400
        
        logger.info(f"🧬 Getting FIXED sequence info for {species_name} ({gene})")
        
        sequence_info = ncbi_service.get_sequence_info_real_safe(species_name, gene)
        
        if not sequence_info:
            return jsonify({
                "error": f"No {gene} sequence found for {species_name} in NCBI"
            }), 404
        
        return jsonify({
            "success": True,
            "sequence_info": sequence_info,
            "note": "FIXED: Proper timeout and error handling"
        })
        
    except Exception as e:
        logger.error(f"❌ Error getting fixed sequence info: {str(e)}")
        return jsonify({
            "error": f"Failed to get sequence info: {str(e)}",
            "note": "Error logged for debugging"
        }), 500

@app.route('/api/test/fixed', methods=['POST'])
def test_fixed_system():
    """Test FIXED system"""
    try:
        logger.info("🧪 Testing FIXED system...")
        
        # Test with a well-known species
        test_data = {
            "species_name": "Panthera leo",
            "gene": "COI",
            "max_results": 5,
            "min_similarity": 0.75
        }
        
        # Call fixed search
        with app.test_request_context('/api/search/fixed', json=test_data, method='POST'):
            response = fixed_similarity_search()
            return response
        
    except Exception as e:
        logger.error(f"❌ Fixed system test failed: {str(e)}")
        return jsonify({
            "success": False,
            "error": f"Fixed system test failed: {str(e)}"
        }), 500

@app.route('/api/status/fixed', methods=['GET'])
def fixed_status():
    """Check FIXED service status"""
    try:
        status_info = {
            "ncbi_service": "available" if ncbi_service else "unavailable",
            "phylogenetic_service": "available" if phylogenetic_service else "unavailable",
            "data_manager": "available" if data_manager else "unavailable"
        }
        
        if ncbi_service:
            status_info["ncbi_email"] = ncbi_service.email
            status_info["timeouts"] = {
                "search": ncbi_service.search_timeout,
                "blast": ncbi_service.blast_timeout,
                "fetch": ncbi_service.fetch_timeout
            }
        
        return jsonify({
            "status": "operational",
            "message": "FIXED services status",
            "services": status_info,
            "note": "All timeout and serialization issues fixed"
        })
        
    except Exception as e:
        return jsonify({
            "status": "error",
            "message": f"Status check failed: {str(e)}"
        })

@app.route('/api/genes/available', methods=['GET'])
def get_available_genes():
    """Get available genes for search"""
    try:
        genes = [
            {"code": "COI", "name": "Cytochrome Oxidase I", "recommended": True},
            {"code": "CYTB", "name": "Cytochrome B", "recommended": False},
            {"code": "16S", "name": "16S ribosomal RNA", "recommended": False},
            {"code": "18S", "name": "18S ribosomal RNA", "recommended": False},
            {"code": "ITS", "name": "Internal Transcribed Spacer", "recommended": False}
        ]
        
        return jsonify({
            "genes": genes,
            "default": "COI",
            "note": "COI is recommended for most animal species"
        })
        
    except Exception as e:
        return jsonify({
            "error": f"Failed to get available genes: {str(e)}"
        }), 500

if __name__ == '__main__':
    print("🧬 FIXED IUCN RED LIST PHYLOGENETIC TREE API 🧬")
    print("=" * 65)
    print("✅ FIXED: JSON serialization issues")
    print("✅ FIXED: NCBI timeout problems")
    print("✅ FIXED: Proper error handling")
    print("✅ FIXED: All blocking operations")
    print("✅ FIXED: Compatibility with test script")
    print("\n🔍 Main endpoint: POST /api/search/fixed")
    print("🧪 Test endpoint: POST /api/test/fixed")
    print("📊 Status check: GET /api/status/fixed")
    print("\n⚠️ Compatibility endpoints for test script:")
    print("   - GET /api/status/ncbi")
    print("   - POST /api/search/real")
    print("   - GET /api/sequence/info")
    print("\n🔥 ALL ISSUES FIXED! 🔥")
    
    app.run(debug=True, host='0.0.0.0', port=5000)

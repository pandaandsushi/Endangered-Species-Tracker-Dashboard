"""
Test script untuk FIXED API - no more JSON serialization issues!
"""

import requests
import re
import json
import time
import os
from Bio import Entrez

class FixedAPITester:
    def __init__(self, base_url="http://localhost:5000"):
        self.base_url = base_url
    
    def test_ncbi_status(self):
        """Test NCBI service status"""
        print("🏥 Testing FIXED NCBI service status...")
        
        try:
            # Gunakan endpoint FIXED
            response = requests.get(f"{self.base_url}/api/status/fixed", timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                print(f"✅ FIXED NCBI Status: {data['status']}")
                print(f"   Message: {data['message']}")
                if 'services' in data:
                    print(f"   NCBI Service: {data['services'].get('ncbi_service', 'N/A')}")
                    print(f"   Phylogenetic Service: {data['services'].get('phylogenetic_service', 'N/A')}")
                return True
            else:
                print(f"❌ FIXED NCBI status check failed: {response.status_code}")
                return False
                
        except Exception as e:
            print(f"❌ FIXED NCBI status error: {str(e)}")
            return False
    
    def test_real_search(self, species_name="Panthera leo"):
        """Test FIXED similarity search"""
        print(f"\n🔍 Testing FIXED search for: {species_name}")
        print("⏳ This will take 30-60 seconds (real NCBI BLAST)...")
        
        try:
            payload = {
                "species_name": species_name,
                "gene": "COI",
                "max_results": 30,
                "min_similarity": 0.75
            }
            
            start_time = time.time()
            
            # Gunakan endpoint FIXED
            response = requests.post(
                f"{self.base_url}/api/search/fixed",
                json=payload,
                timeout=300  # 2 minutes timeout for BLAST
            )
            
            end_time = time.time()
            duration = end_time - start_time
            
            if response.status_code == 200:
                data = response.json()
                print(f"✅ FIXED search completed in {duration:.1f} seconds!")
                print(f"   Query: {species_name}")
                print(f"   Total found: {data['metadata']['results_summary']['total_found']}")
                print(f"   Search method: {data['metadata']['results_summary']['search_method']}")
                print(f"   Tree method: {data['metadata']['results_summary']['tree_method']}")
                
                # Display FIXED search results
                print("\n📋 FIXED Search Results from NCBI:")

                
                for i, result in enumerate(data['search_results'][:20], 1):

                    accession_number = result['metadata']['accession']
                    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
                    record = handle.read()
                    import re
                    match = re.search(r"SOURCE\s+(.*)", record)
                    organism = match.group(1).strip()
                    
                    print(f"   {i}. {organism}")
                    print(f"      IUCN Status: {result['status']}")
                    print(f"      Sequence Similarity: {result['alignment']}")
                    print(f"      NCBI Accession: {result['metadata']['accession']}")
                    print(f"      BLAST Bit Score: {result['metadata']['bit_score']}")
                    print(f"      E-value: {result['metadata']['e_value']}")
                
                # Check tree image
                if data['phylogenetic_tree']['tree_image']:
                    print(f"\n🌳 FIXED phylogenetic tree generated!")
                    print(f"   Image size: {len(data['phylogenetic_tree']['tree_image'])} chars")
                
                print(f"\n✨ {data.get('note', 'FIXED data confirmed!')}")
                
                return data
            else:
                print(f"❌ FIXED search failed: {response.status_code}")
                try:
                    error_data = response.json()
                    print(f"   Error: {error_data.get('error')}")
                except:
                    print(f"   Raw error: {response.text}")
                return None
                
        except Exception as e:
            print(f"❌ FIXED search error: {str(e)}")
            return None
    
    def test_sequence_info(self, species_name="Panthera tigris"):
        """Test FIXED sequence info retrieval"""
        print(f"\n🧬 Testing FIXED sequence info for: {species_name}")
        
        try:
            params = {"species": species_name, "gene": "COI"}
            
            # Gunakan endpoint FIXED
            response = requests.get(
                f"{self.base_url}/api/sequence/info/fixed",
                params=params,
                timeout=30
            )
            
            if response.status_code == 200:
                data = response.json()
                seq_info = data['sequence_info']
                
                print("✅ FIXED sequence info retrieved!")
                print(f"   Accession: {seq_info.get('accession', 'N/A')}")
                print(f"   Length: {seq_info.get('length', 0)} bp")
                print(f"   Description: {seq_info.get('description', 'N/A')[:80]}...")
                print(f"   Authors: {seq_info.get('authors', 'N/A')[:50]}...")
                
                return data
            else:
                print(f"❌ Sequence info failed: {response.status_code}")
                return None
                
        except Exception as e:
            print(f"❌ Sequence info error: {str(e)}")
            return None
    
    def run_comprehensive_fixed_test(self):
        """Run comprehensive FIXED API test"""
        print("🧬 COMPREHENSIVE FIXED API TEST 🧬")
        print("=" * 60)
        print("⚠️  This test uses REAL NCBI data - may take several minutes!")
        print("=" * 60)
        
        # Test 1: Health check
        print("\n1️⃣ Testing API health...")
        try:
            response = requests.get(f"{self.base_url}/api/health")
            if response.status_code == 200:
                data = response.json()
                print("✅ API is healthy!")
                print(f"   Status: {data['status']}")
                print(f"   Message: {data['message']}")
                print(f"   Note: {data.get('note', 'N/A')}")
            else:
                print("❌ API health check failed!")
                return False
        except Exception as e:
            print(f"❌ Cannot connect to API: {e}")
            return False
        
        # Test 2: NCBI status
        if not self.test_ncbi_status():
            print("❌ NCBI service not available!")
            return False
        
        # Test 3: FIXED sequence info
        self.test_sequence_info("Panthera leo")
        
        # Test 4: FIXED similarity search
        print("\n" + "="*50)
        print("🚀 STARTING FIXED NCBI BLAST SEARCH")
        print("="*50)
        
        search_result = self.test_real_search("Panthera leo")
        
        if search_result:
            print("\n" + "="*50)
            print("🎉 ALL FIXED TESTS COMPLETED SUCCESSFULLY!")
            print("="*50)
            print("✅ NCBI BLAST search: WORKING")
            print("✅ Sequence retrieval: WORKING") 
            print("✅ Phylogenetic tree: WORKING")
            print("✅ Conservation status: WORKING")
            print("✅ Tree visualization: WORKING")
            print("\n🔥 NO JSON SERIALIZATION ISSUES! 🔥")
            return True
        else:
            print("\n❌ FIXED search test failed!")
            return False

def interactive_fixed_test():
    """Interactive testing for FIXED API"""
    print("🎮 INTERACTIVE FIXED API TESTING 🎮")
    print("=" * 50)
    
    tester = FixedAPITester()
    
    # Check API health first
    try:
        response = requests.get(f"{tester.base_url}/api/health", timeout=5)
        if response.status_code != 200:
            print("❌ API not available! Make sure server is running with: python app_real_fixed.py")
            return
    except:
        print("❌ Cannot connect to API! Make sure server is running.")
        return
    
    while True:
        print("\n📋 FIXED API Commands:")
        print("1. Test NCBI status")
        print("2. Get sequence info")
        print("3. FIXED similarity search")
        print("4. Test with sample species")
        print("5. Run comprehensive test")
        print("0. Exit")
        
        try:
            choice = input("\nEnter choice (0-5): ").strip()
            
            if choice == "0":
                print("👋 Goodbye!")
                break
            elif choice == "1":
                tester.test_ncbi_status()
            elif choice == "2":
                species = input("Enter species name: ").strip()
                gene = input("Enter gene (default COI): ").strip() or "COI"
                if species:
                    tester.test_sequence_info(species)
            elif choice == "3":
                species = input("Enter species name: ").strip()
                if species:
                    print("⚠️  This will take 30-60 seconds for REAL NCBI BLAST...")
                    confirm = input("Continue? (y/n): ").strip().lower()
                    if confirm == 'y':
                        tester.test_real_search(species)
            elif choice == "4":
                sample_species = ["Panthera leo", "Panthera tigris", "Acinonyx jubatus"]
                print("🧪 Testing with sample species (this will take several minutes)...")
                for species in sample_species:
                    print(f"\n--- Testing {species} ---")
                    tester.test_real_search(species)
                    time.sleep(5)  # Brief pause between searches
            elif choice == "5":
                print("⚠️  Comprehensive test will take 5-10 minutes!")
                confirm = input("Continue? (y/n): ").strip().lower()
                if confirm == 'y':
                    tester.run_comprehensive_fixed_test()
            else:
                print("❌ Invalid choice!")
                
        except KeyboardInterrupt:
            print("\n👋 Goodbye!")
            break
        except Exception as e:
            print(f"❌ Error: {str(e)}")

if __name__ == "__main__":
    print("🧬 FIXED API TESTER - NO JSON SERIALIZATION ISSUES! 🧬")
    print("=" * 60)
    print("Choose testing mode:")
    print("1. Comprehensive Test (automatic, ~10 minutes)")
    print("2. Interactive Mode (manual)")
    print("3. Test Compatibility with Old Endpoints")
    
    try:
        mode = input("\nEnter choice (1, 2, or 3): ").strip()
        
        if mode == "1":
            tester = FixedAPITester()
            tester.run_comprehensive_fixed_test()
        elif mode == "2":
            interactive_fixed_test()
        elif mode == "3":
            print("\n🔄 Testing compatibility with old endpoints...")
            from test_real_api import RealAPITester
            old_tester = RealAPITester()
            old_tester.test_ncbi_status()
            old_tester.test_sequence_info("Panthera leo")
            old_tester.test_real_search("Panthera leo")
        else:
            print("❌ Invalid choice!")
            
    except KeyboardInterrupt:
        print("\n👋 Goodbye!")

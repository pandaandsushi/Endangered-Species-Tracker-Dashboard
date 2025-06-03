
import type React from "react"
import { useState } from "react"
import { Search, Loader2, Download, ZoomIn, ZoomOut, AlertCircle, TreePine, Database, Dna } from "lucide-react"

interface Species {
  name: string
  status: string
  alignment: string
  metadata?: {
    accession: string
    bit_score: number
    e_value: string
    identity: string
    sequence_length: number
  }
}

interface SearchResult {
  success: boolean
  query_species: {
    name: string
    status: string
    accession: string
  }
  search_results: Species[]
  phylogenetic_tree: {
    tree_data: any
    tree_image: string | null
    image_info: {
      file_path: string
      filename: string
      url: string
      base64_data: string
      method: string
      success: boolean
    } | null
  }
  metadata: {
    total_found: number
    search_method: string
    tree_method: string
    gene: string
    results_summary: {
      total_found: number
      search_method: string
      tree_method: string
    }
  }
}

export default function Dashboard() {
  const [query, setQuery] = useState("Panthera leo")
  const [searchResults, setSearchResults] = useState<Species[]>([])
  const [selected, setSelected] = useState<Species | null>(null)
  const [loading, setLoading] = useState(false)
  const [treeImage, setTreeImage] = useState<string | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [searchMetadata, setSearchMetadata] = useState<any>(null)
  const [imageInfo, setImageInfo] = useState<any>(null)
  const [gene, setGene] = useState("COI")
  const [minSimilarity, setMinSimilarity] = useState(0.7)

  // API base URL
  const API_BASE_URL = "http://localhost:5000"

  const performSearch = async () => {
    if (!query.trim()) return

    setLoading(true)
    setError(null)
    setTreeImage(null)
    setImageInfo(null)

    try {
      const response = await fetch(`${API_BASE_URL}/api/search/phylogenetic`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          species_name: query,
          gene: gene,
          max_results: 20,
          min_similarity: minSimilarity,
        }),
      })

      const data: SearchResult = await response.json()

      if (data.success) {
        setSearchResults(data.search_results)
        setSearchMetadata(data.metadata)

        // Handle tree image
        if (data.phylogenetic_tree?.tree_image) {
          setTreeImage(data.phylogenetic_tree.tree_image)
        }

        // Set image info
        if (data.phylogenetic_tree?.image_info) {
          setImageInfo(data.phylogenetic_tree.image_info)
        }

        // Set first result as selected
        if (data.search_results.length > 0) {
          setSelected(data.search_results[0])
        }

        console.log("✅ ML Phylogenetic search completed!")
        console.log("🌳 Tree method:", data.metadata.results_summary.tree_method)
      } 
    } catch (err) {
      setError("Failed to connect to server")
      console.error("Search error:", err)
    } finally {
      setLoading(false)
    }
  }

  const handleKeyPress = (e: React.KeyboardEvent) => {
    if (e.key === "Enter") {
      performSearch()
    }
  }

  const getStatusColor = (status: string) => {
    const statusColors: { [key: string]: string } = {
      EX: "text-black",
      EW: "text-purple-600",
      CR: "text-red-600",
      EN: "text-orange-500",
      VU: "text-yellow-500",
      NT: "text-yellow-400",
      LC: "text-green-500",
      DD: "text-gray-500",
      "Query Species": "text-blue-600",
    }
    return statusColors[status] || "text-gray-500"
  }

  const getMethodBadge = (method: string) => {
    if (method.includes("Maximum_Likelihood")) {
      return <span className="bg-green-100 text-green-800 text-xs px-2 py-1 rounded-full">ML</span>
    } else if (method.includes("UPGMA")) {
      return <span className="bg-blue-100 text-blue-800 text-xs px-2 py-1 rounded-full">UPGMA</span>
    } else {
      return <span className="bg-gray-100 text-gray-800 text-xs px-2 py-1 rounded-full">Distance</span>
    }
  }

  const downloadImage = () => {
    if (treeImage) {
      const link = document.createElement("a")
      link.href = treeImage
      link.download = `phylogenetic_tree_${query.replace(" ", "_")}.png`
      link.click()
    }
  }

  return (
    <div className="min-h-screen bg-gradient-to-b from-pink-100 to-pink-200 p-6 text-gray-900">
      <div className="max-w-7xl mx-auto">
        <h1 className="text-4xl font-bold mb-2 flex items-center gap-3">
          <TreePine className="w-8 h-8 text-green-600" />
          IUCN Red List Phylogenetic Tree
        </h1>
        <p className="text-gray-600 mb-6">Maximum Likelihood phylogenetic analysis with real NCBI data</p>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Search Section */}
          <div className="lg:col-span-1">
            <div className="bg-white rounded-lg shadow-lg p-6">
              <h2 className="font-semibold mb-4 flex items-center gap-2">
                <Search className="w-5 h-5" />
                Search Parameters
              </h2>

              <div className="space-y-4">
                <div>
                  <label className="block font-medium mb-1">Species Name</label>
                  <input
                    placeholder="e.g., Panthera leo"
                    type="text"
                    value={query}
                    onChange={(e) => setQuery(e.target.value)}
                    onKeyPress={handleKeyPress}
                    className="w-full border rounded-md px-3 py-2 focus:outline-none focus:ring-2 focus:ring-pink-300"
                  />
                </div>

                <div>
                  <label className="block font-medium mb-1">Gene</label>
                  <select
                    value={gene}
                    onChange={(e) => setGene(e.target.value)}
                    className="w-full border rounded-md px-3 py-2 focus:outline-none focus:ring-2 focus:ring-pink-300"
                  >
                    <option value="COI">COI (Cytochrome Oxidase I)</option>
                    <option value="CYTB">CYTB (Cytochrome B)</option>
                    <option value="16S">16S rRNA</option>
                    <option value="18S">18S rRNA</option>
                  </select>
                </div>

                <div>
                  <label className="block font-medium mb-1">Min Similarity: {(minSimilarity * 100).toFixed(0)}%</label>
                  <input
                    type="range"
                    min="0.5"
                    max="0.95"
                    step="0.05"
                    value={minSimilarity}
                    onChange={(e) => setMinSimilarity(Number.parseFloat(e.target.value))}
                    className="w-full"
                  />
                </div>

                <button
                  onClick={performSearch}
                  disabled={loading}
                  className="w-full bg-pink-500 text-white py-3 rounded-md hover:bg-pink-600 disabled:opacity-50 flex items-center justify-center gap-2 font-medium"
                >
                  {loading ? (
                    <>
                      <Loader2 className="w-5 h-5 animate-spin" />
                      Analyzing...
                    </>
                  ) : (
                    <>
                      <Dna className="w-5 h-5" />
                      Run Analysis
                    </>
                  )}
                </button>
              </div>

              {error && (
                <div className="mt-4 bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded flex items-center gap-2">
                  <AlertCircle className="w-4 h-4" />
                  <span className="text-sm">{error}</span>
                </div>
              )}

              {searchMetadata && (
                <div className="mt-4 bg-blue-50 border border-blue-200 rounded-lg p-4">
                  <h3 className="font-medium text-blue-800 mb-2 flex items-center gap-2">
                    <Database className="w-4 h-4" />
                    Analysis Results
                  </h3>
                  <div className="text-sm space-y-1">
                    <p>
                      <strong>Species Found:</strong> {searchMetadata.results_summary.total_found}
                    </p>
                    <p>
                      <strong>Gene:</strong> {searchMetadata.gene}
                    </p>
                    <p>
                      <strong>Search Method:</strong> {searchMetadata.results_summary.search_method}
                    </p>
                    <p className="flex items-center gap-2">
                      <strong>Tree Method:</strong>
                      {getMethodBadge(searchMetadata.results_summary.tree_method)}
                      {searchMetadata.results_summary.tree_method}
                    </p>
                  </div>
                </div>
              )}

              {imageInfo && (
                <div className="mt-4 bg-green-50 border border-green-200 rounded-lg p-4">
                  <h3 className="font-medium text-green-800 mb-2">🖼️ Image Generation</h3>
                  <div className="text-sm space-y-1">
                    <p>
                      <strong>Method:</strong> {imageInfo.method}
                    </p>
                    <p>
                      <strong>File:</strong> {imageInfo.filename}
                    </p>
                    <p>
                      <strong>Size:</strong> {(imageInfo.size_bytes / 1024).toFixed(1)} KB
                    </p>
                    <p>
                      <strong>Status:</strong> {imageInfo.success ? "✅ Success" : "❌ Failed"}
                    </p>
                  </div>
                </div>
              )}
            </div>

            {/* Search Results */}
            {searchResults.length > 0 && (
              <div className="mt-6 bg-white rounded-lg shadow-lg p-6">
                <h2 className="font-semibold mb-4">Similar Species</h2>
                <div className="space-y-2 max-h-96 overflow-y-auto">
                  {searchResults.map((species, i) => (
                    <div
                      key={i}
                      className={`p-3 rounded-lg cursor-pointer transition-colors ${
                        selected?.name === species.name
                          ? "bg-pink-100 border-2 border-pink-300"
                          : "bg-gray-50 hover:bg-gray-100 border-2 border-transparent"
                      }`}
                      onClick={() => setSelected(species)}
                    >
                      <div className="flex justify-between items-start">
                        <div>
                          <p className="font-medium italic">{species.name}</p>
                          <p className={`text-sm ${getStatusColor(species.status)}`}>{species.status}</p>
                        </div>
                        <div className="text-right">
                          <p className="text-sm font-medium">{species.alignment}</p>
                          <p className="text-xs text-gray-500">similarity</p>
                        </div>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>

          {/* Phylogenetic Tree Section */}
          <div className="lg:col-span-2">
            <div className="bg-white rounded-lg shadow-lg p-6">
              <div className="flex justify-between items-center mb-4">
                <h2 className="font-semibold flex items-center gap-2">
                  <TreePine className="w-5 h-5" />
                  Phylogenetic Tree
                </h2>
                {treeImage && (
                  <button
                    onClick={downloadImage}
                    className="flex items-center gap-2 text-sm bg-green-500 text-white px-4 py-2 rounded-lg hover:bg-green-600"
                  >
                    <Download className="w-4 h-4" />
                    Download PNG
                  </button>
                )}
              </div>

              <div className="bg-gray-50 rounded-lg p-4 flex justify-center items-center min-h-[500px]">
                {loading ? (
                  <div className="flex flex-col items-center">
                    <Loader2 className="w-12 h-12 animate-spin mb-4 text-pink-500" />
                    <span className="text-lg font-medium">Generating phylogenetic tree...</span>
                    <span className="text-sm text-gray-500 mt-2">
                      Running NCBI BLAST and Maximum Likelihood analysis
                    </span>
                    <span className="text-xs text-gray-400 mt-1">This may take 30-90 seconds</span>
                  </div>
                ) : treeImage ? (
                  <div className="w-full">
                    <img
                      src={treeImage || "/placeholder.svg"}
                      alt="Phylogenetic Tree"
                      className="max-w-full h-auto rounded-lg shadow-sm"
                      onError={(e) => {
                        console.error("Image failed to load:", e)
                        setError("Failed to load tree image")
                      }}
                    />
                  </div>
                ) : (
                  <div className="text-center text-gray-500">
                    <TreePine className="w-16 h-16 mx-auto mb-4 text-gray-300" />
                    <p className="text-lg">No tree generated yet</p>
                    <p className="text-sm">Enter a species name and click "Run Analysis"</p>
                  </div>
                )}
              </div>

              {treeImage && (
                <div className="flex justify-between items-center mt-4">
                  <div className="text-xs space-y-1">
                    <p className="font-semibold mb-2">Conservation Status Legend:</p>
                    <div className="grid grid-cols-2 gap-1">
                      <p>
                        <span className="text-red-600 font-bold">●</span> CR (Critically Endangered)
                      </p>
                      <p>
                        <span className="text-orange-500 font-bold">●</span> EN (Endangered)
                      </p>
                      <p>
                        <span className="text-yellow-500 font-bold">●</span> VU (Vulnerable)
                      </p>
                      <p>
                        <span className="text-yellow-400 font-bold">●</span> NT (Near Threatened)
                      </p>
                      <p>
                        <span className="text-green-500 font-bold">●</span> LC (Least Concern)
                      </p>
                      <p>
                        <span className="text-gray-500 font-bold">●</span> DD (Data Deficient)
                      </p>
                    </div>
                  </div>

                  <div className="flex gap-2">
                    <button className="text-sm bg-gray-200 px-3 py-1 rounded hover:bg-gray-300 flex items-center gap-1">
                      <ZoomIn className="w-4 h-4" />
                      Zoom In
                    </button>
                    <button className="text-sm bg-gray-200 px-3 py-1 rounded hover:bg-gray-300 flex items-center gap-1">
                      <ZoomOut className="w-4 h-4" />
                      Zoom Out
                    </button>
                  </div>
                </div>
              )}
            </div>

            {/* Selected Species Details */}
            {selected && (
              <div className="mt-6 bg-white rounded-lg shadow-lg p-6">
                <h2 className="font-semibold mb-4">Selected Species Details</h2>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                  <div>
                    <h3 className="font-medium text-gray-700 mb-2">Basic Information</h3>
                    <div className="space-y-2">
                      <p>
                        <span className="font-medium">Scientific Name:</span>
                        <span className="italic ml-2">{selected.name}</span>
                      </p>
                      <p>
                        <span className="font-medium">Conservation Status:</span>
                        <span className={`ml-2 font-medium ${getStatusColor(selected.status)}`}>{selected.status}</span>
                      </p>
                      <p>
                        <span className="font-medium">Sequence Similarity:</span>
                        <span className="ml-2 font-medium text-blue-600">{selected.alignment}</span>
                      </p>
                    </div>
                  </div>

                  {selected.metadata && (
                    <div>
                      <h3 className="font-medium text-gray-700 mb-2">Sequence Data</h3>
                      <div className="space-y-2">
                        <p>
                          <span className="font-medium">NCBI Accession:</span>
                          <span className="ml-2 font-mono text-sm">{selected.metadata.accession}</span>
                        </p>
                        <p>
                          <span className="font-medium">Identity:</span>
                          <span className="ml-2">{selected.metadata.identity}</span>
                        </p>
                        <p>
                          <span className="font-medium">E-value:</span>
                          <span className="ml-2 font-mono text-sm">{selected.metadata.e_value}</span>
                        </p>
                        <p>
                          <span className="font-medium">Sequence Length:</span>
                          <span className="ml-2">{selected.metadata.sequence_length} bp</span>
                        </p>
                      </div>
                    </div>
                  )}
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  )
}

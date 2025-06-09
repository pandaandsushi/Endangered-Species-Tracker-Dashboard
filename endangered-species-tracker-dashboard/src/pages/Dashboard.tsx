"use client"

import type React from "react"
import { useState } from "react"
import { Search, Loader2, Download, ZoomIn, ZoomOut, AlertCircle, TreePine, Database, Dna, Menu, X } from "lucide-react"
import "../../styles/dashboard.css"

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
    interactive_html: any
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
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false)
  const [interactiveTreeUrl, setInteractiveTreeUrl] = useState<string | null>(null)

  // API base URL
  const API_BASE_URL = "http://localhost:5000"

  const performSearch = async () => {
    if (!query.trim()) return

    setLoading(true)
    setError(null)
    setTreeImage(null)
    setImageInfo(null)
    setInteractiveTreeUrl(null)

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

        if (data.phylogenetic_tree?.tree_image) {
          setTreeImage(data.phylogenetic_tree.tree_image)
        }

        if (data.phylogenetic_tree?.interactive_html) {
          setInteractiveTreeUrl(`${API_BASE_URL}${data.phylogenetic_tree.interactive_html}`)
        }

        if (data.phylogenetic_tree?.image_info) {
          setImageInfo(data.phylogenetic_tree.image_info)
        }

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
    const statusClasses: { [key: string]: string } = {
      EX: "status-ex",
      EW: "status-ew",
      CR: "status-cr",
      EN: "status-en",
      VU: "status-vu",
      NT: "status-nt",
      LC: "status-lc",
      DD: "status-dd",
      "Query Species": "status-query",
    }
    return statusClasses[status] || "status-dd"
  }

  const getMethodBadge = (method: string) => {
    if (method.includes("Maximum_Likelihood")) {
      return <span className="method-badge method-ml">ML</span>
    } else if (method.includes("UPGMA")) {
      return <span className="method-badge method-upgma">UPGMA</span>
    } else {
      return <span className="method-badge method-distance">Distance</span>
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

  const scrollToSection = (sectionId: string) => {
    const element = document.getElementById(sectionId)
    if (element) {
      element.scrollIntoView({ behavior: "smooth" })
    }
    setMobileMenuOpen(false)
  }

  return (
    <div className="dashboard-container">
      {/* Header */}
      <header className="header">
        <div className="header-container">
          <div className="header-content">
            {/* Logo */}
            <div className="logo-section">
              <h1 className="logo-text">THE HOMINIDS</h1>
              <img src="/LOGOSVG.svg" alt="The Hominids Logo" className="logo-img" />
            </div>

            {/* Desktop Navigation */}
            <nav className="nav-desktop">
              <button onClick={() => scrollToSection("overview")} className="nav-button">
                Overview
              </button>
              <button onClick={() => scrollToSection("phylogenetic-tree")} className="nav-button">
                Phylogenetic Tree
              </button>
            </nav>

            {/* Mobile menu button */}
            <button onClick={() => setMobileMenuOpen(!mobileMenuOpen)} className="mobile-menu-button">
              {mobileMenuOpen ? <X size={24} /> : <Menu size={24} />}
            </button>
          </div>

          {/* Mobile Navigation */}
          <div className={`nav-mobile ${mobileMenuOpen ? "open" : ""}`}>
            <button onClick={() => scrollToSection("overview")} className="nav-button">
              Overview
            </button>
            <button onClick={() => scrollToSection("phylogenetic-tree")} className="nav-button">
              Phylogenetic Tree
            </button>
          </div>
        </div>
      </header>

      {/* Hero Section */}
      <section className="hero-section">
        <div className="hero-container">
          <img
            src="../../public/BGConservation.png"
            alt="Prioritizing Conservation Through the Lens of Evolution"
            className="hero-image"
          />
        </div>
      </section>

      {/* Overview Section */}
      <section id="overview" className="overview-section">
        <div className="overview-container">
          <div className="section-title">
            <h2>
              <span className="decorative-tilde">~</span>
              Overview
              <span className="decorative-tilde">~</span>
            </h2>
          </div>

          <div className="overview-content">
            <p className="overview-text">
              The Earth is currently facing a biodiversity crisis, with over one million species of plants and animals
              estimated to be at risk of extinction in the coming decades, according to the Intergovernmental
              Science-Policy Platform on Biodiversity and Ecosystem Services (IPBES). Species extinction is happening at
              rates 100–1000 times faster than natural background rates, driven by habitat destruction, climate change,
              pollution, and human exploitation.
            </p>

            <p className="overview-text">
              Conservation efforts often focus on protecting individual species—but what if we could also understand how
              these species are genetically connected? By examining their evolutionary history, we can make more
              informed decisions about which species are most irreplaceable, and which ones might benefit from genetic
              rescue strategies like assisted gene flow.
            </p>
          </div>

          {/* Why This Project Matters */}
          <div className="project-matters">
            <div className="decorative-circle circle-top-left"></div>
            <div className="decorative-circle circle-top-right"></div>
            <div className="decorative-circle circle-bottom-left"></div>
            <div className="decorative-circle circle-bottom-right"></div>

            <h3>Why This Project Matters</h3>

            <div className="project-content">
              <div>
                <p className="project-description">
                  This project aims to help researchers, conservationists, and the public visualize the evolutionary
                  relationships among endangered species using real biological data. It combines genetic information
                  from public databases with conservation status data to build phylogenetic trees—diagrams that show how
                  species are related through evolution.
                </p>
                
              </div>

              <div className="project-features">
                <p className="project-description bold">With this platform, you can:</p>
                <div className="feature-item">
                  <span className="feature-bullet">•</span>
                  <span>Explore species listed as threatened, endangered, or critically endangered</span>
                </div>
                <div className="feature-item">
                  <span className="feature-bullet">•</span>
                  <span>See how closely related different species are</span>
                </div>
                <div className="feature-item">
                  <span className="feature-bullet">•</span>
                  <span>Discover which species are evolutionarily unique</span>
                </div>
                <div className="feature-item">
                  <span className="feature-bullet">•</span>
                  <span>Support data-driven decisions in biodiversity conservation</span>
                </div>
              </div>
            </div>
          </div>

          {/* How It Works */}
          <div className="how-it-works">
            <h3>
              How It Works
              <span className="decorative-tilde">~</span>
            </h3>

            <div className="steps-grid">
              <div className="step-card">
                <div className="step-number">1</div>
                <h4 className="step-title">Fetch Threatened Species List</h4>
                <p className="step-description">using IUCN Red List API</p>
              </div>

              <div className="step-card">
                <div className="step-number">2</div>
                <h4 className="step-title">Get DNA Sequences from NCBI GenBank</h4>
                <p className="step-description">Retrieve genetic data for analysis</p>
              </div>

              <div className="step-card">
                <div className="step-number">3</div>
                <h4 className="step-title">Align DNA and Build Phylogenetic Tree</h4>
                <p className="step-description">Create evolutionary relationship diagrams</p>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Phylogenetic Tree Section */}
      <section id="phylogenetic-tree" className="tree-section">
        <div className="tree-container">
          <div className="tree-title">
            <h2>
              <TreePine size={32} />
              IUCN Red List Phylogenetic Tree
            </h2>
            <p className="tree-subtitle">Maximum Likelihood phylogenetic analysis with real NCBI data</p>
          </div>

          <div className="tree-grid">
            {/* Search Section */}
            <div>
              <div className="search-panel">
                <h3 className="panel-title">
                  <Search size={20} />
                  Search Parameters
                </h3>

                <div className="form-group">
                  <label className="form-label">Species Name</label>
                  <input
                    placeholder="e.g., Panthera leo"
                    type="text"
                    value={query}
                    onChange={(e) => setQuery(e.target.value)}
                    onKeyPress={handleKeyPress}
                    className="form-input"
                  />
                </div>

                <div className="form-group">
                  <label className="form-label">Gene</label>
                  <select value={gene} onChange={(e) => setGene(e.target.value)} className="form-select">
                    <option value="COI">COI (Cytochrome Oxidase I)</option>
                    <option value="CYTB">CYTB (Cytochrome B)</option>
                    <option value="16S">16S rRNA</option>
                    <option value="18S">18S rRNA</option>
                  </select>
                </div>

                <div className="form-group">
                  <label className="form-label">Min Similarity: {(minSimilarity * 100).toFixed(0)}%</label>
                  <input
                    type="range"
                    min="0.5"
                    max="0.95"
                    step="0.05"
                    value={minSimilarity}
                    onChange={(e) => setMinSimilarity(Number.parseFloat(e.target.value))}
                    className="form-range"
                  />
                </div>

                <button onClick={performSearch} disabled={loading} className="search-button">
                  {loading ? (
                    <>
                      <Loader2 size={20} className="loading-spinner" />
                      Analyzing...
                    </>
                  ) : (
                    <>
                      <Dna size={20} />
                      Run Analysis
                    </>
                  )}
                </button>

                {error && (
                  <div className="alert alert-error">
                    <AlertCircle size={16} />
                    <span>{error}</span>
                  </div>
                )}

                {searchMetadata && (
                  <div className="alert alert-info">
                    <div>
                      <h4 style={{ display: "flex", alignItems: "center", gap: "0.5rem", marginBottom: "0.5rem" }}>
                        <Database size={16} />
                        Analysis Results
                      </h4>
                      <div style={{ fontSize: "0.875rem" }}>
                        <p>
                          <strong>Species Found:</strong> {searchMetadata.results_summary.total_found}
                        </p>
                        <p>
                          <strong>Gene:</strong> {searchMetadata.gene}
                        </p>
                        <p>
                          <strong>Search Method:</strong> {searchMetadata.results_summary.search_method}
                        </p>
                        <p style={{ display: "flex", alignItems: "center", gap: "0.5rem" }}>
                          <strong>Tree Method:</strong>
                          {getMethodBadge(searchMetadata.results_summary.tree_method)}
                          {searchMetadata.results_summary.tree_method}
                        </p>
                      </div>
                    </div>
                  </div>
                )}

                {imageInfo && (
                  <div className="alert alert-success">
                    <div>
                      <h4 style={{ marginBottom: "0.5rem" }}>🖼️ Image Generation</h4>
                      <div style={{ fontSize: "0.875rem" }}>
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
                  </div>
                )}
              </div>

              {/* Search Results */}
              {searchResults.length > 0 && (
                <div className="results-panel">
                  <h3 className="panel-title">Similar Species</h3>
                  <div className="results-list">
                    {searchResults.map((species, i) => (
                      <div
                        key={i}
                        className={`species-item ${selected?.name === species.name ? "selected" : ""}`}
                        onClick={() => setSelected(species)}
                      >
                        <div className="species-info">
                          <div>
                            <p className="species-name">{species.name}</p>
                            <p className={`species-status ${getStatusColor(species.status)}`}>{species.status}</p>
                          </div>
                          <div className="species-alignment">
                            <p className="alignment-value">{species.alignment}</p>
                            <p className="alignment-label">similarity</p>
                          </div>
                        </div>
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </div>

            {/* Phylogenetic Tree Display */}
            <div>
              <div className="tree-display">
                <div className="tree-header">
                  <h3 className="panel-title">
                    <TreePine size={20} />
                    Phylogenetic Tree {interactiveTreeUrl && !treeImage ? "(Interactive)" : "(Static)"}
                  </h3>
                  <div style={{ display: "flex", gap: "0.5rem" }}>
                    {interactiveTreeUrl && (
                      <button onClick={() => window.open(interactiveTreeUrl, "_blank")} className="download-button">
                        <ZoomIn size={16} />
                        Open in New Tab
                      </button>
                    )}
                    {treeImage && (
                      <button onClick={downloadImage} className="download-button">
                        <Download size={16} />
                        Download PNG
                      </button>
                    )}
                  </div>
                </div>

                <div className="tree-container-display">
                  {loading ? (
                    <div className="loading-state">
                      <Loader2 size={48} className="loading-spinner" />
                      <span className="loading-title">Generating phylogenetic tree...</span>
                      <span className="loading-subtitle">Running NCBI BLAST and Maximum Likelihood analysis</span>
                      <span className="loading-note">This may take 30-90 seconds</span>
                    </div>
                  ) : interactiveTreeUrl ? (
                    <div style={{ width: "100%", height: "600px" }}>
                      <iframe
                        src={interactiveTreeUrl}
                        style={{
                          width: "100%",
                          height: "100%",
                          border: "1px solid #e2e8f0",
                          borderRadius: "8px",
                        }}
                        title="Interactive Phylogenetic Tree"
                        onError={(e) => {
                          console.error("Interactive tree failed to load:", e)
                          setError("Failed to load interactive tree")
                        }}
                      />
                    </div>
                  ) : treeImage ? (
                    <div style={{ width: "100%" }}>
                      <img
                        src={treeImage || "/placeholder.svg"}
                        alt="Phylogenetic Tree"
                        className="tree-image"
                        onError={(e) => {
                          console.error("Image failed to load:", e)
                          setError("Failed to load tree image")
                        }}
                      />
                      {interactiveTreeUrl && (
                        <div style={{ marginTop: "1rem", textAlign: "center" }}>
                          <button onClick={() => setInteractiveTreeUrl(interactiveTreeUrl)} className="zoom-button">
                            View Interactive Tree
                          </button>
                        </div>
                      )}
                    </div>
                  ) : (
                    <div className="empty-state">
                      <TreePine size={64} className="empty-icon" />
                      <p className="empty-title">No tree generated yet</p>
                      <p className="empty-subtitle">Enter a species name and click "Run Analysis"</p>
                    </div>
                  )}
                </div>

                {treeImage && (
                  <div className="tree-controls">
                    <div className="legend">
                      <p className="legend-title">Conservation Status Legend:</p>
                      <div className="legend-grid">
                        <p className="legend-item">
                          <span className="legend-dot status-cr">●</span> CR (Critically Endangered)
                        </p>
                        <p className="legend-item">
                          <span className="legend-dot status-en">●</span> EN (Endangered)
                        </p>
                        <p className="legend-item">
                          <span className="legend-dot status-vu">●</span> VU (Vulnerable)
                        </p>
                        <p className="legend-item">
                          <span className="legend-dot status-nt">●</span> NT (Near Threatened)
                        </p>
                        <p className="legend-item">
                          <span className="legend-dot status-lc">●</span> LC (Least Concern)
                        </p>
                        <p className="legend-item">
                          <span className="legend-dot status-dd">●</span> DD (Data Deficient)
                        </p>
                      </div>
                    </div>
                  </div>
                )}
              </div>

              {/* Selected Species Details */}
              {selected && (
                <div className="species-details">
                  <h3 className="panel-title">Selected Species Details</h3>
                  <div className="details-grid">
                    <div className="details-section">
                      <h4>Basic Information</h4>
                      <div className="details-list">
                        <p className="detail-item">
                          <span className="detail-label">Scientific Name:</span>
                          <span className="detail-value">{selected.name}</span>
                        </p>
                        <p className="detail-item">
                          <span className="detail-label">Conservation Status:</span>
                          <span className={`detail-value status ${getStatusColor(selected.status)}`}>
                            {selected.status}
                          </span>
                        </p>
                        <p className="detail-item">
                          <span className="detail-label">Sequence Similarity:</span>
                          <span className="detail-value similarity">{selected.alignment}</span>
                        </p>
                      </div>
                    </div>

                    {selected.metadata && (
                      <div className="details-section">
                        <h4>Sequence Data</h4>
                        <div className="details-list">
                          <p className="detail-item">
                            <span className="detail-label">NCBI Accession:</span>
                            <span className="detail-value mono">{selected.metadata.accession}</span>
                          </p>
                          <p className="detail-item">
                            <span className="detail-label">Identity:</span>
                            <span className="detail-value">{selected.metadata.identity}</span>
                          </p>
                          <p className="detail-item">
                            <span className="detail-label">E-value:</span>
                            <span className="detail-value mono">{selected.metadata.e_value}</span>
                          </p>
                          <p className="detail-item">
                            <span className="detail-label">Sequence Length:</span>
                            <span className="detail-value">{selected.metadata.sequence_length} bp</span>
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
      </section>
    </div>
  )
}

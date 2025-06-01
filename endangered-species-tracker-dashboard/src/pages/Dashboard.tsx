import React, { useState } from 'react';


const speciesList = [
  { name: 'Panthera leo', status: 'Vulnerable', alignment: '98.2%', origin: 'Southern Africa' },
  { name: 'Panthera tigris', status: 'Endangered', alignment: '98.1%', origin: 'Asia' },
  { name: 'Panthera onca', status: 'Near Threatened', alignment: '97.8%', origin: 'South America' },
  { name: 'Panthera pardus', status: 'Vulnerable', alignment: '98.1%', origin: 'Africa & Asia' }
];

export default function Dashboard() {
  const [query, setQuery] = useState('Panthera');
  const [selected, setSelected] = useState(speciesList[0]);

  const filteredSpecies = speciesList.filter(s => s.name.toLowerCase().includes(query.toLowerCase()));

  return (
    <div className="min-h-screen bg-gradient-to-b from-pink-100 to-pink-200 p-6 text-gray-900">
      <h1 className="text-3xl font-bold mb-4">IUCN Red List Phylogenetic Tree</h1>
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        <div>
          <label className="block font-semibold mb-1">Search</label>
          <input
            placeholder="Search species by name..."
            type="text"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            className="w-full border rounded-md px-4 py-2 mb-4 shadow"
          />
          <h2 className="font-semibold mb-2">Search Result</h2>
          <div className="bg-white rounded-lg shadow p-4 overflow-x-auto">
            <table className="w-full text-left">
              <thead>
                <tr className="border-b">
                  <th>Name</th>
                  <th>Status</th>
                  <th>Alignment</th>
                  <th>Metadata</th>
                </tr>
              </thead>
              <tbody>
                {filteredSpecies.map((s, i) => (
                  <tr
                    key={i}
                    className="hover:bg-pink-50 cursor-pointer"
                    onClick={() => setSelected(s)}
                  >
                    <td>{s.name}</td>
                    <td className="italic">{s.status}</td>
                    <td>{s.alignment}</td>
                    <td className="text-blue-600 underline">Detail</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
        <div>
          <h2 className="font-semibold mb-2">Phylogenetic Tree</h2>
          <div className="bg-white rounded-lg shadow p-4 flex justify-center items-center h-[300px]">
            {/* <img
              src="/phylogenetic-tree.png" 
              alt="Phylogenetic Tree"
              className="max-h-full max-w-full"
            /> */}
          </div>
          <p className="text-sm mt-2 flex justify-end">+ Zoom -</p>
          <div className="text-xs mt-2">
            <p><span className="text-red-600 font-bold">●</span> CR (Critically Endangered)</p>
            <p><span className="text-orange-500 font-bold">●</span> EN (Endangered)</p>
            <p><span className="text-yellow-500 font-bold">●</span> VU (Vulnerable)</p>
            <p><span className="text-green-500 font-bold">●</span> LC (Least Concern)</p>
          </div>
        </div>
      </div>
      <div className="mt-6 bg-white rounded-lg shadow p-4 max-w-md">
        <h2 className="font-semibold mb-2">Selected Species</h2>
        <p><span className="font-medium">Name</span> : <i>{selected.name}</i></p>
        <p><span className="font-medium">Status</span> : <i>{selected.status}</i></p>
        <p><span className="font-medium">Origin of habitat</span> : <i>{selected.origin}</i></p>
      </div>
    </div>
  );
}

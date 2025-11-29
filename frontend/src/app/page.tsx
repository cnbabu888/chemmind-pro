"use client";

import React, { useState } from 'react';
import SciFinderLayout from '../components/SciFinderLayout';
import MoleculeSearch from '../components/MoleculeSearch';
import ResultsGrid from '../components/ResultsGrid';
import { AlertCircle, Box, ShieldAlert } from 'lucide-react';
import { useRouter } from 'next/navigation';
import Molecule3DViewer from '../components/Molecule3DViewer';
import SafetyDashboard from '../components/SafetyDashboard';

import BrainAnalysis from '../components/BrainAnalysis';
import AIChatBar from '../components/AIChatBar';

export default function Home() {
    const router = useRouter();
    const [results, setResults] = useState<any[]>([]);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const [activeTab, setActiveTab] = useState<'substances' | 'reactions'>('substances');
    const [selectedSmiles, setSelectedSmiles] = useState<string | null>(null);

    const handleResultClick = (item: any) => {
        // If it's a retrosynthesis result or has a product SMILES, navigate to retro page
        if (item.product) {
            router.push(`/retrosynthesis?smiles=${encodeURIComponent(item.product)}`);
        }
    };

    const handleSearch = async (query: string) => {
        setIsLoading(true);
        setError(null);
        setResults([]);

        try {
            const API_URL = process.env.NEXT_PUBLIC_API_URL || 'http://127.0.0.1:8000';

            // 1. Resolve Name to SMILES (if needed)
            let smiles = query;
            if (!query.includes('=')) { // Simple heuristic for SMILES vs Name
                try {
                    const searchRes = await fetch(`${API_URL}/api/search?query=${encodeURIComponent(query)}`);
                    const searchData = await searchRes.json();
                    if (searchData.smiles) {
                        smiles = searchData.smiles;
                    }
                } catch (e) {
                    console.warn("Search API failed, assuming input is SMILES or raw query");
                }
            }

            // 2. Perform Retrosynthesis (Mock)
            const retroRes = await fetch(`${API_URL}/api/retrosynthesis`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ smiles }),
            });

            if (!retroRes.ok) throw new Error("Retrosynthesis service unavailable");
            const retroData = await retroRes.json();

            // 3. Perform Forward Prediction (Mock - using split parts of the query if valid)
            // For demo, we just use the query as "reactants" if it has a dot, else we skip
            let predictionData = null;
            if (smiles.includes('.')) {
                const predRes = await fetch(`${API_URL}/api/prediction`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ smiles }),
                });
                if (predRes.ok) predictionData = await predRes.json();
            }

            // Transform data for Grid
            const newResults = [];

            // Add Retrosynthesis Result
            if (retroData && retroData.trees) {
                newResults.push({
                    id: 'retro-1',
                    type: 'Retrosynthesis Route',
                    reactants: 'Precursors (See Tree)',
                    product: smiles,
                    confidence: 0.95,
                    conditions: 'Multi-step synthesis suggested'
                });
            }

            // Add Prediction Result
            if (predictionData) {
                newResults.push({
                    id: 'pred-1',
                    type: predictionData.reaction_type,
                    reactants: predictionData.reactants,
                    product: predictionData.product,
                    confidence: predictionData.confidence,
                    conditions: 'Standard Conditions'
                });
            }

            // Add some dummy results to populate the grid for "SciFinder" feel
            newResults.push({
                id: 'related-1',
                type: 'Literature Precedent',
                reactants: 'Similar Reactant A + B',
                product: 'Analog of ' + query,
                yield: '88%',
                conditions: 'Pd(PPh3)4, K2CO3, Toluene, 100°C'
            });

            setResults(newResults);

            // Auto-select the first valid SMILES for visualization
            if (predictionData && predictionData.product) {
                setSelectedSmiles(predictionData.product);
            } else if (smiles && !smiles.includes('.')) {
                setSelectedSmiles(smiles);
            }

        } catch (err) {
            console.log("Search failed using API URL:", process.env.NEXT_PUBLIC_API_URL || 'http://127.0.0.1:8000');
            console.error("Full error object:", err);
            setError("Failed to fetch results. Please ensure the backend is running.");
        } finally {
            setIsLoading(false);
        }
    };

    return (
        <SciFinderLayout>
            <div className="mb-8">
                <h1 className="text-2xl font-bold text-slate-800 mb-2">Substance & Reaction Search</h1>
                <p className="text-slate-500">Access millions of reactions and substances from the ChemMind database.</p>
            </div>

            <MoleculeSearch onSearch={handleSearch} isLoading={isLoading} />

            {error && (
                <div className="bg-red-50 border border-red-200 rounded-md p-4 mb-6 flex items-center text-red-700">
                    <AlertCircle className="h-5 w-5 mr-2" />
                    {error}
                </div>
            )}

            {results.length > 0 && (
                <div className="space-y-6">
                    <ResultsGrid results={results} onResultClick={handleResultClick} />

                    {/* Placeholder for Visualizations */}
                    {/* Visualizations */}
                    <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                        {/* 3D Viewer */}
                        <div className="bg-white rounded-lg shadow-sm border border-slate-200 overflow-hidden flex flex-col h-96">
                            <div className="px-4 py-3 border-b border-slate-100 bg-slate-50 flex justify-between items-center">
                                <h3 className="font-semibold text-slate-700 flex items-center text-sm">
                                    <Box className="h-4 w-4 mr-2 text-blue-600" /> 3D Structure Visualization
                                </h3>
                                {selectedSmiles && <span className="text-[10px] font-mono text-slate-400 truncate max-w-[150px]">{selectedSmiles}</span>}
                            </div>
                            <div className="flex-1 relative bg-white">
                                {selectedSmiles ? (
                                    <Molecule3DViewer smiles={selectedSmiles} />
                                ) : (
                                    <div className="absolute inset-0 flex items-center justify-center text-slate-400 text-sm">
                                        Select a molecule to view 3D structure
                                    </div>
                                )}
                            </div>
                        </div>

                        {/* Safety Dashboard */}
                        <div className="h-96">
                            {selectedSmiles ? (
                                <SafetyDashboard smiles={selectedSmiles} />
                            ) : (
                                <div className="bg-white rounded-lg shadow-sm border border-slate-200 h-full flex items-center justify-center text-slate-400 text-sm">
                                    <div className="flex flex-col items-center">
                                        <ShieldAlert className="h-8 w-8 mb-2 text-slate-300" />
                                        <span>Safety Profile Unavailable</span>
                                    </div>
                                </div>
                            )}
                        </div>
                    </div>

                    {/* AI Brain Analysis */}
                    <div className="h-96">
                        <BrainAnalysis />
                    </div>
                </div>
            )}

            {!isLoading && results.length === 0 && !error && (
                <div className="text-center py-20">
                    <div className="inline-flex items-center justify-center h-16 w-16 rounded-full bg-slate-100 mb-4">
                        <AlertCircle className="h-8 w-8 text-slate-300" />
                    </div>
                    <h3 className="text-lg font-medium text-slate-900">No results to display</h3>
                    <p className="text-slate-500 mt-1">Try searching for a substance like "Aspirin" or a SMILES string.</p>
                </div>
            )}

            <div className="h-20"></div> {/* Spacer for fixed bottom bar */}
            <AIChatBar />
        </SciFinderLayout>
    );
}

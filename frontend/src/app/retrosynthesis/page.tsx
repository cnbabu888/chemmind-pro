"use client";

import React, { useState, useEffect } from 'react';
import ReactionTree from '../../components/ReactionTree';
import { useSearchParams } from 'next/navigation';

export default function RetrosynthesisPage() {
    const searchParams = useSearchParams();
    const [smiles, setSmiles] = useState('Nc1ccccc1'); // Aniline default
    const [treeData, setTreeData] = useState<any>(null);
    const [loading, setLoading] = useState(false);

    // Auto-run if smiles param exists
    useEffect(() => {
        const paramSmiles = searchParams.get('smiles');
        if (paramSmiles) {
            setSmiles(paramSmiles);
            handlePredict(paramSmiles);
        }
    }, [searchParams]);

    const handlePredict = async (targetSmiles: string = smiles) => {
        setLoading(true);
        try {
            const response = await fetch('http://localhost:8000/api/retrosynthesis', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({ smiles: targetSmiles }),
            });
            const data = await response.json();
            setTreeData(data);
        } catch (error) {
            console.error("Failed to fetch retrosynthesis data", error);
            alert("Error fetching data. Ensure backend is running.");
        } finally {
            setLoading(false);
        }
    };

    return (
        <main className="flex min-h-screen flex-col items-center p-10 bg-white dark:bg-zinc-900">
            <h1 className="text-4xl font-bold mb-8 text-gray-800 dark:text-white">Retrosynthesis Engine</h1>

            <div className="w-full max-w-2xl flex gap-4 mb-8">
                <input
                    type="text"
                    value={smiles}
                    onChange={(e) => setSmiles(e.target.value)}
                    placeholder="Enter Target SMILES"
                    className="flex-1 px-4 py-3 border rounded-xl shadow-sm focus:ring-2 focus:ring-blue-500 outline-none dark:bg-zinc-800 dark:text-white"
                />
                <button
                    onClick={() => handlePredict()}
                    disabled={loading}
                    className="px-6 py-3 bg-blue-600 text-white font-semibold rounded-xl hover:bg-blue-700 transition-colors disabled:opacity-50"
                >
                    {loading ? 'Analyzing...' : 'Plan Route'}
                </button>
            </div>

            <div className="w-full h-[600px] bg-gray-50 dark:bg-zinc-800 rounded-2xl border border-gray-200 dark:border-zinc-700 overflow-hidden shadow-inner">
                {treeData ? (
                    <ReactionTree data={treeData} />
                ) : (
                    <div className="flex items-center justify-center h-full text-gray-400">
                        Enter a molecule to see the synthesis tree
                    </div>
                )}
            </div>
        </main>
    );
}

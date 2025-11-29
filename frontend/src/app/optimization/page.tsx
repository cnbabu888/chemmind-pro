"use client";

import React, { useState } from 'react';

export default function OptimizationPage() {
    const [reactants, setReactants] = useState('B(O)O.Br'); // Suzuki default
    const [result, setResult] = useState<any>(null);
    const [loading, setLoading] = useState(false);

    const handleOptimize = async () => {
        setLoading(true);
        try {
            const response = await fetch('http://localhost:8000/api/optimization', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({ smiles: reactants }),
            });
            const data = await response.json();
            setResult(data);
        } catch (error) {
            console.error("Failed to fetch optimization data", error);
            alert("Error fetching data. Ensure backend is running.");
        } finally {
            setLoading(false);
        }
    };

    return (
        <main className="flex min-h-screen flex-col items-center p-10 bg-white dark:bg-zinc-900">
            <h1 className="text-4xl font-bold mb-8 text-gray-800 dark:text-white">Condition Optimization</h1>

            <div className="w-full max-w-3xl flex flex-col gap-4 mb-8">
                <label className="text-sm text-gray-500 dark:text-gray-400">
                    Enter Reaction Mixture (dot-separated SMILES):
                </label>
                <div className="flex gap-4">
                    <input
                        type="text"
                        value={reactants}
                        onChange={(e) => setReactants(e.target.value)}
                        placeholder="e.g. c1ccccc1B(O)O.c1ccccc1Br"
                        className="flex-1 px-4 py-3 border rounded-xl shadow-sm focus:ring-2 focus:ring-purple-500 outline-none dark:bg-zinc-800 dark:text-white"
                    />
                    <button
                        onClick={handleOptimize}
                        disabled={loading}
                        className="px-6 py-3 bg-purple-600 text-white font-semibold rounded-xl hover:bg-purple-700 transition-colors disabled:opacity-50"
                    >
                        {loading ? 'Optimizing...' : 'Suggest Conditions'}
                    </button>
                </div>
            </div>

            {result && (
                <div className="w-full max-w-5xl">
                    <div className="mb-6 p-4 bg-purple-50 dark:bg-purple-900/20 border border-purple-200 dark:border-purple-800 rounded-xl">
                        <h2 className="text-lg font-semibold text-purple-900 dark:text-purple-200">
                            Detected Reaction: {result.reaction_type}
                        </h2>
                    </div>

                    <div className="grid gap-6">
                        {result.conditions.map((condition: any, idx: number) => (
                            <div key={idx} className="bg-white dark:bg-zinc-800 rounded-xl border border-gray-200 dark:border-zinc-700 shadow-sm overflow-hidden">
                                <div className="bg-gray-50 dark:bg-zinc-700/50 px-6 py-3 border-b border-gray-200 dark:border-zinc-700 flex justify-between items-center">
                                    <span className="font-semibold text-gray-700 dark:text-gray-200">Rank #{condition.rank}</span>
                                    <span className="px-3 py-1 bg-green-100 dark:bg-green-900 text-green-700 dark:text-green-300 text-sm font-bold rounded-full">
                                        Yield: {condition.predicted_yield}
                                    </span>
                                </div>
                                <div className="p-6 grid grid-cols-2 md:grid-cols-4 gap-6">
                                    <div>
                                        <span className="text-xs uppercase tracking-wider text-gray-500 dark:text-gray-400">Catalyst</span>
                                        <p className="font-medium text-gray-900 dark:text-white mt-1">{condition.catalyst}</p>
                                    </div>
                                    <div>
                                        <span className="text-xs uppercase tracking-wider text-gray-500 dark:text-gray-400">Solvent</span>
                                        <p className="font-medium text-gray-900 dark:text-white mt-1">{condition.solvent}</p>
                                    </div>
                                    <div>
                                        <span className="text-xs uppercase tracking-wider text-gray-500 dark:text-gray-400">Base/Reagent</span>
                                        <p className="font-medium text-gray-900 dark:text-white mt-1">{condition.base}</p>
                                    </div>
                                    <div>
                                        <span className="text-xs uppercase tracking-wider text-gray-500 dark:text-gray-400">Temp / Time</span>
                                        <p className="font-medium text-gray-900 dark:text-white mt-1">{condition.temperature} / {condition.time}</p>
                                    </div>
                                </div>
                            </div>
                        ))}
                    </div>
                </div>
            )}
        </main>
    );
}

import React, { useEffect, useState } from 'react';
import { ExternalLink, FlaskConical, Search } from 'lucide-react';

interface SimilarMolecule {
    name: string;
    smiles: string;
    similarity: number;
}

interface SimilarityPanelProps {
    smiles: string;
}

const SimilarityPanel: React.FC<SimilarityPanelProps> = ({ smiles }) => {
    const [results, setResults] = useState<SimilarMolecule[]>([]);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);

    useEffect(() => {
        if (!smiles) return;

        const fetchSimilar = async () => {
            setLoading(true);
            setError(null);
            try {
                const response = await fetch(`${process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000'}/api/molecule/similar?smiles=${encodeURIComponent(smiles)}&threshold=0.2`);
                if (!response.ok) throw new Error('Failed to fetch similar molecules');
                const data = await response.json();
                setResults(data.results || []);
            } catch (err) {
                setError('Could not load similar molecules.');
                console.error(err);
            } finally {
                setLoading(false);
            }
        };

        fetchSimilar();
    }, [smiles]);

    if (loading) return <div className="p-4 text-slate-500 text-sm animate-pulse">Searching chemical database...</div>;
    if (error) return <div className="p-4 text-red-500 text-sm">{error}</div>;
    if (results.length === 0) return <div className="p-4 text-slate-500 text-sm italic">No similar molecules found in database.</div>;

    return (
        <div className="space-y-4">
            <div className="flex items-center justify-between">
                <h4 className="text-xs font-semibold text-slate-500 uppercase flex items-center">
                    <Search className="h-3 w-3 mr-1" /> Structural Analogs
                </h4>
                <span className="text-xs text-slate-400">{results.length} matches</span>
            </div>

            <div className="space-y-3">
                {results.map((mol, idx) => (
                    <div key={idx} className="bg-white border border-slate-200 rounded-lg p-3 hover:shadow-md transition-shadow">
                        <div className="flex justify-between items-start mb-2">
                            <span className="font-medium text-slate-800 text-sm">{mol.name}</span>
                            <span className={`text-xs font-bold px-1.5 py-0.5 rounded ${mol.similarity > 0.7 ? 'bg-green-100 text-green-700' : 'bg-blue-50 text-blue-600'}`}>
                                {Math.round(mol.similarity * 100)}% Sim
                            </span>
                        </div>

                        <div className="flex items-center space-x-3">
                            <div className="bg-slate-50 p-1 rounded border border-slate-100 w-16 h-16 flex items-center justify-center">
                                <img
                                    src={`${process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000'}/api/molecule/image?smiles=${encodeURIComponent(mol.smiles)}`}
                                    alt={mol.name}
                                    className="max-w-full max-h-full object-contain"
                                />
                            </div>
                            <div className="flex-1 space-y-2">
                                <div className="text-[10px] font-mono text-slate-500 truncate w-32" title={mol.smiles}>
                                    {mol.smiles}
                                </div>
                                <div className="flex space-x-1">
                                    <button className="flex-1 bg-blue-50 text-blue-600 text-[10px] py-1 rounded hover:bg-blue-100 transition-colors flex items-center justify-center">
                                        <FlaskConical className="h-3 w-3 mr-1" /> Reactant
                                    </button>
                                    <button className="flex-1 bg-slate-50 text-slate-600 text-[10px] py-1 rounded hover:bg-slate-100 transition-colors flex items-center justify-center">
                                        <ExternalLink className="h-3 w-3 mr-1" /> View
                                    </button>
                                </div>
                            </div>
                        </div>
                    </div>
                ))}
            </div>
        </div>
    );
};

export default SimilarityPanel;

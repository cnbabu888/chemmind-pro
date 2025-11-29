import React, { useEffect, useState } from 'react';
import { ShieldAlert, Skull, Flame, AlertTriangle, Info } from 'lucide-react';

interface SafetyData {
    molecule_name: string;
    ghs_hazards: { code: string; description: string; category: string }[];
    toxicity: {
        ld50_rat_oral: string;
        carcinogenicity: string;
        mutagenicity: string;
    };
    handling: string[];
    ppe: string[];
    signal_word: string;
}

interface SafetyPanelProps {
    smiles: string;
}

const SafetyPanel: React.FC<SafetyPanelProps> = ({ smiles }) => {
    const [data, setData] = useState<SafetyData | null>(null);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);

    useEffect(() => {
        if (!smiles) return;

        const fetchSafety = async () => {
            setLoading(true);
            setError(null);
            try {
                const response = await fetch(`${process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000'}/api/safety`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ smiles })
                });
                if (!response.ok) throw new Error('Failed to fetch safety data');
                const result = await response.json();
                setData(result);
            } catch (err) {
                setError('Could not load safety profile.');
                console.error(err);
            } finally {
                setLoading(false);
            }
        };

        fetchSafety();
    }, [smiles]);

    if (loading) return <div className="p-4 text-slate-500 text-sm animate-pulse">Analyzing safety profile...</div>;
    if (error) return <div className="p-4 text-red-500 text-sm">{error}</div>;
    if (!data) return null;

    const isDanger = data.signal_word === 'DANGER';

    return (
        <div className="space-y-6">
            {/* Header Badge */}
            <div className={`flex items-center justify-between p-3 rounded-md border ${isDanger ? 'bg-red-50 border-red-200' : 'bg-amber-50 border-amber-200'}`}>
                <div className="flex items-center space-x-2">
                    {isDanger ? <Skull className="h-5 w-5 text-red-600" /> : <AlertTriangle className="h-5 w-5 text-amber-600" />}
                    <span className={`font-bold ${isDanger ? 'text-red-700' : 'text-amber-700'}`}>
                        SIGNAL WORD: {data.signal_word}
                    </span>
                </div>
            </div>

            {/* GHS Hazards */}
            <div>
                <h4 className="text-xs font-semibold text-slate-500 uppercase mb-2 flex items-center">
                    <Flame className="h-3 w-3 mr-1" /> GHS Hazards
                </h4>
                <div className="space-y-2">
                    {data.ghs_hazards.map((hazard, idx) => (
                        <div key={idx} className="flex items-start space-x-2 text-sm bg-white p-2 rounded border border-slate-100">
                            <span className="font-mono font-bold text-slate-600 bg-slate-100 px-1 rounded text-xs">{hazard.code}</span>
                            <span className="text-slate-700">{hazard.description}</span>
                        </div>
                    ))}
                </div>
            </div>

            {/* Toxicity Data */}
            <div>
                <h4 className="text-xs font-semibold text-slate-500 uppercase mb-2 flex items-center">
                    <ShieldAlert className="h-3 w-3 mr-1" /> Toxicity Profile (AI Predicted)
                </h4>
                <div className="grid grid-cols-2 gap-2 text-sm">
                    <div className="bg-slate-50 p-2 rounded">
                        <span className="block text-xs text-slate-400">LD50 (Rat, Oral)</span>
                        <span className="font-medium text-slate-700">{data.toxicity.ld50_rat_oral}</span>
                    </div>
                    <div className="bg-slate-50 p-2 rounded">
                        <span className="block text-xs text-slate-400">Carcinogenicity</span>
                        <span className="font-medium text-slate-700">{data.toxicity.carcinogenicity}</span>
                    </div>
                </div>
            </div>

            {/* Handling & PPE */}
            <div>
                <h4 className="text-xs font-semibold text-slate-500 uppercase mb-2 flex items-center">
                    <Info className="h-3 w-3 mr-1" /> Handling & PPE
                </h4>
                <div className="bg-blue-50 p-3 rounded-md border border-blue-100">
                    <div className="flex flex-wrap gap-2 mb-3">
                        {data.ppe.map((item, idx) => (
                            <span key={idx} className="bg-white text-blue-700 px-2 py-1 rounded text-xs font-medium border border-blue-200 shadow-sm">
                                {item}
                            </span>
                        ))}
                    </div>
                    <ul className="list-disc list-inside text-xs text-blue-800 space-y-1">
                        {data.handling.map((instruction, idx) => (
                            <li key={idx}>{instruction}</li>
                        ))}
                    </ul>
                </div>
            </div>
        </div>
    );
};

export default SafetyPanel;

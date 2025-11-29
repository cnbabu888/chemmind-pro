import React, { useEffect, useState } from 'react';
import { AlertTriangle, Shield, Skull, Flame, Activity, Info } from 'lucide-react';

interface SafetyData {
    molecule_name: string;
    ghs_hazards: {
        code: string;
        description: string;
        category: string;
    }[];
    toxicity: {
        ld50_rat_oral: string;
        carcinogenicity: string;
        mutagenicity: string;
    };
    handling: string[];
    ppe: string[];
    signal_word: string;
    structural_alerts?: string[];
}


interface SafetyDashboardProps {
    smiles: string;
}

const SafetyDashboard: React.FC<SafetyDashboardProps> = ({ smiles }) => {
    const [data, setData] = useState<SafetyData | null>(null);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);

    useEffect(() => {
        if (!smiles) return;

        const fetchData = async () => {
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
                console.error(err);
                setError('Could not load safety profile.');
            } finally {
                setLoading(false);
            }
        };

        fetchData();
    }, [smiles]);

    if (loading) return (
        <div className="h-full flex flex-col items-center justify-center text-slate-400 space-y-2">
            <Activity className="h-6 w-6 animate-pulse" />
            <span className="text-xs">Analyzing Safety Profile...</span>
        </div>
    );

    if (error) return (
        <div className="h-full flex items-center justify-center text-red-400 text-xs p-4 text-center border border-red-100 bg-red-50 rounded">
            {error}
        </div>
    );

    if (!data) return null;

    const isDanger = data.signal_word === 'DANGER';

    return (
        <div className="bg-white rounded-lg border border-slate-200 shadow-sm overflow-hidden h-full flex flex-col">
            {/* Header */}
            <div className={`px-4 py-3 border-b ${isDanger ? 'bg-red-50 border-red-100' : 'bg-yellow-50 border-yellow-100'} flex justify-between items-center`}>
                <div className="flex items-center space-x-2">
                    {isDanger ? <Skull className="h-5 w-5 text-red-600" /> : <AlertTriangle className="h-5 w-5 text-yellow-600" />}
                    <div>
                        <h3 className={`font-bold text-sm ${isDanger ? 'text-red-800' : 'text-yellow-800'}`}>
                            {data.molecule_name || 'Unknown Compound'}
                        </h3>
                        <p className="text-[10px] text-slate-500 uppercase tracking-wider font-semibold">
                            Signal Word: <span className={isDanger ? 'text-red-600' : 'text-yellow-600'}>{data.signal_word}</span>
                        </p>
                    </div>
                </div>
                {/* GHS Pictograms (Mock) */}
                <div className="flex space-x-1">
                    {data.ghs_hazards.some(h => h.category.includes('Flammability')) && (
                        <div className="h-8 w-8 border-2 border-red-600 transform rotate-45 flex items-center justify-center bg-white" title="Flammable">
                            <Flame className="h-4 w-4 text-black transform -rotate-45" />
                        </div>
                    )}
                    {isDanger && (
                        <div className="h-8 w-8 border-2 border-red-600 transform rotate-45 flex items-center justify-center bg-white" title="Toxic">
                            <Skull className="h-4 w-4 text-black transform -rotate-45" />
                        </div>
                    )}
                </div>
            </div>

            {/* Structural Alerts Banner */}
            {data.structural_alerts && data.structural_alerts.length > 0 && (
                <div className="bg-orange-50 px-4 py-2 border-b border-orange-100 flex items-start">
                    <AlertTriangle className="h-4 w-4 text-orange-600 mt-0.5 mr-2 flex-shrink-0" />
                    <div>
                        <span className="text-xs font-bold text-orange-800 block">Structural Alerts Detected:</span>
                        <div className="flex flex-wrap gap-1 mt-1">
                            {data.structural_alerts.map((alert, idx) => (
                                <span key={idx} className="text-[10px] bg-orange-100 text-orange-800 px-1.5 py-0.5 rounded border border-orange-200">
                                    {alert}
                                </span>
                            ))}
                        </div>
                    </div>
                </div>
            )}

            <div className="p-4 overflow-y-auto flex-1 space-y-4">
                {/* Hazards */}
                <div>
                    <h4 className="text-xs font-semibold text-slate-700 mb-2 flex items-center">
                        <AlertTriangle className="h-3 w-3 mr-1" /> GHS Hazard Statements
                    </h4>
                    <ul className="space-y-1">
                        {data.ghs_hazards.map((h, idx) => (
                            <li key={idx} className="text-xs text-slate-600 flex items-start">
                                <span className="font-mono text-slate-400 mr-2 text-[10px]">{h.code}</span>
                                {h.description}
                            </li>
                        ))}
                    </ul>
                </div>

                <div className="grid grid-cols-2 gap-4">
                    {/* Toxicity */}
                    <div>
                        <h4 className="text-xs font-semibold text-slate-700 mb-2 flex items-center">
                            <Skull className="h-3 w-3 mr-1" /> Toxicity Data
                        </h4>
                        <div className="space-y-2">
                            <div className="bg-slate-50 p-2 rounded border border-slate-100">
                                <span className="text-[10px] text-slate-400 block">LD50 (Rat, Oral)</span>
                                <span className="text-xs font-medium text-slate-700">{data.toxicity.ld50_rat_oral}</span>
                            </div>
                            <div className="bg-slate-50 p-2 rounded border border-slate-100">
                                <span className="text-[10px] text-slate-400 block">Carcinogenicity</span>
                                <span className="text-xs font-medium text-slate-700">{data.toxicity.carcinogenicity}</span>
                            </div>
                        </div>
                    </div>

                    {/* PPE */}
                    <div>
                        <h4 className="text-xs font-semibold text-slate-700 mb-2 flex items-center">
                            <Shield className="h-3 w-3 mr-1" /> Required PPE
                        </h4>
                        <div className="flex flex-wrap gap-1">
                            {data.ppe.map((item, idx) => (
                                <span key={idx} className="px-2 py-1 bg-blue-50 text-blue-700 text-[10px] rounded-full border border-blue-100 font-medium">
                                    {item}
                                </span>
                            ))}
                        </div>
                    </div>
                </div>

                {/* Handling */}
                <div>
                    <h4 className="text-xs font-semibold text-slate-700 mb-2 flex items-center">
                        <Info className="h-3 w-3 mr-1" /> Handling Precautions
                    </h4>
                    <ul className="list-disc list-inside space-y-1">
                        {data.handling.map((item, idx) => (
                            <li key={idx} className="text-xs text-slate-600 pl-1">
                                {item}
                            </li>
                        ))}
                    </ul>
                </div>
            </div>
        </div>
    );
};

export default SafetyDashboard;

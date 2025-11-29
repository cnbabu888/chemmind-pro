import React from 'react';
import { Shield, AlertTriangle, Flame, Skull, Zap, Info } from 'lucide-react';

interface SafetyPanelProps {
    data: {
        structural_alerts: string[];
        ghs: {
            ghs_pictograms?: string[];
            signal_word?: string;
            hazard_statements?: string[];
            precautionary_summary?: string;
            ppe?: string[];
        };
    };
    onClose: () => void;
}

const SafetyPanel: React.FC<SafetyPanelProps> = ({ data, onClose }) => {
    const getIcon = (name: string) => {
        const n = name.toLowerCase();
        if (n.includes('flame') || n.includes('fire')) return <Flame className="h-6 w-6 text-red-500" />;
        if (n.includes('skull') || n.includes('toxic')) return <Skull className="h-6 w-6 text-purple-600" />;
        if (n.includes('corrosive')) return <Zap className="h-6 w-6 text-yellow-600" />; // Proxy
        if (n.includes('health')) return <Shield className="h-6 w-6 text-blue-600" />;
        return <AlertTriangle className="h-6 w-6 text-orange-500" />;
    };

    return (
        <div className="absolute top-16 right-4 w-80 bg-white shadow-xl border border-slate-300 rounded-md overflow-hidden z-40 font-sans">
            <div className="bg-gradient-to-r from-red-600 to-orange-500 px-3 py-2 flex justify-between items-center text-white">
                <div className="flex items-center space-x-2">
                    <Shield className="h-4 w-4" />
                    <span className="font-bold text-sm">Safety Dashboard</span>
                </div>
                <button onClick={onClose} className="hover:bg-white/20 rounded p-1">
                    <span className="text-xs">✕</span>
                </button>
            </div>

            <div className="p-4 max-h-[80vh] overflow-y-auto">
                {/* Signal Word */}
                {data.ghs.signal_word && (
                    <div className={`text-center font-black text-xl mb-4 uppercase tracking-wider ${data.ghs.signal_word === 'Danger' ? 'text-red-600' : 'text-orange-500'}`}>
                        {data.ghs.signal_word}
                    </div>
                )}

                {/* GHS Pictograms */}
                {data.ghs.ghs_pictograms && data.ghs.ghs_pictograms.length > 0 && (
                    <div className="flex justify-center space-x-4 mb-6">
                        {data.ghs.ghs_pictograms.map((pic, i) => (
                            <div key={i} className="flex flex-col items-center" title={pic}>
                                <div className="border-2 border-red-600 p-2 rounded-lg transform rotate-45 mb-2">
                                    <div className="transform -rotate-45">
                                        {getIcon(pic)}
                                    </div>
                                </div>
                                <span className="text-[10px] text-slate-500 font-semibold">{pic}</span>
                            </div>
                        ))}
                    </div>
                )}

                {/* Structural Alerts */}
                {data.structural_alerts.length > 0 && (
                    <div className="mb-4">
                        <h4 className="text-xs font-bold text-slate-700 uppercase mb-2 border-b border-slate-200 pb-1">Structural Alerts</h4>
                        <ul className="space-y-1">
                            {data.structural_alerts.map((alert, i) => (
                                <li key={i} className="flex items-start space-x-2 text-xs bg-red-50 p-2 rounded border border-red-100 text-red-800">
                                    <AlertTriangle className="h-3 w-3 mt-0.5 flex-shrink-0" />
                                    <span>{alert}</span>
                                </li>
                            ))}
                        </ul>
                    </div>
                )}

                {/* Hazard Statements */}
                {data.ghs.hazard_statements && (
                    <div className="mb-4">
                        <h4 className="text-xs font-bold text-slate-700 uppercase mb-2 border-b border-slate-200 pb-1">Hazards</h4>
                        <ul className="list-disc list-inside text-xs text-slate-600 space-y-1">
                            {data.ghs.hazard_statements.map((h, i) => <li key={i}>{h}</li>)}
                        </ul>
                    </div>
                )}

                {/* Precautionary & PPE */}
                <div className="bg-blue-50 p-3 rounded border border-blue-100 text-xs">
                    <div className="flex items-start space-x-2 mb-2">
                        <Info className="h-4 w-4 text-blue-600 mt-0.5 flex-shrink-0" />
                        <p className="text-blue-800 font-medium">{data.ghs.precautionary_summary}</p>
                    </div>
                    {data.ghs.ppe && (
                        <div className="flex flex-wrap gap-1 mt-2">
                            {data.ghs.ppe.map((item, i) => (
                                <span key={i} className="px-2 py-0.5 bg-white border border-blue-200 rounded-full text-[10px] text-blue-600 font-semibold">
                                    {item}
                                </span>
                            ))}
                        </div>
                    )}
                </div>

                <div className="mt-4 text-[9px] text-slate-400 text-center italic">
                    AI-generated safety data. Consult official SDS.
                </div>
            </div>
        </div>
    );
};

export default SafetyPanel;

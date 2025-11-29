import React from 'react';
import { BookOpen, FlaskConical, Thermometer, Clock, Percent } from 'lucide-react';

interface Condition {
    rank: number;
    catalyst: string;
    solvent: string;
    base?: string;
    temperature: string;
    time: string;
    predicted_yield: string;
    source?: string;
}

interface ConditionDetailProps {
    conditions: Condition[];
    reactionType: string;
}

const ConditionDetail: React.FC<ConditionDetailProps> = ({ conditions, reactionType }) => {
    return (
        <div className="bg-slate-50 p-4 rounded-lg border border-slate-200 mt-2">
            <div className="flex items-center mb-3">
                <FlaskConical className="h-4 w-4 text-blue-600 mr-2" />
                <h4 className="font-semibold text-slate-800 text-sm">Optimized Conditions: {reactionType}</h4>
            </div>

            <div className="overflow-x-auto">
                <table className="w-full text-xs text-left">
                    <thead className="bg-slate-100 text-slate-500 font-semibold border-b border-slate-200">
                        <tr>
                            <th className="px-3 py-2">Rank</th>
                            <th className="px-3 py-2">Catalyst / Reagents</th>
                            <th className="px-3 py-2">Solvent</th>
                            <th className="px-3 py-2">Temp / Time</th>
                            <th className="px-3 py-2">Yield</th>
                            <th className="px-3 py-2">Source</th>
                        </tr>
                    </thead>
                    <tbody className="divide-y divide-slate-200">
                        {conditions.map((cond, idx) => (
                            <tr key={idx} className="bg-white hover:bg-blue-50 transition-colors">
                                <td className="px-3 py-2 font-medium text-slate-700">#{cond.rank}</td>
                                <td className="px-3 py-2">
                                    <div className="font-medium text-slate-900">{cond.catalyst}</div>
                                    {cond.base && <div className="text-slate-500 mt-0.5">+ {cond.base}</div>}
                                </td>
                                <td className="px-3 py-2 text-slate-600">{cond.solvent}</td>
                                <td className="px-3 py-2 text-slate-600">
                                    <div className="flex items-center"><Thermometer className="h-3 w-3 mr-1" /> {cond.temperature}</div>
                                    <div className="flex items-center mt-0.5"><Clock className="h-3 w-3 mr-1" /> {cond.time}</div>
                                </td>
                                <td className="px-3 py-2 font-bold text-green-700">
                                    <div className="flex items-center"><Percent className="h-3 w-3 mr-1" /> {cond.predicted_yield}</div>
                                </td>
                                <td className="px-3 py-2 text-blue-600 italic">
                                    {cond.source ? (
                                        <div className="flex items-center cursor-pointer hover:underline">
                                            <BookOpen className="h-3 w-3 mr-1" />
                                            {cond.source}
                                        </div>
                                    ) : (
                                        <span className="text-slate-400">AI Prediction</span>
                                    )}
                                </td>
                            </tr>
                        ))}
                    </tbody>
                </table>
            </div>
        </div>
    );
};

export default ConditionDetail;

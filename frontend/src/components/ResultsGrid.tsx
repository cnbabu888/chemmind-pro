import React, { useState } from 'react';
import { ArrowRight, Beaker, AlertCircle, CheckCircle2, ChevronDown, ChevronUp } from 'lucide-react';
import ConditionDetail from './ConditionDetail';

interface ResultItem {
    id: string;
    reactants: string;
    product: string;
    yield?: string;
    conditions?: string;
    confidence?: number;
    type: string;
    source?: string; // Added source field
    conditions_data?: any;
}



interface ResultsGridProps {
    results: ResultItem[];
    title?: string;
    onResultClick?: (item: ResultItem) => void;
}

const ResultsGrid: React.FC<ResultsGridProps> = ({ results, title = "Reaction Results", onResultClick }) => {
    const [expandedRowId, setExpandedRowId] = useState<string | null>(null);

    if (!results || results.length === 0) {
        return null;
    }

    const toggleExpand = (e: React.MouseEvent, id: string) => {
        e.stopPropagation();
        setExpandedRowId(expandedRowId === id ? null : id);
    };

    return (
        <div className="bg-white rounded-lg shadow-sm border border-slate-200 overflow-hidden">
            <div className="bg-slate-50 px-6 py-4 border-b border-slate-200 flex justify-between items-center">
                <h3 className="font-semibold text-slate-800">{title} <span className="text-slate-400 font-normal text-sm ml-2">({results.length} found)</span></h3>
                <div className="flex space-x-2">
                    <button className="text-xs font-medium text-[#0056b3] hover:underline">Export PDF</button>
                    <span className="text-slate-300">|</span>
                    <button className="text-xs font-medium text-[#0056b3] hover:underline">Save to Project</button>
                </div>
            </div>

            <div className="overflow-x-auto">
                <table className="w-full text-left text-sm text-slate-600">
                    <thead className="bg-slate-50 text-xs uppercase font-semibold text-slate-500">
                        <tr>
                            <th className="px-4 py-3 border-b border-slate-200 w-10"></th>
                            <th className="px-6 py-3 border-b border-slate-200">Status</th>
                            <th className="px-6 py-3 border-b border-slate-200">Reaction</th>
                            <th className="px-6 py-3 border-b border-slate-200">Type</th>
                            <th className="px-6 py-3 border-b border-slate-200">Conditions</th>
                            <th className="px-6 py-3 border-b border-slate-200">Source</th>
                            <th className="px-6 py-3 border-b border-slate-200 text-right">Confidence / Yield</th>
                        </tr>
                    </thead>
                    <tbody className="divide-y divide-slate-100">
                        {results.map((item) => (
                            <React.Fragment key={item.id}>
                                <tr
                                    className={`hover:bg-slate-50 transition-colors ${onResultClick ? 'cursor-pointer' : ''} ${expandedRowId === item.id ? 'bg-blue-50/30' : ''}`}
                                    onClick={() => onResultClick && onResultClick(item)}
                                >
                                    <td className="px-4 py-4 text-center">
                                        <button
                                            onClick={(e) => toggleExpand(e, item.id)}
                                            className="p-1 hover:bg-slate-200 rounded text-slate-400 hover:text-slate-600 transition-colors"
                                        >
                                            {expandedRowId === item.id ? <ChevronUp className="h-4 w-4" /> : <ChevronDown className="h-4 w-4" />}
                                        </button>
                                    </td>
                                    <td className="px-6 py-4 whitespace-nowrap">
                                        {item.confidence && item.confidence > 0.8 ? (
                                            <div className="flex items-center text-green-600"><CheckCircle2 className="h-4 w-4 mr-1" /> High</div>
                                        ) : (
                                            <div className="flex items-center text-amber-500"><AlertCircle className="h-4 w-4 mr-1" /> Medium</div>
                                        )}
                                    </td>
                                    <td className="px-6 py-4">
                                        <div className="flex items-center space-x-2">
                                            <span className="font-mono bg-slate-100 px-2 py-1 rounded text-xs text-slate-700 max-w-[150px] truncate" title={item.reactants}>{item.reactants}</span>
                                            <ArrowRight className="h-4 w-4 text-slate-400" />
                                            <span className="font-mono bg-slate-100 px-2 py-1 rounded text-xs text-slate-700 max-w-[150px] truncate" title={item.product}>{item.product}</span>
                                        </div>
                                    </td>
                                    <td className="px-6 py-4">
                                        <span className="inline-flex items-center px-2.5 py-0.5 rounded-full text-xs font-medium bg-blue-50 text-blue-700">
                                            {item.type}
                                        </span>
                                    </td>
                                    <td className="px-6 py-4 text-xs text-slate-500">
                                        {item.conditions || "Not specified"}
                                    </td>
                                    <td className="px-6 py-4 text-xs text-slate-500 italic">
                                        {item.source || "AI Prediction"}
                                    </td>
                                    <td className="px-6 py-4 text-right font-medium text-slate-900">
                                        {item.yield || (item.confidence ? `${Math.round(item.confidence * 100)}%` : "-")}
                                    </td>
                                </tr>
                                {expandedRowId === item.id && (
                                    <tr>
                                        <td colSpan={7} className="px-6 pb-4 pt-0 bg-blue-50/30">
                                            {item.conditions_data ? (
                                                <ConditionDetail
                                                    conditions={item.conditions_data.conditions}
                                                    reactionType={item.conditions_data.reaction_type}
                                                />
                                            ) : (
                                                <div className="p-4 text-center text-slate-500 text-sm italic">
                                                    No detailed conditions available for this reaction.
                                                </div>
                                            )}
                                        </td>
                                    </tr>
                                )}
                            </React.Fragment>
                        ))}
                    </tbody>
                </table >
            </div >
        </div >
    );
};

export default ResultsGrid;

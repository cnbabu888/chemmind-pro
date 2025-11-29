import React, { useState, useEffect } from 'react';
import { Search, Upload, FileText, PenTool, Image as ImageIcon } from 'lucide-react';
import StructureEditor from './StructureEditor';

interface MoleculeSearchProps {
    onSearch: (query: string) => void;
    isLoading?: boolean;
}

const MoleculeSearch: React.FC<MoleculeSearchProps> = ({ onSearch, isLoading = false }) => {
    const [query, setQuery] = useState('');
    const [searchType, setSearchType] = useState<'structure' | 'reaction'>('structure');
    const [debouncedQuery, setDebouncedQuery] = useState(query);
    const [isEditorOpen, setIsEditorOpen] = useState(false);
    const [previewSmiles, setPreviewSmiles] = useState<string | null>(null);
    const [previewUrl, setPreviewUrl] = useState<string | null>(null);
    const [isResolving, setIsResolving] = useState(false);

    // Debounce query for preview
    useEffect(() => {
        const handler = setTimeout(() => {
            setDebouncedQuery(query);
        }, 800); // Increased debounce for API calls
        return () => clearTimeout(handler);
    }, [query]);

    // Resolve structure and fetch preview
    useEffect(() => {
        const resolveAndPreview = async () => {
            if (!debouncedQuery || debouncedQuery.length < 2) {
                setPreviewUrl(null);
                setPreviewSmiles(null);
                return;
            }

            let smilesToPreview = debouncedQuery;
            const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'http://127.0.0.1:8000';

            // Simple heuristic: If it doesn't look like SMILES (has no special chars usually found in SMILES but not names, or has spaces), treat as name
            // This is a rough check. Better to try SMILES parsing first or rely on backend.
            // For now: if it has spaces or is purely alphabetic (like "Aspirin"), treat as name.
            const isLikelyName = !debouncedQuery.includes('=') && !debouncedQuery.includes('#') && (/^[a-zA-Z0-9\s\-\(\)]+$/.test(debouncedQuery) || debouncedQuery.includes(' '));

            if (isLikelyName && searchType === 'structure') {
                setIsResolving(true);
                try {
                    const res = await fetch(`${apiUrl}/api/molecule/name_to_structure`, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ name: debouncedQuery })
                    });
                    const data = await res.json();
                    if (data.smiles) {
                        smilesToPreview = data.smiles;
                    }
                } catch (e) {
                    console.error("Name resolution failed", e);
                    setPreviewUrl('error');
                } finally {
                    setIsResolving(false);
                }
            }

            setPreviewSmiles(smilesToPreview);
            const url = `${apiUrl}/api/molecule/image?smiles=${encodeURIComponent(smilesToPreview)}`;
            setPreviewUrl(url);
        };

        resolveAndPreview();
    }, [debouncedQuery, searchType]);

    const handleSubmit = (e: React.FormEvent) => {
        e.preventDefault();
        // Use the resolved SMILES if available, otherwise the query
        if (previewSmiles) {
            onSearch(previewSmiles);
        } else if (query.trim()) {
            onSearch(query);
        }
    };

    const handleEditorSave = (smiles: string) => {
        setQuery(smiles);
        setPreviewSmiles(smiles); // Immediate update
    };

    const handleEditPreview = (e: React.MouseEvent) => {
        e.stopPropagation();
        e.preventDefault();
        if (previewSmiles) {
            setQuery(previewSmiles); // Ensure editor gets the SMILES
            setIsEditorOpen(true);
        }
    };

    return (
        <div className="bg-white rounded-lg shadow-sm border border-slate-200 p-6 mb-6">
            <div className="flex space-x-6 mb-4 border-b border-slate-100 pb-2">
                <button
                    onClick={() => setSearchType('structure')}
                    className={`pb-2 text-sm font-medium transition-colors border-b-2 ${searchType === 'structure' ? 'border-[#0056b3] text-[#0056b3]' : 'border-transparent text-slate-500 hover:text-slate-700'}`}
                >
                    Structure Search
                </button>
                <button
                    onClick={() => setSearchType('reaction')}
                    className={`pb-2 text-sm font-medium transition-colors border-b-2 ${searchType === 'reaction' ? 'border-[#0056b3] text-[#0056b3]' : 'border-transparent text-slate-500 hover:text-slate-700'}`}
                >
                    Reaction Search
                </button>
            </div>

            <div className="flex gap-6">
                <div className="flex-1">
                    <form onSubmit={handleSubmit} className="space-y-4">
                        <div className="relative">
                            <div className="absolute inset-y-0 left-0 pl-3 flex items-center pointer-events-none">
                                <Search className="h-5 w-5 text-slate-400" />
                            </div>
                            <input
                                type="text"
                                value={query}
                                onChange={(e) => setQuery(e.target.value)}
                                className="block w-full pl-10 pr-24 py-3 border border-slate-300 rounded-md leading-5 bg-white placeholder-slate-400 focus:outline-none focus:ring-2 focus:ring-[#0056b3] focus:border-transparent sm:text-sm transition-shadow"
                                placeholder={searchType === 'structure' ? "Enter CAS RN, Name, or SMILES (e.g., 'Aspirin', 'CC(=O)O')..." : "Enter Reaction SMILES..."}
                            />
                            <div className="absolute inset-y-0 right-0 flex items-center pr-2 space-x-1">
                                <button
                                    type="button"
                                    onClick={() => setIsEditorOpen(true)}
                                    className="p-1.5 text-slate-400 hover:text-[#0056b3] rounded-md hover:bg-blue-50 transition-colors"
                                    title="Draw Structure"
                                >
                                    <PenTool className="h-4 w-4" />
                                </button>
                                <button type="button" className="p-1.5 text-slate-400 hover:text-slate-600 rounded-md hover:bg-slate-100">
                                    <Upload className="h-4 w-4" />
                                </button>
                            </div>
                        </div>

                        <div className="flex justify-between items-center">
                            <div className="text-xs text-slate-500">
                                {searchType === 'structure'
                                    ? "Supported formats: SMILES, InChI, CAS RN, Name"
                                    : "Format: Reactants>>Products (e.g., CC(=O)O.OCC>>CC(=O)OCC)"}
                            </div>
                            <button
                                type="submit"
                                disabled={isLoading || !query.trim()}
                                className={`px-6 py-2 rounded-md text-white font-medium text-sm transition-all shadow-sm ${isLoading || !query.trim() ? 'bg-slate-300 cursor-not-allowed' : 'bg-[#0056b3] hover:bg-[#004494] hover:shadow-md'}`}
                            >
                                {isLoading ? 'Searching...' : 'Search'}
                            </button>
                        </div>
                    </form>
                </div>

                {/* Live Preview Panel */}
                <div
                    id="molecule-preview-container"
                    onClick={() => previewSmiles && setIsEditorOpen(true)}
                    className={`w-48 h-32 border border-slate-200 rounded-lg bg-slate-50 flex items-center justify-center overflow-hidden relative group transition-all duration-300 hover:scale-[1.5] hover:z-20 origin-top-right shadow-sm hover:shadow-xl bg-white ${previewSmiles ? 'cursor-pointer' : ''}`}
                >
                    {isResolving ? (
                        <div className="flex flex-col items-center text-slate-400">
                            <div className="animate-spin rounded-full h-5 w-5 border-b-2 border-blue-500 mb-2"></div>
                            <span className="text-[10px]">Resolving...</span>
                        </div>
                    ) : previewUrl === 'error' ? (
                        <div className="flex flex-col items-center text-red-400 p-2 text-center">
                            <span className="text-[10px]">Could not resolve structure.</span>
                        </div>
                    ) : previewUrl ? (
                        <>
                            <img
                                src={previewUrl}
                                alt="Structure Preview"
                                className="max-w-full max-h-full object-contain p-2"
                                onError={(e) => {
                                    (e.target as HTMLImageElement).style.display = 'none';
                                    setPreviewUrl(null); // Hide on error
                                }}
                            />
                            {/* Edit Overlay */}
                            <div className="absolute inset-0 bg-black/0 hover:bg-black/10 transition-colors flex items-center justify-center">
                                <button
                                    className="bg-white text-blue-600 px-3 py-1 rounded-full text-xs font-medium shadow-md hover:bg-blue-50 flex items-center transform scale-90 hover:scale-100 transition-all"
                                >
                                    <PenTool className="h-3 w-3 mr-1" /> Edit
                                </button>
                            </div>
                        </>
                    ) : (
                        <div className="text-center p-4">
                            <ImageIcon className="h-8 w-8 text-slate-300 mx-auto mb-1" />
                            <span className="text-xs text-slate-400">Preview</span>
                        </div>
                    )}
                </div>
            </div>

            {isEditorOpen && (
                <StructureEditor
                    initialSmiles={query}
                    onSave={handleEditorSave}
                    onClose={() => setIsEditorOpen(false)}
                />
            )}
        </div>
    );
};

export default MoleculeSearch;

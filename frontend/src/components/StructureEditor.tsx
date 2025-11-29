import React, { useState, useRef, useEffect } from 'react';
import { X, Save, Maximize2, Minimize2, ChevronDown, Move } from 'lucide-react';
import { Editor } from 'ketcher-react';
import { Ketcher } from 'ketcher-core';
import 'ketcher-react/dist/index.css';
import Draggable from 'react-draggable'; // You might need to install this: npm install react-draggable

interface StructureEditorProps {
    initialSmiles?: string;
    onSave: (smiles: string) => void;
    onClose: () => void;
}

const StructureEditor: React.FC<StructureEditorProps> = ({ initialSmiles = '', onSave, onClose }) => {
    const [ketcher, setKetcher] = useState<Ketcher | null>(null);
    const [isFullScreen, setIsFullScreen] = useState(false);
    const [isLoading, setIsLoading] = useState(true);
    const [nameInput, setNameInput] = useState('');
    const [isResolving, setIsResolving] = useState(false);

    // Analysis Window State
    const [showAnalysis, setShowAnalysis] = useState(true);
    const [analysisData, setAnalysisData] = useState({
        formula: '',
        exactMass: '',
        molWt: '',
        mZ: '',
        elemental: ''
    });

    // Safety timeout for loading
    useEffect(() => {
        const timer = setTimeout(() => setIsLoading(false), 3000);
        return () => clearTimeout(timer);
    }, []);

    const handleInit = (ketcherInstance: Ketcher) => {
        setKetcher(ketcherInstance);
        setIsLoading(false);
        if (initialSmiles) {
            ketcherInstance.setMolecule(initialSmiles);
            updateAnalysis(initialSmiles);
        }

        // Subscribe to changes (mocking subscription by checking on mouse up or key up if possible, or just manual update for now)
        // Ketcher doesn't have a simple "onChange" prop easily accessible here without deeper integration.
        // We will update analysis on specific actions or periodically check? 
        // For now, let's update when "Get Name" or "Transfer" is clicked, or maybe add a "Calculate" button in Analysis window.
        // Better: Use an interval to check if SMILES changed? Or just rely on manual refresh for MVP.
    };

    const updateAnalysis = async (smiles: string) => {
        if (!smiles) return;
        try {
            const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'http://127.0.0.1:8000';
            const response = await fetch(`${apiUrl}/api/molecule/info`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ smiles })
            });
            const data = await response.json();
            if (data.formula) {
                setAnalysisData({
                    formula: data.formula,
                    exactMass: data.exact_mass?.toFixed(4) || '',
                    molWt: data.mol_wt?.toFixed(2) || '',
                    mZ: data.m_z?.toFixed(4) || '',
                    elemental: data.elemental_analysis || '',
                    logp: data.logp?.toFixed(2) || '',
                    tpsa: data.tpsa?.toFixed(2) || ''
                });
            }
        } catch (error) {
            console.error("Analysis failed", error);
        }
    };

    // Poll for changes (simple way to keep analysis live)
    useEffect(() => {
        const interval = setInterval(async () => {
            if (ketcher && showAnalysis) {
                try {
                    const currentSmiles = await ketcher.getSmiles();
                    // Only update if changed? For now just fetch. Optimization: compare with lastSmiles.
                    updateAnalysis(currentSmiles);
                } catch (e) { /* ignore empty structure errors */ }
            }
        }, 2000);
        return () => clearInterval(interval);
    }, [ketcher, showAnalysis]);

    const handleNameLoad = async () => {
        if (!nameInput.trim() || !ketcher) return;
        setIsResolving(true);
        try {
            const response = await fetch('http://localhost:8000/api/molecule/name_to_structure', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ name: nameInput })
            });
            const data = await response.json();
            if (data.smiles) {
                ketcher.setMolecule(data.smiles);
                updateAnalysis(data.smiles);
            } else {
                alert("Could not resolve name to structure.");
            }
        } catch (error) {
            console.error("Error resolving name:", error);
            alert("Error connecting to service.");
        } finally {
            setIsResolving(false);
        }
    };

    const handleStructureToName = async () => {
        if (!ketcher) return;
        try {
            const smiles = await ketcher.getSmiles();
            const res = await fetch('http://localhost:8000/api/molecule/structure_to_name', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ smiles })
            });
            const data = await res.json();
            if (data.name) {
                setNameInput(data.name);
            } else {
                alert("Could not generate name.");
            }
        } catch (e) {
            console.error(e);
            alert("Error generating name.");
        }
    };

    const handleSaveRequest = async () => {
        if (ketcher) {
            try {
                const smiles = await ketcher.getSmiles();
                onSave(smiles);
                onClose();
            } catch (error) {
                console.error("Error getting SMILES:", error);
            }
        }
    };

    const toggleFullScreen = () => setIsFullScreen(!isFullScreen);

    // Menu State
    const [activeMenu, setActiveMenu] = useState<string | null>(null);
    const fileInputRef = useRef<HTMLInputElement>(null);

    // Name <-> Structure State
    const [showNameModal, setShowNameModal] = useState(false);
    const [nameToStructInput, setNameToStructInput] = useState('');

    // NMR Modal State
    const [showNMR, setShowNMR] = useState(false);
    const [nmrData, setNmrData] = useState<any>(null);
    const [nmrType, setNmrType] = useState<'1H' | '13C'>('1H');

    const handleStructureAction = async (action: string) => {
        setActiveMenu(null);
        const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'http://127.0.0.1:8000';

        // Actions that require Ketcher
        if (action !== 'name_to_struct' && !ketcher) {
            return;
        }

        try {
            switch (action) {
                case 'cleanup':
                    await ketcher.layout();
                    break;
                case 'aromatize':
                    await ketcher.aromatize();
                    break;
                case 'dearomatize':
                    await ketcher.dearomatize();
                    break;
                case 'add_h':
                    break;
                case 'name_to_struct':
                    setShowNameModal(true);
                    break;
                case 'struct_to_name':
                    const s = await ketcher.getSmiles();
                    setIsLoading(true);
                    try {
                        const res = await fetch(`${apiUrl}/api/molecule/structure_to_name`, {
                            method: 'POST',
                            headers: { 'Content-Type': 'application/json' },
                            body: JSON.stringify({ smiles: s })
                        });
                        const data = await res.json();
                        if (data.name) {
                            alert(`IUPAC Name: ${data.name}`);
                            // Ideally, we could add this as a text object to the canvas
                            // ketcher.addText(data.name, {x: 0, y: 0}); // Hypothetical
                        } else {
                            alert("Could not generate name.");
                        }
                    } catch (e) {
                        console.error("Name generation failed", e);
                    } finally {
                        setIsLoading(false);
                    }
                    break;
                case 'predict_nmr':
                case 'predict_c13':
                    const smiles = await ketcher.getSmiles();
                    setIsLoading(true);
                    const endpoint = action === 'predict_nmr' ? 'predict_nmr' : 'predict_c13';
                    setNmrType(action === 'predict_nmr' ? '1H' : '13C');
                    try {
                        const res = await fetch(`${apiUrl}/api/molecule/${endpoint}`, {
                            method: 'POST',
                            headers: { 'Content-Type': 'application/json' },
                            body: JSON.stringify({ smiles })
                        });
                        const data = await res.json();
                        if (data.peaks) {
                            setNmrData(data);
                            setShowNMR(true);
                        } else {
                            console.error("No peaks data", data);
                            alert("Could not generate NMR prediction. Please try again.");
                        }
                    } catch (e) {
                        console.error("NMR Prediction failed", e);
                        alert("Error connecting to AI service.");
                    } finally {
                        setIsLoading(false);
                    }
                    break;
                case '3d':
                    const smiles3d = await ketcher.getSmiles();
                    setIsLoading(true);
                    try {
                        const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'http://127.0.0.1:8000';
                        const res = await fetch(`${apiUrl}/api/molecule/3d`, {
                            method: 'POST',
                            headers: { 'Content-Type': 'application/json' },
                            body: JSON.stringify({ smiles: smiles3d })
                        });
                        const data = await res.json();
                        if (data.molblock) {
                            ketcher.setMolecule(data.molblock);
                        }
                    } catch (e) {
                        console.error("3D generation failed", e);
                    } finally {
                        setIsLoading(false);
                    }
                    break;
            }
        } catch (e) {
            console.error("Structure action failed", e);
        }
    };

    // Close menus when clicking outside
    useEffect(() => {
        const handleClickOutside = () => setActiveMenu(null);
        window.addEventListener('click', handleClickOutside);
        return () => window.removeEventListener('click', handleClickOutside);
    }, []);

    const handleMenuClick = (e: React.MouseEvent, menu: string) => {
        e.stopPropagation();
        setActiveMenu(activeMenu === menu ? null : menu);
    };

    const handleFileAction = async (action: string) => {
        if (!ketcher) return;
        switch (action) {
            case 'new':
                ketcher.setMolecule('');
                break;
            case 'open':
                fileInputRef.current?.click();
                break;
            case 'save_mol':
                const mol = await ketcher.getMolfile();
                const blob = new Blob([mol], { type: 'chemical/x-mdl-molfile' });
                const url = URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = 'structure.mol';
                a.click();
                break;
            case 'exit':
                onClose();
                break;
        }
        setActiveMenu(null);
    };

    const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (file && ketcher) {
            const reader = new FileReader();
            reader.onload = (e) => {
                const content = e.target?.result as string;
                ketcher.setMolecule(content);
            };
            reader.readAsText(file);
        }
    };

    const handleViewAction = (action: string) => {
        switch (action) {
            case 'analysis':
                setShowAnalysis(!showAnalysis);
                break;
            case 'zoom_in':
                // Ketcher API might vary, mocking zoom or using internal if available
                // ketcher.setZoom(ketcher.getZoom() + 0.1); // Hypothetical
                break;
            case 'zoom_out':
                break;
        }
        setActiveMenu(null);
    };

    // Menu Bar Items
    const MenuItem = ({ label, id, children }: { label: string, id: string, children?: React.ReactNode }) => (
        <div className="relative">
            <div
                onClick={(e) => handleMenuClick(e, id)}
                className={`px-3 py-1 cursor-pointer hover:bg-[#cce8ff] text-slate-800 text-[11px] select-none ${activeMenu === id ? 'bg-[#cce8ff]' : ''}`}
            >
                {label}
            </div>
            {activeMenu === id && (
                <div className="absolute left-0 top-full bg-[#f0f0f0] border border-[#999] shadow-lg min-w-[150px] z-50 flex flex-col py-1">
                    {children}
                </div>
            )}
        </div>
    );

    const MenuAction = ({ label, onClick, shortcut }: { label: string, onClick: () => void, shortcut?: string }) => (
        <div
            onClick={(e) => { e.stopPropagation(); onClick(); }}
            className="px-4 py-1 hover:bg-[#0078d4] hover:text-white cursor-pointer flex justify-between items-center group"
        >
            <span>{label}</span>
            {shortcut && <span className="text-slate-500 text-[9px] group-hover:text-white ml-4">{shortcut}</span>}
        </div>
    );

    const MenuSeparator = () => <div className="h-[1px] bg-[#dcdcdc] my-1 mx-1"></div>;

    return (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40 backdrop-blur-[2px] p-4 font-sans">
            {/* Hidden File Input */}
            <input
                type="file"
                ref={fileInputRef}
                className="hidden"
                accept=".mol,.rxn,.smi,.sdf"
                onChange={handleFileUpload}
            />

            {/* Window Container - ChemDraw Style */}
            <div className={`bg-[#f0f0f0] shadow-[0_0_20px_rgba(0,0,0,0.3)] flex flex-col border border-[#666] transition-all duration-300 ${isFullScreen ? 'w-full h-full' : 'w-[95vw] h-[90vh]'}`}>

                {/* Title Bar */}
                <div className="flex justify-between items-center px-2 py-1 bg-white border-b border-[#a0a0a0] select-none h-8">
                    <div className="flex items-center space-x-2">
                        <div className="w-4 h-4 bg-[#009e0f] rounded-[2px] flex items-center justify-center text-[9px] text-white font-bold shadow-sm">Cd</div>
                        <h3 className="text-[11px] text-slate-900">ChemMind Professional - [Untitled Document 1]</h3>
                    </div>
                    <div className="flex items-center space-x-1">
                        <button onClick={toggleFullScreen} className="p-1 hover:bg-slate-200 rounded-sm" title={isFullScreen ? "Restore" : "Maximize"}>
                            {isFullScreen ? <Minimize2 className="h-3 w-3 text-slate-600" /> : <Maximize2 className="h-3 w-3 text-slate-600" />}
                        </button>
                        <button onClick={onClose} className="p-1 hover:bg-[#e81123] hover:text-white rounded-sm group" title="Close">
                            <X className="h-3 w-3 text-slate-600 group-hover:text-white" />
                        </button>
                    </div>
                </div>

                {/* Menu Bar */}
                <div className="flex items-center bg-[#f5f6f7] border-b border-[#dadbdc] px-1 h-6 relative">
                    <MenuItem label="File" id="file">
                        <MenuAction label="New Document" onClick={() => handleFileAction('new')} shortcut="Ctrl+N" />
                        <MenuAction label="Open..." onClick={() => handleFileAction('open')} shortcut="Ctrl+O" />
                        <MenuSeparator />
                        <MenuAction label="Save As Molfile..." onClick={() => handleFileAction('save_mol')} shortcut="Ctrl+S" />
                        <MenuSeparator />
                        <MenuAction label="Page Setup..." onClick={() => { }} />
                        <MenuAction label="Print..." onClick={() => { }} shortcut="Ctrl+P" />
                        <MenuSeparator />
                        <MenuAction label="Exit" onClick={() => handleFileAction('exit')} />
                    </MenuItem>
                    <MenuItem label="Edit" id="edit">
                        <MenuAction label="Undo" onClick={() => { }} shortcut="Ctrl+Z" />
                        <MenuAction label="Redo" onClick={() => { }} shortcut="Ctrl+Y" />
                        <MenuSeparator />
                        <MenuAction label="Cut" onClick={() => { }} shortcut="Ctrl+X" />
                        <MenuAction label="Copy" onClick={() => { }} shortcut="Ctrl+C" />
                        <MenuAction label="Paste" onClick={() => { }} shortcut="Ctrl+V" />
                        <MenuSeparator />
                        <MenuAction label="Clear" onClick={() => handleFileAction('new')} />
                    </MenuItem>
                    <MenuItem label="View" id="view">
                        <MenuAction label="Analysis Window" onClick={() => handleViewAction('analysis')} shortcut="F9" />
                        <MenuSeparator />
                        <MenuAction label="Zoom In" onClick={() => handleViewAction('zoom_in')} shortcut="Ctrl++" />
                        <MenuAction label="Zoom Out" onClick={() => handleViewAction('zoom_out')} shortcut="Ctrl+-" />
                    </MenuItem>
                    <MenuItem label="Object" id="object" />
                    <MenuItem label="Structure" id="structure">
                        <MenuAction label="Convert Name to Structure..." onClick={() => handleStructureAction('name_to_struct')} shortcut="Ctrl+Shift+N" />
                        <MenuAction label="Convert Structure to Name" onClick={() => handleStructureAction('struct_to_name')} shortcut="Ctrl+Shift+S" />
                        <MenuSeparator />
                        <MenuAction label="Clean Up Structure" onClick={() => handleStructureAction('cleanup')} shortcut="Ctrl+Shift+K" />
                        <MenuSeparator />
                        <MenuAction label="Check Aromaticity" onClick={() => handleStructureAction('aromatize')} />
                        <MenuAction label="Dearomatize" onClick={() => handleStructureAction('dearomatize')} />
                        <MenuSeparator />
                        <MenuAction label="Predict 1H NMR (AI)" onClick={() => handleStructureAction('predict_nmr')} />
                        <MenuAction label="Predict 13C NMR (AI)" onClick={() => handleStructureAction('predict_c13')} />
                        <MenuAction label="3D Optimization" onClick={() => handleStructureAction('3d')} />
                    </MenuItem>
                    <MenuItem label="Text" id="text" />
                    <MenuItem label="Curves" id="curves" />
                    <MenuItem label="Colors" id="colors" />
                    <MenuItem label="Search" id="search" />
                    <MenuItem label="Window" id="window" />
                    <MenuItem label="Help" id="help" />
                </div>

                {/* Toolbar / Quick Actions */}
                <div className="flex items-center bg-[#f5f6f7] border-b border-[#dadbdc] px-2 py-1 space-x-4 h-8">
                    {/* Name <-> Structure */}
                    <div className="flex items-center space-x-2 bg-white border border-[#bfbfbf] rounded-[2px] px-1 h-6">
                        <input
                            type="text"
                            value={nameInput}
                            onChange={(e) => setNameInput(e.target.value)}
                            onKeyDown={(e) => e.key === 'Enter' && handleNameLoad()}
                            className="text-[11px] outline-none w-48 text-slate-700 placeholder-slate-400 border-none h-full"
                            placeholder="Type name (e.g. Aspirin)..."
                        />
                        <div className="h-4 w-[1px] bg-slate-200 mx-1"></div>
                        <button onClick={handleNameLoad} disabled={isResolving} className="text-[10px] text-slate-600 hover:text-blue-600 px-1 font-medium">
                            {isResolving ? '...' : 'Load'}
                        </button>
                        <button onClick={handleStructureToName} className="text-[10px] text-slate-600 hover:text-blue-600 px-1 font-medium border-l border-slate-200 pl-2">
                            Get Name
                        </button>
                    </div>

                    <div className="flex-1"></div>

                    <button
                        onClick={() => setShowAnalysis(!showAnalysis)}
                        className={`text-[10px] px-2 py-0.5 border rounded-[2px] ${showAnalysis ? 'bg-[#cce8ff] border-[#99d1ff] text-[#005a9e]' : 'bg-white border-[#bfbfbf] text-slate-700'} hover:bg-[#e5f1fb]`}
                    >
                        Analysis Window
                    </button>

                    <button
                        onClick={handleSaveRequest}
                        className="text-[10px] font-semibold bg-[#0078d4] text-white px-3 py-0.5 rounded-[2px] hover:bg-[#106ebe] shadow-sm"
                    >
                        Transfer Structure
                    </button>
                </div>

                {/* Main Content Area - 3 Column Layout */}
                <div className="flex-1 bg-[#808080] relative p-[1px] overflow-hidden flex">

                    {/* Left Toolbar (Templates) */}
                    <div className="w-10 bg-[#f0f0f0] border-r border-[#808080] flex flex-col items-center py-1 space-y-1 z-10">
                        {/* Mock Template Icons */}
                        <div className="w-8 h-8 bg-white border border-[#bfbfbf] hover:bg-[#e5f1fb] hover:border-[#0078d4] cursor-pointer flex items-center justify-center text-[10px] text-slate-600" title="Cyclohexane Ring">⬡</div>
                        <div className="w-8 h-8 bg-white border border-[#bfbfbf] hover:bg-[#e5f1fb] hover:border-[#0078d4] cursor-pointer flex items-center justify-center text-[10px] text-slate-600" title="Benzene Ring">⌬</div>
                        <div className="w-8 h-8 bg-white border border-[#bfbfbf] hover:bg-[#e5f1fb] hover:border-[#0078d4] cursor-pointer flex items-center justify-center text-[10px] text-slate-600" title="Cyclopentane Ring">⬠</div>
                        <div className="w-8 h-1 bg-slate-300 my-1"></div>
                        <div className="w-8 h-8 bg-white border border-[#bfbfbf] hover:bg-[#e5f1fb] hover:border-[#0078d4] cursor-pointer flex items-center justify-center text-[10px] text-slate-600" title="Chain">〰</div>
                    </div>

                    {/* Center Canvas (Ketcher) */}
                    <div className="flex-1 bg-white relative ketcher-container shadow-inner overflow-hidden">
                        <Editor
                            staticResourcesUrl={process.env.PUBLIC_URL ? process.env.PUBLIC_URL + "/ketcher" : "/ketcher"}
                            structServiceProvider={ketcher as any}
                            onInit={handleInit}
                            errorHandler={(err) => console.error(err)}
                        />

                        {/* Floating Analysis Window (Now inside canvas area but floating) */}
                        {showAnalysis && (
                            <Draggable handle=".analysis-handle" bounds="parent">
                                <div className="absolute top-4 right-4 w-64 bg-[#f0f0f0] border border-[#666] shadow-[2px_2px_5px_rgba(0,0,0,0.2)] text-[11px] font-sans z-10">
                                    <div className="analysis-handle flex justify-between items-center px-1 py-[2px] bg-gradient-to-r from-[#0078d4] to-[#50a0e0] text-white cursor-move select-none">
                                        <span className="font-semibold ml-1">Analysis</span>
                                        <button onClick={() => setShowAnalysis(false)} className="hover:bg-red-500/50 rounded-sm p-[1px]">
                                            <X className="h-3 w-3 text-white" />
                                        </button>
                                    </div>
                                    <div className="p-2 bg-white text-[10px] space-y-1">
                                        <div className="flex justify-between"><span className="text-slate-500">Formula:</span> <span className="font-mono font-semibold">{analysisData.formula}</span></div>
                                        <div className="flex justify-between"><span className="text-slate-500">Exact Mass:</span> <span>{analysisData.exactMass}</span></div>
                                        <div className="flex justify-between"><span className="text-slate-500">Mol. Wt.:</span> <span>{analysisData.molWt}</span></div>
                                        <div className="flex justify-between"><span className="text-slate-500">m/z:</span> <span>{analysisData.mZ}</span></div>
                                        <div className="flex justify-between"><span className="text-slate-500">LogP:</span> <span>{analysisData.logp}</span></div>
                                        <div className="flex justify-between"><span className="text-slate-500">tPSA:</span> <span>{analysisData.tpsa}</span></div>
                                        <div className="border-t border-slate-200 pt-1 mt-1">
                                            <span className="text-slate-500 block mb-1">Elemental Analysis:</span>
                                            <span className="font-mono text-[9px] text-slate-700 break-words">{analysisData.elemental}</span>
                                        </div>
                                    </div>
                                </div>
                            </Draggable>
                        )}

                        {/* NMR Prediction Modal */}
                        {showNMR && nmrData && (
                            <Draggable handle=".nmr-handle" bounds="parent">
                                <div className="absolute top-10 left-10 w-80 bg-[#f0f0f0] border border-[#666] shadow-[2px_2px_10px_rgba(0,0,0,0.3)] text-[11px] font-sans z-30">
                                    <div className="nmr-handle flex justify-between items-center px-2 py-1 bg-gradient-to-r from-[#0078d4] to-[#50a0e0] text-white cursor-move select-none">
                                        <span className="font-semibold">{nmrType} NMR Prediction (AI)</span>
                                        <button onClick={() => setShowNMR(false)} className="hover:bg-red-500/50 rounded-sm p-[1px]">
                                            <X className="h-3 w-3 text-white" />
                                        </button>
                                    </div>
                                    <div className="p-2 bg-white max-h-64 overflow-y-auto">
                                        <table className="w-full text-left border-collapse">
                                            <thead>
                                                <tr className="border-b border-slate-300">
                                                    <th className="py-1 font-semibold text-slate-700">Shift (ppm)</th>
                                                    <th className="py-1 font-semibold text-slate-700">Mult.</th>
                                                    <th className="py-1 font-semibold text-slate-700">Count</th>
                                                    <th className="py-1 font-semibold text-slate-700">Assign.</th>
                                                </tr>
                                            </thead>
                                            <tbody>
                                                {nmrData.peaks?.map((peak: any, i: number) => (
                                                    <tr key={i} className="border-b border-slate-100 hover:bg-blue-50">
                                                        <td className="py-1 font-mono text-blue-600">{peak.shift}</td>
                                                        <td className="py-1">{peak.multiplicity}</td>
                                                        <td className="py-1">{peak.count}</td>
                                                        <td className="py-1 text-slate-500 truncate max-w-[80px]" title={peak.assignment}>{peak.assignment}</td>
                                                    </tr>
                                                ))}
                                            </tbody>
                                        </table>
                                        <div className="mt-2 text-[9px] text-slate-400 italic text-center">
                                            * Predictions generated by AI. Verify experimentally.
                                        </div>
                                    </div>
                                </div>
                            </Draggable>
                        )}

                        {/* Name to Structure Modal */}
                        {showNameModal && (
                            <div className="absolute inset-0 flex items-center justify-center bg-black/50 z-50">
                                <div className="bg-white p-4 rounded shadow-lg w-96">
                                    <h3 className="text-sm font-bold mb-2">Convert Name to Structure</h3>
                                    <input
                                        type="text"
                                        className="w-full border border-gray-300 p-1 mb-2 text-sm"
                                        placeholder="Enter chemical name (e.g., Aspirin)"
                                        value={nameToStructInput}
                                        onChange={(e) => setNameToStructInput(e.target.value)}
                                        onKeyDown={(e) => {
                                            if (e.key === 'Enter') {
                                                // Trigger conversion logic (needs to be implemented or extracted)
                                                // For now, let's just use a simple inline handler or call a function
                                                const convert = async () => {
                                                    setIsLoading(true);
                                                    try {
                                                        const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'http://127.0.0.1:8000';
                                                        const res = await fetch(`${apiUrl}/api/molecule/name_to_structure`, {
                                                            method: 'POST',
                                                            headers: { 'Content-Type': 'application/json' },
                                                            body: JSON.stringify({ name: nameToStructInput })
                                                        });
                                                        const data = await res.json();
                                                        if (data.smiles && ketcher) {
                                                            ketcher.setMolecule(data.smiles);
                                                            setShowNameModal(false);
                                                            setNameToStructInput('');
                                                        } else {
                                                            alert("Could not resolve name.");
                                                        }
                                                    } catch (e) {
                                                        console.error(e);
                                                        alert("Error converting name.");
                                                    } finally {
                                                        setIsLoading(false);
                                                    }
                                                };
                                                convert();
                                            }
                                        }}
                                    />
                                    <div className="flex justify-end space-x-2">
                                        <button onClick={() => setShowNameModal(false)} className="px-3 py-1 bg-gray-200 text-xs rounded hover:bg-gray-300">Cancel</button>
                                        <button
                                            onClick={async () => {
                                                setIsLoading(true);
                                                try {
                                                    const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'http://127.0.0.1:8000';
                                                    const res = await fetch(`${apiUrl}/api/molecule/name_to_structure`, {
                                                        method: 'POST',
                                                        headers: { 'Content-Type': 'application/json' },
                                                        body: JSON.stringify({ name: nameToStructInput })
                                                    });
                                                    const data = await res.json();
                                                    if (data.smiles && ketcher) {
                                                        ketcher.setMolecule(data.smiles);
                                                        setShowNameModal(false);
                                                        setNameToStructInput('');
                                                    } else {
                                                        alert("Could not resolve name.");
                                                    }
                                                } catch (e) {
                                                    console.error(e);
                                                    alert("Error converting name.");
                                                } finally {
                                                    setIsLoading(false);
                                                }
                                            }}
                                            className="px-3 py-1 bg-blue-600 text-white text-xs rounded hover:bg-blue-700"
                                        >
                                            Convert
                                        </button>
                                    </div>
                                </div>
                            </div>
                        )}

                        {isLoading && (
                            <div className="absolute inset-0 flex items-center justify-center bg-white/90 z-20">
                                <div className="flex flex-col items-center">
                                    <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-[#0078d4] mb-2"></div>
                                    <span className="text-sm text-slate-600 font-medium">Loading ChemMind Professional...</span>
                                </div>
                            </div>
                        )}
                    </div>

                    {/* Right Panel (Properties) */}
                    <div className="w-48 bg-[#f0f0f0] border-l border-[#808080] flex flex-col z-10">
                        <div className="px-2 py-1 bg-[#e1e1e1] border-b border-[#bfbfbf] text-[11px] font-semibold text-slate-700 flex justify-between items-center">
                            <span>Properties</span>
                            <ChevronDown className="h-3 w-3 text-slate-500" />
                        </div>
                        <div className="p-2 space-y-2 overflow-y-auto flex-1">
                            {/* Mock Properties */}
                            <div className="text-[10px] text-slate-600">
                                <div className="font-medium mb-1 text-slate-800">Atom</div>
                                <div className="grid grid-cols-2 gap-1">
                                    <span>Symbol:</span> <span className="bg-white border border-[#bfbfbf] px-1">C</span>
                                    <span>Charge:</span> <span className="bg-white border border-[#bfbfbf] px-1">0</span>
                                </div>
                            </div>
                            <div className="h-[1px] bg-[#bfbfbf]"></div>
                            <div className="text-[10px] text-slate-600">
                                <div className="font-medium mb-1 text-slate-800">Bond</div>
                                <div className="grid grid-cols-2 gap-1">
                                    <span>Type:</span> <span className="bg-white border border-[#bfbfbf] px-1">Single</span>
                                    <span>Stereo:</span> <span className="bg-white border border-[#bfbfbf] px-1">None</span>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>

                {/* Status Bar */}
                <div className="px-2 py-[2px] bg-[#0078d4] text-white text-[10px] flex justify-between items-center h-5">
                    <div className="flex space-x-4">
                        <span>Ready</span>
                        <span>x: 0.00 y: 0.00</span>
                    </div>
                    <div className="flex space-x-2">
                        <span>Zoom: 100%</span>
                        <span>Page 1 of 1</span>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default StructureEditor;

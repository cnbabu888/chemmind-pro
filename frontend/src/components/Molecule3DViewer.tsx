import React, { useEffect, useRef, useState } from 'react';

interface Molecule3DViewerProps {
    smiles: string;
}

declare global {
    interface Window {
        $3Dmol: any;
    }
}

const Molecule3DViewer: React.FC<Molecule3DViewerProps> = ({ smiles }) => {
    const viewerRef = useRef<HTMLDivElement>(null);
    const [molBlock, setMolBlock] = useState<string | null>(null);
    const [loading, setLoading] = useState(false);
    const [viewer, setViewer] = useState<any>(null);

    // Load 3Dmol.js script
    useEffect(() => {
        if (!window.$3Dmol) {
            const script = document.createElement('script');
            script.src = "https://3dmol.org/build/3Dmol-min.js";
            script.async = true;
            document.body.appendChild(script);
        }
    }, []);

    // Fetch 3D coordinates
    useEffect(() => {
        if (!smiles) return;

        const fetch3D = async () => {
            setLoading(true);
            try {
                const response = await fetch(`${process.env.NEXT_PUBLIC_API_URL || 'http://127.0.0.1:8000'}/api/molecule/3d?smiles=${encodeURIComponent(smiles)}`);
                if (!response.ok) throw new Error('Failed to fetch 3D structure');
                const data = await response.json();
                if (data.mol_block) {
                    setMolBlock(data.mol_block);
                }
            } catch (err) {
                console.error(err);
            } finally {
                setLoading(false);
            }
        };

        fetch3D();
    }, [smiles]);

    // Render 3D structure
    useEffect(() => {
        if (molBlock && viewerRef.current && window.$3Dmol) {
            // Initialize viewer if not exists
            let v = viewer;
            if (!v) {
                const config = { backgroundColor: 'white' };
                v = window.$3Dmol.createViewer(viewerRef.current, config);
                setViewer(v);
            }

            v.clear();
            v.addModel(molBlock, "mol");
            v.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.25 } });
            v.zoomTo();
            v.render();
        }
    }, [molBlock, viewer]);

    if (loading) return <div className="h-64 flex items-center justify-center bg-slate-50 text-slate-400 text-xs">Generating 3D Conformer...</div>;

    return (
        <div className="relative h-64 w-full bg-white rounded border border-slate-200 overflow-hidden">
            <div ref={viewerRef} className="w-full h-full" />
            {!window.$3Dmol && (
                <div className="absolute inset-0 flex items-center justify-center bg-white/80 text-xs text-slate-500">
                    Loading 3D Engine...
                </div>
            )}
            <div className="absolute bottom-2 right-2 text-[10px] text-slate-400">
                Powered by 3Dmol.js & RDKit
            </div>
        </div>
    );
};

export default Molecule3DViewer;

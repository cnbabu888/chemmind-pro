import React, { useCallback, useState } from 'react';
import ReactFlow, {
    MiniMap,
    Controls,
    Background,
    useNodesState,
    useEdgesState,
    addEdge,
    Node,
    Edge,
    Position,
} from 'reactflow';
import 'reactflow/dist/style.css';
import dagre from 'dagre';
import { X, ExternalLink, Beaker, FlaskConical } from 'lucide-react';

import SafetyPanel from './SafetyPanel';
import Molecule3DViewer from './Molecule3DViewer';

const nodeWidth = 180;
const nodeHeight = 50;

const getLayoutedElements = (nodes: Node[], edges: Edge[], direction = 'TB') => {
    const dagreGraph = new dagre.graphlib.Graph();
    dagreGraph.setDefaultEdgeLabel(() => ({}));

    const isHorizontal = direction === 'LR';
    dagreGraph.setGraph({ rankdir: direction });

    nodes.forEach((node) => {
        dagreGraph.setNode(node.id, { width: nodeWidth, height: nodeHeight });
    });

    edges.forEach((edge) => {
        dagreGraph.setEdge(edge.source, edge.target);
    });

    dagre.layout(dagreGraph);

    nodes.forEach((node) => {
        const nodeWithPosition = dagreGraph.node(node.id);
        node.targetPosition = isHorizontal ? Position.Left : Position.Top;
        node.sourcePosition = isHorizontal ? Position.Right : Position.Bottom;

        node.position = {
            x: nodeWithPosition.x - nodeWidth / 2,
            y: nodeWithPosition.y - nodeHeight / 2,
        };

        return node;
    });

    return { nodes, edges };
};

interface ReactionTreeProps {
    data: any;
}

const ReactionTree: React.FC<ReactionTreeProps> = ({ data }) => {
    const [selectedNode, setSelectedNode] = useState<any>(null);
    const [activeTab, setActiveTab] = useState<'info' | '3d' | 'safety'>('info');

    const initialNodes: Node[] = [];
    const initialEdges: Edge[] = [];

    let nodeIdCounter = 0;
    const processNode = (nodeData: any, parentId: string | null = null) => {
        const currentId = `node-${nodeIdCounter++}`;
        const isReaction = nodeData.is_reaction;

        // Custom styling based on type
        const style = isReaction
            ? { background: '#f0f9ff', border: '1px solid #0056b3', borderRadius: '50%', width: 60, height: 60, display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '10px' }
            : { background: '#fff', border: '1px solid #cbd5e1', borderRadius: '8px', padding: '8px', width: 180, fontSize: '12px' };

        const label = isReaction
            ? (nodeData.metadata?.reaction_name || 'Rxn')
            : (nodeData.smiles.length > 20 ? nodeData.smiles.substring(0, 18) + '...' : nodeData.smiles);

        initialNodes.push({
            id: currentId,
            data: { label: label, fullData: nodeData, type: isReaction ? 'reaction' : 'chemical', smiles: nodeData.smiles, metadata: nodeData.metadata, in_stock: nodeData.in_stock },
            position: { x: 0, y: 0 },
            style: style,
            type: isReaction ? 'default' : 'input',
        });

        if (parentId) {
            initialEdges.push({
                id: `edge-${parentId}-${currentId}`,
                source: parentId,
                target: currentId,
                animated: true,
                style: { stroke: '#94a3b8' },
            });
        }

        if (nodeData.children) {
            nodeData.children.forEach((child: any) => processNode(child, currentId));
        }
    };

    if (data && data.trees && data.trees.length > 0) {
        processNode(data.trees[0]);
    }

    const { nodes: layoutedNodes, edges: layoutedEdges } = getLayoutedElements(
        initialNodes,
        initialEdges
    );

    const [nodes, setNodes, onNodesChange] = useNodesState(layoutedNodes);
    const [edges, setEdges, onEdgesChange] = useEdgesState(layoutedEdges);

    const onConnect = useCallback(
        (params: any) => setEdges((eds) => addEdge(params, eds)),
        [setEdges]
    );

    const onNodeClick = (event: React.MouseEvent, node: Node) => {
        setSelectedNode(node.data);
        setActiveTab('info'); // Reset tab on new selection
    };

    return (
        <div className="relative w-full h-[600px] flex border border-slate-200 rounded-lg overflow-hidden bg-slate-50">
            <div className="flex-1 h-full relative">
                <ReactFlow
                    nodes={nodes}
                    edges={edges}
                    onNodesChange={onNodesChange}
                    onEdgesChange={onEdgesChange}
                    onConnect={onConnect}
                    onNodeClick={onNodeClick}
                    fitView
                    attributionPosition="bottom-left"
                >
                    <Controls />
                    <MiniMap />
                    <Background gap={12} size={1} color="#e2e8f0" />
                </ReactFlow>
            </div>

            {/* Side Details Panel */}
            {selectedNode && (
                <div className="w-96 bg-white border-l border-slate-200 shadow-xl overflow-y-auto flex flex-col absolute right-0 top-0 bottom-0 z-10">
                    <div className="p-4 border-b border-slate-100 flex justify-between items-center bg-slate-50">
                        <h3 className="font-semibold text-slate-800">
                            {selectedNode.type === 'chemical' ? 'Molecule Details' : 'Reaction Details'}
                        </h3>
                        <button onClick={() => setSelectedNode(null)} className="text-slate-400 hover:text-slate-600">
                            <X className="h-4 w-4" />
                        </button>
                    </div>

                    {selectedNode.type === 'chemical' ? (
                        <div className="flex flex-col h-full">
                            {/* Tabs */}
                            <div className="flex border-b border-slate-200">
                                <button
                                    onClick={() => setActiveTab('info')}
                                    className={`flex-1 py-2 text-xs font-medium ${activeTab === 'info' ? 'text-blue-600 border-b-2 border-blue-600' : 'text-slate-500 hover:text-slate-700'}`}
                                >
                                    Info
                                </button>
                                <button
                                    onClick={() => setActiveTab('3d')}
                                    className={`flex-1 py-2 text-xs font-medium ${activeTab === '3d' ? 'text-blue-600 border-b-2 border-blue-600' : 'text-slate-500 hover:text-slate-700'}`}
                                >
                                    3D View
                                </button>
                                <button
                                    onClick={() => setActiveTab('safety')}
                                    className={`flex-1 py-2 text-xs font-medium ${activeTab === 'safety' ? 'text-blue-600 border-b-2 border-blue-600' : 'text-slate-500 hover:text-slate-700'}`}
                                >
                                    Safety
                                </button>
                            </div>

                            <div className="p-4 flex-1 overflow-y-auto">
                                {activeTab === 'info' && (
                                    <div className="space-y-4">
                                        <div className="bg-white border border-slate-200 rounded p-2 flex justify-center">
                                            <img
                                                src={`${process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000'}/api/molecule/image?smiles=${encodeURIComponent(selectedNode.smiles)}`}
                                                alt="Structure"
                                                className="max-h-40 object-contain"
                                            />
                                        </div>
                                        <div>
                                            <label className="text-xs font-bold text-slate-400 uppercase">SMILES</label>
                                            <p className="text-xs font-mono bg-slate-50 p-2 rounded border border-slate-100 break-all mt-1">
                                                {selectedNode.smiles}
                                            </p>
                                        </div>
                                        <div className="flex space-x-2 pt-2">
                                            <button className="flex-1 bg-blue-600 text-white text-xs py-2 rounded hover:bg-blue-700 transition-colors">
                                                Use as Reactant
                                            </button>
                                            <button className="flex-1 bg-white border border-slate-300 text-slate-700 text-xs py-2 rounded hover:bg-slate-50 transition-colors">
                                                Search Similar
                                            </button>
                                        </div>
                                        {selectedNode.in_stock !== undefined && (
                                            <div className={`text-xs font-medium px-2 py-1 rounded inline-block ${selectedNode.in_stock ? 'bg-green-100 text-green-700' : 'bg-amber-100 text-amber-700'}`}>
                                                {selectedNode.in_stock ? 'Commercially Available' : 'Not In Stock'}
                                            </div>
                                        )}
                                    </div>
                                )}

                                {activeTab === '3d' && (
                                    <div className="space-y-4">
                                        <Molecule3DViewer smiles={selectedNode.smiles} />
                                        <p className="text-xs text-slate-500 text-center">
                                            Interactive 3D visualization generated by RDKit & 3Dmol.js
                                        </p>
                                    </div>
                                )}

                                {activeTab === 'safety' && (
                                    <SafetyPanel smiles={selectedNode.smiles} />
                                )}
                            </div>
                        </div>
                    ) : (
                        <div className="p-4 space-y-4">
                            {/* Reaction Details */}
                            <div className="bg-blue-50 p-3 rounded-lg border border-blue-100">
                                <div className="text-xs text-blue-600 font-medium uppercase mb-1">Reaction Name</div>
                                <div className="text-sm font-semibold text-slate-800">{selectedNode.metadata?.reaction_name || 'Unknown'}</div>
                            </div>
                            <div>
                                <label className="text-xs font-bold text-slate-400 uppercase">Confidence Score</label>
                                <div className="mt-1 w-full bg-slate-200 rounded-full h-2">
                                    <div
                                        className="bg-green-500 h-2 rounded-full"
                                        style={{ width: `${(selectedNode.metadata?.template_score || 0) * 100}%` }}
                                    ></div>
                                </div>
                                <p className="text-xs text-right text-slate-500 mt-1">
                                    {Math.round((selectedNode.metadata?.template_score || 0) * 100)}%
                                </p>
                            </div>
                        </div>
                    )}
                </div>
            )}
        </div>
    );
};

export default ReactionTree;

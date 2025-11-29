"use client";

import React, { useState, useRef } from 'react';
import jsPDF from 'jspdf';
import html2canvas from 'html2canvas';
import { Download } from 'lucide-react';

export default function PredictionPage() {
    const [reactants, setReactants] = useState('CC(=O)O.Nc1ccccc1'); // Acetic acid + Aniline
    const [result, setResult] = useState<any>(null);
    const [loading, setLoading] = useState(false);
    const [productImage, setProductImage] = useState<string | null>(null);
    const reportRef = useRef<HTMLDivElement>(null);

    const handlePredict = async () => {
        setLoading(true);
        setResult(null);
        setProductImage(null);
        try {
            // 1. Get Prediction
            const response = await fetch('http://localhost:8000/api/prediction', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ smiles: reactants }),
            });
            const data = await response.json();
            setResult(data);

            // 2. Get Image for Product
            if (data.product) {
                const imgResponse = await fetch(`http://localhost:8000/api/molecule/image?smiles=${encodeURIComponent(data.product)}`);
                if (imgResponse.ok) {
                    const blob = await imgResponse.blob();
                    setProductImage(URL.createObjectURL(blob));
                }
            }

        } catch (error) {
            console.error("Failed to fetch prediction", error);
            alert("Error fetching prediction. Ensure backend is running.");
        } finally {
            setLoading(false);
        }
    };

    const handleExportPDF = async () => {
        if (!reportRef.current) return;

        try {
            const canvas = await html2canvas(reportRef.current);
            const imgData = canvas.toDataURL('image/png');
            const pdf = new jsPDF('p', 'mm', 'a4');
            const pdfWidth = pdf.internal.pageSize.getWidth();
            const pdfHeight = (canvas.height * pdfWidth) / canvas.width;

            pdf.addImage(imgData, 'PNG', 0, 0, pdfWidth, pdfHeight);
            pdf.save('reaction_report.pdf');
        } catch (error) {
            console.error("Failed to export PDF", error);
            alert("Failed to generate PDF.");
        }
    };

    return (
        <main className="flex min-h-screen flex-col items-center p-10 bg-white dark:bg-zinc-900">
            <div className="w-full max-w-4xl flex justify-between items-center mb-8">
                <h1 className="text-4xl font-bold text-gray-800 dark:text-white">Forward Prediction</h1>
                {result && (
                    <button
                        onClick={handleExportPDF}
                        className="flex items-center gap-2 px-4 py-2 bg-blue-100 dark:bg-blue-900 text-blue-700 dark:text-blue-300 rounded-lg hover:bg-blue-200 dark:hover:bg-blue-800 transition-colors"
                    >
                        <Download className="w-4 h-4" />
                        Export Report
                    </button>
                )}
            </div>

            <div className="w-full max-w-2xl flex flex-col gap-4 mb-8">
                <label className="text-sm text-gray-500 dark:text-gray-400">
                    Enter Reactants (dot-separated SMILES):
                </label>
                <div className="flex gap-4">
                    <input
                        type="text"
                        value={reactants}
                        onChange={(e) => setReactants(e.target.value)}
                        placeholder="e.g. CC(=O)O.Nc1ccccc1"
                        className="flex-1 px-4 py-3 border rounded-xl shadow-sm focus:ring-2 focus:ring-green-500 outline-none dark:bg-zinc-800 dark:text-white"
                    />
                    <button
                        onClick={handlePredict}
                        disabled={loading}
                        className="px-6 py-3 bg-green-600 text-white font-semibold rounded-xl hover:bg-green-700 transition-colors disabled:opacity-50"
                    >
                        {loading ? 'Predicting...' : 'Run Reaction'}
                    </button>
                </div>
            </div>

            {result && (
                <div ref={reportRef} className="w-full max-w-4xl grid grid-cols-1 md:grid-cols-2 gap-8 bg-white dark:bg-zinc-900 p-4">
                    {/* Result Details */}
                    <div className="p-6 bg-gray-50 dark:bg-zinc-800 rounded-2xl border border-gray-200 dark:border-zinc-700">
                        <h2 className="text-xl font-semibold mb-4 text-gray-800 dark:text-white">Prediction Result</h2>
                        <div className="space-y-3">
                            <div>
                                <span className="text-sm font-medium text-gray-500">Reaction Type:</span>
                                <p className="text-lg text-gray-900 dark:text-gray-100">{result.reaction_type}</p>
                            </div>
                            <div>
                                <span className="text-sm font-medium text-gray-500">Confidence:</span>
                                <div className="flex items-center gap-2">
                                    <div className="w-full bg-gray-200 rounded-full h-2.5 dark:bg-gray-700 max-w-[200px]">
                                        <div className="bg-green-600 h-2.5 rounded-full" style={{ width: `${result.confidence * 100}%` }}></div>
                                    </div>
                                    <span className="text-sm text-gray-700 dark:text-gray-300">{(result.confidence * 100).toFixed(0)}%</span>
                                </div>
                            </div>
                            <div>
                                <span className="text-sm font-medium text-gray-500">Product SMILES:</span>
                                <p className="font-mono text-sm bg-white dark:bg-zinc-900 p-2 rounded border border-gray-200 dark:border-zinc-700 break-all">
                                    {result.product}
                                </p>
                            </div>
                        </div>
                    </div>

                    {/* Product Visualization */}
                    <div className="flex items-center justify-center p-6 bg-white dark:bg-zinc-800 rounded-2xl border border-gray-200 dark:border-zinc-700 h-[400px]">
                        {productImage ? (
                            <img src={productImage} alt="Predicted Product" className="max-w-full max-h-full object-contain" />
                        ) : (
                            <span className="text-gray-400">Product structure</span>
                        )}
                    </div>
                </div>
            )}
        </main>
    );
}

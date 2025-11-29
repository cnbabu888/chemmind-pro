import React, { useState } from 'react';
import { Brain, Send, Sparkles, Bot } from 'lucide-react';

const BrainAnalysis: React.FC = () => {
    const [prompt, setPrompt] = useState('');
    const [provider, setProvider] = useState<'gemini' | 'openai'>('gemini');
    const [response, setResponse] = useState<string | null>(null);
    const [loading, setLoading] = useState(false);

    const handleAnalyze = async () => {
        if (!prompt.trim()) return;

        setLoading(true);
        setResponse(null);
        try {
            const res = await fetch(`${process.env.NEXT_PUBLIC_API_URL || 'http://127.0.0.1:8000'}/api/analyze`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ prompt, provider })
            });
            const data = await res.json();
            setResponse(data.response);
        } catch (error) {
            console.error(error);
            setResponse("Error analyzing prompt. Please try again.");
        } finally {
            setLoading(false);
        }
    };

    return (
        <div className="bg-white rounded-lg shadow-sm border border-slate-200 overflow-hidden flex flex-col h-full">
            <div className="px-4 py-3 border-b border-slate-100 bg-slate-50 flex justify-between items-center">
                <h3 className="font-semibold text-slate-700 flex items-center text-sm">
                    <Brain className="h-4 w-4 mr-2 text-purple-600" /> AI Research Assistant
                </h3>
                <div className="flex bg-white rounded-md border border-slate-200 p-0.5">
                    <button
                        onClick={() => setProvider('gemini')}
                        className={`px-3 py-1 text-xs font-medium rounded-sm transition-colors ${provider === 'gemini' ? 'bg-purple-100 text-purple-700' : 'text-slate-500 hover:text-slate-700'}`}
                    >
                        Gemini Pro
                    </button>
                    <button
                        onClick={() => setProvider('openai')}
                        className={`px-3 py-1 text-xs font-medium rounded-sm transition-colors ${provider === 'openai' ? 'bg-green-100 text-green-700' : 'text-slate-500 hover:text-slate-700'}`}
                    >
                        ChatGPT-4
                    </button>
                </div>
            </div>

            <div className="flex-1 p-4 flex flex-col overflow-hidden">
                <div className="flex-1 overflow-y-auto mb-4 bg-slate-50 rounded border border-slate-100 p-3 text-sm text-slate-700 whitespace-pre-wrap">
                    {response ? (
                        <div className="animate-in fade-in duration-300">
                            <div className="flex items-center mb-2 text-xs font-semibold text-slate-500">
                                {provider === 'gemini' ? <Sparkles className="h-3 w-3 mr-1" /> : <Bot className="h-3 w-3 mr-1" />}
                                {provider === 'gemini' ? 'Gemini Response' : 'ChatGPT Response'}
                            </div>
                            {response}
                        </div>
                    ) : (
                        <div className="h-full flex flex-col items-center justify-center text-slate-400">
                            <Brain className="h-8 w-8 mb-2 opacity-20" />
                            <p className="text-center text-xs">
                                Ask complex questions about reaction mechanisms, <br />
                                retrosynthesis strategies, or spectral analysis.
                            </p>
                        </div>
                    )}
                </div>

                <div className="relative">
                    <textarea
                        value={prompt}
                        onChange={(e) => setPrompt(e.target.value)}
                        placeholder="Ask the AI assistant (e.g., 'Explain the mechanism of Fischer Esterification')..."
                        className="w-full pl-3 pr-12 py-2 text-sm border border-slate-300 rounded-md focus:outline-none focus:ring-2 focus:ring-purple-500 focus:border-transparent resize-none h-20"
                        onKeyDown={(e) => {
                            if (e.key === 'Enter' && !e.shiftKey) {
                                e.preventDefault();
                                handleAnalyze();
                            }
                        }}
                    />
                    <button
                        onClick={handleAnalyze}
                        disabled={loading || !prompt.trim()}
                        className="absolute bottom-2 right-2 p-1.5 bg-purple-600 text-white rounded hover:bg-purple-700 disabled:opacity-50 transition-colors"
                    >
                        {loading ? <div className="h-4 w-4 border-2 border-white/30 border-t-white rounded-full animate-spin" /> : <Send className="h-4 w-4" />}
                    </button>
                </div>
            </div>
        </div>
    );
};

export default BrainAnalysis;

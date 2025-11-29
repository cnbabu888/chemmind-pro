import React, { useState, useRef, useEffect } from 'react';
import { Send, Sparkles, Bot, Zap, Search, ChevronUp, ChevronDown, X, User } from 'lucide-react';

interface AIChatBarProps {
    onAnalyze?: (prompt: string, provider: string, deepSearch: boolean) => void;
}

interface Message {
    role: 'user' | 'assistant';
    content: string;
    provider?: string;
    deepSearch?: boolean;
}

const AIChatBar: React.FC<AIChatBarProps> = () => {
    const [messages, setMessages] = useState<Message[]>([]);
    const [prompt, setPrompt] = useState('');
    const [provider, setProvider] = useState<'gemini' | 'openai' | 'claude' | 'grok'>('gemini');
    const [deepSearch, setDeepSearch] = useState(false);
    const [loading, setLoading] = useState(false);
    const [isExpanded, setIsExpanded] = useState(false);
    const messagesEndRef = useRef<HTMLDivElement>(null);

    const scrollToBottom = () => {
        messagesEndRef.current?.scrollIntoView({ behavior: "smooth" });
    };

    useEffect(() => {
        scrollToBottom();
    }, [messages, isExpanded]);

    const handleAnalyze = async () => {
        if (!prompt.trim()) return;

        const userMessage: Message = { role: 'user', content: prompt };
        setMessages(prev => [...prev, userMessage]);
        setPrompt('');
        setLoading(true);
        setIsExpanded(true);

        try {
            const res = await fetch(`${process.env.NEXT_PUBLIC_API_URL || 'http://127.0.0.1:8000'}/api/analyze`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ prompt: userMessage.content, provider, deep_search: deepSearch })
            });
            const data = await res.json();

            const aiMessage: Message = {
                role: 'assistant',
                content: data.response,
                provider,
                deepSearch
            };
            setMessages(prev => [...prev, aiMessage]);
        } catch (error) {
            console.error(error);
            const errorMessage: Message = { role: 'assistant', content: "Error analyzing prompt. Please try again.", provider };
            setMessages(prev => [...prev, errorMessage]);
        } finally {
            setLoading(false);
        }
    };

    const getProviderIcon = (p?: string) => {
        switch (p) {
            case 'gemini': return <Sparkles className="h-4 w-4 text-blue-500" />;
            case 'openai': return <Bot className="h-4 w-4 text-green-500" />;
            case 'claude': return <Zap className="h-4 w-4 text-orange-500" />;
            case 'grok': return <Search className="h-4 w-4 text-slate-800" />;
            default: return <Sparkles className="h-4 w-4" />;
        }
    };

    return (
        <div className={`fixed bottom-0 left-0 right-0 bg-white border-t border-slate-200 shadow-[0_-4px_6px_-1px_rgba(0,0,0,0.1)] transition-all duration-300 z-50 flex flex-col ${isExpanded ? 'h-[60vh]' : 'h-auto'}`}>

            {/* Chat Area (Collapsible) */}
            {isExpanded && (
                <div className="flex-1 overflow-y-auto p-4 bg-slate-50">
                    <div className="max-w-4xl mx-auto space-y-6">
                        <div className="flex justify-between items-center mb-2 sticky top-0 bg-slate-50 z-10 py-2 border-b border-slate-200">
                            <h3 className="text-xs font-semibold text-slate-500 uppercase tracking-wider">Conversation History</h3>
                            <button onClick={() => setIsExpanded(false)} className="text-slate-400 hover:text-slate-600">
                                <ChevronDown className="h-5 w-5" />
                            </button>
                        </div>

                        {messages.map((msg, idx) => (
                            <div key={idx} className={`flex ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`}>
                                <div className={`flex max-w-[80%] ${msg.role === 'user' ? 'flex-row-reverse' : 'flex-row'} gap-3`}>

                                    {/* Avatar/Icon */}
                                    <div className={`flex-shrink-0 h-8 w-8 rounded-full flex items-center justify-center ${msg.role === 'user' ? 'bg-slate-200' : 'bg-white border border-slate-200'}`}>
                                        {msg.role === 'user' ? <User className="h-5 w-5 text-slate-500" /> : getProviderIcon(msg.provider)}
                                    </div>

                                    {/* Message Bubble */}
                                    <div className={`p-4 rounded-2xl text-sm whitespace-pre-wrap shadow-sm ${msg.role === 'user'
                                            ? 'bg-blue-600 text-white rounded-tr-none'
                                            : 'bg-white border border-slate-200 text-slate-800 rounded-tl-none'
                                        }`}>
                                        {msg.role === 'assistant' && msg.deepSearch && (
                                            <div className="mb-2">
                                                <span className="text-[10px] bg-blue-50 text-blue-600 px-2 py-0.5 rounded-full font-medium border border-blue-100">Deep Search Analysis</span>
                                            </div>
                                        )}
                                        {msg.content}
                                    </div>
                                </div>
                            </div>
                        ))}

                        {loading && (
                            <div className="flex justify-start">
                                <div className="flex max-w-[80%] gap-3">
                                    <div className="flex-shrink-0 h-8 w-8 rounded-full bg-white border border-slate-200 flex items-center justify-center">
                                        {getProviderIcon(provider)}
                                    </div>
                                    <div className="bg-white border border-slate-200 p-4 rounded-2xl rounded-tl-none shadow-sm flex items-center space-x-2">
                                        <div className="w-2 h-2 bg-slate-400 rounded-full animate-bounce" style={{ animationDelay: '0ms' }} />
                                        <div className="w-2 h-2 bg-slate-400 rounded-full animate-bounce" style={{ animationDelay: '150ms' }} />
                                        <div className="w-2 h-2 bg-slate-400 rounded-full animate-bounce" style={{ animationDelay: '300ms' }} />
                                    </div>
                                </div>
                            </div>
                        )}
                        <div ref={messagesEndRef} />
                    </div>
                </div>
            )}

            {/* Input Bar */}
            <div className="bg-white p-4 border-t border-slate-100">
                <div className="max-w-4xl mx-auto flex items-center gap-4">

                    {/* Expand Toggle (if messages exist but collapsed) */}
                    {!isExpanded && messages.length > 0 && (
                        <button onClick={() => setIsExpanded(true)} className="p-2 text-slate-500 hover:bg-slate-100 rounded-full relative">
                            <ChevronUp className="h-5 w-5" />
                            <span className="absolute top-0 right-0 h-2 w-2 bg-red-500 rounded-full"></span>
                        </button>
                    )}

                    {/* Model Selector */}
                    <div className="relative group">
                        <select
                            value={provider}
                            onChange={(e) => setProvider(e.target.value as any)}
                            className="appearance-none bg-slate-50 border border-slate-300 text-slate-700 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block w-full p-2.5 pl-9 pr-8 cursor-pointer"
                        >
                            <option value="gemini">Gemini</option>
                            <option value="openai">ChatGPT</option>
                            <option value="claude">Claude</option>
                            <option value="grok">Grok</option>
                        </select>
                        <div className="absolute inset-y-0 left-0 flex items-center pl-2.5 pointer-events-none">
                            {getProviderIcon(provider)}
                        </div>
                    </div>

                    {/* Deep Search Toggle */}
                    <label className="flex items-center cursor-pointer">
                        <input
                            type="checkbox"
                            checked={deepSearch}
                            onChange={(e) => setDeepSearch(e.target.checked)}
                            className="sr-only peer"
                        />
                        <div className="relative w-9 h-5 bg-slate-200 peer-focus:outline-none peer-focus:ring-4 peer-focus:ring-blue-300 rounded-full peer peer-checked:after:translate-x-full rtl:peer-checked:after:-translate-x-full peer-checked:after:border-white after:content-[''] after:absolute after:top-[2px] after:start-[2px] after:bg-white after:border-gray-300 after:border after:rounded-full after:h-4 after:w-4 after:transition-all peer-checked:bg-blue-600"></div>
                        <span className="ms-2 text-sm font-medium text-slate-700">Deep Search</span>
                    </label>

                    {/* Input Field */}
                    <div className="flex-1 relative">
                        <input
                            type="text"
                            value={prompt}
                            onChange={(e) => setPrompt(e.target.value)}
                            onKeyDown={(e) => e.key === 'Enter' && handleAnalyze()}
                            placeholder="Ask ChemMind..."
                            className="block w-full p-3 pl-4 text-sm text-slate-900 border border-slate-300 rounded-lg bg-slate-50 focus:ring-blue-500 focus:border-blue-500 shadow-inner"
                        />
                        <button
                            onClick={handleAnalyze}
                            disabled={!prompt.trim() || loading}
                            className="absolute right-2 bottom-1.5 bg-blue-600 hover:bg-blue-700 text-white font-medium rounded-lg text-sm px-3 py-1.5 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
                        >
                            <Send className="h-4 w-4" />
                        </button>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default AIChatBar;

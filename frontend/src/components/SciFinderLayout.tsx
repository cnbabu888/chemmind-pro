import React, { useState } from 'react';
import { Search, FlaskConical, History, Menu, Bell, User, Settings } from 'lucide-react';

interface SciFinderLayoutProps {
    children: React.ReactNode;
}

const SciFinderLayout: React.FC<SciFinderLayoutProps> = ({ children }) => {
    const [isSidebarOpen, setIsSidebarOpen] = useState(true);

    return (
        <div className="min-h-screen flex flex-col bg-slate-50 text-slate-900 font-sans">
            {/* Top Navigation Bar */}
            <header className="bg-[#0056b3] text-white shadow-md z-50">
                <div className="max-w-7xl mx-auto px-4 h-16 flex items-center justify-between">
                    <div className="flex items-center space-x-6">
                        <div className="flex items-center space-x-2 cursor-pointer">
                            <FlaskConical className="h-8 w-8" />
                            <span className="text-xl font-bold tracking-tight">ChemMind</span>
                        </div>
                        <nav className="hidden md:flex space-x-1">
                            <button className="px-3 py-2 rounded hover:bg-white/10 text-sm font-medium transition-colors">Explore</button>
                            <button className="px-3 py-2 rounded bg-white/20 text-sm font-medium transition-colors">Substances</button>
                            <button className="px-3 py-2 rounded hover:bg-white/10 text-sm font-medium transition-colors">Reactions</button>
                            <button className="px-3 py-2 rounded hover:bg-white/10 text-sm font-medium transition-colors">Retrosynthesis</button>
                        </nav>
                    </div>

                    <div className="flex items-center space-x-4">
                        <div className="relative">
                            <input
                                type="text"
                                placeholder="Quick Search..."
                                className="bg-white/10 border border-white/20 rounded-full px-4 py-1.5 text-sm text-white placeholder-white/60 focus:outline-none focus:bg-white/20 w-64 transition-all"
                            />
                            <Search className="h-4 w-4 absolute right-3 top-2 text-white/60" />
                        </div>
                        <button className="p-2 hover:bg-white/10 rounded-full"><Bell className="h-5 w-5" /></button>
                        <button className="p-2 hover:bg-white/10 rounded-full"><Settings className="h-5 w-5" /></button>
                        <div className="h-8 w-8 rounded-full bg-white/20 flex items-center justify-center cursor-pointer hover:bg-white/30">
                            <User className="h-5 w-5" />
                        </div>
                    </div>
                </div>
            </header>

            <div className="flex flex-1 overflow-hidden">
                {/* Sidebar */}
                <aside
                    className={`bg-white border-r border-slate-200 transition-all duration-300 ease-in-out flex flex-col ${isSidebarOpen ? 'w-64' : 'w-16'}`}
                >
                    <div className="p-4 border-b border-slate-100 flex justify-between items-center">
                        {isSidebarOpen && <span className="font-semibold text-slate-700 text-sm">History</span>}
                        <button onClick={() => setIsSidebarOpen(!isSidebarOpen)} className="p-1.5 hover:bg-slate-100 rounded text-slate-500">
                            <Menu className="h-4 w-4" />
                        </button>
                    </div>

                    <div className="flex-1 overflow-y-auto p-2 space-y-1">
                        {/* Mock History Items */}
                        {[1, 2, 3].map((i) => (
                            <div key={i} className="group flex items-center space-x-3 p-2 hover:bg-slate-50 rounded cursor-pointer border border-transparent hover:border-slate-200">
                                <History className="h-4 w-4 text-slate-400 group-hover:text-[#0056b3]" />
                                {isSidebarOpen && (
                                    <div className="flex flex-col">
                                        <span className="text-sm text-slate-700 font-medium truncate w-40">Reaction #{100 + i}</span>
                                        <span className="text-xs text-slate-400">2 mins ago</span>
                                    </div>
                                )}
                            </div>
                        ))}
                    </div>
                </aside>

                {/* Main Content Area */}
                <main className="flex-1 overflow-y-auto bg-[#f8fafc] p-6">
                    <div className="max-w-7xl mx-auto">
                        {children}
                    </div>
                </main>
            </div>
        </div>
    );
};

export default SciFinderLayout;

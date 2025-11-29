"use client";

import React from 'react';
import ChatWindow from '../../components/ChatWindow';

export default function ChatPage() {
    return (
        <main className="flex min-h-screen flex-col items-center p-10 bg-white dark:bg-zinc-900">
            <h1 className="text-4xl font-bold mb-2 text-gray-800 dark:text-white">AI Chat Mode</h1>
            <p className="text-gray-500 dark:text-gray-400 mb-8">Your conversational chemistry expert.</p>

            <ChatWindow />
        </main>
    );
}

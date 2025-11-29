import type { Metadata } from "next";
import { Inter } from "next/font/google";
import "./globals.css";
import "./ketcher.css";
import Sidebar from "../components/Sidebar";
import ThemeToggle from "../components/ThemeToggle";

const inter = Inter({ subsets: ["latin"] });

export const metadata: Metadata = {
    title: "ChemMind - AI Chemistry Assistant",
    description: "Next-generation AI-powered chemistry intelligence platform",
};

export default function RootLayout({
    children,
}: Readonly<{
    children: React.ReactNode;
}>) {
    return (
        <html lang="en">
            <body className={`${inter.className} bg-gray-50 dark:bg-black`}>
                <div className="flex">
                    <Sidebar />
                    <div className="flex-1 ml-64">
                        <ThemeToggle />
                        {children}
                    </div>
                </div>
            </body>
        </html>
    );
}

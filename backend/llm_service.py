import os
import google.generativeai as genai
from typing import Optional
try:
    from openai import OpenAI
except ImportError:
    OpenAI = None

class LLMService:
    def __init__(self):
        # Gemini Setup
        self.gemini_key = os.environ.get("GEMINI_API_KEY")
        self.gemini_model = None
        if self.gemini_key:
            genai.configure(api_key=self.gemini_key)
            self.gemini_model = genai.GenerativeModel('gemini-pro-latest')
        else:
            print("WARNING: GEMINI_API_KEY not found.")

        # OpenAI Setup
        self.openai_key = os.environ.get("OPENAI_API_KEY")
        self.openai_client = None
        if self.openai_key and OpenAI:
            self.openai_client = OpenAI(api_key=self.openai_key)
        elif not OpenAI:
            print("WARNING: openai package not installed.")
        else:
            print("WARNING: OPENAI_API_KEY not found.")

    def generate_response(self, prompt: str, provider: str = "gemini", deep_search: bool = False) -> Optional[str]:
        # ChemMind Persona & System Prompt
        system_instruction = """
        You are ChemMind, an advanced AI research assistant for chemists.
        
        CRITICAL OUTPUT RULES:
        1. IF asked about synthesis, reaction conditions, or yields:
           - You MUST provide a DIRECT, numbered list of the best known conditions first.
           - Format: "1. [Reagents], [Solvent], [Temp/Time] -> [Yield] (Ref: [Citation])"
           - Do NOT start with "Here is how you make..." or conversational filler.
           - List top 1-3 methods max, unless asked for more.
        
        2. IF asked for "details", "mechanism", "explain", OR if "Deep Search" is active:
           - Provide the list of conditions first (as above).
           - THEN provide a detailed explanation, mechanism, and context.
           
        3. IF asked about properties or general questions:
           - Be concise and data-driven.
        """

        # Construct final prompt based on mode
        if deep_search:
            final_prompt = f"{system_instruction}\n\nUSER REQUEST (DEEP SEARCH ACTIVE): {prompt}\n\nProvide a comprehensive analysis including mechanism and context after the conditions list."
        else:
            final_prompt = f"{system_instruction}\n\nUSER REQUEST: {prompt}\n\nPrioritize direct data/conditions. Keep explanations brief unless asked."

        if provider == "openai":
            return self._generate_openai(final_prompt)
        elif provider == "claude":
            return self._generate_claude(final_prompt)
        elif provider == "grok":
            return self._generate_grok(final_prompt)
        return self._generate_gemini(final_prompt)

    def _generate_claude(self, prompt: str) -> str:
        return "Claude integration coming soon! (Placeholder response)"

    def _generate_grok(self, prompt: str) -> str:
        return "Grok integration coming soon! (Placeholder response)"

    def _generate_gemini(self, prompt: str) -> Optional[str]:
        if not self.gemini_model:
            return "Error: Gemini API key not configured."
        try:
            response = self.gemini_model.generate_content(prompt)
            return response.text
        except Exception as e:
            print(f"Gemini Error: {e}")
            return f"Error calling Gemini: {str(e)}"

    def _generate_openai(self, prompt: str) -> Optional[str]:
        if not self.openai_client:
            return "Error: OpenAI API key not configured or package missing."
        try:
            response = self.openai_client.chat.completions.create(
                model="gpt-4",
                messages=[{"role": "user", "content": prompt}]
            )
            return response.choices[0].message.content
        except Exception as e:
            print(f"OpenAI Error: {e}")
            return f"Error calling OpenAI: {str(e)}"

# Global instance
llm_service = LLMService()

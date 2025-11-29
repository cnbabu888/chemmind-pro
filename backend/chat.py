import random
from typing import List, Dict

class ChatEngine:
    def __init__(self):
        pass

    def get_response(self, message: str, history: List[Dict[str, str]]) -> str:
        """
        Generates a response to the user's message.
        """
        message = message.lower()
        
        if "mechanism" in message:
            return "The mechanism likely involves nucleophilic attack of the amine on the carbonyl carbon, followed by elimination of water to form the amide bond."
        
        if "yield" in message:
            return "Based on similar reactions in the literature, you can expect a yield between 75% and 85% under standard conditions."
        
        if "safety" in message:
            return "Warning: Aniline is toxic by inhalation and skin contact. Acetic anhydride is corrosive. Ensure proper ventilation and use PPE."
        
        if "hello" in message or "hi" in message:
            return "Hello! I am your AI Chemistry Assistant. How can I help you with your synthesis today?"

        return "That's an interesting question. As an AI, I can help you analyze reaction pathways, predict products, or suggest optimization conditions. Could you be more specific?"

# Global instance
chat_engine = ChatEngine()

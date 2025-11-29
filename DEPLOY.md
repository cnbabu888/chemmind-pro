# Deployment Instructions for ChemMind

ChemMind is containerized using Docker and Docker Compose for easy deployment.

## Prerequisites
-   [Docker](https://docs.docker.com/get-docker/) installed.
-   [Docker Compose](https://docs.docker.com/compose/install/) installed.
-   A valid [Google Gemini API Key](https://ai.google.dev/).

## Setup & Run

1.  **Navigate to the project root**:
    ```bash
    cd chemmind
    ```

2.  **Set your API Key**:
    You can export it as an environment variable:
    ```bash
    export GEMINI_API_KEY="your_actual_api_key_here"
    ```
    *Alternatively, you can create a `.env` file in the root directory:*
    ```env
    GEMINI_API_KEY=your_actual_api_key_here
    ```

3.  **Build and Run**:
    ```bash
    docker-compose up --build
    ```
    *This may take a few minutes the first time as it builds the frontend and backend images.*

4.  **Access the Application**:
    Open your browser and navigate to:
    [http://localhost:3000](http://localhost:3000)

## Architecture
-   **Frontend**: Next.js app running on port `3000`.
-   **Backend**: FastAPI app running on port `8000`.
-   **Communication**: The frontend communicates with the backend via `http://localhost:8000` (browser-side).

## Troubleshooting
-   **"AI service unavailable"**: Ensure `GEMINI_API_KEY` is set correctly in the backend container. Check logs with `docker-compose logs backend`.
-   **CORS Errors**: If accessing from a different machine, update `NEXT_PUBLIC_API_URL` in `docker-compose.yml` and CORS settings in `backend/main.py`.

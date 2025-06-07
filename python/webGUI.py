# Backend: FastAPI server
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
import numpy as np
import plotly.graph_objects as go
import json

app = FastAPI()

# Enable CORS for frontend communication
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/plot-data")
def get_plot_data():
    # Simulating computational data (replace this with your actual logic)
    x = np.linspace(0, 10, 100)
    y = np.sin(x)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name='Sine Wave'))
    
    return json.loads(fig.to_json())  # Send Plotly JSON to frontend

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)

# Frontend: React app with Plotly
# import React, { useEffect, useState } from "react";
# import Plot from "react-plotly.js";
# import axios from "axios";

# const PlotComponent = () => {
#     const [plotData, setPlotData] = useState(null);

#     useEffect(() => {
#         axios.get("http://localhost:8000/plot-data")
#             .then(response => setPlotData(response.data))
#             .catch(error => console.error("Error fetching plot data:", error));
#     }, []);

#     return (
#         <div>
#             <h1>Plotly Graph</h1>
#             {plotData && <Plot data={plotData.data} layout={plotData.layout} />}
#         </div>
#     );
# };

# export default PlotComponent;

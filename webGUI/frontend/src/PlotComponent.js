import React, { useEffect, useState } from "react";
import Plot from "react-plotly.js";
import axios from "axios";

const PlotComponent = () => {
    const [plotData, setPlotData] = useState(null);

    useEffect(() => {
        axios.get("http://localhost:8000/plot-data")
            .then(response => setPlotData(response.data))
            .catch(error => console.error("Error fetching plot data:", error));
    }, []);

    return (
        <div>
            <h1>Plotly Graph</h1>
            {plotData && <Plot data={plotData.data} layout={plotData.layout} />}
        </div>
    );
};

export default PlotComponent;

import {
    processIvDataInput
} from "./ivDataInput.js"

let inData = null;

function main() {

    let socket = io();

    let voltage, current, workFunction, voltageMult, currentMult, workFunctionMult, result;

    loadInitEventListeners();

    const ctx = document.getElementById('myChart').getContext('2d');

    const indicatorOptions = {

        radius: 4,
        borderWidth: 1,
        borderColor: 'red',
        backgroundColor: 'transparent'

    };

    const getLabelAndValue = Chart.controllers.line.prototype.getLabelAndValue;

    Chart.controllers.line.prototype.getLabelAndValue = function (index) {

        if (index === -1) {

            const meta = this.getMeta();
            const pt = meta._pt;
            const vScale = meta.vScale;
            const xScale = meta.xScale;

            return {

                label: '',
                value: "x: " + Math.round(xScale.getValueForPixel(pt.x) * 10000) / 10000 +
                    ", y: " + Math.round(vScale.getValueForPixel(pt.y) * 10000) / 10000

            };

        }

        return getLabelAndValue.call(this, index);

    }

    // The interaction mode
    Chart.Interaction.modes.interpolate = function (chart, e, option) {

        const x = e.x;
        const items = [];
        const metas = chart.getSortedVisibleDatasetMetas();

        for (let i = 0; i < metas.length; i++) {

            const meta = metas[i];
            const pt = meta.dataset.interpolate({
                x
            }, "x");

            if (pt) {

                const element = new Chart.elements.PointElement({
                    ...pt,
                    options: {
                        ...indicatorOptions
                    }
                });

                meta._pt = element;

                items.push({

                    element,
                    index: -1,
                    datasetIndex: meta.index

                });

            } else {

                meta._pt = null;

            }
        }

        return items;

    };

    // Plugin to draw the indicators
    Chart.register({

        id: 'indicators',

        afterDraw(chart) {

            const metas = chart.getSortedVisibleDatasetMetas();

            for (let i = 0; i < metas.length; i++) {

                const meta = metas[i];

                if (meta._pt) {

                    meta._pt.draw(chart.ctx);

                }

            }

        },

        afterEvent(chart, args) {

            if (args.event.type === 'mouseout') {

                const metas = chart.getSortedVisibleDatasetMetas();

                for (let i = 0; i < metas.length; i++) {

                    metas[i]._pt = null;

                }

                args.changed = true;

            }

        }
    })

    const data = [{
            x: 0.10207040421489234,
            y: 5617.000000000002
        }, {
            x: 0.10559087321846734,
            y: 4214.878587327415
        }, {
            x: 0.10911134222204238,
            y: 3164.722642398102
        },
        {
            x: 0.1126318112256174,
            y: 2377.5417614456496
        }, {
            x: 0.11615228022919243,
            y: 1787.0286368310246
        }, {
            x: 0.11967274923276744,
            y: 1343.742827965096
        }, {
            x: 0.21120494332571796,
            y: 0.8305729972325626
        }
    ]

    let regressionData;

    let myChart = new Chart(ctx, {

        data: {

            datasets: [

                {
                    label: "Input data points",
                    fill: false,
                    lineTension: 0.1,
                    borderColor: "rgb(145, 64, 2)",
                    backgroundColor: "rgb(145, 64, 2)",
                    type: 'scatter',
                    data: regressionData,
                    interpolate: false,
                    borderWidth: 3

                },

                {
                    label: "Fitted line",
                    fill: false,
                    lineTension: 0,
                    borderColor: "rgb(0, 79, 32)",
                    backgroundColor: "rgb(0, 79, 32)",
                    data: data,
                    type: 'scatter',
                    showLine: true,
                    pointRadius: 0,
                    pointHoverRadius: 0,
                    interpolate: true,
                    borderWidth: 5

                }
            ]
        },

        options: {

            fill: false,
            lineTension: 0.1,

            interaction: {
                intersect: false,
                mode: 'interpolate',
                axis: "x"

            },

            scales: {
                y: {

                    title: {

                        display: true,
                        text: "Current",
                        font:{
                            size: 18
                        }

                    },

                    type: "logarithmic",
                    position: "bottom",

                    ticks: {

                        callback: function (value, index, ticks) {

                            if (value === 1000000) {

                                return "1 [MA]"

                            }

                            if (value === 100000) {

                                return "100 [kA]"

                            }

                            if (value === 10000) {

                                return "10 [kA]"

                            }

                            if (value === 1000) {

                                return "1 [kA]"

                            }

                            if (value === 100) {

                                return "100 [A]"

                            }

                            if (value === 10) {

                                return "10 [A]"

                            }

                            if (value === 1) {

                                return "1 [A]"

                            }

                            if (value === 0.1) {

                                return "100 [mA]"

                            }

                            if (value === 0.01) {

                                return "10 [mA]"

                            }

                            if (value === 0.001) {

                                return "1 [mA]"

                            }

                            return null
                        }
                    }

                },

                x: {

                    title: {

                        text: "1 / Local Field, (V/nm)⁻¹]",
                        display: true,
                        font:{
                            size: 18
                        }
                    },

                    type: "linear",
                    position: "bottom",
                    ticks: {

                        callback: function (value, index, ticks) {

                            if (String(value).length > 4) {
                                value = Math.round(value * 100) / 100;
                            }

                            return value;
                        }

                    }

                },
            },

            plugins: {

                tooltip: {

                    mode: 'interpolate',
                    intersect: false,
                    enabled: true

                },

                title: {

                    display: true,
                    text: "Fitted parameters: Radius: 50.000 nm, β: 0.017306 nm⁻¹, σAeff: 3.4500 nm²",
                    font: {
                        size: 24
                    }

                },

                zoom: {

                    zoom: {

                        wheel: {
                            enabled: true
                        },
                        pinch: {
                            enabled: false
                        },
                        drag: {
                            enabled: true,
                            modifierKey: 'ctrl'
                        },
                        mode: 'xy'

                    },

                    pan: {

                        enabled: true,
                        mode: 'xy'

                    },

                    limits: {

                        x: {
                            min: 0.1,
                            max: 1
                        },
                        y: {
                            min: -100,
                            max: 5e4
                        }

                    }
                },

                hover: {

                    intercept: false

                }
            }
        },

        tooltips: {

            mode: "interpolate",
            intersect: true

        },

        plugins: [
            //tooltipLine
        ]
    });

    function checkValidity() {

        if (errorDivs.length > 0) {

            errorDivs.forEach(div => {
                div.remove();
            })

            errorDivs = [];

        }

        voltage = document.getElementById("voltage_in").value;
        current = document.getElementById("current_in").value;
        workFunction = document.getElementById("work_function_in").value;

        voltageMult = document.getElementById("voltage_mult_in").value;
        currentMult = document.getElementById("current_mult_in").value;
        workFunctionMult = document.getElementById("work_function_mult_in").value;


        if (voltage == "") voltage = "2.413e+02, 2.511e+02, 2.622e+02, 2.706e+02, 2.803e+02, 2.915e+02, 2.999e+02, 3.096e+02, 3.208e+02, 3.305e+02, 3.403e+02, 3.515e+02, 3.612e+02, 3.710e+02, 3.808e+02, 3.891e+02, 4.003e+02, 4.100e+02, 4.184e+02, 4.338e+02, 4.491e+02, 4.644e+02, 4.826e+02, 4.993e+02";
        if (current == "") current = "8.719e-01, 1.582e+00, 2.967e+00, 5.038e+00, 8.555e+00, 1.406e+01, 2.309e+01, 3.670e+01, 5.643e+01, 8.678e+01, 1.249e+02, 1.797e+02, 2.502e+02, 3.371e+02, 4.540e+02, 6.116e+02, 7.459e+02, 9.720e+02, 1.267e+03, 1.764e+03, 2.376e+03, 3.096e+03, 4.310e+03, 5.617e+03";
        if (workFunction == "") workFunction = "4.5";

        let _voltage = processIvDataInput(voltage);
        let _current = processIvDataInput(current);

        let _workFunction = workFunction;

        for (let i = 0; i < _voltage.length; i++) {
            _voltage[i] = _voltage[i] * voltageMult;
        }

        for (let i = 0; i < _current.length; i++) {
            _current[i] = _current[i] * currentMult;
        }

        _workFunction = _workFunction * workFunctionMult;

        let canCompute = true;

        if (_voltage.length != _current.length) {

            if (_voltage.length > _current.length) {

                raiseInputError("2000");
                canCompute = false;

            } else if (_current.length > _voltage.length) {

                canCompute = false;
                raiseInputError("2001");

            }


        } else if ((_current.length < 3 || _voltage.length < 3) && (canCompute == true)) {

            canCompute = false;
            raiseInputError("2005");


        } else {

            if (canCompute == true) {

                _voltage.forEach(voltage => {

                    if ((voltage < 0 || voltage > 100000) && (canCompute == true)) {

                        canCompute = false;
                        raiseInputError("2002");

                    }


                });

                _current.forEach(current => {

                    if ((current < 0 || current > 10000) && (canCompute == true)) {

                        canCompute = false;
                        raiseInputError("2003");

                    };

                });

                if ((parseFloat(_workFunction) < 0.5 || parseFloat(_workFunction) > 10) && (canCompute == true)) {

                    raiseInputError("2004");
                    canCompute = false;

                };

                let data = [_voltage, _current, _workFunction];

                if (canCompute) {
                    socket.emit('calculateIv', data);
                    $('#loadingModal').modal('show');
                }


            }

        }
    }

    socket.on('calculatedData', (data) => {

        updateGraph(data);
        $('#loadingModal').modal('hide');
        console.log(data);

    })

    socket.on('logServerSideError', (arg) => {

        alert("Got an error from server while computing data. Check console for more info.");
        console.log("Error server-side: ");
        console.log(arg);

    })

    function updateGraph(dict) {

        const rad = Math.round(dict.Radius * 1000) / 1000;
        const beta = Math.round(dict.beta * 100000) / 100000;
        const sigmaAeff = Math.round(dict.sigma_Aeff * 1000) / 1000;

        const xData = dict.xplot_line;
        const yData = dict.yplot_line;

        const xReg = dict.xplot_mrk;
        const yReg = dict.yplot_mrk;

        updateTitle();
        updateAxes();
        updatePoints();
        updateRegressionLine();

        function updateRegressionLine() {

            const regressionPoints = createPoints(xReg, yReg);

            let dataSet = myChart.data.datasets[0];

            regressionPoints.forEach(point => {
                dataSet.data.push(point);
            })

            myChart.update();

        }

        function updateTitle() {

            if (inVoltageMode) {
                myChart.options.plugins.title.text = "Fitted parameters: Radius: " + rad + " nm, β: " + beta + " nm⁻¹, σAeff: " + sigmaAeff + " nm²";
                return;
            } else {
                myChart.options.plugins.title.text = "Fitted parameters: Radius: " + rad + " nm, β: " + beta + ", σAeff: " + sigmaAeff + " nm²";
            }

        }

        function updateAxes() {

        }

        function updatePoints() {

            myChart.data.labels.pop();

            myChart.data.datasets.forEach((dataset) => {

                dataset.data = null;

            });

            myChart.update();

            let dataSet = myChart.data.datasets[1];

            const points = createPoints(xData, yData)

            points.forEach((point) => {

                dataSet.data.push(point);

            })

            myChart.update();
        }

        function createPoints(xS, yS) {

            let res = [];

            for (let i = 0; i < xS.length; i++) {

                res.push({
                    x: xS[i],
                    y: yS[i]
                });

            }

            return res;
        }

    }

    function loadInitEventListeners(){

        let enterButton = document.getElementById("enterButton");
        helpButton = document.getElementById("helpIcon");


        helpButton.addEventListener("click", addHelpDiv);
        enterButton.addEventListener("click", checkValidity);
    
    }

}

let errorDivs = [];
let inVoltageMode = true;
let helpButton, docDiv;

main();
loadEventListeners();

function downloadData(){

    if(inData != null){ 

       download(inData, "data", "json")

    } else{

        raiseInputError("2008");

    }

    function download(data, filename) {

        var file = new Blob([JSON.stringify(data)], {type : 'application/json'});

        if (window.navigator.msSaveOrOpenBlob) // IE10+

            window.navigator.msSaveOrOpenBlob(file, filename);

        else { // Others

            var a = document.createElement("a"),
                    url = URL.createObjectURL(file);
                    
            a.href = url;
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            
            setTimeout(function() {
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);  
            }, 0); 
        }
    }
    
}

function loadEventListeners() {

    let voltageUnitDiv = document.getElementById("voltage_mult_in");

    const opt1 = ["[V]", "[kV]"];
    const opt2 = ["[V/m]", "[kV/m]"];

    document.getElementById('voltageSelectDiv').addEventListener('change', function () {

        if (voltageSelectDiv.value == 1) {

            voltageUnitDiv.options[0].text = opt1[0];
            voltageUnitDiv.options[1].text = opt1[1];
            inVoltageMode = true;

        }

        if (voltageSelectDiv.value == 2) {

            voltageUnitDiv.options[0].text = opt2[0];
            voltageUnitDiv.options[1].text = opt2[1];
            inVoltageMode = false;

        }


    });

    let canvasButton = document.getElementById("canvasButton");

    canvasButton.addEventListener('click', function(){
        downloadData();
    })

}

export function raiseInputError(id) {

    switch (id) {

        case "2000":
            addErrorDiv("Can not create graph as in voltage data there are more points than in current data");
            break;
        case "2001":
            addErrorDiv("Can not create graph as in current data there are more points than in voltage data");
            break;
        case "2002":
            addErrorDiv("One of the voltage values is out of bounds 0 < x < 100'000 V");
            break;
        case "2003":
            addErrorDiv("One of the current values is out of bounds 0 < x < 10'000 A");
            break;
        case "2004":
            addErrorDiv("Work function value is out of bounds 0.5 < x < 10 eV");
            break;
        case "2005":
            addErrorDiv("One must enter at least 3 points for voltage and current data");
            break;
        case "2006":
            addErrorDiv("One of the input lines has no separator between values! Check console for more info");
            break;
        case "2007":
            addErrorDiv("One of the input lines has a data of unknown type! Check console for more info");
            break;
        case "2008":
            addErrorDiv("Can not download data file as there is no data to download. (Did you click on Enter button?)");
            break;
        default:
            addErrorDiv("Unknown error");
    }

    function addErrorDiv(message) {

        console.log(message);

        const template = `
        <section class="container showcase alert alert-danger alert-dismissible my-auto animated bounceInLeft">
            <strong class="mx-2">Error!</strong> ${message}!
            <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
        </section>
        `;

        let relativeDiv = document.getElementById("myChart");
        let errorDiv = document.createElement("div");
        errorDiv.innerHTML = template;

        let section = errorDiv.children[0];

        errorDivs.push(section)

        insertAfter(section, relativeDiv);

    }

}


function addHelpDiv(){

    const title = "How to use:"
    const message = "Find general tips and info below";

    let helpButton = document.getElementById("helpIcon");

    const template1 = `
    <section class="alert alert-success alert-dismissible d-flex align-items-center fade show animated bounceInLeft">
        <strong class="mx-2">${title}</strong> ${message}!
        <button type="button" class="btn-close" data-bs-dismiss="alert" id="howToUseBtn"></button>
    </section>    
    
    `;

    const template2 = `
    <section>
    <h2>What does this program do?</h2><br>
<p>
    This program is a web interface for GETELEC that can be used to predict the current created by an emitter.
</p><br>
    <h2>Controls</h2><br>
<p>
    Zoom in and out with <code>mouse wheel</code>, or hold <code>CTRL</code> and press <code>left click</code>
    to select the area to zoom into with a mouse. 
</p><br>
    <h2>Data input</h2><br>
<p>
    Start by selecting appropriate units. Program allows user to choose between <code>Voltage</code>
    and <code>Applied Field</code>, that could be measured in Volts or kiloVolts. The same logic applies
    to <code>Current</code> and <code>Work Function</code> fields. Please note that <code>Work Function</code>
    should be measured in electronvolts.
</p><br>
<p>
    The minimum amount of voltage-current input points for the calculations to work is 3. The upper and lower
    bounds for the values are the following:
    <ul>
    <li><code>Voltage</code>: 0 - 100 kV </li>
    <li><code>Applied Field</code>: 0 - 100 kV/m </li>
    <li><code>Current</code>: 0 - 10 kA </li>
    <li><code>Work Function</code>: 0.5 - 10 eV </li>
    
    </ul>
</p><br>
<p> The program supports different formats of input, but the recommended way is to use 
    <code>comma separated</code> values,
    as in the following example: </p>

    <ul><li>291.5, 299.9, 309.6, 320.8</li></ul>

<p> One could enter the values using the <code>standard view</code>: </p> 

    <ul><li>2.915e+02, 2.999e+02, 3.096e+02, 3.208e+02</li></ul>

    If program does not understand your input, a message will appear on the screen.
    </p><br>
    <h2>Data output</h2><br>
<p>
    Save the chart <code>.png</code> by <code>right click</code> with mouse and then press <code>save image as</code>.

    In order to save the raw computed data, press the <code>download button</code>.

    The saved file consists of following JSON content:
    <ul>
    <li><code>xplot_mrk</code>: List of x values for input data</li>
    <li><code>yplot_mrk</code>: List of y values for input data</li>
    <li><code>xplot_line</code>: List of x values for fitted line</li>
    <li><code>yplot_line</code>: List of y values for fitted line</li>
    <li><code>β</code>: Field enhancement factor coefficient</li>
    <li><code>Radius</code>: Local radius of the emitter tip [nm]</li>
    <li><code>sigma_Aeff</code>: Effective surface area of Emitter [nm²]</li>
    <li><code>xAxisUnit</code>: Units of the x axis, default 1 / (Local Field [V/nm]</li>
    <li><code>yAxisUnit</code>: Units of the y axis, default Current [Amps]</li>
    </ul>
</p><br>
    <h2>Data processing</h2><br>
<p>
    Upon successful submitance of data, a graph should appear. <code>Hover</code> over with <code>mouse</code> to see point coordinates.
    Exact set of point data can be seen in <code>console</code> in developer tools (<code>F12</code>) .
    All the calculations are done on the server. This page simply displays generated values. For more info on
    how the back-end works, check out the <code>source code</code> and full <code>documentation</code>
    <a href="https://github.com/AndKyr/GETELEC">here</a>.
</p>




</section>
    `

    let relativeDiv = document.getElementById("canvasButton");
    let errorDiv = document.createElement("section");
    docDiv = document.createElement("section");
    docDiv.innerHTML = template2;

    errorDiv.innerHTML = template1;

    let section1 = errorDiv.children[0];
    let section2 = docDiv;

    insertAfter(section1, relativeDiv);
    insertAfter(section2, section1);


    helpButton.style.display = "none";

    document.getElementById("howToUseBtn").addEventListener("click", addHelpIcon);
    
}

function addHelpIcon(){

    helpButton.style.display = "block";
    docDiv.innerHTML = "";

}

function insertAfter(newNode, existingNode) {

    existingNode.parentNode.insertBefore(newNode, existingNode.nextSibling);

}
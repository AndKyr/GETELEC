let inData = null;
let graphData = null;

let errorDivs = [];
let inVoltageMode = true;
let helpButton, docDiv;

main();


function main() {

    //Connection with node.js server

    let socket = io();

    let voltage, current, workFunction, voltageMult, currentMult, workFunctionMult, result;

    loadEventListeners();

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
                    //", y: " + Math.round(vScale.getValueForPixel(pt.y) * 10000) / 10000
                    ", y: " + vScale.getValueForPixel(pt.y).toExponential(4)

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
    

    //Dummy data for visuals
    const data = [{x: 0.09747321398292647, y: 0.00005653999999999998}, {x: 0.10166524343024846, y: 0.0000354075119705509}, {x: 0.10585727287757045, y: 0.000022046283752204346}, {x: 0.11004930232489243, y: 0.000013643221928263145}, {x: 0.11424133177221442, y: 0.000008388666961697198}, {x: 0.11843336121953642, y: 0.000005122791752667756}, {x: 0.12262539066685842, y: 0.000003106062691860675}, {x: 0.1268174201141804, y: 0.0000018691917259290598}]

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

                            if (value === 0.0001) {

                                return "0.1 [mA]"

                            }

                            if (value === 0.00001) {

                                return "0.01 [mA]"

                            }

                            if (value === 0.000001) {

                                return "1000 [nA]"

                            } 

                            if (value === 0.0000001) {

                                return "100 [nA]"

                            }

                            if (value === 0.00000001) {

                                return "10 [nA]"

                            }

                            if (value === 0.000000001) {

                                return "1 [nA]"

                            }

                            if (value === 0.0000000001) {

                                return "100 [pA]"

                            }

                            if (value === 0.00000000001) {

                                return "10 [pA]"

                            }

                            if (value === 0.000000000001) {

                                return "1 [pA]"

                            }

                            if (value === 0.0000000000001) {

                                return "100 [fA]"

                            }

                            if (value === 0.00000000000001) {

                                return "10 [fA]"

                            }

                            if (value === 0.000000000000001) {

                                return "1 [fA]"

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
                            min: 0,
                            max: 3
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


    //General function that check for input validity and then sends packet to server
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


        if (voltage == "") voltage = "6666.66666667, 6528.49740933, 6395.93908629, 6268.65671642, 6146.34146341, 6028.70813397, 5915.49295775, 5806.4516129, 5701.35746606, 5600, 5502.18340611, 5407.72532189, 5316.4556962, 5228.21576763, 5142.85714286, 5060.24096386, 4980.23715415, 4902.72373541, 4827.5862069 , 4754.71698113, 4684.01486989, 4615.38461538, 4548.73646209, 4483.98576512, 4421.05263158, 4359.8615917 , 4300.34129693, 4242.42424242, 4186.04651163, 4131.14754098, 4077.66990291, 4025.55910543, 3925.23364486, 3876.92307692, 3829.78723404, 3783.78378378, 3738.87240356, 3695.01466276, 3652.17391304, 3610.31518625, 3569.40509915, 3529.41176471, 3490.30470914, 3452.05479452, 3414.63414634, 3378.01608579, 3342.17506631, 3307.08661417, 3272.72727273, 3239.07455013, 3206.10687023, 3173.80352645, 3142.1446384 , 3111.11111111, 3080.68459658, 3050.84745763, 3021.58273381, 2992.87410926, 2964.70588235, 2937.06293706, 2909.93071594, 2883.29519451, 2857.14285714";
        if (current == "") current = "5.65354016e-05, 4.28464283e-05, 3.24541251e-05, 2.45683389e-05, 1.85875533e-05, 1.40539879e-05, 1.06193664e-05, 8.01882052e-06, 6.05098032e-06, 4.56287757e-06, 3.43825950e-06, 2.58892706e-06, 1.94793392e-06, 1.46451644e-06, 1.10020423e-06, 8.25856631e-07, 6.19415049e-07, 4.64192603e-07, 3.47574447e-07, 2.60030670e-07, 1.94366991e-07, 1.45156296e-07, 1.08307623e-07, 8.07395611e-08, 6.01329743e-08, 4.47437594e-08, 3.32614033e-08, 2.47019853e-08, 1.83274191e-08, 1.35845100e-08, 1.00590041e-08, 7.44096972e-09, 4.05929377e-09, 2.99355710e-09, 2.20530977e-09, 1.62290148e-09, 1.19302882e-09, 8.76074733e-10, 6.42626597e-10, 4.70867981e-10, 3.44633799e-10, 2.51960496e-10, 1.83999483e-10, 1.34216889e-10, 9.77911792e-11, 7.11688901e-11, 5.17338317e-11, 3.75620267e-11, 2.72401384e-11, 1.97311118e-11, 1.42748585e-11, 1.03149269e-11, 7.44441693e-12, 5.36613434e-12, 3.86326821e-12, 2.77783297e-12, 1.99486455e-12, 1.43077440e-12, 1.02488643e-12, 7.33204336e-13, 5.23857148e-13, 3.73798141e-13, 2.66375053e-13";
        if (workFunction == "") workFunction = "4.5";

        let _voltage = processDataInput(voltage);
        let _current = processDataInput(current);

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

                addErrorDiv("Can not create graph as in voltage data there are more points than in current data");
                canCompute = false;

            } else if (_current.length > _voltage.length) {

                canCompute = false;
                addErrorDiv("Can not create graph as in current data there are more points than in voltage data");

            }


        } else if ((_current.length < 3 || _voltage.length < 3) && (canCompute == true)) {

            canCompute = false;
            addErrorDiv("One must enter at least 3 points for voltage and current data");


        } else {

            if (canCompute == true) {

                _voltage.forEach(voltage => {

                    if ((voltage < 0 || voltage > 100000) && (canCompute == true)) {

                        canCompute = false;
                        addErrorDiv("One of the voltage values is out of bounds 0 < x < 100'000 V");

                    }


                });

                _current.forEach(current => {

                    if ((current < 0 || current > 10000) && (canCompute == true)) {

                        canCompute = false;
                        addErrorDiv("One of the current values is out of bounds 0 < x < 10'000 A");

                    };

                });

                if ((parseFloat(_workFunction) < 0.5 || parseFloat(_workFunction) > 10) && (canCompute == true)) {

                    addErrorDiv("Work function value is out of bounds 0.5 < x < 10 eV");
                    canCompute = false;

                };

                let calculateR = document.getElementById("chooseR").checked;
                if(calculateR == true) calculateR = 1;
                if(calculateR == false) calculateR = 0;

                inData = [_voltage, _current, _workFunction, calculateR];

                if (canCompute) {
                    socket.emit('calculateIv', inData);
                    $('#loadingModal').modal('show');
                }


            }

        }
    };

    //Event - on receive data from server

    socket.on('calculatedData', (data) => {

        graphData = data;
        updateGraph(data);
        $('#loadingModal').modal('hide');
        console.log(data);

    });

    //Event - error message from server received

    socket.on('logServerSideError', (arg) => {

        $('#loadingModal').modal('hide');

        alert("Got an error from server while computing data. Check console for more info.");
        console.log("Error server-side: ");
        console.log(arg);

    });

    //Function that draws data on existing graphs

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
                if(isNaN(rad) == false) myChart.options.plugins.title.text = "Fitted parameters: Radius: " + dict.Radius.toExponential(4) + " nm, β (F/V): " + dict.beta.toExponential(4) + " nm⁻¹, σAeff: " + dict.sigma_Aeff.toExponential(6) + " nm²";
                if(isNaN(rad) == true) myChart.options.plugins.title.text = "Fitted parameters: β (F/V): " + dict.beta.toExponential(4) + " nm⁻¹, σAeff: " + dict.sigma_Aeff.toExponential(6) + " nm²";
            } else {
                myChart.options.plugins.title.text = "Fitted parameters: Radius: " + dict.Radius.toExponential(4) + " nm, β (F/V): " + beta.toExponential(4) + ", σAeff: " + dict.sigma_Aeff.toExponential(6) + " nm²";
            }
            
            myChart.options.plugins.title.text += dict.orthodoxyMessage

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

    function loadEventListeners(){

        //click event on the enter data button

        let enterButton = document.getElementById("enterButton");

        //click event to get help text

        let helpButton = document.getElementById("helpIcon");

        try {helpButton.addEventListener("click", addHelpDiv);} catch (e) {console.log(e)};
        try {enterButton.addEventListener("click", checkValidity);} catch (e) {console.log(e)};


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
    
        let downloadButton = document.getElementById('downloadDataBtn');
        downloadButton.addEventListener("click", downloadData);
    
    }

}

//
// Below are functions used to display messages/convert data
//

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

    let relativeDiv = document.getElementById("downloadButton");
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

const typeDict = {
    0: "String",
    1: "Array",
    99: "Unknown"
}

const separatorDict = {
    0: ".",
    1: ",",
    2: " ",
    99: "Unknown"
}

const multDict = {
    0: "*",
    1: "e",
    99: "Unknown"
}

const commaDict = {
    0: ",",
    1: ".",
    99: "Unknown"
}

function processDataInput(input) {

    let sp = "";
    let cm = "";
    let data = new Input(input);

    if (data.type == 0) data.processString();
    if (data.type == 1) data.processArray();
    if (data.type == 99) data.processUnknownInput();

    switch (data.separator) {
        case 0:
            sp = ".";
            break;

        case 1:
            sp = ",";
            break;

        case 2:
            sp = " ";
            break;

        default:
            //console.log("Separator error");
            break;
    }

    switch (data.comma) {
        case 0:
            cm = ",";
            break;

        case 1:
            cm = ".";
            break;

        default:
            //console.log("Comma error");
            break;
    }

    let _dsb;
    let dsb = [];

    if(data.type == 0){
        _dsb = data.data;

    } else {
        _dsb = beautifyResult(data.data.split(sp));
    }

    if (cm == ",") {
        for (let i = 0; i < _dsb.length; i++) {
            _dsb[i] = _dsb[i].replace(cm, ".");
        }
    }

    if(_dsb.separator == 99){

        dsb.push(parseFloat(_dsb));

    } else {

        _dsb = _dsb.split(separatorDict[data.separator]);
        beautifyResult(_dsb);

        _dsb.forEach(el => {
            dsb.push(parseFloat(el));
        });
        
    }

    return dsb;
}

class Input {

    constructor(input) {

        this.data = input;
        this.type = getTypeOfInput(input);
        this.separator = getSeparator(input);
        this.mult = getMult(input);
        this.comma = getComma(input);

    }

    processString() {

        if (this.separator == 99) {
            // addErrorDiv("One of the input lines has no separator between values! Check console for more info");
            // console.log("The following line has no separator: ");
            // console.log(this);
            // return;
        }

        if (this.comma == 99) {
            /*            console.log("The following line has no commas: ");
                        console.log(this);
                        return;*/
        }

        //this.data = beautifyResult(this.data.split(separatorDict.separator));
        //console.log(this);

        //  divide string to floats?
        //  here?
        //
        //

    }

    processArray() {

        this.data = this.data.slice(1, -1);
        this.processString();

    }

    processUnknownInput() {

        addErrorDiv("One of the input lines has a data of unknown type! Check console for more info");
        console.log("Unknown data: ");
        console.log(this);

    }

    logProperties() {

        console.log(this);

    }

}

function getTypeOfInput(data) {

    if (data[0] == "[") {
        return 1;
    } else
    if (!isNaN(data[0])) {
        return 0;
    } else if(
        data[0] == "-"){
            return 0;
        } else {
        return 99;
    }

}

function getMult(data) {

    if (data.includes("*", 1)) {
        return 0;
    }

    if (data.includes("e", 1)) {
        return 1;
    }

    return 99;

}

function getComma(data) {

    if (getSeparator(data) == 1) {
        if (data.includes(".")) {
            return 1;
        }
    } else if (getSeparator(data) == 0) {
        if (data.includes(",")) {
            return 0;
        }
    } else if (getSeparator(data) == 2) {
        if (data.includes(",")) {
            return 0;
        } else if (data.includes(".")) {
            return 1;
        }
    }


    return 99;
}

function getSeparator(data) {

    let mult = getMult(data);

    if (mult != 99) {

        if (mult == 0) {

            mult = "*";
            let pos = data.indexOf(mult);

            for (let i = pos; i > 0; i--) {
                if (isNaN(data[i])) {
                    if (data[i] == ".") {
                        return 0
                    }
                    if (data[i] == ",") {
                        return 1
                    }
                }
            }
        } else if (mult == 1) {

            mult = "e";
            let pos = data.indexOf(mult);

            for (let i = pos; i < data.length; i++) {
                if (isNaN(data[i])) {
                    if (data[i] == ".") {
                        return 0
                    }
                    if (data[i] == ",") {
                        return 1
                    }
                }
            }
        }

    } else if (mult == 99) {

        let spacePos = data.indexOf(" ", 2);

        if(spacePos != -1){

            if (data[spacePos - 1] == ".") {

                return 0;
    
            } else if (data[spacePos - 1] == ",") {
    
                return 1;
    
            } else {
    
                if (!isNaN(data[spacePos - 1])) {
    
                    return 2;
    
                }
    
                return 99;
            }
     
        } else{

            return 1;

        }

    } else {

        return 99;

    }

}

function beautifyResult(data) {

    for(let i = 0; i < data.length; i++){
        if(endsWith(data[i], " ")){
            data[i] = data[i].slice(0, -1);
        }

        if(data[i][0] == " "){
            data[i] = data[i].slice(1);
        }

    }

    return data;
}

function endsWith(str, suffix) {

    return str.indexOf(suffix, str.length - suffix.length) !== -1;

}

function insertAfter(newNode, existingNode) {

    existingNode.parentNode.insertBefore(newNode, existingNode.nextSibling);

}

//Download data function

function downloadData() {

    const blob = new Blob([csvmaker(graphData)], { type: 'text/csv' });
 
    const url = window.URL.createObjectURL(blob);
 
    const a = document.createElement('a');
 
    a.setAttribute('href', url);
 
    a.setAttribute('download', 'data.csv');

    a.click();
}
 
const csvmaker = function (data) {
 
    // Empty array for storing the values
    let csvRows = [];
    
    const headers = Object.keys(data);

    for(let key of headers){

        const line = data[key]
        csvRows.push(key + "," + line);

    }
 
    return csvRows.join('\n')
}

// Lower and upper parameters for frontend. If user enters out of bounds number, error msg should appear.

const bounds = {

    field: {min: 0.5, max: 20},
    radius: {min: 0.5, max: 995},

    workFunction: {min: 0, max: 101},
    temperature: {min: 9, max: 3001}

}

function main(){

    // Connection with node.js socket.io server

    let socket = io();


    //local global variables for main function

    let materialType, field, radius, workFunction, temperature,
        fieldMult, radiusMult, workFunctionMult, temperatureMult,
        calculateNH, calculateES, calculateEC, gammaMetal, gammaSemi,
        ec, ef, eg, me, mp, _field, _radius, _workFunction, _temperature,
        _ec, _ef, _eg, _me, _mp, _gammaMetal, _gammaSemi, sweepParam,
        data, chart1, chart2, chart3, _data
    

    // boolean used to load charts only when users inputs data

    let chartsLoaded = false;

    //Loads event listener when user loads the page

    loadInitEventListeners();

    //Big function that aims to check if data can be sent to server

    function checkValidity(){

        //Removes old errors so new can be displayed without stacking

        function removeErrorDivs(){

            if (errorDivs.length > 0) {

                errorDivs.forEach(div => {

                    div.remove();

                })

                errorDivs = [];
    
            }
        }

        //Get values that user typed into input fields. Loads all fields. Syncs them with main global vars.

        function getValuesFromUserInputFields(){

            materialType = document.getElementById("pickMaterialType").value;
            sweepParam = document.getElementById("pickChangingVar").value;
            field = document.getElementById("field_in").value;

            radius = document.getElementById("radius_in").value;
            workFunction = document.getElementById("wf_in").value;
            temperature = document.getElementById("temperature_in").value;
    
            fieldMult = document.getElementById("field_mult_in").value;
            radiusMult = document.getElementById("radius_mult_in").value;

            //gammaMult = document.getElementById("gamma_mult_in").value;
            workFunctionMult = document.getElementById("wf_mult_in").value;
            temperatureMult = document.getElementById("temperature_mult_in").value;
    
            gammaMetal = document.getElementById("gammaMetalParam").value;
            gammaSemi = document.getElementById("gammaSemiParam").value
    
            ec = document.getElementById("ecParam").value;
            ef = document.getElementById("efParam").value;
            eg = document.getElementById("egParam").value;

            me = document.getElementById("meParam").value;
            mp = document.getElementById("mpParam").value;
    
            calculateEC = document.getElementById("chooseEC").checked;
            calculateES = document.getElementById("chooseES").checked;
            calculateNH = document.getElementById("chooseNH").checked;
    
        }

        //Applies default values if fields are empty, proccesses inputs

        function processUserInputFields(){

            if(field == "") field = "10";
            if(radius == "") radius = "50";

            if(gammaMetal == "") gammaMetal = "10";
            if(gammaSemi == "") gammaSemi = "10";

            if(workFunction == "") workFunction = "4.5";
            if(temperature == "") temperature = "300";

            if(ec == "") ec = "4.05";
            if(ef == "") ef = "4.61";
            if(eg == "") eg = "1.12";

            if(me == "") me = "0.98";
            if(mp == "") mp = "0.5";

            _field = processDataInput(field);
            _radius = processDataInput(radius);

            _gammaMetal = processDataInput(gammaMetal);
            _gammaSemi = processDataInput(gammaSemi);

            _workFunction = processDataInput(workFunction);
            _temperature = processDataInput(temperature);

            _ec = processDataInput(ec);
            _ef = processDataInput(ef);
            _eg = processDataInput(eg);

            _me = processDataInput(me);
            _mp = processDataInput(mp);

        }

        //Binarize the booleans so its smaller packets to send

        function getWhatToCompute(){

            if(calculateEC){

                calculateEC = 1;

            } else {

                calculateEC = 0;

            }
    
            if(calculateES){

                calculateES = 1;

            } else {

                calculateES = 0;

            }
    
            if(calculateNH){

                calculateNH = 1;

            } else {
                
                calculateNH = 0;

            }
                
        }

        //BUG temperature should probably be not multiplied
        //This function updates all the values considering picked multipliers (i.e Units) in the fields prompts.

        function updateValuesUnits(){

            for(let i = 0; i < _field.length; i++){

                _field[i] = _field[i] * fieldMult;

            }
    
            for(let i = 0; i < _radius.length; i++){

                _radius[i] = _radius[i] * radiusMult;

            }
    
            for(let i = 0; i < _workFunction.length; i++){

                _workFunction[i] = _workFunction[i] * workFunctionMult;

            }

            for(let i = 0; i < _temperature.length; i++){

                _temperature[i] = _temperature[i] * temperatureMult;

            }
    
        }

        //This function converts sweepparam to an integer, for smaller packets
        function convertSweepParamToNumber(){

            if(typeof(sweepParam) == 'string'){

                switch(sweepParam){

                    case "1":

                        sweepParam = 2;

                        break;
        
                    case "2":

                        sweepParam = 3;

                        break;
        
                    case "4":

                        sweepParam = 4;

                        break;
                    
                    case "5":

                        sweepParam = 5;

                        break;
        
                }

            }

        }

        //This function converts material type to an integer, for smaller packets

        function convertMaterialTypeToNumber(){

            if(typeof(materialType) == 'string'){

                switch(materialType){

                    case "1":
                        
                        materialType = 1;

                        break;

                    case "2":

                        materialType = 2;

                        break;

                }

            }

        }

        //Main logic spot. Checks if data can be sent to server

        function checkForCanCompute(){

            if(calculateEC == false && calculateES == false && calculateNH == false){

                addErrorDiv("Please select what do You want to compute, ie. Emitted Current, Nottigham Heat or Electron spectrum");

                return false;

            }


            if(data[sweepParam].length >= 3){

                for(let i = 2; i < 5; i++){

                    if(i != sweepParam){

                        if(data[i].length > 1){

                            addErrorDiv("Only sweeping parameter accepts multiple values");

                            return false;

                        } else { 

                            let pointsCount = document.getElementById('numberOfPointsInBounds').value;

                            if(pointsCount != undefined){
                                
                                if(pointsCount > 151){

                                    addErrorDiv('Please do not enter more than 151 points for sweeping parameter');
                                    return false;
    
                                } else {
    
                                    return true; 
    
                                }    

                            } else return true;
                        
                        }

                    }

                }

            } else {


                addErrorDiv("Sweeping parameter must have at least 3 values");

                return false;

            }

            return false;

        }

        //Below are functions that are ran once main is called

        removeErrorDivs();
    
        getValuesFromUserInputFields();

        processUserInputFields();

        getWhatToCompute();

        updateValuesUnits();

        convertSweepParamToNumber();

        convertMaterialTypeToNumber();

        //data object that is to be sent to server, can be jsoned etc.

        data = [[materialType], [sweepParam], _field,
             _radius, _workFunction, _temperature,
            _ec, _ef, _eg, _gammaMetal, _gammaSemi, _me, _mp,
            [calculateEC], [calculateNH], [calculateES]];
        
    
        //if all conditions are met, send data to server. Displays loading modal

        if(checkForCanCompute()){

            updateSliderBounds(data[`${sweepParam}`].length);

            socket.emit('calculateEmission', data);

            console.log(data);

            $('#loadingModal').modal('show');

        }

    }

    function loadInitEventListeners(){

        let enterButton = document.getElementById("enterButton");
        let pickChangingVarGenerationMethodDiv = document.getElementById("autoGenerateValuesToggle");

        let generateValuesBtn = document.getElementById("generateValuesButton");
        let pickMaterialTypeDiv = document.getElementById("pickMaterialType");

        let advancedModeToggleDiv = document.getElementById("advancedModeToggle");
        let advancedModeToggleMainDiv = document.getElementById("advancedModeToggleMainDiv");

        let pickChangingVarDiv = document.getElementById("pickChangingVar");
        let preselectSemiPropertiesDiv = document.getElementById("preselectSemiProperties");

        let slider = document.getElementById('myRange');

        //advancedModeToggleMainDiv.hidden = true;

        enterButton.addEventListener("click", checkValidity);
        pickChangingVarGenerationMethodDiv.addEventListener("change", updateAutoGenerateValuesDiv);

        generateValuesBtn.addEventListener("click", autoGenerateValues);
        pickMaterialTypeDiv.addEventListener("change", updateWFName);

        advancedModeToggleDiv.addEventListener("change", updateAdvancedMethods);
        pickChangingVarDiv.addEventListener("change", updatePropertiesPresets);
        
        preselectSemiPropertiesDiv.addEventListener("change", updatePreselectSemiProperties);
        enterButton.addEventListener("click", checkValidity);

        updatePropertiesPresets();
        
    }

    function loadCharts(){

        const ctx1 = document.getElementById("currentChart");
        const ctx2 = document.getElementById("heatChart");
        const ctx3 = document.getElementById("spectrumChart");


        //Data used for debugging, can see it as first data before anything loads

        const dddata = [{

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


        //Emitted current chart object

        chart1 = new Chart(ctx1, {

            type: 'line',

            data: {

                datasets: [{

                    label: 'emitted current',

                    data: dddata,

                    backgroundColor: [

                        'rgba(255, 99, 132, 0.2)'

                    ],

                    borderColor: [

                        'rgba(255, 99, 132, 1)'

                    ],

                    borderWidth: 1,

                    interpolate: true

                }]

            },

            options: {

                fill: false,

                lineTension: 0.1,

                scales: {

                    y: {
    
                        title: {
    
                            display: true,

                            text: "Emitted current density",

                            font:{

                                size: 18

                            }
    
                        },
    
                        type: "logarithmic",

                        position: "bottom",

                        ticks: {
    
                            callback: function (value, index, ticks) {
    
                                if (value === 1000000) {
    
                                    return "1 [MA/nm^2]"
    
                                }
    
                                if (value === 100000) {
    
                                    return "100 [kA/nm^2]"
    
                                }
    
                                if (value === 10000) {
    
                                    return "10 [kA/nm^2]"
    
                                }
    
                                if (value === 1000) {
    
                                    return "1 [kA/nm^2]"
    
                                }
    
                                if (value === 100) {
    
                                    return "100 [A/nm^2]"
    
                                }
    
                                if (value === 10) {
    
                                    return "10 [A/nm^2]"
    
                                }
    
                                if (value === 1) {
    
                                    return "1 [A/nm^2]"
    
                                }
    
                                if (value === 0.1) {
    
                                    return "100 [mA/nm^2]"
    
                                }
    
                                if (value === 0.01) {
    
                                    return "10 [mA/nm^2]"
    
                                }
    
                                if (value === 0.001) {
    
                                    return "1 [mA/nm^2]"
    
                                }

                                if (value === 0.0001) {
    
                                    return "0.1 [mA/nm^2]"
    
                                }

                                if (value === 0.00001) {
    
                                    return "0.01 [mA/nm^2]"
    
                                }

                                if (value === 0.000001) {
    
                                    return "1000 [nA/nm^2]"
    
                                } 

                                if (value === 0.0000001) {
    
                                    return "100 [nA/nm^2]"
    
                                }

                                if (value === 0.00000001) {
    
                                    return "10 [nA/nm^2]"
    
                                }

                                if (value === 0.000000001) {
    
                                    return "1 [nA/nm^2]"
    
                                }

                                if (value === 0.0000000001) {
    
                                    return "100 [pA/nm^2]"
    
                                }

                                if (value === 0.00000000001) {
    
                                    return "10 [pA/nm^2]"
    
                                }

                                if (value === 0.000000000001) {
    
                                    return "1 [pA/nm^2]"
    
                                }

                                if (value === 0.0000000000001) {
    
                                    return "100 [fA/nm^2]"
    
                                }

                                if (value === 0.00000000000001) {
    
                                    return "10 [fA/nm^2]"
    
                                }

                                if (value === 0.000000000000001) {
    
                                    return "1 [fA/nm^2]"
    
                                }
    
                                return null
                            }
                        }
    
                    },
    
                    x: {
    
                        title: {
    
                            text: "Sweeping parameter",

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

                        callbacks: {
                            
                            label: function(context) {
                                let label = context.dataset.label || '';
        
                                if (label) {
                                    label += ': ';
                                }
                                if (context.parsed.y !== null) {
                                    label += context.parsed.y.toExponential(3);
                                }
                                return label;
                            }
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

                                min: -1,
                                max: 3000

                             },
                             y: {

                                 min: -10,
                                 max: 5e4

                             }
    
                            }
                    },

                    title: {

                        display: true,

                        text: "Chart 1",

                        font: {

                            size: 24

                        }
    
                    }

                }

            }

        });

        //Nottingham heat chart object

        chart2 = new Chart(ctx2, {

            type: 'line',
            
            data: {

                datasets: [{

                    label: 'positive Nottingham heat',

                    data: dddata,

                    backgroundColor: [

                        'rgba(255, 99, 132, 0.2)'

                    ],

                    borderColor: [

                        'rgba(255, 99, 132, 1)'

                    ],

                    borderWidth: 1,

                    interpolate: true

                },

                {

                    label: 'negative Nottingham heat',

                    data: dddata,

                    backgroundColor: [

                        'rgba(0, 32, 255, 0.2)'

                    ],

                    borderColor: [

                        'rgba(0, 32, 255, 1)'

                    ],

                    borderWidth: 1,

                    interpolate: true
                }
                
            ]

            },

            options: {

                spanGaps: true,

                fill: false,

                lineTension: 0.1,

                scales: {

                    y: {
    
                        title: {
    
                            display: true,

                            text: "Nottigham heat",

                            font:{

                                size: 18

                            }
    
                        },
    
                        type: "logarithmic",

                        position: "bottom",

                        ticks: {
    
                            callback: function (value, index, ticks) {
    
                                if (value === 1000000) {
    
                                    return "1 [MW/nm^2]"
    
                                }
    
                                if (value === 100000) {
    
                                    return "100 [kW/nm^2]"
    
                                }
    
                                if (value === 10000) {
    
                                    return "10 [kW/nm^2]"
    
                                }
    
                                if (value === 1000) {
    
                                    return "1 [kW/nm^2]"
    
                                }
    
                                if (value === 100) {
    
                                    return "100 [W/nm^2]"
    
                                }
    
                                if (value === 10) {
    
                                    return "10 [W/nm^2]"
    
                                }
    
                                if (value === 1) {
    
                                    return "1 [W/nm^2]"
    
                                }
    
                                if (value === 0.1) {
    
                                    return "100 [mW/nm^2]"
    
                                }
    
                                if (value === 0.01) {
    
                                    return "10 [mW/nm^2]"
    
                                }
    
                                if (value === 0.001) {
    
                                    return "1 [mW/nm^2]"
    
                                }

                                if (value === 0.0001) {
    
                                    return "0.1 [mW/nm^2]"
    
                                }

                                if (value === 0.00001) {
    
                                    return "0.01 [mW/nm^2]"
    
                                }

                                if (value === 0.000001) {
    
                                    return "1000 [nW/nm^2]"
    
                                } 

                                if (value === 0.0000001) {
    
                                    return "100 [nW/nm^2]"
    
                                }

                                if (value === 0.00000001) {
    
                                    return "10 [nW/nm^2]"
    
                                }

                                if (value === 0.000000001) {
    
                                    return "1 [nW/nm^2]"
    
                                }

                                if (value === 0.0000000001) {
    
                                    return "100 [pW/nm^2]"
    
                                }

                                if (value === 0.00000000001) {
    
                                    return "10 [pW/nm^2]"
    
                                }

                                if (value === 0.000000000001) {
    
                                    return "1 [pW/nm^2]"
    
                                }

                                if (value === 0.0000000000001) {
    
                                    return "100 [fW/nm^2]"
    
                                }

                                if (value === 0.00000000000001) {
    
                                    return "10 [fW/nm^2]"
    
                                }

                                if (value === 0.000000000000001) {
    
                                    return "1 [fW/nm^2]"
    
                                }

                                return null
                            }
                        }
    
                    },
    
                    x: {
    
                        title: {
    
                            text: "Sweeping parameter",

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

                        callbacks: {
                            
                            label: function(context) {
                                let label = context.dataset.label || '';
        
                                if (label) {
                                    label += ': ';
                                }
                                if (context.parsed.y !== null) {
                                    label += context.parsed.y.toExponential(3);
                                }
                                return label;
                            }
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

                                min: -1,
                                max: 3000
                                
                             },

                             y: {

                                min: -10,
                                max: 5e4

                             }
    
                        }
                    },

                    title: {

                        display: true,

                        text: "Chart 2",

                        font: {

                            size: 24

                        }
    
                    }

                }

            }

        });

        //Electron spectrum chart object

        chart3 = new Chart(ctx3, {

            type: 'line',

            data: {

                datasets: [{

                    label: 'electron spectrum',

                    data: dddata,

                    backgroundColor: [

                        'rgba(255, 99, 132, 0.2)'

                    ],

                    borderColor: [

                        'rgba(255, 99, 132, 1)'

                    ],

                    borderWidth: 1,

                    interpolate: true

                }

                ]

            },

            options: {

                fill: false,

                lineTension: 0.1,

                scales: {

                    y: {
    
                        title: {
    
                            display: true,

                            text: "Electron count",

                            font:{

                                size: 18

                            }
    
                        },
    
                        type: "logarithmic",

                        position: "bottom",

                        ticks: {
    
                            callback: function (value, index, ticks) {
    
                                if (value === 1000000) {
    
                                    return "1 [M]"
    
                                }
    
                                if (value === 100000) {
    
                                    return "100 [k]"
    
                                }
    
                                if (value === 10000) {
    
                                    return "10 [k]"
    
                                }
    
                                if (value === 1000) {
    
                                    return "1 [k]"
    
                                }
    
                                if (value === 100) {
    
                                    return "100 []"
    
                                }
    
                                if (value === 10) {
    
                                    return "10 []"
    
                                }
    
                                if (value === 1) {
    
                                    return "1 []"
    
                                }
    
                                if (value === 0.1) {
    
                                    return "100 [m]"
    
                                }
    
                                if (value === 0.01) {
    
                                    return "10 [m]"
    
                                }
    
                                if (value === 0.001) {
    
                                    return "1 [m]"
    
                                }

                                if (value === 0.0001) {
    
                                    return "0.1 [m]"
    
                                }

                                if (value === 0.00001) {
    
                                    return "0.01 [m]"
    
                                }

                                if (value === 0.000001) {
    
                                    return "1000 [n]"
    
                                } 

                                if (value === 0.0000001) {
    
                                    return "100 [n]"
    
                                }

                                if (value === 0.00000001) {
    
                                    return "10 [n]"
    
                                }

                                if (value === 0.000000001) {
    
                                    return "1 [n]"
    
                                }

                                if (value === 0.0000000001) {
    
                                    return "100 [p]"
    
                                }

                                if (value === 0.00000000001) {
    
                                    return "10 [p]"
    
                                }

                                if (value === 0.000000000001) {
    
                                    return "1 [p]"
    
                                }

                                if (value === 0.0000000000001) {
    
                                    return "100 [f]"
    
                                }

                                if (value === 0.00000000000001) {
    
                                    return "10 [f]"
    
                                }

                                if (value === 0.000000000000001) {
    
                                    return "1 [f]"
    
                                }
    
                                return null
                            }
                        }
    
                    },
    
                    x: {
    
                        title: {
    
                            text: "Energy space, [eV]",

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

                animation: {

                    duration: 0

                },

                plugins: {

                    tooltip: {

                        callbacks: {
                            
                            label: function(context) {
                                let label = context.dataset.label || '';
        
                                if (label) {
                                    label += ': ';
                                }
                                if (context.parsed.y !== null) {
                                    label += context.parsed.y.toExponential(3);
                                }
                                return label;
                            }
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

                                min: -1,
                                max: 3000

                             },
                             y: {

                                 min: -10,
                                 max: 5e4

                             }
    
                            }
                    },

                    title: {

                        display: true,

                        text: "Chart 3",

                        font: {

                            size: 24

                        }
    
                    }

                }

            }

        });

    }

    //Listens for data from server. Once received, updates graphs with given data. Hides loading modal

    socket.on("calculatedEmission", (data) =>{

        if(!chartsLoaded){

            loadCharts();

            _data = data;

            let slider = document.getElementById('myRange');

            slider.addEventListener('input', updateESGraph);
            //slider.addEventListener('input', readSliderValue);

            chartsLoaded = true;

        }

        console.log(data);
        console.log(`Size of data: ` + memorySizeOf(data));

        updateGraphs();

        $('#loadingModal').modal('hide');

    })

    socket.on('logServerSideError', (arg) => {

        $('#loadingModal').modal('hide');

        alert("Got an error from server while computing data. Check console for more info.");
        console.log("Error server-side: ");
        console.log(arg);

    })
    //Big function that updates graphs with given data.

    function updateGraphs(update = 'all'){
        
        data = _data;
    
        let _sweepParam = data.sweepParam;
        let _materialType = data.materialType;
        let _sweepValues;
    
        let _field = data.field;
        let _radius = data.radius;
    
        let _wf = data.work_function;
        let _temp = data.temperature;
    
        let data1 = data.metalEC;
        let data2 = data.metalNH;
    
        let data4 = data.semiEC;
    
        let data5 = data.semiNH;
    
        let data6 = data.metalESelcount;
        let data6e = data.metalESenergy;
    
        //For ease, copies the array values of parameters to a new array, so that its easier to access right data.
    
        updateSweepValues();
    
        if(update == 'all'){
    
            updateECGraph();
            updateNHGraph();
            updateESGraph();
    
        } else if(update == 'es'){
    
            updateESGraph();
    
        }
    
        let value = document.getElementById('myRangeValue');
        let slider = document.getElementById('myRange');
        value.value = `${data.sweepParam}: ${_sweepValues[slider.value]}`;
    
        function updateECGraph(){
    
            if((data1 === undefined || data1.length == 0) && (data4 === undefined || data4.length == 0)){
    
                return;
    
            } else {
    
                let points = [];
    
                if(_materialType == "1"){
    
                    for(let i = 0; i < _sweepValues.length; i++){
    
                        points.push({x: _sweepValues[i], y: data1[i]});
    
                    }
    
                } else {
                    
                    for(let i = 0; i < _sweepValues.length; i++){
    
                        points.push({x: _sweepValues[i], y: data4[i]});
    
                    }
    
                }
    
                chart1.resetZoom();
    
                updateSweepParameterBounds(chart1, data.sweepParam);
    
                updateTitle(chart1, data);
    
                updatePoints(chart1, points);    
                
            }
    
        }
    
        function updateESGraph(){
    
    
            if((data6e === undefined || data6e.length == 0)){
    
                return;
    
            } else {
    
                let sliderValue = document.getElementById('myRange').value;
    
                let points = [];
    
                if(_materialType == "1"){
    
                    for(let i = 0; i < data6e[0].length; i++){
    
                        points.push({x: data6e[sliderValue][i], y: data6[sliderValue][i]});
    
                    }
    
                } else {
    
                    alert("Semiconductor electron spectrum calculations are not yet implemented");
    
                }
    
                chart3.options.plugins.zoom.limits.x = {min: -10, max: 5};
    
                chart3.resetZoom();
    
                chart3.update();
    
                updateTitle(chart3, data);
    
                updatePoints(chart3, points);    
    
            }
    
        }
    
        function updateNHGraph(){
    
            if((data2 === undefined || data2.length == 1) && (data5 === undefined || data5.length == 0)){
    
                return;
    
            } else {
    
                let points = [];
    
                if(_materialType == "1"){
    
                    for(let i = 0; i < _sweepValues.length; i++){
    
                        points.push({x: _sweepValues[i], y: data2[i]})
    
                    }
    
                } else {
                    
                    for(let i = 0; i < _sweepValues.length; i++){
    
                        points.push({x: _sweepValues[i], y: data5[i]})
    
                    }
    
                }
    
                chart2.resetZoom();
    
                updateSweepParameterBounds(chart2, data.sweepParam);
    
                updateTitle(chart2, data);
    
                updatePoints(chart2, points);    
    
            }
    
        }
    
        function updateTitle(chart, data){
    
            if(data.materialType == "1"){
    
                switch(data.sweepParam){
    
                    case "field":
                        
                        chart.options.plugins.title.text = `Radius: ${data.radius[0]} nm, Work Function: ${data.work_function[0]} eV, Temperature: ${data.temperature[0]} K` 
    
                        break;
    
                    case "radius":
    
                        chart.options.plugins.title.text = `Field: ${data.field[0]} V/nm, Work Function: ${data.work_function[0]} eV, Temperature: ${data.temperature[0]} K` 
    
                        break;
    
                    case "wf":
    
                        chart.options.plugins.title.text = `Field: ${data.field[0]} V/nm, Radius: ${data.radius[0]} nm, Temperature: ${data.temperature[0]} K` 
    
                        break;
    
                    case "temp":
    
                        chart.options.plugins.title.text = `Field: ${data.field[0]} V/nm, Radius: ${data.radius[0]} nm, Work Function: ${data.work_function[0]} eV` 
    
                        break;
    
    
                }
    
            } else  {
                
                switch(data.sweepParam){
    
                    case "field":
    
                        chart.options.plugins.title.text = [`Radius: ${data.radius[0]} nm, Work Function: ${data.work_function[0]} eV, Temperature: ${data.temperature[0]} K, Ec: ${data.ec} eV`, `Ef: ${data.ef} eV, Eg: ${data.eg} eV, Me: ${data.me} m/me, Mp: ${data.mp} m/me`]; 
                        
                        break;
    
                    case "radius":
    
                        chart.options.plugins.title.text = [`Field: ${data.field[0]} V/nm, Work Function: ${data.work_function[0]} eV, Temperature: ${data.temperature[0]} K, Ec: ${data.ec} eV`, `Ef: ${data.ef} eV, Eg: ${data.eg} eV, Me: ${data.me} m/me, Mp: ${data.mp} m/me`] 
                        
                        break;
    
                    case "wf":
    
                        chart.options.plugins.title.text = [`Field: ${data.field[0]} V/nm, Radius: ${data.radius[0]} nm, Temperature: ${data.temperature[0]} K, Ec: ${data.ec} eV`, `Ef: ${data.ef} eV, Eg: ${data.eg} eV, Me: ${data.me} m/me, Mp: ${data.mp} m/me`]
                        
                        break;
    
                    case "temp":
    
                        chart.options.plugins.title.text = [`Field: ${data.field[0]} V/nm, Radius: ${data.radius[0]} nm, Work Function: ${data.work_function[0]} eV, Ec: ${data.ec} eV`, `Ef: ${data.ef} eV, Eg: ${data.eg} eV, Me: ${data.me} m/me, Mp: ${data.mp} m/me`]
                        
                        break;
                            
                }
    
            }
    
            chart.update();
    
        }
    
        function updateSweepParameterBounds(chart, sweepParam){
    
            switch(sweepParam){
    
                case "field":
                    
                    chart.options.scales.x.title.text = "Field, [V/nm]"
    
                    chart.options.plugins.zoom.limits.x = {min: bounds.field.min, max: bounds.field.max};
    
                    chart.update();
    
                    break;
    
                case "radius":
    
                    chart.options.scales.x.title.text = "Radius, [nm]"
    
                    chart.options.plugins.zoom.limits.x = {min: bounds.radius.min, max: bounds.radius.max};
    
                    chart.update();
    
                    break;
    
                case "wf":
    
                    chart.options.scales.x.title.text = "Work Function, [eV]"
    
                    chart.options.plugins.zoom.limits.x = {min: bounds.workFunction.min, max: bounds.workFunction.max};
    
                    chart.update();
    
                    break;
                    
                case "temp":
    
                    chart.options.scales.x.title.text = "Temperature, [K]"
    
                    chart.options.plugins.zoom.limits.x = {min: bounds.temperature.min, max: bounds.temperature.max};
    
                    chart.update();
    
                    break;
    
            }
        }
    
        function updatePoints(chart, points){
    
            chart.data.datasets.forEach((dataset) =>{
    
                dataset.data = null;
                
            });
    
            chart.update();
    
            if(chart.data.datasets.length > 1){
    
                let posDataSet = chart.data.datasets[0];
                let negDataSet = chart.data.datasets[1];
    
                points.forEach((point) => {
    
                    if(point.y > 0){
    
                        posDataSet.data.push(point);
                        negDataSet.data.push({x: point.x, y: null});
    
                    } else{
    
                        negDataSet.data.push({x: point.x, y: -point.y});
                        posDataSet.data.push({x: point.x, y: null});
    
                    }
    
                });
    
            } else{
    
                let dataSet = chart.data.datasets[0];
    
                points.forEach((point) => {
    
                    dataSet.data.push(point);
    
                })
    
            }
    
            chart.update();
        }
         
        function updateSweepValues(){
    
            switch(_sweepParam){
    
                case "field":
    
                    _sweepValues = _field;
    
                    break;
    
                case "radius":
    
                    _sweepValues = _radius;
    
                    break;
    
                case "wf":
    
                    _sweepValues = _wf;
    
                    break;
    
                case "temp":
    
                    _sweepValues = _temp;
    
                    break;
    
            }
    
        }
    
    }
    
    
    function updateESGraph(){
    
        updateGraphs('es');
    
    }
    

}

//Global var, used to store all error message divs for quick access

let errorDivs = [];

//Main function call

main();

//Used to log slider value

function readSliderValue(){

    let slider = document.getElementById('myRange');
    console.log(slider.value);

}

//Update amount of points for slider
function updateSliderBounds(points){

    let slider = document.getElementById('myRange');
    slider.max = `${points-1}`;

}

//Function that is used to predefine properties for predetermined materials for the user.
//Its goal is to give a set of materials and apply their values to corresponding fields.

function updatePreselectSemiProperties(){

    let preselectSemiPropertiesDiv = document.getElementById("preselectSemiProperties");

    let ec = document.getElementById("ecParam");
    let ef = document.getElementById("efParam");

    let eg = document.getElementById("egParam");
    let gamma = document.getElementById("gammaSemiParam");

    let me = document.getElementById("meParam");
    let mp = document.getElementById("mpParam");

    const matPropDict = {

        //[name, energygap300K, me, mp, ef, ec, eg]
        //Taken from Physics of Semiconductor devices 2nd edition book, S. M. Sze

        1: ["Si", 1.12, 0.98, 0.49, 4.61, 4.05, 5.17],
        2: ["Ge", 0.66, 1.64, 0.28, 4.33, 4, 4.66],
        3: ["GaAs", 1.42, 0.067, 0.082, 4.78, 4.07, 5.49]

    }

    switch(preselectSemiPropertiesDiv.value){

        case "1":
            
            setValues(1);

            break;

        case "2":

            setValues(2);

            break;

        case "3":

            setValues(3);

            break;
    }

    function setValues(num){

        eg.value = matPropDict[num][1];
        me.value = matPropDict[num][2];

        mp.value = matPropDict[num][3];
        ef.value = matPropDict[num][4];

        ec.value = matPropDict[num][5];
        gamma.value = 10;
        

    }

}

//This function is used to predefine default values for some input fields if user does not enter anything
//Also is an example for user of what data they can enter

function updatePropertiesPresets(){

    const fieldLong = "4.5, 7, 9.5, 13, 15.5";
    const fieldShort = "10";

    const radiusLong = "50, 60, 70, 80, 90";
    const radiusShort = "50";

    const wfLong = "4, 5, 6, 7, 8";
    const wfShort = "4.5";

    const tempLong = "273.15, 298.15, 300, 350, 500";
    const tempShort = "300";

    let fieldDiv = document.getElementById("field_in");
    let radiusDiv = document.getElementById("radius_in");

    let wfDiv = document.getElementById("wf_in");
    let tempDiv = document.getElementById("temperature_in");

    let pickChangingVarDiv = document.getElementById("pickChangingVar");

    switch(pickChangingVarDiv.value){

        case "1":

            fieldDiv.placeholder = fieldLong;
            radiusDiv.placeholder = radiusShort;

            wfDiv.placeholder = wfShort;
            tempDiv.placeholder = tempShort;

            break;

        case "2":

            fieldDiv.placeholder = fieldShort;
            radiusDiv.placeholder = radiusLong;

            wfDiv.placeholder = wfShort;
            tempDiv.placeholder = tempShort;

            break;

        case "4":

            fieldDiv.placeholder = fieldShort;
            radiusDiv.placeholder = radiusShort;

            wfDiv.placeholder = wfLong;
            tempDiv.placeholder = tempShort;

            break;

        case "5":

            fieldDiv.placeholder = fieldShort;
            radiusDiv.placeholder = radiusShort;

            wfDiv.placeholder = wfShort;
            tempDiv.placeholder = tempLong;

            break;

    }


}

//Hide or display the material properties window

function updateAdvancedMethods(){

    let advancedModeToggleDiv = document.getElementById("advancedModeToggle");
    let materialSelectDiv = document.getElementById("pickMaterialType");

    let advancedParametersSemi = document.getElementById("advancedParametersSemi");
    let advancedParametersMetal = document.getElementById("advancedParametersMetal");


    if(advancedModeToggleDiv.checked == false){

        advancedParametersMetal.hidden = true;
        advancedParametersSemi.hidden = true;

    } else if(advancedModeToggleDiv.checked == true){

        if(materialSelectDiv.value == 1){

            advancedParametersMetal.hidden = false;
            advancedParametersSemi.hidden = true;

        } else if(materialSelectDiv.value == 2){

            advancedParametersMetal.hidden = true;
            advancedParametersSemi.hidden = false;

        }

    }

}


//Function used to generate n values from a to b. Checks for good bounds inside. Returns array of numbers!

function autoGenerateValues(){

    //TEMPORARY. BETTER INPUT TRANSFORMING LATER

    let lowerBound, upperBound, count

    let varNum = document.getElementById("pickChangingVar").value;

    try{

        lowerBound = parseFloat(document.getElementById("lowerBound").value);

        upperBound = parseFloat(document.getElementById("upperBound").value);

        count = parseInt(document.getElementById("numberOfPointsInBounds").value);

    } catch (e){

        
        console.log(e);

    }

    let canCompute = true;

    canCompute = checkForGoodBounds();

    function checkForGoodBounds(){

        if (errorDivs.length > 0) {

            errorDivs.forEach(div => {

                div.remove();

            })

            errorDivs = [];

        }

        let result = true;

        switch(varNum){

            case "1":

                if(lowerBound <= bounds.field.min || upperBound >= bounds.field.max){

                    addErrorDiv(`Field values must be in range of ${bounds.field.min} < x < ${bounds.field.max} V/nm`);

                    result = false;

                }

                break;

            case "2":

                if(lowerBound <= bounds.radius.min || upperBound >= bounds.radius.max){

                    addErrorDiv(`Radius values must be in range of ${bounds.radius.min} < x < ${bounds.radius.max} nm`);

                    result = false;

                }

                break;

            case "4":

                if(lowerBound <= bounds.workFunction.min || upperBound >= bounds.workFunction.max){

                    addErrorDiv(`Work Function values must be in range of ${bounds.workFunction.min} < x < ${bounds.workFunction.max} eV`);

                    result = false;

                }

                break;

            case "5":

                if(lowerBound <= bounds.temperature.min || upperBound >= bounds.temperature.max){

                    addErrorDiv(`Temperature values must be in range of ${bounds.temperature.min} < x < ${bounds.temperature.max} K`);

                    result = false;

                }

                break;  

        }

        return result;

    }

    if(canCompute){

        const varTypeDict = {

            1: "Field",
            2: "Radius",
            3: "Gamma",
            4: "Work_Function",
            5: "Temperature"

        }
    
        const dTot = upperBound - lowerBound;
        const dSt = dTot / (count - 1);
    
        let result = [];
    
        for(let i = 0; i < count; i++){
    
            if(i == 0){
    
                result.push((lowerBound + i * dSt).toFixed(2));
    
            } else {
    
                result.push(" " + (lowerBound + i * dSt).toFixed(2));
    
    
            }
        }
    
        switch(varNum){

            case "1":

                document.getElementById("field_in").value = result;

                break;

            case "2":

                document.getElementById("radius_in").value = result;

                break;

            case "3":

                document.getElementById("gamma_in").value = result;

                break;

            case "4":

                document.getElementById("wf_in").value = result;

                break;

            case "5":

                document.getElementById("temperature_in").value = result;

                break;

        }
    
    }

    //updateSliderBounds(count);

}

//In case of semiconductors Work Function makes no sence, so remove it when user picks semi, and display when metals are picked

function updateWFName(){

    let pickMaterialTypeDiv = document.getElementById("pickMaterialType");
    let pickChangingVarDiv = document.getElementById("pickChangingVar");

    let wfMain = document.getElementById("wf_main");
    let wfSelectId = document.getElementById("wfSelectId");

    let advancedModeToggleMainDiv = document.getElementById("advancedModeToggleMainDiv");
    let advancedModeToggle = document.getElementById("advancedModeToggle");

    let advancedModeSemi = document.getElementById("advancedParametersSemi");
    let advancedModeMetal = document.getElementById("advancedParametersMetal");



    if(pickMaterialTypeDiv.value == 1){

        //wfName.textContent = "Work Function";
        //pickChangingVarDiv.options[3].textContent = "Work Function";
        //advancedModeToggleMainDiv.hidden = true;

        wfMain.hidden = false;
        wfSelectId.hidden = false;

        if(advancedModeToggle.checked == true){

            advancedModeMetal.hidden = false;
            advancedModeSemi.hidden = true;

        } else {

            advancedModeMetal.hidden = true;
            advancedModeSemi.hidden = true;

        }

    } else if(pickMaterialTypeDiv.value == 2){

        //wfName.textContent = "Work Function";
        //pickChangingVarDiv.options[3].textContent = "Work Function";
        wfMain.hidden = true;
        wfSelectId.hidden = true;
        
        if(advancedModeToggle.checked == true){

            advancedModeMetal.hidden = true;
            advancedModeSemi.hidden = false;

        }   else {

            advancedModeMetal.hidden = true;
            advancedModeSemi.hidden = true;

        }

    }

}

//Hides or displays the window where user can define parameters to generate an array of equally distanced numbers (auto-gen values for sweep param)

function updateAutoGenerateValuesDiv(){

    let autoGenerateValuesDiv = document.getElementById("autoGeneratedInput");
    let pickChangingVarGenerationMethodDiv = document.getElementById("autoGenerateValuesToggle");

    if(pickChangingVarGenerationMethodDiv.checked == true){

        autoGenerateValuesDiv.hidden = false;

    } else if(pickChangingVarGenerationMethodDiv.checked == false){

        autoGenerateValuesDiv.hidden = true;

    }

}

//A general function used to display error divs on the page. Takes string as argument. Pushes div to errorDivs!

function addErrorDiv(message) {

    console.log(message);

    const template = `
    <section class="container showcase alert alert-danger alert-dismissible my-auto animated bounceInLeft">
        <strong class="mx-2">Error!</strong> ${message}!
        <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
    </section>
    `;

    let relativeDiv = document.getElementById("TopBarView2");

    let errorDiv = document.createElement("div");

    errorDiv.innerHTML = template;

    let section = errorDiv.children[0];

    errorDivs.push(section)

    insertAfter(section, relativeDiv);

}

//Used to store types as integers for smaller packets

const typeDict = {

    0: "String",
    1: "Array",
    99: "Unknown"

}

//Used to store separator as integers for smaller packets. A separator is something that separates two numbers, not a comma/dot in a single number!

const separatorDict = {

    0: ".",
    1: ",",
    2: " ",
    99: "Unknown"
    
}

//Used to store multiplies as integers for smaller packets. 

const multDict = {

    0: "*",
    1: "e",
    99: "Unknown"

}

//Used to store commas as integers for smaller packets. A comma is something that can be found inside a single number, and separates the whole part and fraction.

const commaDict = {

    0: ",",
    1: ".",
    99: "Unknown"

}

//Purpose of code below is to proccess user input and understand whats inside it

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

//A class for input, used for better OOP

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
            // raiseInputError("2006");
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

    //Removes "[" and "]"

    processArray() {

        this.data = this.data.slice(1, -1);

        this.processString();

    }

    //If cant determine type, log an error

    processUnknownInput() {

        addErrorDiv("One of the input lines has a data of unknown type! Check console for more info");

        console.log("Unknown data: ");

        console.log(this);

    }

    //Log instance of class

    logProperties() {

        console.log(this);

    }

}

//returns type of data considering dataTypeDict

function getTypeOfInput(data) {

    if (data[0] == "[") {

        return 1;

    } else if (!isNaN(data[0])) {
        
        return 0;

    } else if(data[0] == "-"){

            return 0;

        } else {

        return 99;

    }

}

//returns type of multiplier considering dataMultDict

function getMult(data) {

    if (data.includes("*", 1)) {

        return 0;

    }

    if (data.includes("e", 1)) {

        return 1;

    }

    return 99;

}

//returnc type of comma in a number considering dataCommaDict

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

//returns type of separator between different numbers considering dataSeparatorDict

function getSeparator(data) {

    let mult = getMult(data);

    if (mult != 99) {

        if (mult == 0) {

            mult = "*";

            let pos = data.indexOf(mult);

            for (let i = pos; i > 0; i--) {

                if (isNaN(data[i])) {

                    if (data[i] == ".") {

                        return 0;

                    }

                    if (data[i] == ",") {

                        return 1;

                    }

                }

            }

        } else if (mult == 1) {

            mult = "e";

            let pos = data.indexOf(mult);

            for (let i = pos; i < data.length; i++) {

                if (isNaN(data[i])) {

                    if (data[i] == ".") {

                        return 0;

                    }

                    if (data[i] == ",") {

                        return 1;
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

//Basically this functions removes spaces

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

//Checks if string ends with a given string

function endsWith(str, suffix) {

    return str.indexOf(suffix, str.length - suffix.length) !== -1;

}

//Inserts one HTML div right after the other

function insertAfter(newNode, existingNode) {

    existingNode.parentNode.insertBefore(newNode, existingNode.nextSibling);

}

function memorySizeOf(obj) {

    let bytes = 0;

    function sizeOf(obj) {

        if(obj !== null && obj !== undefined) {

            switch(typeof obj) {

                case 'number':

                    bytes += 8;
                    break;

                case 'string':

                    bytes += obj.length * 2;
                    break;

                case 'boolean':

                    bytes += 4;
                    break;

                case 'object':

                    var objClass = Object.prototype.toString.call(obj).slice(8, -1);

                    if(objClass === 'Object' || objClass === 'Array') {

                        for(var key in obj) {

                            if(!obj.hasOwnProperty(key)) continue;

                            sizeOf(obj[key]);
                        }

                    } else {
                    
                      bytes += obj.toString().length * 2;
                
                   }
                
                   break;
            }

        }

        return bytes;
    };

    function formatByteSize(bytes) {

        if(bytes < 1024) return bytes + " bytes";

        else if(bytes < 1048576) return(bytes / 1024).toFixed(3) + " KiB";

        else if(bytes < 1073741824) return(bytes / 1048576).toFixed(3) + " MiB";

        else return(bytes / 1073741824).toFixed(3) + " GiB";

    };

    return formatByteSize(sizeOf(obj));

};

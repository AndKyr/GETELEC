function main(){

    let socket = io();

    let materialType, field, radius, workFunction, temperature,
        fieldMult, radiusMult, workFunctionMult, temperatureMult,
        calculateNH, calculateES, calculateEC, gammaMetal, gammaSemi,
        ec, ef, eg, me, mp, _field, _radius, _workFunction, _temperature,
        _ec, _ef, _eg, _me, _mp, _gammaMetal, _gammaSemi, sweepParam,
        data

    loadInitEventListeners();

    function checkValidity(){

        function removeErrorDivs(){

            if (errorDivs.length > 0) {

                errorDivs.forEach(div => {
                    div.remove();
                })
    
                errorDivs = [];
    
            }
        }

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

        function processUserInputFields(){

            if(field == "") field = "250";
            if(radius == "") radius = "50";
            if(gammaMetal == "") gammaMetal = "10";
            if(gammaSemi == "") gammaSemi = "10";
            if(workFunction == "") workFunction = "4.5";
            if(temperature == "") temperature = "300";
            if(ec == "") ec = "1.12";
            if(ef == "") ef = "-0.6";
            if(eg == "") eg = "1.14";
            if(me == "") me = "0.33";
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

        function checkForCanCompute(){

            if(data[sweepParam].length >= 3){

                for(let i = 2; i < 5; i++){

                    if(i != sweepParam){

                        if(data[i].length > 1){

                            raiseInputError("3003");
                            return false;

                        } else{ return true; }
                    }

                }

            } else {


                raiseInputError("3002");
                return false;

            }

            return false;

        }

        removeErrorDivs();
    
        getValuesFromUserInputFields();

        processUserInputFields();

        getWhatToCompute();

        updateValuesUnits();

        convertSweepParamToNumber();

        convertMaterialTypeToNumber();

        data = [[materialType], [sweepParam], _field, _radius, _workFunction, _temperature,
            _ec, _ef, _eg, _gammaMetal, _gammaSemi, _me, _mp,
            [calculateEC], [calculateNH], [calculateES]];
        
        if(checkForCanCompute()){

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


        const chart1 = new Chart(ctx1, {
            type: 'line',
            data: {

                labels: ['Red', 'Blue', 'Yellow', 'Green', 'Purple', 'Orange'],
                datasets: [{
                    label: '# of Votes',
                    data: [12, 19, 3, 5, 2, 3],
                    backgroundColor: [
                        'rgba(255, 99, 132, 0.2)',
                        'rgba(54, 162, 235, 0.2)',
                        'rgba(255, 206, 86, 0.2)',
                        'rgba(75, 192, 192, 0.2)',
                        'rgba(153, 102, 255, 0.2)',
                        'rgba(255, 159, 64, 0.2)'
                    ],
                    borderColor: [
                        'rgba(255, 99, 132, 1)',
                        'rgba(54, 162, 235, 1)',
                        'rgba(255, 206, 86, 1)',
                        'rgba(75, 192, 192, 1)',
                        'rgba(153, 102, 255, 1)',
                        'rgba(255, 159, 64, 1)'
                    ],
                    borderWidth: 1
                }]

            },
            options: {

            }
        });

        const chart2 = new Chart(ctx2, {
            type: 'line',
            data: {

                labels: ['Red', 'Blue', 'Yellow', 'Green', 'Purple', 'Orange'],
                datasets: [{
                    label: '# of Votes',
                    data: [12, 19, 3, 5, 2, 3],
                    backgroundColor: [
                        'rgba(255, 99, 132, 0.2)',
                        'rgba(54, 162, 235, 0.2)',
                        'rgba(255, 206, 86, 0.2)',
                        'rgba(75, 192, 192, 0.2)',
                        'rgba(153, 102, 255, 0.2)',
                        'rgba(255, 159, 64, 0.2)'
                    ],
                    borderColor: [
                        'rgba(255, 99, 132, 1)',
                        'rgba(54, 162, 235, 1)',
                        'rgba(255, 206, 86, 1)',
                        'rgba(75, 192, 192, 1)',
                        'rgba(153, 102, 255, 1)',
                        'rgba(255, 159, 64, 1)'
                    ],
                    borderWidth: 1
                }]

            },
            options:{

            }

        });

        const chart3 = new Chart(ctx3, {
            type: 'line',
            data: {

                labels: ['Red', 'Blue', 'Yellow', 'Green', 'Purple', 'Orange'],
                datasets: [{
                    label: '# of Votes',
                    data: [12, 19, 3, 5, 2, 3],
                    backgroundColor: [
                        'rgba(255, 99, 132, 0.2)',
                        'rgba(54, 162, 235, 0.2)',
                        'rgba(255, 206, 86, 0.2)',
                        'rgba(75, 192, 192, 0.2)',
                        'rgba(153, 102, 255, 0.2)',
                        'rgba(255, 159, 64, 0.2)'
                    ],
                    borderColor: [
                        'rgba(255, 99, 132, 1)',
                        'rgba(54, 162, 235, 1)',
                        'rgba(255, 206, 86, 1)',
                        'rgba(75, 192, 192, 1)',
                        'rgba(153, 102, 255, 1)',
                        'rgba(255, 159, 64, 1)'
                    ],
                    borderWidth: 1
                }]

            },
            options:{

            }
        })

    }

    loadCharts();

    socket.on("calculatedEmission", (data) =>{

        updateGraphs(data);
        console.log(data);
        $('#loadingModal').modal('hide');

    })

    function updateGraphs(data){

        let _sweepParam = data.sweepParam;
        let _materialType = data.materialType;
        let _sweepValues;

        let _field = data.field;
        let _radius = data.radius;
        let _wf = data.wf;
        let _temp = data.temp;

        let data1 = data.metalEC;
        let data2 = data.metalNH;
        let data3 = data.metalES;
        let data4 = data.semiECl;
        let data5 = data.semiNH;
        let data6 = data.semiES;

        updateSweepValues();

        updateECGraph();

        updateNHGraph();

        updateESGraph();

        function updateECGraph(){

            if(_materialType == "1"){


            }

        }

        function updateESGraph(){

            return;

        }

        function updateNHGraph(){

            return;

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
}


let errorDivs = [];
main();

function updatePreselectSemiProperties(){

    let preselectSemiPropertiesDiv = document.getElementById("preselectSemiProperties");

    let ec = document.getElementById("ecParam");
    let ef = document.getElementById("efParam");
    let eg = document.getElementById("egParam");
    let gamma = document.getElementById("gammaSemiParam");
    let me = document.getElementById("meParam");
    let mp = document.getElementById("mpParam");

    const matPropDict = {

        //[name, energygap273K, me, mp, dielectricconstant, energygap300K, fermilevel, ec]
        1: ["Ge", 0.67, 0.2, 0.3, 16, 0.66, -0.27],
        2: ["GaAs", 1.39, 0.072, 0.5, 13, 1.43, -0.1],
        3: ["GaSb", 0.67, 0.047, 0.5, 15, 0.68, -0.1],
        4: ["InSb", 0.16, 0.013, 0.6, 18, 0.17, -0.07],
        5: ["InAs", 0.33, 0.02, 0.4, 14.5, 0.36, -0.15],
        6: ["Si", 1.14, 0.33, 0.5, 12, 1.11, -0.6]

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

        case "4":

            setValues(4);
            break;

        case "5":

            setValues(5);
            break;

        case "6":

            setValues(6);
            break;
        
    }

    function setValues(num){

        eg.value = matPropDict[num][1];
        me.value = matPropDict[num][2];
        mp.value = matPropDict[num][3];
        ef.value = matPropDict[num][6];

    }

}

function updatePropertiesPresets(){

    const fieldLong = "2.413e+02, 2.511e+02, 2.622e+02, 2.706e+02, 2.803e+02";
    const fieldShort = "3.000e+02";
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
        return;

    }

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

        } else{

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

function updateWFName(){

    let pickMaterialTypeDiv = document.getElementById("pickMaterialType");
    let pickChangingVarDiv = document.getElementById("pickChangingVar");
    let wfName = document.getElementById("wf_name");
    let advancedModeToggleMainDiv = document.getElementById("advancedModeToggleMainDiv");
    let advancedModeToggle = document.getElementById("advancedModeToggle");
    let advancedModeSemi = document.getElementById("advancedParametersSemi");
    let advancedModeMetal = document.getElementById("advancedParametersMetal");



    if(pickMaterialTypeDiv.value == 1){

        //wfName.textContent = "Work Function";
        //pickChangingVarDiv.options[3].textContent = "Work Function";
        //advancedModeToggleMainDiv.hidden = true;

        if(advancedModeToggle.checked == true){

            advancedModeMetal.hidden = false;
            advancedModeSemi.hidden = true;

        } else{

            advancedModeMetal.hidden = true;
            advancedModeSemi.hidden = true;

        }

    } else if(pickMaterialTypeDiv.value == 2){

        //wfName.textContent = "Work Function";
        //pickChangingVarDiv.options[3].textContent = "Work Function";
        
        if(advancedModeToggle.checked == true){

            advancedModeMetal.hidden = true;
            advancedModeSemi.hidden = false;

        }   else {

            advancedModeMetal.hidden = true;
            advancedModeSemi.hidden = true;

        }

    }

}

function updateAutoGenerateValuesDiv(){

    let autoGenerateValuesDiv = document.getElementById("autoGeneratedInput");
    let pickChangingVarGenerationMethodDiv = document.getElementById("autoGenerateValuesToggle");


    if(pickChangingVarGenerationMethodDiv.checked == true){

        autoGenerateValuesDiv.hidden = false;

    } else if(pickChangingVarGenerationMethodDiv.checked == false){

        autoGenerateValuesDiv.hidden = true;

    }

}

function raiseInputError(id) {

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
        case "3000":
            addErrorDiv("Can not create graph as all inputs must have equal amount of data points");
            break;
        case "3001":
            addErrorDiv("All input fields must have at least 3 values each");
            break;
        case "3002":
            addErrorDiv("Sweeping parameter must have at least 3 values");
            break;
        case "3003":
            addErrorDiv("Only sweeping parameter accepts multiple values");
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
    
        let relativeDiv = document.getElementById("TopBarView2");
        let errorDiv = document.createElement("div");
        errorDiv.innerHTML = template;
    
        let section = errorDiv.children[0];
    
        errorDivs.push(section)
    
        insertAfter(section, relativeDiv);
    
    }

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

    processArray() {

        this.data = this.data.slice(1, -1);
        this.processString();

    }

    processUnknownInput() {

        raiseInputError("2007");
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
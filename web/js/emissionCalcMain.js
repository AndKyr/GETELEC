function main(){

    let socket = io();

    let materialType, field, radius, gamma, workFunction, temperature,
        fieldMult, radiusMult, gammaMult, workFunctionMult, temperatureMult,
        calculateNH, calculateES, calculateEC

    loadInitEventListeners();

    function checkValidity(){

        if (errorDivs.length > 0) {

            errorDivs.forEach(div => {
                div.remove();
            })

            errorDivs = [];

        }

        materialType = document.getElementById("pickMaterialType").value;
        field = document.getElementById("field_in").value;
        radius = document.getElementById("radius_in").value;
        gamma = document.getElementById("gamma_in").value;
        workFunction = document.getElementById("wf_in").value;
        temperature = document.getElementById("temperature_in").value;

        fieldMult = document.getElementById("field_mult_in").value;
        radiusMult = document.getElementById("radius_mult_in").value;
        //gammaMult = document.getElementById("gamma_mult_in").value;
        workFunctionMult = document.getElementById("wf_mult_in").value;
        temperatureMult = document.getElementById("temperature_mult_in").value;

        calculateEC = document.getElementById("chooseEC").checked;
        calculateES = document.getElementById("chooseES").checked;
        calculateNH = document.getElementById("chooseNH").checked;

        if(field == "") field = "1, 5, 20"
        if(radius == "") radius = "15, 25, 50"
        if(gamma == "") gamma = "1, 10, 100"
        if(workFunction == "") workFunction = "10, 10, 10"
        if(temperature == "") temperature = "299.99, 300, 300.01"

        let _field = processDataInput(field);
        let _radius = processDataInput(radius);
        let _workFunction = processDataInput(workFunction);
        let _gamma = processDataInput(gamma);
        let _temperature = processDataInput(temperature);

        if(calculateEC){
            calculateEC = "1";
        } else {
            calculateEC = "0";
        }

        if(calculateES){
            calculateES = "1";
        } else {
            calculateES = "0";
        }

        if(calculateNH){
            calculateNH = "1";
        } else {
            calculateNH = "0";
        }

        for(let i = 0; i < _field.length; i++){
            _field[i] = _field[i] * fieldMult;
        }

        for(let i = 0; i < _radius.length; i++){
            _radius[i] = _radius[i] * radiusMult;
        }

        for(let i = 0; i < _workFunction.length; i++){
            _workFunction[i] = _workFunction[i] * workFunctionMult;
        }

        // for(let i = 0; i < _gamma.length; i++){
        //     _gamma[i] = _gamma[i] * gammaMult;
        // }

        for(let i = 0; i < _temperature.length; i++){
            _temperature[i] = _temperature[i] * temperatureMult;
        }

        if(materialType == "metal") materialType == 0;
        if(materialType == "semiconductor") materialType == 1;

        let canCompute = true;

        //THIS IS BROKEN
        if( (_field.length == _radius.length_ == _workFunction.length == _gamma.length == _temperature.length)){

            canCompute = false;

        }
        
        else if(_field.length < 3 || _radius.length < 3 || _workFunction.length < 3 || _gamma.length < 3 || _temperature.length < 3){


            canCompute = false;

        } else {

            if(canCompute == true){

                //HERE WILL CHECK FOR VALID INPUT RANGE

                let data = [materialType, _field, _radius, _gamma, _workFunction, _temperature, calculateNH, calculateES, calculateEC];
                
                socket.emit('calculateEmission', data);
                // $('#loadingModal').modal('show');

            }
        }
        
    }

    function loadInitEventListeners(){

        let enterButton = document.getElementById("enterButton");

        enterButton.addEventListener("click", checkValidity);
        
    }
}

let errorDivs = [];
main();


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

    let _dsb = beautifyResult(data.data.split(sp));
    let dsb = [];

    if (cm != ".") {
        for (let i = 0; i < _dsb.length; i++) {
            _dsb[i] = _dsb[i].replace(cm, ".");
        }
    }

    _dsb.forEach(el => {
        dsb.push(parseFloat(el));
    });

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
            raiseInputError("2006");
            console.log("The following line has no separator: ");
            console.log(this);
            return;
        }

        if (this.comma == 99) {
            /*            console.log("The following line has no commas: ");
                        console.log(this);
                        return;*/
        }

        let separator = this.separator;

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
            return 0;
        }
    } else if (getSeparator(data) == 0) {
        if (data.includes(",")) {
            return 1;
        }
    } else if (getSeparator(data) == 2) {
        if (data.includes(",")) {
            return 1
        } else if (data.includes(".")) {
            return 0
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

        if (data[spacePos - 1] == ".") {
            return 0;
        } else if (data[spacePos - 1] == ",") {
            return 1
        } else {
            if (!isNaN(data[spacePos - 1])) {
                return 2
            }
            return 99;
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




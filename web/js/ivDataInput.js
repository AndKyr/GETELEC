import {
    raiseInputError
} from "./ivCalculationsMain.js"

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

export function processIvDataInput(input) {

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

    if (endsWith(data, " ")) {
        data.slice(0, -1);
    }

    if (data[0] == " ") {
        data.slice(1);
    }

    return data;
}

function endsWith(str, suffix) {

    return str.indexOf(suffix, str.length - suffix.length) !== -1;

}
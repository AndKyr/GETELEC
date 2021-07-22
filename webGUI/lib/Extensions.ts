




declare global {
    
    
    
    interface Array<T> {
    
        // @ts-ignore
        removeElementAtIndex(index: number);
    
        // @ts-ignore
        removeElement(element: T);
    
        // @ts-ignore
        insertElementAtIndex(index: number, element: T);
        
        contains(element: T): boolean;
        
        groupedBy(keyFunction: (item: T) => any): { [key: string]: Array<T> } & Object;
        
        copy(): Array<T>;
        
        arrayByRepeating(numberOfRepetitions: number): Array<T>;
        arrayByTrimmingToLengthIfLonger(maxLength: number): Array<T>;
        
        anyMatch(predicate: (value: T, index: number, obj: T[]) => boolean): boolean
        noneMatch(predicate: (value: T, index: number, obj: T[]) => boolean): boolean
        
        allMatch(predicate: (value: T, index: number, obj: T[]) => boolean): boolean
    
        // @ts-ignore
        readonly firstElement: T;
        // @ts-ignore
        readonly lastElement: T;
        readonly summedValue: T;
        
        max(): number;
        min(): number;
        
    }
    
    
    
    
    
    interface String {
    
        // @ts-ignore
        contains(string): boolean;
        
        readonly numericalValue: number;
        isAString: boolean;
        
    }
    
    
    
    interface Number {
        
        isANumber: boolean;
        
        readonly integerValue: number;
        
    }
    
    
    
    interface Date {
        
        readonly dateString: string;
        
    }
    
    
    interface Object {
        
        forEachValue(callbackFunction: (value: any, key: string, index: number) => void): void;
        
        readonly allValues: Array<any>;
        
        readonly allKeys: string[];
        
    }
    
    
}

export {};



const YES = true;
const NO = false;







Array.prototype.removeElementAtIndex = function (this: any[], index: number) {
    
    if (index >= 0 && index < this.length) {
        
        this.splice(index, 1);
        
    }
    
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "removeElementAtIndex", { enumerable: false });






// @ts-ignore
Array.prototype.removeElement = function (this: any[], element) {
    
    this.removeElementAtIndex(this.indexOf(element));
    
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "removeElement", { enumerable: false });






if ("insertElementAtIndex" in Array.prototype == NO) {
    
    (Array.prototype as any).insertElementAtIndex = function (this: Array<any>, index: number, element: any) {
        
        if (index >= 0 && index <= this.length) {
            
            this.splice(index, 0, element)
            
        }
        
    }
    
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "insertElementAtIndex", { enumerable: false });





if ("replaceElementAtIndex" in Array.prototype == NO) {
    
    (Array.prototype as any).replaceElementAtIndex = function (this: Array<any>, index: number, element: any) {
        
        this.removeElementAtIndex(index)
        this.insertElementAtIndex(index, element)
        
    }
    
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "replaceElementAtIndex", { enumerable: false });





// @ts-ignore
Array.prototype.contains = function (this: any[], element) {
    
    var result = (this.indexOf(element) != -1);
    return result;
    
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "contains", { enumerable: false });



// @ts-ignore
Array.prototype.groupedBy = function (this: any[], funcProp) {
    return this.reduce(function (acc, val) {
        (acc[funcProp(val)] = acc[funcProp(val)] || []).push(val);
        return acc;
    }, {});
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "groupedBy", { enumerable: false });


// @ts-ignore
if (!Array.prototype.flat) {
    
    Object.defineProperty(Array.prototype, 'flat', {
        value: function(depth = 1) {
            // @ts-ignore
            return this.reduce(function (flat, toFlatten) {
                // @ts-ignore
                return flat.concat((Array.isArray(toFlatten) && (depth>1)) ? toFlatten.flat(depth-1) : toFlatten);
            }, []);
        }
    })
    
}

Array.prototype.max = function() {
    return Math.max.apply(null, this);
};

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "max", { enumerable: false });

Array.prototype.min = function() {
    return Math.min.apply(null, this);
};

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "min", { enumerable: false });








Object.defineProperty(Array.prototype, "firstElement", {
    
    get: function firstElement(this: Array<any>) {
        var result = this[0];
        return result;
    },
    enumerable: false
    
});



Object.defineProperty(Array.prototype, "lastElement", {
    
    get: function lastElement(this: any[]) {
        var result = this[this.length - 1];
        return result;
    },
    enumerable: false
    
});




(Array.prototype as any).copy = function (this: any[]) {
    
    var result = this.slice(0);
    return result;
    
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "copy", { enumerable: false });



Array.prototype.arrayByRepeating = function (this: any[], numberOfRepetitions: number) {
    // @ts-ignore
    var result = [];
    for (var i = 0; i < numberOfRepetitions; i++) {
        this.forEach(function (element, index, array) {
            result.push(element);
        })
    }
    // @ts-ignore
    return result;
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "arrayByRepeating", { enumerable: false });




Array.prototype.arrayByTrimmingToLengthIfLonger = function (this: any[], maxLength: number) {
    var result = [];
    for (var i = 0; i < maxLength && i < this.length; i++) {
        result.push(this[i]);
    }
    return result;
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "arrayByTrimmingToLengthIfLonger", { enumerable: false });



Object.defineProperty(Array.prototype, "summedValue", {
    
    get: function summedValue(this: any[]) {
        var result = this.reduce(function (a, b) { return a + b; }, 0);
        return result;
    },
    enumerable: false
    
});





Array.prototype.anyMatch = function (this: any[], functionToCall: (value: any, index: number, array: any[]) => boolean) {
    
    var result = (this.findIndex(functionToCall) > -1);
    
    return result;
    
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "anyMatch", { enumerable: false });



Array.prototype.noneMatch = function (this: any[], functionToCall: (value: any, index: number, array: any[]) => boolean) {
    
    var result = (this.findIndex(functionToCall) == -1);
    
    return result;
    
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "noneMatch", { enumerable: false });


Array.prototype.allMatch = function (this: any[], functionToCall: (value: any, index: number, array: any[]) => boolean) {
    
    function reversedFunction(value: any, index: number, array: any[]) {
        return !functionToCall(value, index, array);
    }
    
    var result = (this.findIndex(reversedFunction) == -1);
    
    return result;
    
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "allMatch", { enumerable: false });





// Object.prototype.forEach = function (this: Object, callbackFunction: (value: any) => void) {
//     var keys = Object.keys(this);
//     keys.forEach(function (key, index, array) {
//         callbackFunction(this[key]);
//     }.bind(this))
// }








if ("contains" in String.prototype == NO) {
    
    // @ts-ignore
    (String.prototype as any).contains = function (this: String, string) {
        
        var result = (this.indexOf(string) != -1);
        return result;
        
    }
    
    // Hide method from for-in loops
    Object.defineProperty(String.prototype, "contains", { enumerable: false });
    
}



if ("numericalValue" in String.prototype == NO) {
    
    Object.defineProperty(String.prototype, "numericalValue", {
        
        get: function numericalValue(this: string) {
            var result = Number(this);
            return result;
        },
        enumerable: false
        
    });
    
}




if ("isAString" in String.prototype == NO) {
    
    (String.prototype as any).isAString = YES;
    
}




if ("isANumber" in Number.prototype == NO) {
    
    (Number.prototype as any).isANumber = YES;
    
}



if ("integerValue" in Number.prototype == NO) {
    
    Object.defineProperty(Number.prototype, "integerValue", {
        
        get: function (this: number) {
            var result = parseInt("" + (Math.round(this) + 0.5));
            return result;
        },
        enumerable: false
        
    });
    
}





if ("dateString" in Date.prototype == NO) {
    
    Object.defineProperty(Date.prototype, "dateString", {
        
        get: function dateString(this: Date) {
            
            var result = ("0" + this.getDate()).slice(-2) + "-" + ("0" + (this.getMonth() + 1)).slice(-2) + "-" + this.getFullYear() + " " + ("0" + this.getHours()).slice(-2) + ":" + ("0" + this.getMinutes()).slice(-2);
            
            return result;
            
        },
        enumerable: false
        
    });
    
    
    
}

if ("forEachValue" in Object.prototype == NO) {
    
    (Object.prototype as any).forEachValue = function (this: Object, callbackFunction: (value: any, key: string, index: number) => void) {
        
        if (this) {
            
            var keys = Object.keys(this);
            // @ts-ignore
            keys.forEach(function (key, index, array) {
                
                // @ts-ignore
                callbackFunction(this[key], key, index);
                
            }.bind(this))
            
        }
        
    }
    
    // Hide method from for-in loops
    Object.defineProperty(Object.prototype, "forEachValue", { enumerable: false });
    
}


if ("allValues" in Object.prototype == NO) {
    
    Object.defineProperty(Object.prototype, "allValues", {
        
        get: function (this: Object) {
            // @ts-ignore
            var values = [];
            this.forEachValue(function (value: any) {
                values.push(value);
            })
            // @ts-ignore
            return values;
        },
        enumerable: false
        
    });
    
}


if ("allKeys" in Object.prototype == NO) {
    
    Object.defineProperty(Object.prototype, "allKeys", {
        
        get: function (this: Object) {
            var values = Object.keys(this);
            return values;
        },
        enumerable: false
        
    });
    
}























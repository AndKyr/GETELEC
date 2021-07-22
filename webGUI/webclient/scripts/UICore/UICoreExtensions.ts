if ("removeElementAtIndex" in Array.prototype == NO) {
    
    (Array.prototype as any).removeElementAtIndex = function (this: Array<any>, index: number) {
        
        if (index >= 0 && index < this.length) {
            
            this.splice(index, 1)
            
        }
        
    }
    
}


interface Array<T> {
    
    removeElementAtIndex(index: number);
    
}


if ("removeElement" in Array.prototype == NO) {
    
    (Array.prototype as any).removeElement = function (this: Array<any>, element) {
        
        this.removeElementAtIndex(this.indexOf(element))
        
    }
    
}


interface Array<T> {
    
    removeElement(element: T);
    
}


if ("insertElementAtIndex" in Array.prototype == NO) {
    
    (Array.prototype as any).insertElementAtIndex = function (this: Array<any>, index: number, element: any) {
        
        if (index >= 0 && index <= this.length) {
            
            this.splice(index, 0, element)
            
        }
        
    }
    
}


interface Array<T> {
    
    insertElementAtIndex(index: number, element: T);
    
}


if ("replaceElementAtIndex" in Array.prototype == NO) {
    
    (Array.prototype as any).replaceElementAtIndex = function (this: Array<any>, index: number, element: any) {
        
        this.removeElementAtIndex(index)
        this.insertElementAtIndex(index, element)
        
    }
    
}


interface Array<T> {
    
    replaceElementAtIndex(index: number, element: T);
    
}


if ("contains" in Array.prototype == NO) {
    
    (Array.prototype as any).contains = function (this: Array<any>, element) {
        
        const result = (this.indexOf(element) != -1)
        return result
        
    }
    
}

if ("containsAny" in Array.prototype == NO) {
    
    (Array.prototype as any).containsAny = function (this: Array<any>, elements: any[]) {
        
        const result = this.anyMatch(function (element, index, array) {
            return elements.contains(element)
        })
        
        return result
        
    }
    
}


interface Array<T> {
    
    contains(element: T): boolean;
    
    containsAny(element: T[]): boolean;
    
}

Array.prototype.max = function() {
    return Math.max.apply(null, this);
};

Array.prototype.min = function() {
    return Math.min.apply(null, this);
};

interface Array<T> {
    
    readonly summedValue: T;
    
    max(): number;
    min(): number;
    
    
}


if ("anyMatch" in Array.prototype == NO) {
    
    (Array.prototype as any).anyMatch = function (
        this: Array<any>,
        functionToCall: (value: any, index: number, array: any[]) => boolean
    ) {
        
        const result = (this.findIndex(functionToCall) > -1)
        
        return result
        
    }
    
}

if ("noneMatch" in Array.prototype == NO) {
    
    (Array.prototype as any).noneMatch = function (
        this: Array<any>,
        functionToCall: (value: any, index: number, array: any[]) => boolean
    ) {
        
        const result = (this.findIndex(functionToCall) == -1)
        
        return result
        
    }
    
}

if ("allMatch" in Array.prototype == NO) {
    
    (Array.prototype as any).allMatch = function (
        this: Array<any>,
        functionToCall: (value: any, index: number, array: any[]) => boolean
    ) {
        
        function reversedFunction(value: any, index: number, array: any[]) {
            return !functionToCall(value, index, array)
        }
        
        const result = (this.findIndex(reversedFunction) == -1)
        
        return result
        
    }
    
}


interface Array<T> {
    
    anyMatch(predicate: (value: T, index: number, obj: T[]) => boolean): boolean
    
    noneMatch(predicate: (value: T, index: number, obj: T[]) => boolean): boolean
    
    allMatch(predicate: (value: T, index: number, obj: T[]) => boolean): boolean
    
}


if ("groupedBy" in Array.prototype == NO) {
    
    Array.prototype.groupedBy = function (this: Array<any>, funcProp) {
        return this.reduce(function (acc, val) {
            (acc[funcProp(val)] = acc[funcProp(val)] || []).push(val)
            return acc
        }, {})
    }
    
}


interface Array<T> {
    
    groupedBy(keyFunction: (item: T) => any): { [key: string]: Array<T> };
    
}


if ("firstElement" in Array.prototype == NO) {
    Object.defineProperty(Array.prototype, "firstElement", {
        get: function firstElement(this: Array<any>) {
            const result = this[0]
            return result
        },
        set: function (this: Array<any>, element: any) {
            if (this.length == 0) {
                this.push(element)
                return
            }
            this[0] = element
        }
    })
}

if ("lastElement" in Array.prototype == NO) {
    Object.defineProperty(Array.prototype, "lastElement", {
        get: function lastElement(this: Array<any>) {
            const result = this[this.length - 1]
            return result
        },
        set: function (this: Array<any>, element: any) {
            if (this.length == 0) {
                this.push(element)
                return
            }
            this[this.length - 1] = element
        }
    })
}

if ("everyElement" in Array.prototype == NO) {
    
    Object.defineProperty(Array.prototype, "everyElement", {
        
        get: function everyElement(this: Array<any>) {
            
            var valueKeys = []
            
            const targetFunction = (objects) => {
    
                return this.map((element, index, array) => {
        
                    var elementFunction = (UIObject.valueForKeyPath(valueKeys.join("."), element) as Function).bind(
                        element,
                        objects
                    )
        
                    return elementFunction()
        
                });
    
            }
            
            const result = new Proxy(
                targetFunction,
                {
                    
                    get: (target, key, receiver) => {
    
                        if (key == "UI_elementValues") {
        
                            return this.map((element, index, array) => UIObject.valueForKeyPath(
                                valueKeys.join("."),
                                element
                            ));
                            
                        }
                        
                        valueKeys.push(key)
                        
                        return result
                        
                    },
                    set: (target, key, value, receiver) => {
                        
                        valueKeys.push(key)
                        
                        this.forEach((element, index, array) => {
                            
                            UIObject.setValueForKeyPath(valueKeys.join("."), value, element, YES)
                            
                        })
                        
                        return true
                        
                    }
                    
                }
            )
            
            return result
            
        },
        set: function (this: Array<any>, element: any) {
            
            for (var i = 0; i < this.length; ++i) {
                
                this[i] = element
                
            }
            
        }
        
    })
    
}


type UIEveryElementItem<T> = {
    
    [P in keyof T]: UIEveryElementItem<T[P]>
    
} & {
    
    UI_elementValues?: T[];
    
} & T

interface Array<T> {
    
    firstElement: T;
    lastElement: T;
    
    everyElement: UIEveryElementItem<T>;
    
}


if ("copy" in Array.prototype == NO) {
    
    (Array.prototype as any).copy = function (this: Array<any>) {
        
        const result = this.slice(0)
        return result
        
    }
    
}


interface Array<T> {
    
    copy(): Array<T>;
    
}


if ("arrayByRepeating" in Array.prototype == NO) {
    
    (Array.prototype as any).arrayByRepeating = function (this: Array<any>, numberOfRepetitions: number) {
        const result = []
        for (var i = 0; i < numberOfRepetitions; i++) {
            this.forEach(function (element, index, array) {
                result.push(element)
            })
        }
        return result
    }
    
}


interface Array<T> {
    
    arrayByRepeating(numberOfRepetitions: number): Array<T>;
    
}


if ("arrayByTrimmingToLengthIfLonger" in Array.prototype == NO) {
    (Array.prototype as any).arrayByTrimmingToLengthIfLonger = function (this: Array<any>, maxLength: number) {
        const result = []
        for (var i = 0; i < maxLength && i < this.length; i++) {
            result.push(this[i])
        }
        return result
    }
}


interface Array<T> {
    
    arrayByTrimmingToLengthIfLonger(maxLength: number): Array<T>;
    
}


if ("summedValue" in Array.prototype == NO) {
    
    Object.defineProperty(Array.prototype, "summedValue", {
        get: function summedValue(this: Array<any>) {
            const result = this.reduce(function (a, b) {
                return a + b
            }, 0)
            return result
        }
    })
    
}


interface Array<T> {
    
    readonly summedValue: T;
    
}


// Warn if overriding existing method
if ("isEqualToArray" in Array.prototype == YES) {
    console.warn(
        "Overriding existing Array.prototype.isEqualToArray. Possible causes: New API defines the method, there's a framework conflict or you've got double inclusions in your code.")
}
// attach the .equals method to Array's prototype to call it on any array
Array.prototype.isEqualToArray = function (array: any[], keyPath?: string) {
    
    // if the other array is a falsy value, return
    if (!array) {
        return false
    }
    
    // compare lengths - can save a lot of time 
    if (this.length != array.length) {
        return false
    }
    
    var i = 0
    const l = this.length
    for (; i < l; i++) {
        
        // Check if we have nested arrays
        if (this[i] instanceof Array && array[i] instanceof Array && !keyPath) {
            
            // recurse into the nested arrays
            if (!this[i].isEqualToArray(array[i])) {
                return false
            }
            
        }
        else if (keyPath && UIObject.valueForKeyPath(keyPath, this[i]) != UIObject.valueForKeyPath(keyPath, array[i])) {
            
            return false
            
        }
        else if (this[i] != array[i]) {
            
            // Warning - two different object instances will never be equal: {x:20} != {x:20}
            return false
            
        }
        
    }
    
    return true
    
}

// Hide method from for-in loops
Object.defineProperty(Array.prototype, "isEqualToArray", { enumerable: false })


interface Array<T> {
    
    isEqualToArray(array: Array<T>, keyPath?: string): boolean;
    
}



if ("forEach" in Object.prototype == NO) {
    
    (Object.prototype as any).forEach = function (this: Object, callbackFunction: (value: any, key: string) => void) {
        const keys = Object.keys(this)
        keys.forEach(function (key, index, array) {
            callbackFunction(this[key], key)
        }.bind(this))
    }
    
    // Hide method from for-in loops
    Object.defineProperty(Object.prototype, "forEach", { enumerable: false })
    
}


interface Object {
    
    forEach(callbackFunction: (value: any, key: string) => void): void;
    
}


if ("allValues" in Object.prototype == NO) {
    Object.defineProperty(Object.prototype, "allValues", {
        get: function (this: Object) {
            const values = []
            this.forEach(function (value: any) {
                values.push(value)
            })
            return values
        }
    })
}


interface Object {
    
    readonly allValues: Array<any>;
    
}


if ("allKeys" in Object.prototype == NO) {
    Object.defineProperty(Object.prototype, "allKeys", {
        get: function (this: Object) {
            const values = Object.keys(this)
            return values
        }
    })
}


interface Object {
    
    readonly allKeys: string[];
    
}


if ("objectByCopyingValuesRecursivelyFromObject" in Object.prototype == NO) {
    
    (Object.prototype as any).objectByCopyingValuesRecursivelyFromObject = function (this: Object, object: any) {
        
        
        function isAnObject(item: any) {
            return (item && typeof item === "object" && !Array.isArray(item))
        }
        
        function mergeRecursively(target: any, source: any) {
            
            const output = Object.assign({}, target)
            
            if (isAnObject(target) && isAnObject(source)) {
                
                Object.keys(source).forEach(function (key) {
                    
                    if (isAnObject(source[key])) {
                        
                        // if (!(key in target)) {
                        
                        //     Object.assign(output, { [key]: source[key] });
                        
                        // }
                        // else {
                        
                        output[key] = mergeRecursively(target[key], source[key])
                        
                        //}
                        
                    }
                    else {
                        
                        Object.assign(output, { [key]: source[key] })
                        
                    }
                    
                })
                
            }
            
            return output
            
        }
        
        const result = mergeRecursively(this, object)
        
        return result
        
    }
    
    // Hide method from for-in loops
    Object.defineProperty(Object.prototype, "objectByCopyingValuesRecursivelyFromObject", { enumerable: false })
    
}

if ("copyValuesRecursivelyFromObject" in Object.prototype == NO) {
    
    (Object.prototype as any).copyValuesRecursivelyFromObject = function (this: Object, object: any) {
        
        
        function isAnObject(item: any) {
            return (item && typeof item === "object" && !Array.isArray(item))
        }
        
        function mergeRecursively(target: any, source: any) {
            
            const output = target;
            
            if (isAnObject(target) && isAnObject(source)) {
                
                Object.keys(source).forEach(function (key) {
                    
                    if (isAnObject(source[key])) {
                        
                        // if (!(key in target)) {
                        
                        //     Object.assign(output, { [key]: source[key] });
                        
                        // }
                        // else {
                        
                        output[key] = mergeRecursively(target[key], source[key])
                        
                        //}
                        
                    }
                    else {
                        
                        Object.assign(output, { [key]: source[key] })
                        
                    }
                    
                })
                
            }
            
            return output
            
        }
        
        const result = mergeRecursively(this, object)
        
        return result
        
    }
    
    // Hide method from for-in loops
    Object.defineProperty(Object.prototype, "copyValuesRecursivelyFromObject", { enumerable: false })
    
}


interface Object {
    
    // objectByCopyingValuesRecursivelyFromObject<T>(object: T): this & T;
    copyValuesRecursivelyFromObject(object: Partial<this>): this;
    
}

interface UIObject {
   
   // objectByCopyingValuesRecursivelyFromObject<T>(object: T): this & T;
    copyValuesRecursivelyFromObject(object: RecursivePartial<this>): this;
    
}



if ("contains" in String.prototype == NO) {
    
    (String.prototype as any).contains = function (this: String, string) {
        
        const result = (this.indexOf(string) != -1)
        return result
        
    }
    
}


interface String {
    
    contains(string): boolean;
    
}


if ("capitalizedString" in String.prototype == NO) {
    Object.defineProperty(Object.prototype, "capitalizedString", {
        get: function (this: String) {
            const result = this.charAt(0).toUpperCase() + this.slice(1).toLowerCase()
            return result
        }
    })
}


interface String {
    
    readonly capitalizedString: string;
    
}


if ("numericalValue" in String.prototype == NO) {
    Object.defineProperty(String.prototype, "numericalValue", {
        get: function numericalValue(this: string) {
            const result = Number(this)
            return result
        }
    })
}


interface String {
    
    readonly numericalValue: number;
    
}


if ("isAString" in String.prototype == NO) {
    
    (String.prototype as any).isAString = YES
    
}


interface String {
    
    isAString: boolean;
    
}


if ("isANumber" in Number.prototype == NO) {
    
    (Number.prototype as any).isANumber = YES
    
}


interface Number {
    
    isANumber: boolean;
    
}





if ("integerValue" in Number.prototype == NO) {
    Object.defineProperty(Number.prototype, "integerValue", {
        get: function (this: number) {
            const result = parseInt("" + (Math.round(this) + 0.5))
            return result
        }
    })
}


interface Number {
    
    readonly integerValue: number;
    
}



class PrimitiveNumber {
    
    // @ts-ignore
    static [Symbol.hasInstance](x) {
        return
    }
    
}


if ("integerValue" in Boolean.prototype == NO) {
    
    Object.defineProperty(Boolean.prototype, "integerValue", {
        get: function (this: boolean) {
            
            if (this == true) {
                
                return 1
                
            }
            
            return 0
            
        }
    })
    
}


interface Boolean {
    
    readonly integerValue: number;
    
}


if ("dateString" in Date.prototype == NO) {
    
    Object.defineProperty(Date.prototype, "dateString", {
        get: function dateString(this: Date) {
            
            const result = ("0" + this.getDate()).slice(-2) + "-" + ("0" + (this.getMonth() + 1)).slice(-2) + "-" +
                this.getFullYear() + " " + ("0" + this.getHours()).slice(-2) + ":" +
                ("0" + this.getMinutes()).slice(-2)
            
            return result
            
        }
    })
    
    
    
}


interface Date {
    
    readonly dateString: string;
    
}









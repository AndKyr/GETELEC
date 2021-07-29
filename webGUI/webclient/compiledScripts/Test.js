var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
function NilFunction() {
    return nil;
}
var nil = new Proxy(Object.assign(NilFunction, { "class": nil, "className": "Nil" }), {
    get(target, name) {
        if (name == Symbol.toPrimitive) {
            return function (hint) {
                if (hint == "number") {
                    return 0;
                }
                if (hint == "string") {
                    return "";
                }
                return false;
            };
        }
        if (name == "toString") {
            return function toString() {
                return "";
            };
        }
        return NilFunction();
    },
    set(target, name, value) {
        return NilFunction();
    }
});
function wrapInNil(object) {
    var result = FIRST_OR_NIL(object);
    if (object instanceof Object && !(object instanceof Function)) {
        result = new Proxy(object, {
            get(target, name) {
                if (name == "wrapped_nil_target") {
                    return target;
                }
                const value = Reflect.get(target, name);
                if (typeof value === "object") {
                    return wrapInNil(value);
                }
                if (IS_NOT_LIKE_NULL(value)) {
                    return value;
                }
                return nil;
            }
        });
    }
    return result;
}
const YES = true;
const NO = false;
function IS(object) {
    if (object && object !== nil) {
        return YES;
    }
    return NO;
    //return (object != nil && object);
}
function IS_NOT(object) {
    return !IS(object);
}
function IS_DEFINED(object) {
    if (object != undefined) {
        return YES;
    }
    return NO;
}
function IS_UNDEFINED(object) {
    return !IS_DEFINED(object);
}
function IS_NIL(object) {
    if (object === nil) {
        return YES;
    }
    return NO;
}
function IS_NOT_NIL(object) {
    return !IS_NIL(object);
}
function IS_LIKE_NULL(object) {
    return (IS_UNDEFINED(object) || IS_NIL(object) || object == null);
}
function IS_NOT_LIKE_NULL(object) {
    return !IS_LIKE_NULL(object);
}
function IS_AN_EMAIL_ADDRESS(email) {
    const re = /\S+@\S+\.\S+/;
    return re.test(email);
}
function FIRST_OR_NIL(...objects) {
    const result = objects.find(function (object, index, array) {
        return IS(object);
    });
    return result || nil;
}
function FIRST(...objects) {
    const result = objects.find(function (object, index, array) {
        return IS(object);
    });
    return result || IF(IS_DEFINED(objects.lastElement))(RETURNER(objects.lastElement))();
}
function MAKE_ID(randomPartLength = 15) {
    var result = "";
    const characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    for (var i = 0; i < randomPartLength; i++) {
        result = result + characters.charAt(Math.floor(Math.random() * characters.length));
    }
    result = result + Date.now();
    return result;
}
function RETURNER(value) {
    return function (...objects) {
        return value;
    };
}
function IF(value) {
    var thenFunction = nil;
    var elseFunction = nil;
    const result = function (functionToCall) {
        thenFunction = functionToCall;
        return result.evaluateConditions;
    };
    result.evaluateConditions = function () {
        if (IS(value)) {
            return thenFunction();
        }
        return elseFunction();
    };
    result.evaluateConditions.ELSE_IF = function (otherValue) {
        const functionResult = IF(otherValue);
        elseFunction = functionResult.evaluateConditions;
        const functionResultEvaluateConditionsFunction = function () {
            return result.evaluateConditions();
        };
        functionResultEvaluateConditionsFunction.ELSE_IF = functionResult.evaluateConditions.ELSE_IF;
        functionResultEvaluateConditionsFunction.ELSE = functionResult.evaluateConditions.ELSE;
        functionResult.evaluateConditions = functionResultEvaluateConditionsFunction;
        return functionResult;
    };
    result.evaluateConditions.ELSE = function (functionToCall) {
        elseFunction = functionToCall;
        return result.evaluateConditions();
    };
    return result;
}
// @ts-ignore
if (!window.AutoLayout) {
    // @ts-ignore
    window.AutoLayout = nil;
}
class UIObject {
    constructor() {
        this._class = UIObject;
        this.superclass = nil.class;
    }
    get class() {
        return this.constructor;
    }
    get superclass() {
        return this.constructor.superclass;
    }
    set superclass(superclass) {
        this.constructor.superclass = superclass;
    }
    static wrapObject(object) {
        if (IS_NOT(object)) {
            return nil;
        }
        if (object instanceof UIObject) {
            return object;
        }
        const result = Object.assign(new UIObject(), object);
        return result;
    }
    isKindOfClass(classObject) {
        if (this.isMemberOfClass(classObject)) {
            return YES;
        }
        for (var superclassObject = this.superclass; IS(superclassObject); superclassObject = superclassObject.superclass) {
            if (superclassObject == classObject) {
                return YES;
            }
        }
        return NO;
    }
    isMemberOfClass(classObject) {
        return (this.class == classObject);
    }
    valueForKey(key) {
        return this[key];
    }
    valueForKeyPath(keyPath) {
        return UIObject.valueForKeyPath(keyPath, this);
    }
    static valueForKeyPath(keyPath, object) {
        if (IS_NOT(keyPath)) {
            return object;
        }
        const keys = keyPath.split(".");
        var currentObject = object;
        for (var i = 0; i < keys.length; i++) {
            const key = keys[i];
            if (key.substring(0, 2) == "[]") {
                // This next object will be an array and the rest of the keys need to be run for each of the elements
                currentObject = currentObject[key.substring(2)];
                // CurrentObject is now an array
                const remainingKeyPath = keys.slice(i + 1).join(".");
                const currentArray = currentObject;
                currentObject = currentArray.map(function (subObject, index, array) {
                    const result = UIObject.valueForKeyPath(remainingKeyPath, subObject);
                    return result;
                });
                break;
            }
            currentObject = currentObject[key];
            if (IS_NOT(currentObject)) {
                currentObject = nil;
            }
        }
        return currentObject;
    }
    setValueForKeyPath(keyPath, value, createPath = YES) {
        return UIObject.setValueForKeyPath(keyPath, value, this, createPath);
    }
    static setValueForKeyPath(keyPath, value, currentObject, createPath) {
        const keys = keyPath.split(".");
        var didSetValue = NO;
        keys.forEach(function (key, index, array) {
            if (index == array.length - 1 && IS_NOT_LIKE_NULL(currentObject)) {
                currentObject[key] = value;
                didSetValue = YES;
                return;
            }
            else if (IS_NOT(currentObject)) {
                return;
            }
            const currentObjectValue = currentObject[key];
            if (IS_LIKE_NULL(currentObjectValue) && createPath) {
                currentObject[key] = {};
            }
            currentObject = currentObject[key];
        });
        return didSetValue;
    }
    performFunctionWithSelf(functionToPerform) {
        return functionToPerform(this);
    }
    performFunctionWithDelay(delay, functionToCall) {
        new UITimer(delay, NO, functionToCall);
    }
}
///<reference path="UIObject.ts"/>
class UIColor extends UIObject {
    constructor(stringValue) {
        super();
        this.stringValue = stringValue;
        this._class = UIColor;
        this.superclass = UIObject;
    }
    toString() {
        return this.stringValue;
    }
    static get redColor() {
        return new UIColor("red");
    }
    static get blueColor() {
        return new UIColor("blue");
    }
    static get greenColor() {
        return new UIColor("green");
    }
    static get yellowColor() {
        return new UIColor("yellow");
    }
    static get blackColor() {
        return new UIColor("black");
    }
    static get brownColor() {
        return new UIColor("brown");
    }
    static get whiteColor() {
        return new UIColor("white");
    }
    static get greyColor() {
        return new UIColor("grey");
    }
    static get lightGreyColor() {
        return new UIColor("lightgrey");
    }
    static get transparentColor() {
        return new UIColor("transparent");
    }
    static get undefinedColor() {
        return new UIColor("");
    }
    static get nilColor() {
        return new UIColor("");
    }
    static nameToHex(name) {
        return {
            "aliceblue": "#f0f8ff",
            "antiquewhite": "#faebd7",
            "aqua": "#00ffff",
            "aquamarine": "#7fffd4",
            "azure": "#f0ffff",
            "beige": "#f5f5dc",
            "bisque": "#ffe4c4",
            "black": "#000000",
            "blanchedalmond": "#ffebcd",
            "blue": "#0000ff",
            "blueviolet": "#8a2be2",
            "brown": "#a52a2a",
            "burlywood": "#deb887",
            "cadetblue": "#5f9ea0",
            "chartreuse": "#7fff00",
            "chocolate": "#d2691e",
            "coral": "#ff7f50",
            "cornflowerblue": "#6495ed",
            "cornsilk": "#fff8dc",
            "crimson": "#dc143c",
            "cyan": "#00ffff",
            "darkblue": "#00008b",
            "darkcyan": "#008b8b",
            "darkgoldenrod": "#b8860b",
            "darkgray": "#a9a9a9",
            "darkgreen": "#006400",
            "darkkhaki": "#bdb76b",
            "darkmagenta": "#8b008b",
            "darkolivegreen": "#556b2f",
            "darkorange": "#ff8c00",
            "darkorchid": "#9932cc",
            "darkred": "#8b0000",
            "darksalmon": "#e9967a",
            "darkseagreen": "#8fbc8f",
            "darkslateblue": "#483d8b",
            "darkslategray": "#2f4f4f",
            "darkturquoise": "#00ced1",
            "darkviolet": "#9400d3",
            "deeppink": "#ff1493",
            "deepskyblue": "#00bfff",
            "dimgray": "#696969",
            "dodgerblue": "#1e90ff",
            "firebrick": "#b22222",
            "floralwhite": "#fffaf0",
            "forestgreen": "#228b22",
            "fuchsia": "#ff00ff",
            "gainsboro": "#dcdcdc",
            "ghostwhite": "#f8f8ff",
            "gold": "#ffd700",
            "goldenrod": "#daa520",
            "gray": "#808080",
            "green": "#008000",
            "greenyellow": "#adff2f",
            "honeydew": "#f0fff0",
            "hotpink": "#ff69b4",
            "indianred ": "#cd5c5c",
            "indigo": "#4b0082",
            "ivory": "#fffff0",
            "khaki": "#f0e68c",
            "lavender": "#e6e6fa",
            "lavenderblush": "#fff0f5",
            "lawngreen": "#7cfc00",
            "lemonchiffon": "#fffacd",
            "lightblue": "#add8e6",
            "lightcoral": "#f08080",
            "lightcyan": "#e0ffff",
            "lightgoldenrodyellow": "#fafad2",
            "lightgrey": "#d3d3d3",
            "lightgreen": "#90ee90",
            "lightpink": "#ffb6c1",
            "lightsalmon": "#ffa07a",
            "lightseagreen": "#20b2aa",
            "lightskyblue": "#87cefa",
            "lightslategray": "#778899",
            "lightsteelblue": "#b0c4de",
            "lightyellow": "#ffffe0",
            "lime": "#00ff00",
            "limegreen": "#32cd32",
            "linen": "#faf0e6",
            "magenta": "#ff00ff",
            "maroon": "#800000",
            "mediumaquamarine": "#66cdaa",
            "mediumblue": "#0000cd",
            "mediumorchid": "#ba55d3",
            "mediumpurple": "#9370d8",
            "mediumseagreen": "#3cb371",
            "mediumslateblue": "#7b68ee",
            "mediumspringgreen": "#00fa9a",
            "mediumturquoise": "#48d1cc",
            "mediumvioletred": "#c71585",
            "midnightblue": "#191970",
            "mintcream": "#f5fffa",
            "mistyrose": "#ffe4e1",
            "moccasin": "#ffe4b5",
            "navajowhite": "#ffdead",
            "navy": "#000080",
            "oldlace": "#fdf5e6",
            "olive": "#808000",
            "olivedrab": "#6b8e23",
            "orange": "#ffa500",
            "orangered": "#ff4500",
            "orchid": "#da70d6",
            "palegoldenrod": "#eee8aa",
            "palegreen": "#98fb98",
            "paleturquoise": "#afeeee",
            "palevioletred": "#d87093",
            "papayawhip": "#ffefd5",
            "peachpuff": "#ffdab9",
            "peru": "#cd853f",
            "pink": "#ffc0cb",
            "plum": "#dda0dd",
            "powderblue": "#b0e0e6",
            "purple": "#800080",
            "red": "#ff0000",
            "rosybrown": "#bc8f8f",
            "royalblue": "#4169e1",
            "saddlebrown": "#8b4513",
            "salmon": "#fa8072",
            "sandybrown": "#f4a460",
            "seagreen": "#2e8b57",
            "seashell": "#fff5ee",
            "sienna": "#a0522d",
            "silver": "#c0c0c0",
            "skyblue": "#87ceeb",
            "slateblue": "#6a5acd",
            "slategray": "#708090",
            "snow": "#fffafa",
            "springgreen": "#00ff7f",
            "steelblue": "#4682b4",
            "tan": "#d2b48c",
            "teal": "#008080",
            "thistle": "#d8bfd8",
            "tomato": "#ff6347",
            "turquoise": "#40e0d0",
            "violet": "#ee82ee",
            "wheat": "#f5deb3",
            "white": "#ffffff",
            "whitesmoke": "#f5f5f5",
            "yellow": "#ffff00",
            "yellowgreen": "#9acd32"
        }[name.toLowerCase()];
    }
    static hexToDescriptor(c) {
        if (c[0] === "#") {
            c = c.substr(1);
        }
        const r = parseInt(c.slice(0, 2), 16);
        const g = parseInt(c.slice(2, 4), 16);
        const b = parseInt(c.slice(4, 6), 16);
        const a = parseInt(c.slice(6, 8), 16);
        const result = { "red": r, "green": g, "blue": b, "alpha": a };
        return result;
        //return 'rgb(' + r + ',' + g + ',' + b + ')';
    }
    static rgbToDescriptor(colorString) {
        if (colorString.startsWith("rgba(")) {
            colorString = colorString.slice(5, colorString.length - 1);
        }
        if (colorString.startsWith("rgb(")) {
            colorString = colorString.slice(4, colorString.length - 1) + ", 0";
        }
        const components = colorString.split(",");
        const result = {
            "red": Number(components[0]),
            "green": Number(components[1]),
            "blue": Number(components[2]),
            "alpha": Number(components[3])
        };
        return result;
    }
    get colorDescriptor() {
        var descriptor;
        const colorHEXFromName = UIColor.nameToHex(this.stringValue);
        if (this.stringValue.startsWith("rgb")) {
            descriptor = UIColor.rgbToDescriptor(this.stringValue);
        }
        else if (colorHEXFromName) {
            descriptor = UIColor.hexToDescriptor(colorHEXFromName);
        }
        else {
            descriptor = UIColor.hexToDescriptor(this.stringValue);
        }
        return descriptor;
    }
    colorWithRed(red) {
        const descriptor = this.colorDescriptor;
        const result = new UIColor("rgba(" + red + "," + descriptor.green + "," + descriptor.blue + "," +
            descriptor.alpha + ")");
        return result;
    }
    colorWithGreen(green) {
        const descriptor = this.colorDescriptor;
        const result = new UIColor("rgba(" + descriptor.red + "," + green + "," + descriptor.blue + "," +
            descriptor.alpha + ")");
        return result;
    }
    colorWithBlue(blue) {
        const descriptor = this.colorDescriptor;
        const result = new UIColor("rgba(" + descriptor.red + "," + descriptor.green + "," + blue + "," +
            descriptor.alpha + ")");
        return result;
    }
    colorWithAlpha(alpha) {
        const descriptor = this.colorDescriptor;
        const result = new UIColor("rgba(" + descriptor.red + "," + descriptor.green + "," + descriptor.blue + "," +
            alpha + ")");
        return result;
    }
    static colorWithRGBA(red, green, blue, alpha = 1) {
        const result = new UIColor("rgba(" + red + "," + green + "," + blue + "," + alpha + ")");
        return result;
    }
    static colorWithDescriptor(descriptor) {
        const result = new UIColor("rgba(" + descriptor.red.toFixed(0) + "," + descriptor.green.toFixed(0) + "," +
            descriptor.blue.toFixed(0) + "," + this.defaultAlphaToOne(descriptor.alpha) + ")");
        return result;
    }
    static defaultAlphaToOne(value = 1) {
        if (value != value) {
            value = 1;
        }
        return value;
    }
    colorByMultiplyingRGB(multiplier) {
        const descriptor = this.colorDescriptor;
        descriptor.red = descriptor.red * multiplier;
        descriptor.green = descriptor.green * multiplier;
        descriptor.blue = descriptor.blue * multiplier;
        const result = UIColor.colorWithDescriptor(descriptor);
        return result;
    }
}
/// <reference path="./UIObject.ts" />
class UIPoint extends UIObject {
    constructor(x, y) {
        super();
        this.x = x;
        this.y = y;
        this._class = UIPoint;
        this.superclass = UIObject;
    }
    copy() {
        return new UIPoint(this.x, this.y);
    }
    isEqualTo(point) {
        const result = (this.x == point.x && this.y == point.y);
        return result;
    }
    scale(zoom) {
        const x = this.x;
        const y = this.y;
        this.x = x * zoom;
        this.y = y * zoom;
        return this;
    }
    add(v) {
        this.x = this.x + v.x;
        this.y = this.y + v.y;
        return this;
    }
    subtract(v) {
        this.x = this.x - v.x;
        this.y = this.y - v.y;
        return this;
    }
    to(b) {
        const a = this;
        const ab = b.copy().add(a.copy().scale(-1));
        return ab;
    }
    pointWithX(x) {
        const result = this.copy();
        result.x = x;
        return result;
    }
    pointWithY(y) {
        const result = this.copy();
        result.y = y;
        return result;
    }
    pointByAddingX(x) {
        return this.pointWithX(this.x + x);
    }
    pointByAddingY(y) {
        return this.pointWithY(this.y + y);
    }
    get length() {
        var result = this.x * this.x + this.y * this.y;
        result = Math.sqrt(result);
        return result;
    }
    didChange(b) {
        // Callback to be set by delegate
    }
}
/// <reference path="./UIPoint.ts" />
class UIRectangle extends UIObject {
    constructor(x = 0, y = 0, height = 0, width = 0) {
        super();
        this._class = UIRectangle;
        this.superclass = UIObject;
        this.min = new UIPoint(Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY);
        this.max = new UIPoint(Number.NEGATIVE_INFINITY, Number.NEGATIVE_INFINITY);
        this.min.didChange = this.rectanglePointDidChange;
        this.max.didChange = this.rectanglePointDidChange;
        this._isBeingUpdated = NO;
        this.min = new UIPoint(x, y);
        this.max = new UIPoint(x + width, y + height);
        if (IS_NIL(height)) {
            this.max.y = height;
        }
        if (IS_NIL(width)) {
            this.max.x = width;
        }
    }
    copy() {
        const result = new UIRectangle(this.x, this.y, this.height, this.width);
        return result;
    }
    isEqualTo(rectangle) {
        const result = (IS(rectangle) && this.min.isEqualTo(rectangle.min) && this.max.isEqualTo(rectangle.max));
        return result;
    }
    static zero() {
        const result = new UIRectangle(0, 0, 0, 0);
        return result;
    }
    containsPoint(point) {
        return this.min.x <= point.x && this.min.y <= point.y &&
            point.x <= this.max.x && point.y <= this.max.y;
    }
    updateByAddingPoint(point) {
        if (!point) {
            point = new UIPoint(0, 0);
        }
        this.beginUpdates();
        const min = this.min.copy();
        if (min.x === nil) {
            min.x = this.max.x;
        }
        if (min.y === nil) {
            min.y = this.max.y;
        }
        const max = this.max.copy();
        if (max.x === nil) {
            max.x = this.min.x;
        }
        if (max.y === nil) {
            max.y = this.min.y;
        }
        this.min.x = Math.min(min.x, point.x);
        this.min.y = Math.min(min.y, point.y);
        this.max.x = Math.max(max.x, point.x);
        this.max.y = Math.max(max.y, point.y);
        this.finishUpdates();
    }
    get height() {
        if (this.max.y === nil) {
            return nil;
        }
        return this.max.y - this.min.y;
    }
    set height(height) {
        this.max.y = this.min.y + height;
    }
    get width() {
        if (this.max.x === nil) {
            return nil;
        }
        return this.max.x - this.min.x;
    }
    set width(width) {
        this.max.x = this.min.x + width;
    }
    get x() {
        return this.min.x;
    }
    set x(x) {
        this.beginUpdates();
        const width = this.width;
        this.min.x = x;
        this.max.x = this.min.x + width;
        this.finishUpdates();
    }
    get y() {
        return this.min.y;
    }
    set y(y) {
        this.beginUpdates();
        const height = this.height;
        this.min.y = y;
        this.max.y = this.min.y + height;
        this.finishUpdates();
    }
    get topLeft() {
        return this.min.copy();
    }
    get topRight() {
        return new UIPoint(this.max.x, this.y);
    }
    get bottomLeft() {
        return new UIPoint(this.x, this.max.y);
    }
    get bottomRight() {
        return this.max.copy();
    }
    get center() {
        const result = this.min.copy().add(this.min.to(this.max).scale(0.5));
        return result;
    }
    set center(center) {
        const offset = this.center.to(center);
        this.offsetByPoint(offset);
    }
    offsetByPoint(offset) {
        this.min.add(offset);
        this.max.add(offset);
        return this;
    }
    concatenateWithRectangle(rectangle) {
        this.updateByAddingPoint(rectangle.bottomRight);
        this.updateByAddingPoint(rectangle.topLeft);
        return this;
    }
    intersectionRectangleWithRectangle(rectangle) {
        const result = this.copy();
        result.beginUpdates();
        const min = result.min;
        if (min.x === nil) {
            min.x = rectangle.max.x - Math.min(result.width, rectangle.width);
        }
        if (min.y === nil) {
            min.y = rectangle.max.y - Math.min(result.height, rectangle.height);
        }
        const max = result.max;
        if (max.x === nil) {
            max.x = rectangle.min.x + Math.min(result.width, rectangle.width);
        }
        if (max.y === nil) {
            max.y = rectangle.min.y + Math.min(result.height, rectangle.height);
        }
        result.min.x = Math.max(result.min.x, rectangle.min.x);
        result.min.y = Math.max(result.min.y, rectangle.min.y);
        result.max.x = Math.min(result.max.x, rectangle.max.x);
        result.max.y = Math.min(result.max.y, rectangle.max.y);
        if (result.height < 0) {
            const averageY = (this.center.y + rectangle.center.y) * 0.5;
            result.min.y = averageY;
            result.max.y = averageY;
        }
        if (result.width < 0) {
            const averageX = (this.center.x + rectangle.center.x) * 0.5;
            result.min.x = averageX;
            result.max.x = averageX;
        }
        result.finishUpdates();
        return result;
    }
    get area() {
        const result = this.height * this.width;
        return result;
    }
    intersectsWithRectangle(rectangle) {
        return (this.intersectionRectangleWithRectangle(rectangle).area != 0);
    }
    // add some space around the rectangle
    rectangleWithInsets(left, right, bottom, top) {
        const result = this.copy();
        result.min.x = this.min.x + left;
        result.max.x = this.max.x - right;
        result.min.y = this.min.y + top;
        result.max.y = this.max.y - bottom;
        return result;
    }
    rectangleWithInset(inset) {
        const result = this.rectangleWithInsets(inset, inset, inset, inset);
        return result;
    }
    rectangleWithHeight(height, centeredOnPosition = nil) {
        if (isNaN(centeredOnPosition)) {
            centeredOnPosition = nil;
        }
        const result = this.copy();
        result.height = height;
        if (centeredOnPosition != nil) {
            const change = height - this.height;
            result.offsetByPoint(new UIPoint(0, change * centeredOnPosition).scale(-1));
        }
        return result;
    }
    rectangleWithWidth(width, centeredOnPosition = nil) {
        if (isNaN(centeredOnPosition)) {
            centeredOnPosition = nil;
        }
        const result = this.copy();
        result.width = width;
        if (centeredOnPosition != nil) {
            const change = width - this.width;
            result.offsetByPoint(new UIPoint(change * centeredOnPosition, 0).scale(-1));
        }
        return result;
    }
    rectangleWithHeightRelativeToWidth(heightRatio = 1, centeredOnPosition = nil) {
        const result = this.rectangleWithHeight(this.width * heightRatio, centeredOnPosition);
        return result;
    }
    rectangleWithWidthRelativeToHeight(widthRatio = 1, centeredOnPosition = nil) {
        const result = this.rectangleWithWidth(this.height * widthRatio, centeredOnPosition);
        return result;
    }
    rectangleWithX(x, centeredOnPosition = 0) {
        const result = this.copy();
        result.x = x - result.width * centeredOnPosition;
        return result;
    }
    rectangleWithY(y, centeredOnPosition = 0) {
        const result = this.copy();
        result.y = y - result.height * centeredOnPosition;
        return result;
    }
    rectangleByAddingX(x) {
        const result = this.copy();
        result.x = this.x + x;
        return result;
    }
    rectangleByAddingY(y) {
        const result = this.copy();
        result.y = this.y + y;
        return result;
    }
    rectanglesBySplittingWidth(weights, paddings = 0, absoluteWidths = nil) {
        if (IS_NIL(paddings)) {
            paddings = 1;
        }
        if (!(paddings instanceof Array)) {
            paddings = [paddings].arrayByRepeating(weights.length - 1);
        }
        paddings = paddings.arrayByTrimmingToLengthIfLonger(weights.length - 1);
        if (!(absoluteWidths instanceof Array) && IS_NOT_NIL(absoluteWidths)) {
            absoluteWidths = [absoluteWidths].arrayByRepeating(weights.length);
        }
        const result = [];
        const sumOfWeights = weights.reduce(function (a, b, index) {
            if (IS_NOT_NIL(absoluteWidths[index])) {
                b = 0;
            }
            return a + b;
        }, 0);
        const sumOfPaddings = paddings.summedValue;
        const sumOfAbsoluteWidths = absoluteWidths.summedValue;
        const totalRelativeWidth = this.width - sumOfPaddings - sumOfAbsoluteWidths;
        var previousCellMaxX = this.x;
        for (var i = 0; i < weights.length; i++) {
            var resultWidth;
            if (IS_NOT_NIL(absoluteWidths[i])) {
                resultWidth = absoluteWidths[i] || 0;
            }
            else {
                resultWidth = totalRelativeWidth * (weights[i] / sumOfWeights);
            }
            const rectangle = this.rectangleWithWidth(resultWidth);
            var padding = 0;
            if (paddings.length > i && paddings[i]) {
                padding = paddings[i];
            }
            rectangle.x = previousCellMaxX;
            previousCellMaxX = rectangle.max.x + padding;
            //rectangle = rectangle.rectangleWithInsets(0, padding, 0, 0);
            result.push(rectangle);
        }
        return result;
    }
    rectanglesBySplittingHeight(weights, paddings = 0, absoluteHeights = nil) {
        if (IS_NIL(paddings)) {
            paddings = 1;
        }
        if (!(paddings instanceof Array)) {
            paddings = [paddings].arrayByRepeating(weights.length - 1);
        }
        paddings = paddings.arrayByTrimmingToLengthIfLonger(weights.length - 1);
        if (!(absoluteHeights instanceof Array) && IS_NOT_NIL(absoluteHeights)) {
            absoluteHeights = [absoluteHeights].arrayByRepeating(weights.length);
        }
        const result = [];
        const sumOfWeights = weights.reduce(function (a, b, index) {
            if (IS_NOT_NIL(absoluteHeights[index])) {
                b = 0;
            }
            return a + b;
        }, 0);
        const sumOfPaddings = paddings.summedValue;
        const sumOfAbsoluteHeights = absoluteHeights.summedValue;
        const totalRelativeHeight = this.height - sumOfPaddings - sumOfAbsoluteHeights;
        var previousCellMaxY = this.y;
        for (var i = 0; i < weights.length; i++) {
            var resultHeight;
            if (IS_NOT_NIL(absoluteHeights[i])) {
                resultHeight = absoluteHeights[i] || 0;
            }
            else {
                resultHeight = totalRelativeHeight * (weights[i] / sumOfWeights);
            }
            const rectangle = this.rectangleWithHeight(resultHeight);
            var padding = 0;
            if (paddings.length > i && paddings[i]) {
                padding = paddings[i];
            }
            rectangle.y = previousCellMaxY;
            previousCellMaxY = rectangle.max.y + padding;
            //rectangle = rectangle.rectangleWithInsets(0, 0, padding, 0);
            result.push(rectangle);
        }
        return result;
    }
    rectanglesByEquallySplittingWidth(numberOfFrames, padding = 0) {
        const result = [];
        const totalPadding = padding * (numberOfFrames - 1);
        const resultWidth = (this.width - totalPadding) / numberOfFrames;
        for (var i = 0; i < numberOfFrames; i++) {
            const rectangle = this.rectangleWithWidth(resultWidth, i / (numberOfFrames - 1));
            result.push(rectangle);
        }
        return result;
    }
    rectanglesByEquallySplittingHeight(numberOfFrames, padding = 0) {
        const result = [];
        const totalPadding = padding * (numberOfFrames - 1);
        const resultHeight = (this.height - totalPadding) / numberOfFrames;
        for (var i = 0; i < numberOfFrames; i++) {
            const rectangle = this.rectangleWithHeight(resultHeight, i / (numberOfFrames - 1));
            result.push(rectangle);
        }
        return result;
    }
    distributeViewsAlongWidth(views, weights = 1, paddings, absoluteWidths) {
        if (!(weights instanceof Array)) {
            weights = [weights].arrayByRepeating(views.length);
        }
        const frames = this.rectanglesBySplittingWidth(weights, paddings, absoluteWidths);
        frames.forEach(function (frame, index, array) {
            views[index].frame = frame;
        });
        return this;
    }
    distributeViewsAlongHeight(views, weights = 1, paddings, absoluteHeights) {
        if (!(weights instanceof Array)) {
            weights = [weights].arrayByRepeating(views.length);
        }
        const frames = this.rectanglesBySplittingHeight(weights, paddings, absoluteHeights);
        frames.forEach(function (frame, index, array) {
            views[index].frame = frame;
        });
        return this;
    }
    distributeViewsEquallyAlongWidth(views, padding) {
        const frames = this.rectanglesByEquallySplittingWidth(views.length, padding);
        frames.forEach(function (frame, index, array) {
            views[index].frame = frame;
        });
        return this;
    }
    distributeViewsEquallyAlongHeight(views, padding) {
        const frames = this.rectanglesByEquallySplittingHeight(views.length, padding);
        frames.forEach(function (frame, index, array) {
            views[index].frame = frame;
        });
        return this;
    }
    rectangleForNextRow(padding = 0, height = this.height) {
        const result = this.rectangleWithY(this.max.y + padding);
        if (height != this.height) {
            result.height = height;
        }
        return result;
    }
    rectangleForNextColumn(padding = 0, width = this.width) {
        const result = this.rectangleWithX(this.max.x + padding);
        if (width != this.width) {
            result.width = width;
        }
        return result;
    }
    rectangleForPreviousRow(padding = 0) {
        const result = this.rectangleWithY(this.min.y - this.height - padding);
        return result;
    }
    rectangleForPreviousColumn(padding = 0) {
        const result = this.rectangleWithX(this.min.x - this.width - padding);
        return result;
    }
    // Bounding box
    static boundingBoxForPoints(points) {
        const result = new UIRectangle();
        for (var i = 0; i < points.length; i++) {
            result.updateByAddingPoint(points[i]);
        }
        return result;
    }
    beginUpdates() {
        this._isBeingUpdated = YES;
    }
    finishUpdates() {
        this._isBeingUpdated = NO;
        this.didChange();
    }
    didChange() {
        // Callback to be set by delegate
    }
    _rectanglePointDidChange() {
        if (!this._isBeingUpdated) {
            this.didChange();
        }
    }
}
///// <reference path="../autolayout/src/AutoLayout.js" />
/// <reference path="./UIObject.ts" />
/// <reference path="./UIColor.ts" />
/// <reference path="./UIRectangle.ts" />
class UIView extends UIObject {
    constructor(elementID = ("UIView" +
        UIView.nextIndex), viewHTMLElement = null, elementType = null, initViewData) {
        super();
        this._nativeSelectionEnabled = YES;
        this._enabled = YES;
        this._backgroundColor = UIColor.transparentColor;
        this._localizedTextObject = nil;
        this._controlEventTargets = {}; //{ "PointerDown": Function[]; "PointerMove": Function[]; "PointerLeave": Function[]; "PointerEnter": Function[]; "PointerUpInside": Function[]; "PointerUp": Function[]; "PointerHover": Function[]; };
        this._viewControllerLayoutFunction = nil;
        this._isHidden = NO;
        this.pausesPointerEvents = NO;
        this.stopsPointerEventPropagation = YES;
        this._pointerDragThreshold = 2;
        this.ignoresTouches = NO;
        this.ignoresMouse = NO;
        this.forceIntrinsicSizeZero = NO;
        this.controlEvent = UIView.controlEvent;
        this._class = UIView;
        this.superclass = UIObject;
        // Instance variables
        UIView._UIViewIndex = UIView.nextIndex;
        this._UIViewIndex = UIView._UIViewIndex;
        this._styleClasses = [];
        // Object.defineProperty(this, "styleClasses", { get: this.styleClasses, set: this.setStyleClasses });
        // Object.defineProperty(this, "styleClassName", { get: this.styleClassName });
        this._initViewHTMLElement(elementID, viewHTMLElement, elementType);
        this.subviews = [];
        this.superview = nil;
        // Object.defineProperty(this, "elementID", { get: this.elementID });
        // Object.defineProperty(this, "constraints", { get: this.constraints, set: this.setConstraints });
        this._constraints = [];
        this._updateLayoutFunction = nil;
        //Object.defineProperty(this, "backgroundColor", { get: this.backgroundColor, set: this.setBackgroundColor });
        //this.backgroundColor = "transparent";
        // Object.defineProperty(this, "alpha", { get: this.alpha, set: this.setAlpha });
        // Object.defineProperty(this, "frame", { get: this.frame, set: this.setFrame });
        // Object.defineProperty(this, "bounds", { get: this.bounds, set: this.setBounds });
        // Object.defineProperty(this, "userInteractionEnabled", { get: this.userInteractionEnabled, set: this.setUserInteractionEnabled });
        // this._controlEventTargets = {
        //     "PointerDown": [],
        //     "PointerMove": [],
        //     "PointerLeave": [],
        //     "PointerEnter": [],
        //     "PointerUpInside": [],
        //     "PointerUp": [],
        //     "PointerHover": []
        // }
        this._didLayoutSubviewsDelegateFunction = function () {
        };
        this._frameTransform = "";
        this.initView(this.viewHTMLElement.id, this.viewHTMLElement, initViewData);
        this._initViewCSSSelectorsIfNeeded();
        this._loadUIEvents();
        this.setNeedsLayout();
    }
    static get nextIndex() {
        return UIView._UIViewIndex + 1;
    }
    static get pageHeight() {
        const body = document.body;
        const html = document.documentElement;
        const height = Math.max(body.scrollHeight, body.offsetHeight, html.clientHeight, html.scrollHeight, html.offsetHeight);
        return height;
    }
    static get pageWidth() {
        const body = document.body;
        const html = document.documentElement;
        const width = Math.max(body.scrollWidth, body.offsetWidth, html.clientWidth, html.scrollWidth, html.offsetWidth);
        return width;
    }
    initView(elementID, viewHTMLElement, initViewData) {
    }
    centerInContainer() {
        this.style.left = "50%";
        this.style.top = "50%";
        this.style.transform = "translateX(-50%) translateY(-50%)";
    }
    centerXInContainer() {
        this.style.left = "50%";
        this.style.transform = "translateX(-50%)";
    }
    centerYInContainer() {
        this.style.top = "50%";
        this.style.transform = "translateY(-50%)";
    }
    _initViewHTMLElement(elementID, viewHTMLElement, elementType = "div") {
        if (!IS(elementType)) {
            elementType = "div";
        }
        if (!IS(viewHTMLElement)) {
            this._viewHTMLElement = this.createElement(elementID, elementType);
            this.style.position = "absolute";
            this.style.margin = "0";
        }
        else {
            this._viewHTMLElement = viewHTMLElement;
        }
        if (IS(elementID)) {
            this.viewHTMLElement.id = elementID;
        }
        this.viewHTMLElement.obeyAutolayout = YES;
        this.viewHTMLElement.UIView = this;
        this.addStyleClass(this.styleClassName);
    }
    set nativeSelectionEnabled(selectable) {
        this._nativeSelectionEnabled = selectable;
        if (!selectable) {
            this.style.cssText = this.style.cssText +
                " -webkit-touch-callout: none; -webkit-user-select: none; -khtml-user-select: none; -moz-user-select: none; -ms-user-select: none; user-select: none;";
        }
        else {
            this.style.cssText = this.style.cssText +
                " -webkit-touch-callout: text; -webkit-user-select: text; -khtml-user-select: text; -moz-user-select: text; -ms-user-select: text; user-select: text;";
        }
    }
    get nativeSelectionEnabled() {
        return this._nativeSelectionEnabled;
    }
    get styleClassName() {
        const result = "UICore_UIView_" + this.class.name;
        return result;
    }
    _initViewCSSSelectorsIfNeeded() {
        if (!this.class._areViewCSSSelectorsInitialized) {
            this.initViewStyleSelectors();
            this.class._areViewCSSSelectorsInitialized = YES;
        }
    }
    initViewStyleSelectors() {
        // Override this in a subclass
    }
    initStyleSelector(selector, style) {
        const styleRules = UIView.getStyleRules(selector);
        if (!styleRules) {
            UIView.createStyleSelector(selector, style);
        }
    }
    createElement(elementID, elementType) {
        var result = document.getElementById(elementID);
        if (!result) {
            result = document.createElement(elementType);
        }
        return result;
    }
    get viewHTMLElement() {
        return this._viewHTMLElement;
    }
    get elementID() {
        return this.viewHTMLElement.id;
    }
    setInnerHTML(key, defaultString, parameters) {
        this._innerHTMLKey = key;
        this._defaultInnerHTML = defaultString;
        this._parameters = parameters;
        const languageName = UICore.languageService.currentLanguageKey;
        const result = UICore.languageService.stringForKey(key, languageName, defaultString, parameters);
        this.innerHTML = result;
    }
    _setInnerHTMLFromKeyIfPossible() {
        if (this._innerHTMLKey && this._defaultInnerHTML) {
            this.setInnerHTML(this._innerHTMLKey, this._defaultInnerHTML, this._parameters);
        }
    }
    _setInnerHTMLFromLocalizedTextObjectIfPossible() {
        if (IS(this._localizedTextObject)) {
            this.innerHTML = UICore.languageService.stringForCurrentLanguage(this._localizedTextObject);
        }
    }
    get localizedTextObject() {
        return this._localizedTextObject;
    }
    set localizedTextObject(localizedTextObject) {
        this._localizedTextObject = localizedTextObject;
        this._setInnerHTMLFromLocalizedTextObjectIfPossible();
    }
    get innerHTML() {
        return this.viewHTMLElement.innerHTML;
    }
    set innerHTML(innerHTML) {
        if (this.innerHTML != innerHTML) {
            this.viewHTMLElement.innerHTML = FIRST(innerHTML, "");
        }
    }
    set hoverText(hoverText) {
        this.viewHTMLElement.setAttribute("title", hoverText);
    }
    get hoverText() {
        return this.viewHTMLElement.getAttribute("title");
    }
    get scrollSize() {
        const result = new UIRectangle(0, 0, this.viewHTMLElement.scrollHeight, this.viewHTMLElement.scrollWidth);
        return result;
    }
    get dialogView() {
        if (!IS(this.superview)) {
            return nil;
        }
        if (!(this instanceof UIDialogView)) {
            return this.superview.dialogView;
        }
        return this;
    }
    get rootView() {
        if (IS(this.superview)) {
            return this.superview.rootView;
        }
        return this;
    }
    set enabled(enabled) {
        this._enabled = enabled;
        this.updateContentForCurrentEnabledState();
    }
    get enabled() {
        return this._enabled;
    }
    updateContentForCurrentEnabledState() {
        this.hidden = !this.enabled;
        this.userInteractionEnabled = this.enabled;
    }
    get tabIndex() {
        return Number(this.viewHTMLElement.getAttribute("tabindex"));
    }
    set tabIndex(index) {
        this.viewHTMLElement.setAttribute("tabindex", "" + index);
    }
    get styleClasses() {
        return this._styleClasses;
    }
    set styleClasses(styleClasses) {
        this._styleClasses = styleClasses;
    }
    hasStyleClass(styleClass) {
        // This is for performance reasons
        if (!IS(styleClass)) {
            return NO;
        }
        const index = this.styleClasses.indexOf(styleClass);
        if (index > -1) {
            return YES;
        }
        return NO;
    }
    addStyleClass(styleClass) {
        if (!IS(styleClass)) {
            return;
        }
        if (!this.hasStyleClass(styleClass)) {
            this._styleClasses.push(styleClass);
        }
    }
    removeStyleClass(styleClass) {
        // This is for performance reasons
        if (!IS(styleClass)) {
            return;
        }
        const index = this.styleClasses.indexOf(styleClass);
        if (index > -1) {
            this.styleClasses.splice(index, 1);
        }
    }
    static findViewWithElementID(elementID) {
        const viewHTMLElement = document.getElementById(elementID);
        if (IS_NOT(viewHTMLElement)) {
            return nil;
        }
        // @ts-ignore
        const result = viewHTMLElement.UIView;
        return result;
    }
    static createStyleSelector(selector, style) {
        return;
        // @ts-ignore
        if (!document.styleSheets) {
            return;
        }
        if (document.getElementsByTagName("head").length == 0) {
            return;
        }
        var styleSheet;
        var mediaType;
        if (document.styleSheets.length > 0) {
            for (var i = 0, l = document.styleSheets.length; i < l; i++) {
                if (document.styleSheets[i].disabled) {
                    continue;
                }
                const media = document.styleSheets[i].media;
                mediaType = typeof media;
                if (mediaType === "string") {
                    if (media === "" || (media.indexOf("screen") !== -1)) {
                        styleSheet = document.styleSheets[i];
                    }
                }
                else if (mediaType == "object") {
                    if (media.mediaText === "" || (media.mediaText.indexOf("screen") !== -1)) {
                        styleSheet = document.styleSheets[i];
                    }
                }
                if (typeof styleSheet !== "undefined") {
                    break;
                }
            }
        }
        if (typeof styleSheet === "undefined") {
            const styleSheetElement = document.createElement("style");
            styleSheetElement.type = "text/css";
            document.getElementsByTagName("head")[0].appendChild(styleSheetElement);
            for (i = 0; i < document.styleSheets.length; i++) {
                if (document.styleSheets[i].disabled) {
                    continue;
                }
                styleSheet = document.styleSheets[i];
            }
            mediaType = typeof styleSheet.media;
        }
        if (mediaType === "string") {
            for (var i = 0, l = styleSheet.rules.length; i < l; i++) {
                if (styleSheet.rules[i].selectorText && styleSheet.rules[i].selectorText.toLowerCase() ==
                    selector.toLowerCase()) {
                    styleSheet.rules[i].style.cssText = style;
                    return;
                }
            }
            styleSheet.addRule(selector, style);
        }
        else if (mediaType === "object") {
            var styleSheetLength = 0;
            try {
                styleSheetLength = (styleSheet.cssRules) ? styleSheet.cssRules.length : 0;
            }
            catch (error) {
            }
            for (var i = 0; i < styleSheetLength; i++) {
                if (styleSheet.cssRules[i].selectorText && styleSheet.cssRules[i].selectorText.toLowerCase() ==
                    selector.toLowerCase()) {
                    styleSheet.cssRules[i].style.cssText = style;
                    return;
                }
            }
            styleSheet.insertRule(selector + "{" + style + "}", styleSheetLength);
        }
    }
    static getStyleRules(selector) {
        var selector = selector.toLowerCase();
        for (var i = 0; i < document.styleSheets.length; i++) {
            const styleSheet = document.styleSheets[i];
            var styleRules;
            try {
                styleRules = styleSheet.cssRules ? styleSheet.cssRules : styleSheet.rules;
            }
            catch (error) {
            }
            return styleRules;
        }
    }
    get style() {
        return this.viewHTMLElement.style;
    }
    get computedStyle() {
        return getComputedStyle(this.viewHTMLElement);
    }
    get hidden() {
        return this._isHidden;
    }
    set hidden(v) {
        this._isHidden = v;
        if (this._isHidden) {
            this.style.visibility = "hidden";
        }
        else {
            this.style.visibility = "visible";
        }
    }
    static set pageScale(scale) {
        UIView._pageScale = scale;
        const zoom = scale;
        const width = 100 / zoom;
        const viewHTMLElement = UICore.main.rootViewController.view.viewHTMLElement;
        viewHTMLElement.style.transformOrigin = "left top";
        viewHTMLElement.style.transform = "scale(" + zoom + ")";
        viewHTMLElement.style.width = width + "%";
    }
    static get pageScale() {
        return UIView._pageScale;
    }
    calculateAndSetViewFrame() {
        // Use this method to calculate the frame for the view itself
        // This can be used when adding subviews to existing views like buttons
    }
    get frame() {
        // var result = new UIRectangle(1 * this.viewHTMLElement.offsetLeft, 1 * this.viewHTMLElement.offsetTop, 1 * this.viewHTMLElement.offsetHeight, 1 * this.viewHTMLElement.offsetWidth);
        // result.zIndex = 1 * this.style.zIndex;
        var result = this._frame;
        if (!result) {
            result = new UIRectangle(1 * this.viewHTMLElement.offsetLeft, 1 * this.viewHTMLElement.offsetTop, 1 *
                this.viewHTMLElement.offsetHeight, 1 * this.viewHTMLElement.offsetWidth);
            result.zIndex = 0;
        }
        return result.copy();
    }
    set frame(rectangle) {
        if (IS(rectangle)) {
            this.setFrame(rectangle);
        }
    }
    setFrame(rectangle, zIndex = 0, performUncheckedLayout = NO) {
        const frame = this._frame || new UIRectangle(nil, nil, nil, nil);
        if (zIndex != undefined) {
            rectangle.zIndex = zIndex;
        }
        this._frame = rectangle;
        // This is useless because frames are copied
        // frame.didChange = function () {
        //     // Do nothing
        // }
        // rectangle.didChange = function () {
        //     this.frame = rectangle;
        // }.bind(this);
        if (frame && frame.isEqualTo(rectangle) && !performUncheckedLayout) {
            return;
        }
        UIView._setAbsoluteSizeAndPosition(this.viewHTMLElement, rectangle.topLeft.x, rectangle.topLeft.y, rectangle.width, rectangle.height, rectangle.zIndex);
        if (frame.height != rectangle.height || frame.width != rectangle.width || performUncheckedLayout) {
            this.setNeedsLayout();
            this.boundsDidChange();
            //this.layoutSubviews();
        }
    }
    get bounds() {
        var result;
        // if (IS_NOT(this._frame) && this.style.height == "" && this.style.width  == "" && this.style.left == "" && this.style.right == "" && this.style.bottom == "" && this.style.top == "") {
        //     result = new UIRectangle(0, 0, 0, 0)
        // }
        // else
        if (IS_NOT(this._frame)) {
            result = new UIRectangle(0, 0, 1 * this.viewHTMLElement.offsetHeight, 1 * this.viewHTMLElement.offsetWidth);
        }
        else {
            result = this.frame.copy();
            result.x = 0;
            result.y = 0;
        }
        return result;
    }
    set bounds(rectangle) {
        const frame = this.frame;
        this.frame = new UIRectangle(frame.topLeft.x, frame.topLeft.y, rectangle.height, rectangle.width);
    }
    boundsDidChange() {
    }
    setPosition(left = nil, right = nil, bottom = nil, top = nil, height = nil, width = nil) {
        const previousBounds = this.bounds;
        this.setStyleProperty("left", left);
        this.setStyleProperty("right", right);
        this.setStyleProperty("bottom", bottom);
        this.setStyleProperty("top", top);
        this.setStyleProperty("height", height);
        this.setStyleProperty("width", width);
        const bounds = this.bounds;
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout();
            this.boundsDidChange();
        }
    }
    setSizes(height, width) {
        const previousBounds = this.bounds;
        this.setStyleProperty("height", height);
        this.setStyleProperty("width", width);
        const bounds = this.bounds;
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout();
            this.boundsDidChange();
        }
    }
    setMinSizes(height, width) {
        const previousBounds = this.bounds;
        this.setStyleProperty("minHeight", height);
        this.setStyleProperty("minWidth", width);
        const bounds = this.bounds;
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout();
            this.boundsDidChange();
        }
    }
    setMaxSizes(height, width) {
        const previousBounds = this.bounds;
        this.setStyleProperty("maxHeight", height);
        this.setStyleProperty("maxWidth", width);
        const bounds = this.bounds;
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout();
            this.boundsDidChange();
        }
    }
    setMargin(margin) {
        const previousBounds = this.bounds;
        this.setStyleProperty("margin", margin);
        const bounds = this.bounds;
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout();
            this.boundsDidChange();
        }
    }
    setMargins(left, right, bottom, top) {
        const previousBounds = this.bounds;
        this.setStyleProperty("marginLeft", left);
        this.setStyleProperty("marginRight", right);
        this.setStyleProperty("marginBottom", bottom);
        this.setStyleProperty("marginTop", top);
        const bounds = this.bounds;
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout();
            this.boundsDidChange();
        }
    }
    setPadding(padding) {
        const previousBounds = this.bounds;
        this.setStyleProperty("padding", padding);
        const bounds = this.bounds;
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout();
            this.boundsDidChange();
        }
    }
    setPaddings(left, right, bottom, top) {
        const previousBounds = this.bounds;
        this.setStyleProperty("paddingLeft", left);
        this.setStyleProperty("paddingRight", right);
        this.setStyleProperty("paddingBottom", bottom);
        this.setStyleProperty("paddingTop", top);
        const bounds = this.bounds;
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout();
            this.boundsDidChange();
        }
    }
    setBorder(radius = nil, width = 1, color = UIColor.blackColor, style = "solid") {
        this.setStyleProperty("borderStyle", style);
        this.setStyleProperty("borderRadius", radius);
        this.setStyleProperty("borderColor", color.stringValue);
        this.setStyleProperty("borderWidth", width);
    }
    setStyleProperty(propertyName, value) {
        try {
            if (IS_NIL(value)) {
                return;
            }
            if (IS_DEFINED(value) && value.isANumber) {
                value = "" + value.integerValue + "px";
            }
            this.style[propertyName] = value;
        }
        catch (exception) {
            console.log(exception);
        }
    }
    get userInteractionEnabled() {
        const result = (this.style.pointerEvents != "none");
        return result;
    }
    set userInteractionEnabled(userInteractionEnabled) {
        if (userInteractionEnabled) {
            this.style.pointerEvents = "";
        }
        else {
            this.style.pointerEvents = "none";
        }
    }
    get backgroundColor() {
        return this._backgroundColor;
    }
    set backgroundColor(backgroundColor) {
        this._backgroundColor = backgroundColor;
        this.style.backgroundColor = backgroundColor.stringValue;
    }
    get alpha() {
        return 1 * this.style.opacity;
    }
    set alpha(alpha) {
        this.style.opacity = "" + alpha;
    }
    static animateViewOrViewsWithDurationDelayAndFunction(viewOrViews, duration, delay, timingStyle = "cubic-bezier(0.25,0.1,0.25,1)", transformFunction, transitioncompletionFunction) {
        function callTransitioncompletionFunction() {
            (transitioncompletionFunction || nil)();
            viewOrViews.forEach(function (view, index, array) {
                view.animationDidFinish();
            });
        }
        if (IS_FIREFOX) {
            // Firefox does not fire the transition completion event properly
            new UIObject().performFunctionWithDelay(delay + duration, callTransitioncompletionFunction);
        }
        if (!(viewOrViews instanceof Array)) {
            viewOrViews = [viewOrViews];
        }
        const transitionStyles = [];
        const transitionDurations = [];
        const transitionDelays = [];
        const transitionTimings = [];
        for (var i = 0; i < viewOrViews.length; i++) {
            var view = viewOrViews[i];
            if (view.viewHTMLElement) {
                view = view.viewHTMLElement;
            }
            view.addEventListener("transitionend", transitionDidFinish, true);
            transitionStyles.push(view.style.transition);
            transitionDurations.push(view.style.transitionDuration);
            transitionDelays.push(view.style.transitionDelay);
            transitionTimings.push(view.style.transitionTimingFunction);
            view.style.transition = "all";
            view.style.transitionDuration = "" + duration + "s";
            view.style.transitionDelay = "" + delay + "s";
            view.style.transitionTimingFunction = timingStyle;
        }
        transformFunction();
        const transitionObject = {
            "finishImmediately": finishTransitionImmediately,
            "didFinish": transitionDidFinishManually,
            "views": viewOrViews,
            "registrationTime": Date.now()
        };
        function finishTransitionImmediately() {
            for (var i = 0; i < viewOrViews.length; i++) {
                var view = viewOrViews[i];
                if (view.viewHTMLElement) {
                    view = view.viewHTMLElement;
                }
                view.style.transition = "all";
                view.style.transitionDuration = "" + duration + "s";
                view.style.transitionDelay = "" + delay + "s";
                view.style.transition = transitionStyles[i];
                view.style.transitionDuration = transitionDurations[i];
                view.style.transitionDelay = transitionDelays[i];
                view.style.transitionTimingFunction = transitionTimings[i];
            }
        }
        function transitionDidFinish(event) {
            var view = event.srcElement;
            if (!view) {
                return;
            }
            if (view.viewHTMLElement) {
                view = view.viewHTMLElement;
            }
            view.style.transition = transitionStyles[i];
            view.style.transitionDuration = transitionDurations[i];
            view.style.transitionDelay = transitionDelays[i];
            view.style.transitionTimingFunction = transitionTimings[i];
            callTransitioncompletionFunction();
            view.removeEventListener("transitionend", transitionDidFinish, true);
        }
        function transitionDidFinishManually() {
            for (var i = 0; i < viewOrViews.length; i++) {
                var view = viewOrViews[i];
                if (view.viewHTMLElement) {
                    view = view.viewHTMLElement;
                }
                view.style.transition = transitionStyles[i];
                view.style.transitionDuration = transitionDurations[i];
                view.style.transitionDelay = transitionDelays[i];
                view.style.transitionTimingFunction = transitionTimings[i];
                view.removeEventListener("transitionend", transitionDidFinish, true);
            }
        }
        return transitionObject;
    }
    animationDidFinish() {
    }
    static _setAbsoluteSizeAndPosition(element, left, top, width, height, zIndex = 0) {
        // if (!UIView._transformAttribute) {
        //     UIView._transformAttribute = (('transform' in document.documentElement.style) ? 'transform' : undefined);
        //     UIView._transformAttribute = UIView._transformAttribute || (('-webkit-transform' in document.documentElement.style) ? '-webkit-transform' : 'undefined');
        //     UIView._transformAttribute = UIView._transformAttribute || (('-moz-transform' in document.documentElement.style) ? '-moz-transform' : 'undefined');
        //     UIView._transformAttribute = UIView._transformAttribute || (('-ms-transform' in document.documentElement.style) ? '-ms-transform' : 'undefined');
        //     UIView._transformAttribute = UIView._transformAttribute || (('-o-transform' in document.documentElement.style) ? '-o-transform' : 'undefined');
        // }
        if (!IS(element) || !element.obeyAutolayout && !element.getAttribute("obeyAutolayout")) {
            return;
        }
        if (element.id == "mainView") {
            var asd = 1;
        }
        if (IS(height)) {
            height = height.integerValue + "px";
        }
        if (IS(width)) {
            width = width.integerValue + "px";
        }
        var str = element.style.cssText;
        const frameTransform = UIView._transformAttribute + ": translate3d(" + (1 * left).integerValue + "px, " +
            (1 * top).integerValue + "px, " + zIndex.integerValue + "px)";
        if (element.UIView) {
            str = str + frameTransform + ";";
        }
        else {
            element.UIView._frameTransform = frameTransform;
        }
        if (height == nil) {
            str = str + " height: unset;";
        }
        else {
            str = str + " height:" + height + ";";
        }
        if (width == nil) {
            str = str + " width: unset;";
        }
        else {
            str = str + " width:" + width + ";";
        }
        if (element.id == "mainView") {
            var asd = 1;
        }
        element.style.cssText = element.style.cssText + str;
    }
    static performAutoLayout(parentElement, visualFormatArray, constraintsArray) {
        const view = new AutoLayout.View();
        if (IS(visualFormatArray) && IS(visualFormatArray.length)) {
            view.addConstraints(AutoLayout.VisualFormat.parse(visualFormatArray, { extended: true }));
        }
        if (IS(constraintsArray) && IS(constraintsArray.length)) {
            view.addConstraints(constraintsArray);
        }
        const elements = {};
        for (var key in view.subViews) {
            if (!view.subViews.hasOwnProperty(key)) {
                continue;
            }
            var element = nil;
            try {
                element = parentElement.querySelector("#" + key);
            }
            catch (error) {
                //console.log("Error occurred " + error);
            }
            if (element && !element.obeyAutolayout && !element.getAttribute("obeyAutolayout")) {
            }
            else if (element) {
                element.className += element.className ? " abs" : "abs";
                elements[key] = element;
            }
        }
        var parentUIView = nil;
        if (parentElement.UIView) {
            parentUIView = parentElement.UIView;
        }
        const updateLayout = function () {
            view.setSize(parentElement ? parentElement.clientWidth : window.innerWidth, parentElement ? parentElement.clientHeight : window.innerHeight);
            for (key in view.subViews) {
                if (!view.subViews.hasOwnProperty(key)) {
                    continue;
                }
                const subView = view.subViews[key];
                if (elements[key]) {
                    UIView._setAbsoluteSizeAndPosition(elements[key], subView.left, subView.top, subView.width, subView.height);
                }
            }
            parentUIView.didLayoutSubviews();
        };
        updateLayout();
        return updateLayout;
    }
    static runFunctionBeforeNextFrame(step) {
        if (IS_SAFARI) {
            // This creates a microtask
            Promise.resolve().then(step);
        }
        else {
            window.requestAnimationFrame(step);
        }
    }
    static scheduleLayoutViewsIfNeeded() {
        UIView.runFunctionBeforeNextFrame(UIView.layoutViewsIfNeeded);
    }
    static layoutViewsIfNeeded() {
        for (var i = 0; i < UIView._viewsToLayout.length; i++) {
            const view = UIView._viewsToLayout[i];
            view.layoutIfNeeded();
        }
        UIView._viewsToLayout = [];
    }
    setNeedsLayout() {
        if (this._shouldLayout) {
            return;
        }
        this._shouldLayout = YES;
        // Register view for layout before next frame
        UIView._viewsToLayout.push(this);
        if (UIView._viewsToLayout.length == 1) {
            UIView.scheduleLayoutViewsIfNeeded();
        }
    }
    get needsLayout() {
        return this._shouldLayout;
    }
    layoutIfNeeded() {
        if (!this._shouldLayout) {
            return;
        }
        this._shouldLayout = NO;
        try {
            this.layoutSubviews();
        }
        catch (exception) {
            console.log(exception);
        }
    }
    layoutSubviews() {
        this._shouldLayout = NO;
        // Autolayout
        //window.removeEventListener('resize', this._updateLayoutFunction);
        if (this.constraints.length) {
            this._updateLayoutFunction = UIView.performAutoLayout(this.viewHTMLElement, null, this.constraints);
        }
        //this._updateLayoutFunction = this.layoutSubviews.bind(this);
        //window.addEventListener('resize', this._updateLayoutFunction);
        this._viewControllerLayoutFunction();
        this.applyClassesAndStyles();
        for (var i = 0; i < this.subviews.length; i++) {
            const subview = this.subviews[i];
            subview.calculateAndSetViewFrame();
            //subview.layoutSubviews();
        }
        this.didLayoutSubviews();
    }
    applyClassesAndStyles() {
        //var classesString = "";
        for (var i = 0; i < this.styleClasses.length; i++) {
            const styleClass = this.styleClasses[i];
            if (styleClass) {
                this.viewHTMLElement.classList.add(styleClass);
            }
            //classesString = classesString + " " + styleClass;
        }
        //this.viewHTMLElement.className = classesString;
    }
    didLayoutSubviews() {
        this._didLayoutSubviewsDelegateFunction();
    }
    get constraints() {
        return this._constraints;
    }
    set constraints(constraints) {
        this._constraints = constraints;
    }
    addConstraint(constraint) {
        this.constraints.push(constraint);
    }
    addConstraintsWithVisualFormat(visualFormatArray) {
        this.constraints = this.constraints.concat(AutoLayout.VisualFormat.parse(visualFormatArray, { extended: true }));
    }
    static constraintWithView(view, attribute, relation, toView, toAttribute, multiplier, constant, priority) {
        var UIViewObject = nil;
        var viewID = null;
        if (view) {
            if (view.isKindOfClass && view.isKindOfClass(UIView)) {
                UIViewObject = view;
                view = view.viewHTMLElement;
            }
            viewID = view.id;
        }
        var toUIViewObject = nil;
        var toViewID = null;
        if (toView) {
            if (toView.isKindOfClass && view.isKindOfClass(UIView)) {
                toUIViewObject = toView;
                toView = toView.viewHTMLElement;
            }
            toViewID = toView.id;
        }
        const constraint = {
            view1: viewID,
            attr1: attribute,
            relation: relation,
            view2: toViewID,
            attr2: toAttribute,
            multiplier: multiplier,
            constant: constant,
            priority: priority
        };
        return constraint;
    }
    subviewWithID(viewID) {
        var resultHTMLElement = nil;
        try {
            resultHTMLElement = this.viewHTMLElement.querySelector("#" + viewID);
        }
        catch (error) {
        }
        if (resultHTMLElement && resultHTMLElement.UIView) {
            return resultHTMLElement.UIView;
        }
        return nil;
    }
    rectangleContainingSubviews() {
        const center = this.bounds.center;
        var result = new UIRectangle(center.x, center.y, 0, 0);
        for (var i = 0; i < this.subviews.length; i++) {
            const subview = this.subviews[i];
            var frame = subview.frame;
            const rectangleContainingSubviews = subview.rectangleContainingSubviews();
            frame = frame.concatenateWithRectangle(rectangleContainingSubviews);
            result = result.concatenateWithRectangle(frame);
        }
        return result;
    }
    hasSubview(view) {
        // This is for performance reasons
        if (!IS(view)) {
            return NO;
        }
        for (var i = 0; i < this.subviews.length; i++) {
            const subview = this.subviews[i];
            if (subview == view) {
                return YES;
            }
        }
        return NO;
    }
    get viewBelowThisView() {
        const result = (this.viewHTMLElement.previousElementSibling || {}).UIView;
        return result;
    }
    get viewAboveThisView() {
        const result = (this.viewHTMLElement.nextElementSibling || {}).UIView;
        return result;
    }
    addSubview(view, aboveView) {
        if (!this.hasSubview(view) && IS(view)) {
            view.willMoveToSuperview(this);
            if (IS(aboveView)) {
                this.viewHTMLElement.insertBefore(view.viewHTMLElement, aboveView.viewHTMLElement.nextSibling);
                this.subviews.insertElementAtIndex(this.subviews.indexOf(aboveView), view);
            }
            else {
                this.viewHTMLElement.appendChild(view.viewHTMLElement);
                this.subviews.push(view);
            }
            view.didMoveToSuperview(this);
            if (this.superview && this.isMemberOfViewTree) {
                view.broadcastEventInSubtree({
                    name: UIView.broadcastEventName.AddedToViewTree,
                    parameters: nil
                });
            }
            this.setNeedsLayout();
        }
    }
    addSubviews(views) {
        views.forEach(function (view, index, array) {
            this.addSubview(view);
        }, this);
    }
    moveToBottomOfSuperview() {
        if (IS(this.superview)) {
            const bottomView = this.superview.subviews.firstElement;
            if (bottomView == this) {
                return;
            }
            this.superview.subviews.removeElement(this);
            this.superview.subviews.insertElementAtIndex(0, this);
            this.superview.viewHTMLElement.insertBefore(this.viewHTMLElement, bottomView.viewHTMLElement);
        }
    }
    moveToTopOfSuperview() {
        if (IS(this.superview)) {
            const topView = this.superview.subviews.lastElement;
            if (topView == this) {
                return;
            }
            this.superview.subviews.removeElement(this);
            this.superview.subviews.push(this);
            this.superview.viewHTMLElement.appendChild(this.viewHTMLElement);
        }
    }
    removeFromSuperview() {
        if (IS(this.superview)) {
            this.forEachViewInSubtree(function (view) {
                view.blur();
            });
            const index = this.superview.subviews.indexOf(this);
            if (index > -1) {
                this.superview.subviews.splice(index, 1);
                this.superview.viewHTMLElement.removeChild(this.viewHTMLElement);
                this.superview = nil;
                this.broadcastEventInSubtree({
                    name: UIView.broadcastEventName.RemovedFromViewTree,
                    parameters: nil
                });
            }
        }
    }
    willAppear() {
    }
    willMoveToSuperview(superview) {
        this._setInnerHTMLFromKeyIfPossible();
        this._setInnerHTMLFromLocalizedTextObjectIfPossible();
    }
    didMoveToSuperview(superview) {
        this.superview = superview;
    }
    wasAddedToViewTree() {
    }
    wasRemovedFromViewTree() {
    }
    get isMemberOfViewTree() {
        var element = this.viewHTMLElement;
        for (var i = 0; element; i = i) {
            if (element.parentElement && element.parentElement == document.body) {
                return YES;
            }
            element = element.parentElement;
        }
        return NO;
    }
    get allSuperviews() {
        const result = [];
        var view = this;
        for (var i = 0; IS(view); i = i) {
            result.push(view);
            view = view.superview;
        }
        return result;
    }
    setNeedsLayoutOnAllSuperviews() {
        this.allSuperviews.reverse().forEach(function (view, index, array) {
            view.setNeedsLayout();
        });
    }
    setNeedsLayoutUpToRootView() {
        this.setNeedsLayoutOnAllSuperviews();
        this.setNeedsLayout();
    }
    focus() {
        this.viewHTMLElement.focus();
    }
    blur() {
        this.viewHTMLElement.blur();
    }
    _loadUIEvents() {
        //this.viewHTMLElement = nil;
        const isTouchEventClassDefined = NO || window.TouchEvent;
        const pauseEvent = (event, forced = NO) => {
            if (this.pausesPointerEvents || forced) {
                if (event.stopPropagation) {
                    event.stopPropagation();
                }
                if (event.preventDefault) {
                    event.preventDefault();
                }
                event.cancelBubble = true;
                event.returnValue = false;
                return false;
            }
            if (event.stopPropagation && this.stopsPointerEventPropagation) {
                event.stopPropagation();
            }
        };
        const onMouseDown = (event) => {
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                ((this.ignoresMouse || (IS(this._touchEventTime) && (Date.now() - this._touchEventTime) > 500)) &&
                    event instanceof MouseEvent)) {
                return;
            }
            this.sendControlEventForKey(UIView.controlEvent.PointerDown, event);
            this._isPointerInside = YES;
            this._isPointerValid = YES;
            this._initialPointerPosition = new UIPoint(event.clientX, event.clientY);
            if (isTouchEventClassDefined && event instanceof TouchEvent) {
                this._touchEventTime = Date.now();
                this._initialPointerPosition = new UIPoint(event.touches[0].clientX, event.touches[0].clientY);
                if (event.touches.length > 1) {
                    onTouchCancel(event);
                    return;
                }
            }
            else {
                this._touchEventTime = nil;
                pauseEvent(event);
            }
            this._hasPointerDragged = NO;
        };
        const onTouchStart = onMouseDown;
        const onmouseup = (event) => {
            if (!this._isPointerValid) {
                return;
            }
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return;
            }
            if (this._isPointerInside) {
                onPointerUpInside(event);
                if (!this._hasPointerDragged) {
                    this.sendControlEventForKey(UIView.controlEvent.PointerTap, event);
                }
            }
            // This has to be sent after the more specific event so that UIButton can ignore it when not highlighted
            this.sendControlEventForKey(UIView.controlEvent.PointerUp, event);
            pauseEvent(event);
        };
        const onTouchEnd = onmouseup;
        // function onMouseEnter(event) {
        //     this.sendControlEventForKey(UIView.controlEvent.PointerEnter, event);
        //     this._isPointerInside = YES;
        //     pauseEvent(event);
        // }
        const onmouseout = (event) => {
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return;
            }
            this.sendControlEventForKey(UIView.controlEvent.PointerLeave, event);
            this._isPointerInside = NO;
            pauseEvent(event);
        };
        const onTouchLeave = onmouseout;
        var onTouchCancel = function (event) {
            if (!this._isPointerValid) {
                return;
            }
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return;
            }
            this._isPointerValid = NO;
            this.sendControlEventForKey(UIView.controlEvent.PointerCancel, event);
        }.bind(this);
        const onmouseover = (event) => {
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return;
            }
            this.sendControlEventForKey(UIView.controlEvent.PointerHover, event);
            this._isPointerInside = YES;
            this._isPointerValid = YES;
            pauseEvent(event);
        };
        const onMouseMove = (event) => {
            if (!this._isPointerValid) {
                return;
            }
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return;
            }
            if (IS_NOT(this._initialPointerPosition)) {
                this._initialPointerPosition = new UIPoint(event.clientX, event.clientY);
            }
            if (new UIPoint(event.clientX, event.clientY).to(this._initialPointerPosition).length >
                this._pointerDragThreshold) {
                this._hasPointerDragged = YES;
            }
            this.sendControlEventForKey(UIView.controlEvent.PointerMove, event);
            pauseEvent(event);
        };
        const onTouchMove = function (event) {
            if (!this._isPointerValid) {
                return;
            }
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return;
            }
            if (event.touches.length > 1) {
                onTouchZoom(event);
                return;
            }
            const touch = event.touches[0];
            if (new UIPoint(touch.clientX, touch.clientY).to(this._initialPointerPosition).length >
                this._pointerDragThreshold) {
                this._hasPointerDragged = YES;
            }
            if (this._isPointerInside && this.viewHTMLElement !=
                document.elementFromPoint(touch.clientX, touch.clientY)) {
                this._isPointerInside = NO;
                this.sendControlEventForKey(UIView.controlEvent.PointerLeave, event);
            }
            this.sendControlEventForKey(UIView.controlEvent.PointerMove, event);
            //pauseEvent(event);
        };
        var onTouchZoom = function onTouchZoom(event) {
            this.sendControlEventForKey(UIView.controlEvent.MultipleTouches, event);
        }.bind(this);
        var onPointerUpInside = (event) => {
            pauseEvent(event);
            this.sendControlEventForKey(UIView.controlEvent.PointerUpInside, event);
        };
        function eventKeyIsEnter(event) {
            if (event.keyCode !== 13) {
                return NO;
            }
            return YES;
        }
        function eventKeyIsTab(event) {
            if (event.keyCode !== 9) {
                return NO;
            }
            return YES;
        }
        function eventKeyIsEsc(event) {
            var result = false;
            if ("key" in event) {
                result = (event.key == "Escape" || event.key == "Esc");
            }
            else {
                result = (event.keyCode == 27);
            }
            return result;
        }
        function eventKeyIsLeft(event) {
            if (event.keyCode != "37") {
                return NO;
            }
            return YES;
        }
        function eventKeyIsRight(event) {
            if (event.keyCode != "39") {
                return NO;
            }
            return YES;
        }
        function eventKeyIsDown(event) {
            if (event.keyCode != "40") {
                return NO;
            }
            return YES;
        }
        function eventKeyIsUp(event) {
            if (event.keyCode != "38") {
                return NO;
            }
            return YES;
        }
        const onKeyDown = function (event) {
            if (eventKeyIsEnter(event)) {
                this.sendControlEventForKey(UIView.controlEvent.EnterDown, event);
            }
            if (eventKeyIsEsc(event)) {
                this.sendControlEventForKey(UIView.controlEvent.EscDown, event);
            }
            if (eventKeyIsTab(event) && this._controlEventTargets.TabDown && this._controlEventTargets.TabDown.length) {
                this.sendControlEventForKey(UIView.controlEvent.TabDown, event);
                pauseEvent(event, YES);
            }
            if (eventKeyIsLeft(event)) {
                this.sendControlEventForKey(UIView.controlEvent.LeftArrowDown, event);
            }
            if (eventKeyIsRight(event)) {
                this.sendControlEventForKey(UIView.controlEvent.RightArrowDown, event);
            }
            if (eventKeyIsDown(event)) {
                this.sendControlEventForKey(UIView.controlEvent.DownArrowDown, event);
            }
            if (eventKeyIsUp(event)) {
                this.sendControlEventForKey(UIView.controlEvent.UpArrowDown, event);
            }
        }.bind(this);
        const onKeyUp = function (event) {
            if (eventKeyIsEnter(event)) {
                this.sendControlEventForKey(UIView.controlEvent.EnterUp, event);
            }
        }.bind(this);
        const onfocus = function (event) {
            this.sendControlEventForKey(UIView.controlEvent.Focus, event);
        }.bind(this);
        const onblur = function (event) {
            this.sendControlEventForKey(UIView.controlEvent.Blur, event);
        }.bind(this);
        // Mouse and touch start events
        this._viewHTMLElement.onmousedown = onMouseDown.bind(this);
        this._viewHTMLElement.ontouchstart = onTouchStart.bind(this);
        // this.viewHTMLElement.addEventListener("mousedown", onMouseDown.bind(this), false)
        // this.viewHTMLElement.addEventListener('touchstart', onTouchStart.bind(this), false)
        // //this.viewHTMLElement.addEventListener("mouseenter", onMouseEnter.bind(this), false);
        // Mouse and touch move events
        this._viewHTMLElement.onmousemove = onMouseMove.bind(this);
        this._viewHTMLElement.ontouchmove = onTouchMove.bind(this);
        // this.viewHTMLElement.addEventListener("mousemove", onMouseMove.bind(this), false)
        // this.viewHTMLElement.addEventListener('touchmove', onTouchMove.bind(this), false)
        //this.viewHTMLElement.addEventListener("mousewheel", onmousewheel.bind(this), false)
        this._viewHTMLElement.onmouseover = onmouseover.bind(this);
        // this.viewHTMLElement.addEventListener("mouseover", onmouseover.bind(this), false)
        // Mouse and touch end events
        this._viewHTMLElement.onmouseup = onmouseup.bind(this);
        this._viewHTMLElement.ontouchend = onTouchEnd.bind(this);
        this._viewHTMLElement.ontouchcancel = onTouchCancel.bind(this);
        // this.viewHTMLElement.addEventListener("mouseup", onmouseup.bind(this), false)
        // this.viewHTMLElement.addEventListener('touchend', onTouchEnd.bind(this), false)
        // this.viewHTMLElement.addEventListener('touchcancel', onTouchCancel.bind(this), false)
        this._viewHTMLElement.onmouseout = onmouseout.bind(this);
        // this.viewHTMLElement.addEventListener("mouseout", onmouseout.bind(this), false)
        this._viewHTMLElement.addEventListener("touchleave", onTouchLeave.bind(this), false);
        // this.viewHTMLElement.onkeydown = onkeydown
        // this.viewHTMLElement.onkeyup = onkeyup
        this._viewHTMLElement.addEventListener("keydown", onKeyDown, false);
        this._viewHTMLElement.addEventListener("keyup", onKeyUp, false);
        // Focus events
        this._viewHTMLElement.onfocus = onfocus;
        this._viewHTMLElement.onblur = onblur;
        // this.viewHTMLElement.addEventListener("focus", onfocus, true)
        // this.viewHTMLElement.addEventListener("blur", onblur, true)
    }
    get addControlEventTarget() {
        const eventKeys = [];
        const result = new Proxy(this.constructor.controlEvent, {
            get: (target, key, receiver) => {
                eventKeys.push(key);
                return result;
            },
            set: (target, key, value, receiver) => {
                eventKeys.push(key);
                this.addTargetForControlEvents(eventKeys, value);
                return true;
            }
        });
        return result;
    }
    addTargetForControlEvents(eventKeys, targetFunction) {
        eventKeys.forEach(function (key, index, array) {
            this.addTargetForControlEvent(key, targetFunction);
        }, this);
    }
    addTargetForControlEvent(eventKey, targetFunction) {
        var targets = this._controlEventTargets[eventKey];
        if (!targets) {
            // @ts-ignore
            targets = [];
            this._controlEventTargets[eventKey] = targets;
        }
        if (targets.indexOf(targetFunction) == -1) {
            targets.push(targetFunction);
        }
    }
    removeTargetForControlEvent(eventKey, targetFunction) {
        const targets = this._controlEventTargets[eventKey];
        if (!targets) {
            return;
        }
        const index = targets.indexOf(targetFunction);
        if (index != -1) {
            targets.splice(index, 1);
        }
    }
    removeTargetForControlEvents(eventKeys, targetFunction) {
        eventKeys.forEach(function (key, index, array) {
            this.removeTargetForControlEvent(key, targetFunction);
        }, this);
    }
    sendControlEventForKey(eventKey, nativeEvent) {
        var targets = this._controlEventTargets[eventKey];
        if (!targets) {
            return;
        }
        targets = targets.copy();
        for (var i = 0; i < targets.length; i++) {
            const target = targets[i];
            target(this, nativeEvent);
        }
    }
    broadcastEventInSubtree(event) {
        this.forEachViewInSubtree(function (view) {
            view.didReceiveBroadcastEvent(event);
            if (view._didReceiveBroadcastEventDelegateFunction) {
                view._didReceiveBroadcastEventDelegateFunction(event);
            }
        });
    }
    didReceiveBroadcastEvent(event) {
        if (event.name == UIView.broadcastEventName.PageDidScroll) {
            this._isPointerValid = NO;
        }
        if (event.name == UIView.broadcastEventName.AddedToViewTree) {
            this.wasAddedToViewTree();
        }
        if (event.name == UIView.broadcastEventName.RemovedFromViewTree) {
            this.wasRemovedFromViewTree();
        }
        if (event.name == UIView.broadcastEventName.LanguageChanged || event.name ==
            UIView.broadcastEventName.AddedToViewTree) {
            this._setInnerHTMLFromKeyIfPossible();
            this._setInnerHTMLFromLocalizedTextObjectIfPossible();
        }
    }
    forEachViewInSubtree(functionToCall) {
        functionToCall(this);
        this.subviews.forEach(function (subview, index, array) {
            subview.forEachViewInSubtree(functionToCall);
        });
    }
    rectangleInView(rectangle, view) {
        if (!view.isMemberOfViewTree || !this.isMemberOfViewTree) {
            return nil;
        }
        const viewClientRectangle = view.viewHTMLElement.getBoundingClientRect();
        const viewLocation = new UIPoint(viewClientRectangle.left, viewClientRectangle.top);
        const selfClientRectangle = this.viewHTMLElement.getBoundingClientRect();
        const selfLocation = new UIPoint(selfClientRectangle.left, selfClientRectangle.top);
        const offsetPoint = selfLocation.subtract(viewLocation);
        return rectangle.copy().offsetByPoint(offsetPoint);
    }
    rectangleFromView(rectangle, view) {
        return view.rectangleInView(rectangle, this);
    }
    intrinsicContentSizeWithConstraints(constrainingHeight = 0, constrainingWidth = 0) {
        // This works but is slow
        const result = new UIRectangle(0, 0, 0, 0);
        if (this.rootView.forceIntrinsicSizeZero) {
            return result;
        }
        var temporarilyInViewTree = NO;
        var nodeAboveThisView;
        if (!this.isMemberOfViewTree) {
            document.body.appendChild(this.viewHTMLElement);
            temporarilyInViewTree = YES;
            nodeAboveThisView = this.viewHTMLElement.nextSibling;
        }
        const height = this.style.height;
        const width = this.style.width;
        this.style.height = "" + constrainingHeight;
        this.style.width = "" + constrainingWidth;
        const left = this.style.left;
        const right = this.style.right;
        const bottom = this.style.bottom;
        const top = this.style.top;
        this.style.left = "";
        this.style.right = "";
        this.style.bottom = "";
        this.style.top = "";
        const resultHeight = this.viewHTMLElement.scrollHeight;
        const whiteSpace = this.style.whiteSpace;
        this.style.whiteSpace = "nowrap";
        const resultWidth = this.viewHTMLElement.scrollWidth;
        this.style.whiteSpace = whiteSpace;
        this.style.height = height;
        this.style.width = width;
        this.style.left = left;
        this.style.right = right;
        this.style.bottom = bottom;
        this.style.top = top;
        if (temporarilyInViewTree) {
            document.body.removeChild(this.viewHTMLElement);
            if (this.superview) {
                if (nodeAboveThisView) {
                    this.superview.viewHTMLElement.insertBefore(this.viewHTMLElement, nodeAboveThisView);
                }
                else {
                    this.superview.viewHTMLElement.appendChild(this.viewHTMLElement);
                }
            }
        }
        result.height = resultHeight;
        result.width = resultWidth;
        return result;
    }
    intrinsicContentWidth(constrainingHeight = 0) {
        const result = this.intrinsicContentSizeWithConstraints(constrainingHeight).width;
        return result;
    }
    intrinsicContentHeight(constrainingWidth = 0) {
        const result = this.intrinsicContentSizeWithConstraints(undefined, constrainingWidth).height;
        return result;
    }
    intrinsicContentSize() {
        return nil;
    }
}
UIView._UIViewIndex = -1;
UIView._viewsToLayout = [];
UIView._pageScale = 1;
UIView._transformAttribute = (("transform" in document.documentElement.style) ? "transform" : undefined) ||
    (("-webkit-transform" in document.documentElement.style) ? "-webkit-transform" : "undefined") ||
    (("-moz-transform" in document.documentElement.style) ? "-moz-transform" : "undefined") ||
    (("-ms-transform" in document.documentElement.style) ? "-ms-transform" : "undefined") ||
    (("-o-transform" in document.documentElement.style) ? "-o-transform" : "undefined");
UIView.constraintAttribute = {
    "left": AutoLayout.Attribute.LEFT,
    "right": AutoLayout.Attribute.RIGHT,
    "bottom": AutoLayout.Attribute.BOTTOM,
    "top": AutoLayout.Attribute.TOP,
    "centerX": AutoLayout.Attribute.CENTERX,
    "centerY": AutoLayout.Attribute.CENTERY,
    "height": AutoLayout.Attribute.HEIGHT,
    "width": AutoLayout.Attribute.WIDTH,
    "zIndex": AutoLayout.Attribute.ZINDEX,
    // Not sure what these are for
    "constant": AutoLayout.Attribute.NOTANATTRIBUTE,
    "variable": AutoLayout.Attribute.VARIABLE
};
UIView.constraintRelation = {
    "equal": AutoLayout.Relation.EQU,
    "lessThanOrEqual": AutoLayout.Relation.LEQ,
    "greaterThanOrEqual": AutoLayout.Relation.GEQ
};
UIView.controlEvent = {
    "PointerDown": "PointerDown",
    "PointerMove": "PointerMove",
    "PointerLeave": "PointerLeave",
    "PointerEnter": "PointerEnter",
    "PointerUpInside": "PointerUpInside",
    "PointerTap": "PointerTap",
    "PointerUp": "PointerUp",
    "MultipleTouches": "PointerZoom",
    "PointerCancel": "PointerCancel",
    "PointerHover": "PointerHover",
    "EnterDown": "EnterDown",
    "EnterUp": "EnterUp",
    "EscDown": "EscDown",
    "TabDown": "TabDown",
    "LeftArrowDown": "LeftArrowDown",
    "RightArrowDown": "RightArrowDown",
    "DownArrowDown": "DownArrowDown",
    "UpArrowDown": "UpArrowDown",
    "Focus": "Focus",
    "Blur": "Blur"
};
UIView.broadcastEventName = {
    "LanguageChanged": "LanguageChanged",
    "RemovedFromViewTree": "RemovedFromViewTree",
    "AddedToViewTree": "AddedToViewTree",
    "PageDidScroll": "PageDidScroll"
};
/// <reference path="./UICore/UIView.ts" />
class BottomBarView extends UIView {
    //label: UITextView;
    constructor(elementID) {
        super(elementID);
        this._class = BottomBarView;
        this.superclass = UIView;
    }
    initView(elementID, viewHTMLElement) {
        super.initView(elementID, viewHTMLElement);
        this.backgroundColor = UIColor.colorWithRGBA(50, 50, 50);
        this.setInnerHTML("bottomBarContent", "bottomBarContent");
        // this.label = new UITextView("BottomBarLabel");
        // this.label.setText("bottomBarLabelText","Bottom bar for contact information and such");
        // this.label.textColor = UIColor.whiteColor;
        // this.label.textAlignment = UITextView.textAlignment.center;
        // this.addSubview(this.label);
    }
    layoutSubviews() {
        super.layoutSubviews();
        //this.label.centerInContainer();
    }
}
if ("removeElementAtIndex" in Array.prototype == NO) {
    Array.prototype.removeElementAtIndex = function (index) {
        if (index >= 0 && index < this.length) {
            this.splice(index, 1);
        }
    };
}
if ("removeElement" in Array.prototype == NO) {
    Array.prototype.removeElement = function (element) {
        this.removeElementAtIndex(this.indexOf(element));
    };
}
if ("insertElementAtIndex" in Array.prototype == NO) {
    Array.prototype.insertElementAtIndex = function (index, element) {
        if (index >= 0 && index <= this.length) {
            this.splice(index, 0, element);
        }
    };
}
if ("replaceElementAtIndex" in Array.prototype == NO) {
    Array.prototype.replaceElementAtIndex = function (index, element) {
        this.removeElementAtIndex(index);
        this.insertElementAtIndex(index, element);
    };
}
if ("contains" in Array.prototype == NO) {
    Array.prototype.contains = function (element) {
        const result = (this.indexOf(element) != -1);
        return result;
    };
}
if ("containsAny" in Array.prototype == NO) {
    Array.prototype.containsAny = function (elements) {
        const result = this.anyMatch(function (element, index, array) {
            return elements.contains(element);
        });
        return result;
    };
}
Array.prototype.max = function () {
    return Math.max.apply(null, this);
};
Array.prototype.min = function () {
    return Math.min.apply(null, this);
};
if ("anyMatch" in Array.prototype == NO) {
    Array.prototype.anyMatch = function (functionToCall) {
        const result = (this.findIndex(functionToCall) > -1);
        return result;
    };
}
if ("noneMatch" in Array.prototype == NO) {
    Array.prototype.noneMatch = function (functionToCall) {
        const result = (this.findIndex(functionToCall) == -1);
        return result;
    };
}
if ("allMatch" in Array.prototype == NO) {
    Array.prototype.allMatch = function (functionToCall) {
        function reversedFunction(value, index, array) {
            return !functionToCall(value, index, array);
        }
        const result = (this.findIndex(reversedFunction) == -1);
        return result;
    };
}
if ("groupedBy" in Array.prototype == NO) {
    Array.prototype.groupedBy = function (funcProp) {
        return this.reduce(function (acc, val) {
            (acc[funcProp(val)] = acc[funcProp(val)] || []).push(val);
            return acc;
        }, {});
    };
}
if ("firstElement" in Array.prototype == NO) {
    Object.defineProperty(Array.prototype, "firstElement", {
        get: function firstElement() {
            const result = this[0];
            return result;
        },
        set: function (element) {
            if (this.length == 0) {
                this.push(element);
                return;
            }
            this[0] = element;
        }
    });
}
if ("lastElement" in Array.prototype == NO) {
    Object.defineProperty(Array.prototype, "lastElement", {
        get: function lastElement() {
            const result = this[this.length - 1];
            return result;
        },
        set: function (element) {
            if (this.length == 0) {
                this.push(element);
                return;
            }
            this[this.length - 1] = element;
        }
    });
}
if ("everyElement" in Array.prototype == NO) {
    Object.defineProperty(Array.prototype, "everyElement", {
        get: function everyElement() {
            var valueKeys = [];
            const targetFunction = (objects) => {
                return this.map((element, index, array) => {
                    var elementFunction = UIObject.valueForKeyPath(valueKeys.join("."), element).bind(element, objects);
                    return elementFunction();
                });
            };
            const result = new Proxy(targetFunction, {
                get: (target, key, receiver) => {
                    if (key == "UI_elementValues") {
                        return this.map((element, index, array) => UIObject.valueForKeyPath(valueKeys.join("."), element));
                    }
                    valueKeys.push(key);
                    return result;
                },
                set: (target, key, value, receiver) => {
                    valueKeys.push(key);
                    this.forEach((element, index, array) => {
                        UIObject.setValueForKeyPath(valueKeys.join("."), value, element, YES);
                    });
                    return true;
                }
            });
            return result;
        },
        set: function (element) {
            for (var i = 0; i < this.length; ++i) {
                this[i] = element;
            }
        }
    });
}
if ("copy" in Array.prototype == NO) {
    Array.prototype.copy = function () {
        const result = this.slice(0);
        return result;
    };
}
if ("arrayByRepeating" in Array.prototype == NO) {
    Array.prototype.arrayByRepeating = function (numberOfRepetitions) {
        const result = [];
        for (var i = 0; i < numberOfRepetitions; i++) {
            this.forEach(function (element, index, array) {
                result.push(element);
            });
        }
        return result;
    };
}
if ("arrayByTrimmingToLengthIfLonger" in Array.prototype == NO) {
    Array.prototype.arrayByTrimmingToLengthIfLonger = function (maxLength) {
        const result = [];
        for (var i = 0; i < maxLength && i < this.length; i++) {
            result.push(this[i]);
        }
        return result;
    };
}
if ("summedValue" in Array.prototype == NO) {
    Object.defineProperty(Array.prototype, "summedValue", {
        get: function summedValue() {
            const result = this.reduce(function (a, b) {
                return a + b;
            }, 0);
            return result;
        }
    });
}
// Warn if overriding existing method
if ("isEqualToArray" in Array.prototype == YES) {
    console.warn("Overriding existing Array.prototype.isEqualToArray. Possible causes: New API defines the method, there's a framework conflict or you've got double inclusions in your code.");
}
// attach the .equals method to Array's prototype to call it on any array
Array.prototype.isEqualToArray = function (array, keyPath) {
    // if the other array is a falsy value, return
    if (!array) {
        return false;
    }
    // compare lengths - can save a lot of time 
    if (this.length != array.length) {
        return false;
    }
    var i = 0;
    const l = this.length;
    for (; i < l; i++) {
        // Check if we have nested arrays
        if (this[i] instanceof Array && array[i] instanceof Array && !keyPath) {
            // recurse into the nested arrays
            if (!this[i].isEqualToArray(array[i])) {
                return false;
            }
        }
        else if (keyPath && UIObject.valueForKeyPath(keyPath, this[i]) != UIObject.valueForKeyPath(keyPath, array[i])) {
            return false;
        }
        else if (this[i] != array[i]) {
            // Warning - two different object instances will never be equal: {x:20} != {x:20}
            return false;
        }
    }
    return true;
};
// Hide method from for-in loops
Object.defineProperty(Array.prototype, "isEqualToArray", { enumerable: false });
if ("forEach" in Object.prototype == NO) {
    Object.prototype.forEach = function (callbackFunction) {
        const keys = Object.keys(this);
        keys.forEach(function (key, index, array) {
            callbackFunction(this[key], key);
        }.bind(this));
    };
    // Hide method from for-in loops
    Object.defineProperty(Object.prototype, "forEach", { enumerable: false });
}
if ("allValues" in Object.prototype == NO) {
    Object.defineProperty(Object.prototype, "allValues", {
        get: function () {
            const values = [];
            this.forEach(function (value) {
                values.push(value);
            });
            return values;
        }
    });
}
if ("allKeys" in Object.prototype == NO) {
    Object.defineProperty(Object.prototype, "allKeys", {
        get: function () {
            const values = Object.keys(this);
            return values;
        }
    });
}
if ("objectByCopyingValuesRecursivelyFromObject" in Object.prototype == NO) {
    Object.prototype.objectByCopyingValuesRecursivelyFromObject = function (object) {
        function isAnObject(item) {
            return (item && typeof item === "object" && !Array.isArray(item));
        }
        function mergeRecursively(target, source) {
            const output = Object.assign({}, target);
            if (isAnObject(target) && isAnObject(source)) {
                Object.keys(source).forEach(function (key) {
                    if (isAnObject(source[key])) {
                        // if (!(key in target)) {
                        //     Object.assign(output, { [key]: source[key] });
                        // }
                        // else {
                        output[key] = mergeRecursively(target[key], source[key]);
                        //}
                    }
                    else {
                        Object.assign(output, { [key]: source[key] });
                    }
                });
            }
            return output;
        }
        const result = mergeRecursively(this, object);
        return result;
    };
    // Hide method from for-in loops
    Object.defineProperty(Object.prototype, "objectByCopyingValuesRecursivelyFromObject", { enumerable: false });
}
if ("copyValuesRecursivelyFromObject" in Object.prototype == NO) {
    Object.prototype.copyValuesRecursivelyFromObject = function (object) {
        function isAnObject(item) {
            return (item && typeof item === "object" && !Array.isArray(item));
        }
        function mergeRecursively(target, source) {
            const output = target;
            if (isAnObject(target) && isAnObject(source)) {
                Object.keys(source).forEach(function (key) {
                    if (isAnObject(source[key])) {
                        // if (!(key in target)) {
                        //     Object.assign(output, { [key]: source[key] });
                        // }
                        // else {
                        output[key] = mergeRecursively(target[key], source[key]);
                        //}
                    }
                    else {
                        Object.assign(output, { [key]: source[key] });
                    }
                });
            }
            return output;
        }
        const result = mergeRecursively(this, object);
        return result;
    };
    // Hide method from for-in loops
    Object.defineProperty(Object.prototype, "copyValuesRecursivelyFromObject", { enumerable: false });
}
if ("contains" in String.prototype == NO) {
    String.prototype.contains = function (string) {
        const result = (this.indexOf(string) != -1);
        return result;
    };
}
if ("capitalizedString" in String.prototype == NO) {
    Object.defineProperty(Object.prototype, "capitalizedString", {
        get: function () {
            const result = this.charAt(0).toUpperCase() + this.slice(1).toLowerCase();
            return result;
        }
    });
}
if ("numericalValue" in String.prototype == NO) {
    Object.defineProperty(String.prototype, "numericalValue", {
        get: function numericalValue() {
            const result = Number(this);
            return result;
        }
    });
}
if ("isAString" in String.prototype == NO) {
    String.prototype.isAString = YES;
}
if ("isANumber" in Number.prototype == NO) {
    Number.prototype.isANumber = YES;
}
if ("integerValue" in Number.prototype == NO) {
    Object.defineProperty(Number.prototype, "integerValue", {
        get: function () {
            const result = parseInt("" + (Math.round(this) + 0.5));
            return result;
        }
    });
}
class PrimitiveNumber {
    // @ts-ignore
    static [Symbol.hasInstance](x) {
        return;
    }
}
if ("integerValue" in Boolean.prototype == NO) {
    Object.defineProperty(Boolean.prototype, "integerValue", {
        get: function () {
            if (this == true) {
                return 1;
            }
            return 0;
        }
    });
}
if ("dateString" in Date.prototype == NO) {
    Object.defineProperty(Date.prototype, "dateString", {
        get: function dateString() {
            const result = ("0" + this.getDate()).slice(-2) + "-" + ("0" + (this.getMonth() + 1)).slice(-2) + "-" +
                this.getFullYear() + " " + ("0" + this.getHours()).slice(-2) + ":" +
                ("0" + this.getMinutes()).slice(-2);
            return result;
        }
    });
}
/// <reference path="../UICore/UIObject.ts" />
/// <reference path="../UICore/UICoreExtensions.ts" />
class CBSocketCallbackHolder extends UIObject {
    constructor(socketClient, previousCallbackHolder) {
        super();
        this.messageDescriptors = {};
        this.handlers = {};
        this.onetimeHandlers = {};
        this.keysForIdentifiers = {};
        this.isValid = YES;
        this._verifiedResponseHashesDictionary = {};
        this._socketClient = socketClient;
        if (IS(previousCallbackHolder)) {
            this.handlers = previousCallbackHolder.handlers;
            this._verifiedResponseHashesDictionary = previousCallbackHolder._verifiedResponseHashesDictionary;
        }
    }
    triggerDisconnectHandlers() {
        this.messageDescriptors.forEach(function (descriptor, key) {
            if (descriptor.mainResponseReceived) {
                descriptor.completionFunction(CBSocketClient.disconnectionMessage, nil);
            }
        });
    }
    registerHandler(key, handlerFunction) {
        if (!this.handlers[key]) {
            this.handlers[key] = [];
        }
        this.handlers[key].push(handlerFunction);
    }
    registerOnetimeHandler(key, handlerFunction) {
        if (!this.onetimeHandlers[key]) {
            this.onetimeHandlers[key] = [];
        }
        this.onetimeHandlers[key].push(handlerFunction);
    }
    get storedResponseHashesDictionary() {
        if (IS_NOT(this._storedResponseHashesDictionary)) {
            this._storedResponseHashesDictionary = JSON.parse(localStorage["CBSocketResponseHashesDictionary"] || "{}");
        }
        return this._storedResponseHashesDictionary;
    }
    storedResponseHashObjectForKey(requestKey, requestDataHash) {
        const localStorageKey = this.keyForRequestKeyAndRequestDataHash(requestKey, requestDataHash);
        const hashObject = this.storedResponseHashesDictionary[localStorageKey];
        const result = FIRST(hashObject, {});
        return result;
    }
    storedResponseForKey(requestKey, requestDataHash) {
        const localStorageKey = this.keyForRequestKeyAndRequestDataHash(requestKey, requestDataHash);
        const storedObject = JSON.parse(localStorage[localStorageKey] || "{}");
        return storedObject.responseMessageData;
    }
    keyForRequestKeyAndRequestDataHash(requestKey, requestDataHash) {
        const result = "_CBSCH_LS_key_" + requestKey + "_" + requestDataHash;
        return result;
    }
    storeResponse(requestKey, requestDataHash, responseMessage, responseDataHash) {
        if (!responseMessage.canBeStoredAsResponse ||
            (IS_NOT(responseMessage.messageData) && IS_NOT(responseMessage.messageDataHash))) {
            return;
        }
        const localStorageKey = this.keyForRequestKeyAndRequestDataHash(requestKey, requestDataHash);
        var validityDate;
        if (responseMessage.responseValidityDuration) {
            validityDate = Date.now() + responseMessage.responseValidityDuration;
        }
        const storedResponseHashesDictionary = this.storedResponseHashesDictionary;
        storedResponseHashesDictionary[localStorageKey] = {
            hash: responseDataHash,
            validityDate: validityDate
        };
        this.saveInLocalStorage(localStorageKey, {
            responseMessageData: responseMessage.messageData,
            responseHash: responseDataHash
        });
        this.saveStoredResponseHashesDictionary(storedResponseHashesDictionary);
    }
    saveStoredResponseHashesDictionary(storedResponseHashesDictionary) {
        this.saveInLocalStorage("CBSocketResponseHashesDictionary", storedResponseHashesDictionary);
    }
    saveInLocalStorage(key, object) {
        const stringToSave = JSON.stringify(object);
        if (stringToSave != localStorage[key]) {
            localStorage[key] = stringToSave;
        }
    }
    socketShouldSendMessage(key, message, completionPolicy, completionFunction) {
        var result = YES;
        var triggerStoredResponseImmediately = NO;
        // @ts-ignore
        const messageDataHash = objectHash(message.messageData || nil);
        const descriptorKey = "socketMessageDescriptor_" + key + messageDataHash;
        this.messageDescriptors[descriptorKey] = (this.messageDescriptors[descriptorKey] || []);
        const hashObject = this.storedResponseHashObjectForKey(key, messageDataHash);
        message.storedResponseHash = hashObject.hash;
        if (completionPolicy == CBSocketClient.completionPolicy.first) {
            const descriptorsForKey = (this.messageDescriptors[descriptorKey] || []);
            const matchingDescriptor = descriptorsForKey.find(function (descriptor, index, array) {
                return (descriptor.messageDataHash == messageDataHash);
            });
            if (matchingDescriptor) {
                result = NO;
            }
        }
        if (completionPolicy == CBSocketClient.completionPolicy.storedOrFirst) {
            const descriptorsForKey = (this.messageDescriptors[descriptorKey] || []);
            const matchingDescriptor = descriptorsForKey.find(function (descriptor, index, array) {
                return (descriptor.messageDataHash == messageDataHash);
            });
            const storedResponse = IS(message.storedResponseHash);
            if (matchingDescriptor ||
                (storedResponse && this._verifiedResponseHashesDictionary[message.storedResponseHash])) {
                result = NO;
                triggerStoredResponseImmediately = YES;
            }
        }
        if (completionPolicy == CBSocketClient.completionPolicy.firstOnly) {
            const descriptorsForKey = (this.messageDescriptors[descriptorKey] || []);
            const matchingDescriptor = descriptorsForKey.find(function (descriptor, index, array) {
                return (descriptor.messageDataHash == messageDataHash);
            });
            if (matchingDescriptor) {
                return NO;
            }
        }
        if (hashObject && hashObject.hash && hashObject.validityDate && message.storedResponseHash &&
            this._verifiedResponseHashesDictionary[message.storedResponseHash] && hashObject.validityDate >
            Date.now()) {
            result = NO;
            triggerStoredResponseImmediately = YES;
        }
        if (IS(completionFunction)) {
            this.messageDescriptors[descriptorKey].push({
                key: key,
                message: {
                    identifier: message.identifier,
                    inResponseToIdentifier: message.inResponseToIdentifier,
                    keepWaitingForResponses: message.keepWaitingForResponses
                },
                sentAtTime: Date.now(),
                //completionTriggered: NO,
                messageDataHash: messageDataHash,
                mainResponseReceived: NO,
                anyMainResponseReceived: NO,
                completionPolicy: completionPolicy,
                completionFunction: completionFunction
            });
            this.keysForIdentifiers[message.identifier] = descriptorKey;
        }
        if (triggerStoredResponseImmediately) {
            this.socketDidReceiveMessageForKey(CBSocketClient.responseMessageKey, {
                identifier: nil,
                messageData: nil,
                inResponseToIdentifier: message.identifier,
                useStoredResponse: YES
            }, nil);
        }
        return result;
    }
    static defaultMultipleMessagecompletionFunction(responseMessages, callcompletionFunctions) {
        callcompletionFunctions();
    }
    socketWillSendMultipleMessage(messageToSend, completionFunction = CBSocketCallbackHolder.defaultMultipleMessagecompletionFunction) {
        const key = CBSocketClient.multipleMessageKey;
        // @ts-ignore
        const messageDataHash = objectHash(messageToSend.messageData || nil);
        const descriptorKey = "socketMessageDescriptor_" + key + messageDataHash;
        this.messageDescriptors[descriptorKey] = (this.messageDescriptors[descriptorKey] || []);
        messageToSend.storedResponseHash = this.storedResponseHashObjectForKey(key, messageDataHash).hash;
        this.messageDescriptors[descriptorKey].push({
            key: key,
            message: {
                identifier: messageToSend.identifier,
                inResponseToIdentifier: messageToSend.inResponseToIdentifier,
                keepWaitingForResponses: messageToSend.keepWaitingForResponses
            },
            sentAtTime: Date.now(),
            //completionTriggered: NO,
            messageDataHash: messageDataHash,
            mainResponseReceived: NO,
            anyMainResponseReceived: NO,
            completionPolicy: CBSocketClient.completionPolicy.directOnly,
            completionFunction: function (responseMessage, respondWithMessage) {
                completionFunction(responseMessage.map(function (messageObject, index, array) {
                    return messageObject.message.messageData;
                }), function () {
                    //console.log("Received multiple message response with length of " + responseMessage.length + ".");
                    // Call all completion functions
                    responseMessage.forEach(function (messageObject, index, array) {
                        this._socketClient.didReceiveMessageForKey(messageObject.key, messageObject.message);
                    }.bind(this));
                }.bind(this));
            }.bind(this)
        });
        this.keysForIdentifiers[messageToSend.identifier] = descriptorKey;
    }
    socketDidReceiveMessageForKey(key, message, sendResponseFunction) {
        if (!this.isValid) {
            return;
        }
        // Call static handlers
        if (this.handlers[key]) {
            this.handlers[key].forEach(function (handler, index, array) {
                handler(message.messageData, sendResponseFunction);
            }.bind(this));
        }
        if (this.onetimeHandlers[key]) {
            this.onetimeHandlers[key].forEach(function (handler, index, array) {
                handler(message.messageData, sendResponseFunction);
            }.bind(this));
            delete this.onetimeHandlers[key];
        }
        // Temporary response handlers are evaluated here
        if (message.inResponseToIdentifier &&
            (CBSocketClient.responseMessageKey == key || CBSocketClient.multipleMessageKey == key)) {
            // Find descriptors for the key of the message that is being responded to
            const descriptorKey = this.keysForIdentifiers[message.inResponseToIdentifier];
            const descriptorsForKey = (this.messageDescriptors[descriptorKey] || []);
            // Find response data hash to check for differences
            const responseDataHash = message.messageDataHash;
            // Remove identifier from dictionary
            if (!message.keepWaitingForResponses) {
                delete this.keysForIdentifiers[message.inResponseToIdentifier];
                delete this.messageDescriptors[descriptorKey];
            }
            // Function to call completion function
            const callCompletionFunction = (descriptor, storedResponseCondition = NO) => {
                var messageData = message.messageData;
                if (message.useStoredResponse && storedResponseCondition) {
                    messageData = this.storedResponseForKey(descriptor.key, descriptor.messageDataHash);
                    const responseHash = this.storedResponseHashObjectForKey(descriptor.key, descriptor.messageDataHash).hash;
                    const localStorageKey = this.keyForRequestKeyAndRequestDataHash(descriptor.key, descriptor.messageDataHash);
                    if (message.responseValidityDuration && this.storedResponseHashesDictionary[localStorageKey]) {
                        this.storedResponseHashesDictionary[localStorageKey].validityDate = Date.now() +
                            message.responseValidityDuration;
                        this.saveStoredResponseHashesDictionary(this.storedResponseHashesDictionary);
                    }
                    this._verifiedResponseHashesDictionary[responseHash] = YES;
                    console.log("Using stored response.");
                }
                // Call completionFunction and set response data hash
                descriptor.completionFunction(messageData, sendResponseFunction);
                descriptor.responseDataHash = responseDataHash;
            };
            descriptorsForKey.copy().forEach(function (descriptor, index, array) {
                if ((descriptor.completionPolicy == CBSocketClient.completionPolicy.directOnly &&
                    descriptor.message.identifier == message.inResponseToIdentifier) || descriptor.completionPolicy ==
                    CBSocketClient.completionPolicy.first || descriptor.completionPolicy ==
                    CBSocketClient.completionPolicy.firstOnly || descriptor.completionPolicy ==
                    CBSocketClient.completionPolicy.storedOrFirst) {
                    // Calling completion function and removing descriptor
                    if (!message.keepWaitingForResponses) {
                        this.storeResponse(descriptor.key, descriptor.messageDataHash, message, responseDataHash);
                        descriptorsForKey.removeElement(descriptor);
                        sendResponseFunction.respondingToMainResponse = YES;
                    }
                    callCompletionFunction(descriptor, !message.keepWaitingForResponses);
                }
                else if (descriptor.completionPolicy == CBSocketClient.completionPolicy.all) {
                    // Calling completion function
                    callCompletionFunction(descriptor, !message.keepWaitingForResponses);
                    // Marking descriptor as having been responded to
                    if (!message.keepWaitingForResponses) {
                        if (message.inResponseToIdentifier == descriptor.message.identifier) {
                            sendResponseFunction.respondingToMainResponse = YES;
                            descriptor.mainResponseReceived = YES;
                            descriptorsForKey.removeElement(descriptor);
                        }
                        descriptor.anyMainResponseReceived = YES;
                    }
                }
                else if (descriptor.completionPolicy == CBSocketClient.completionPolicy.allDifferent) {
                    // Calling completionFunction if messageData is different from previous
                    if (descriptor.responseDataHash != responseDataHash) {
                        callCompletionFunction(descriptor, !message.keepWaitingForResponses);
                    }
                    // Marking descriptor as having been responded to
                    if (!message.keepWaitingForResponses) {
                        if (message.inResponseToIdentifier == descriptor.message.identifier) {
                            sendResponseFunction.respondingToMainResponse = YES;
                            descriptor.mainResponseReceived = YES;
                            descriptorsForKey.removeElement(descriptor);
                        }
                        descriptor.anyMainResponseReceived = YES;
                    }
                }
                else if (descriptor.completionPolicy == CBSocketClient.completionPolicy.last &&
                    descriptor.message.identifier == message.inResponseToIdentifier) {
                    if (!message.keepWaitingForResponses) {
                        // Marking descriptor as having been responded to
                        descriptor.mainResponseReceived = YES;
                        descriptor.anyMainResponseReceived = YES;
                        sendResponseFunction.respondingToMainResponse = YES;
                    }
                    else {
                        descriptor.completionFunction(message.messageData, sendResponseFunction);
                    }
                }
                else if (descriptor.completionPolicy == CBSocketClient.completionPolicy.firstAndLast ||
                    descriptor.completionPolicy == CBSocketClient.completionPolicy.firstAndLastIfDifferent) {
                    if (!message.keepWaitingForResponses) {
                        // Only calling completionFunction once as a first response call
                        if (!descriptor.anyMainResponseReceived) {
                            callCompletionFunction(descriptor, !message.keepWaitingForResponses);
                        }
                        // Marking descriptor as having been responded to
                        if (descriptor.message.identifier == message.inResponseToIdentifier) {
                            descriptor.mainResponseReceived = YES;
                            sendResponseFunction.respondingToMainResponse = YES;
                        }
                        descriptor.anyMainResponseReceived = YES;
                    }
                    else if (descriptor.message.identifier == message.inResponseToIdentifier &&
                        message.keepWaitingForResponses) {
                        descriptor.completionFunction(message.messageData, sendResponseFunction);
                    }
                }
            }.bind(this));
            // Last message completion policies
            const allResponsesReceived = descriptorsForKey.allMatch(function (descriptorObject, index, array) {
                return descriptorObject.mainResponseReceived;
            });
            descriptorsForKey.copy().forEach(function (descriptor, index, array) {
                if ((descriptor.completionPolicy == CBSocketClient.completionPolicy.last ||
                    descriptor.completionPolicy == CBSocketClient.completionPolicy.firstAndLast) &&
                    allResponsesReceived && !message.keepWaitingForResponses) {
                    // Calling completionFunction
                    callCompletionFunction(descriptor, !message.keepWaitingForResponses);
                    // Cleaning up
                    descriptorsForKey.removeElement(descriptor);
                }
                else if (descriptor.completionPolicy == CBSocketClient.completionPolicy.firstAndLastIfDifferent &&
                    allResponsesReceived && !message.keepWaitingForResponses) {
                    // Calling completionFunction if needed
                    if (descriptor.responseDataHash != responseDataHash) {
                        callCompletionFunction(descriptor, !message.keepWaitingForResponses);
                    }
                    // Cleaning up
                    descriptorsForKey.removeElement(descriptor);
                }
            }.bind(this));
        }
    }
}
/// <reference path="../socket.io-client.dev.d.ts" />
/// <reference path="../UICore/UIObject.ts" />
/// <reference path="../UICore/UICoreExtensions.ts" />
/// <reference path="../CBCore.ts" />
/// <reference path="./CBSocketCallbackHolder.ts" />
function IS_SOCKET_ERROR(object) {
    const result = (IS(object) && object._isCBSocketErrorMessage);
    return result;
}
function IS_NOT_SOCKET_ERROR(object) {
    return !IS_SOCKET_ERROR(object);
}
class CBSocketClient extends UIObject {
    constructor(core) {
        super();
        this._socket = io();
        this._isConnectionEstablished = NO;
        this._collectMessagesToSendLater = NO;
        this._messagesToBeSent = [];
        this._subscribedKeys = {};
        this._callbackHolder = new CBSocketCallbackHolder(this);
        this._class = CBSocketClient;
        this.superclass = UIObject;
        this._core = core;
        // @ts-ignore
        this.socket.on("connect", function (socket) {
            console.log("Socket.io connected to server. clientID = " + socket + ", socketID = " + socket);
            var instanceIdentifier = localStorage.getItem("InstanceIdentifier");
            if (IS_NOT(instanceIdentifier)) {
                instanceIdentifier = MAKE_ID();
                localStorage.setItem("InstanceIdentifier", instanceIdentifier);
            }
            const handshakeMessage = {
                accessToken: null,
                instanceIdentifier: instanceIdentifier
            };
            this.socket.emit("CBSocketHandshakeInitMessage", {
                identifier: MAKE_ID(),
                messageData: handshakeMessage
            });
        }.bind(this));
        this.socket.on("CBSocketHandshakeResponseMessage", function (message) {
            this._isConnectionEstablished = message.messageData.accepted;
            if (!message.messageData.accepted) {
                console.log("SocketIO connection failed.");
            }
            else {
                console.log("SocketIO connection handshake completed.");
                this._callbackHolder = new CBSocketCallbackHolder(this, this._callbackHolder);
                this.sendUnsentMessages();
            }
        }.bind(this));
        // @ts-ignore
        this.socket.on("disconnect", function (socket) {
            console.log("Socket.io disconnected from server. clientID = " + socket + ".");
            this._isConnectionEstablished = NO;
            this._callbackHolder.isValid = NO;
            this._callbackHolder.triggerDisconnectHandlers();
        }.bind(this));
        this.socket.on("CBPerformReconnect", function (message) {
            console.log("Performing socket reconnection.");
            core.reloadSocketConnection();
            if (message) {
                alert(message);
            }
        }.bind(this));
        this._socket.on(CBSocketClient.responseMessageKey, function (message) {
            this.didReceiveMessageForKey(CBSocketClient.responseMessageKey, message);
        }.bind(this));
        this._socket.on(CBSocketClient.multipleMessageKey, (message) => {
            console.log("Received " + message.messageData.length + " messages.");
            this.didReceiveMessageForKey(CBSocketClient.multipleMessageKey, message);
        });
    }
    get socket() {
        return this._socket;
    }
    cancelUnsentMessages(messagesToCancel) {
        this._messagesToBeSent = this._messagesToBeSent.filter(function (messageObject, index, array) {
            return !messagesToCancel.contains(messageObject);
        }.bind(this));
    }
    sendUnsentMessages(receiveResponsesTogether = NO, completion) {
        if (!this._isConnectionEstablished || this._collectMessagesToSendLater) {
            return;
        }
        const groupedMessages = [];
        const didSendFunctions = [];
        this._messagesToBeSent.copy().forEach(function (messageToBeSentObject, index, array) {
            if (this._isConnectionEstablished) {
                var message = messageToBeSentObject.message;
                if (IS_NOT(message)) {
                    message = "";
                }
                const identifier = MAKE_ID();
                const completion = messageToBeSentObject.completion;
                const messageObject = {
                    messageData: message,
                    identifier: identifier,
                    keepWaitingForResponses: messageToBeSentObject.keepWaitingForResponses,
                    inResponseToIdentifier: messageToBeSentObject.inResponseToMessage.identifier
                };
                const shouldSendMessage = this._callbackHolder.socketShouldSendMessage(messageToBeSentObject.key, messageObject, messageToBeSentObject.completionPolicy, completion);
                if (shouldSendMessage) {
                    groupedMessages.push({
                        key: messageToBeSentObject.key,
                        message: messageObject
                    });
                }
                didSendFunctions.push(messageToBeSentObject.didSendFunction);
            }
        }.bind(this));
        this._messagesToBeSent = [];
        if (IS_NOT(groupedMessages.length)) {
            return;
        }
        if (groupedMessages.length == 1) {
            console.log("sending 1 unsent message.");
        }
        else {
            console.log("Sending " + groupedMessages.length + " unsent messages.");
        }
        const messageObject = {
            messageData: groupedMessages,
            identifier: MAKE_ID(),
            shouldGroupResponses: receiveResponsesTogether
        };
        //if (receiveResponsesTogether) {
        this._callbackHolder.socketWillSendMultipleMessage(messageObject, completion);
        //}
        this.socket.emit(CBSocketClient.multipleMessageKey, messageObject);
        didSendFunctions.forEach(function (didSendFunction, index, array) {
            didSendFunction();
        });
    }
    sendUserBoundMessageForKeyWithPolicy(key, message, completionPolicy, completion) {
        this._sendMessageForKey(key, message, undefined, NO, completionPolicy, YES, nil, completion);
    }
    sendUserBoundMessageForKey(key, message, completion) {
        this._sendMessageForKey(key, message, undefined, NO, undefined, YES, nil, completion);
    }
    sendMessageForKeyWithPolicy(key, message, completionPolicy, completion) {
        this._sendMessageForKey(key, message, undefined, NO, completionPolicy, NO, nil, completion);
    }
    sendMessageForKey(key, message, completion) {
        this._sendMessageForKey(key, message, undefined, NO, undefined, NO, nil, completion);
    }
    resultForMessageForKey(key, message, completionPolicy, isUserBound = NO, intermediateResultMessageFunction = nil) {
        const result = new Promise((resolve, reject) => {
            this._sendMessageForKey(key, message, undefined, NO, completionPolicy, isUserBound, nil, (responseMessage, respondWithMessage) => {
                if (respondWithMessage.respondingToMainResponse) {
                    resolve({
                        responseMessage: responseMessage,
                        result: IF(IS_NOT_SOCKET_ERROR(responseMessage))(() => responseMessage)
                            .ELSE(RETURNER(undefined)),
                        errorResult: IF(IS_SOCKET_ERROR(responseMessage))(() => responseMessage)
                            .ELSE(RETURNER(undefined)),
                        respondWithMessage: respondWithMessage
                    });
                    return;
                }
                intermediateResultMessageFunction(responseMessage, respondWithMessage);
            });
        });
        return result;
    }
    _sendMessageForKey(key, message, inResponseToMessage = {}, keepMessageConnectionOpen = NO, completionPolicy = CBSocketClient.completionPolicy.directOnly, isUserBound = NO, didSendFunction = nil, completion = nil) {
        if (IS_NIL(message)) {
            message = "";
        }
        if (this._isConnectionEstablished && !this._collectMessagesToSendLater) {
            const identifier = MAKE_ID();
            const messageObject = {
                messageData: message,
                identifier: identifier,
                keepWaitingForResponses: keepMessageConnectionOpen,
                inResponseToIdentifier: inResponseToMessage.identifier
            };
            const shouldSendMessage = this._callbackHolder.socketShouldSendMessage(key, messageObject, completionPolicy, completion);
            if (shouldSendMessage) {
                this.socket.emit(key, messageObject);
            }
            didSendFunction();
        }
        else {
            this._messagesToBeSent.push({
                key: key,
                message: message,
                inResponseToMessage: inResponseToMessage,
                keepWaitingForResponses: keepMessageConnectionOpen,
                completionPolicy: completionPolicy,
                isBoundToUserWithID: null,
                didSendFunction: didSendFunction,
                completion: completion
            });
            return this._messagesToBeSent.lastElement;
        }
    }
    sendMessagesAsGroup(functionToCall) {
        const collectMessagesToSendLater = this._collectMessagesToSendLater;
        this._collectMessagesToSendLater = YES;
        functionToCall();
        this._collectMessagesToSendLater = collectMessagesToSendLater;
        this.sendUnsentMessages();
    }
    sendAndReceiveMessagesAsGroup(functionToCall, completion) {
        const collectMessagesToSendLater = this._collectMessagesToSendLater;
        this._collectMessagesToSendLater = YES;
        functionToCall();
        this._collectMessagesToSendLater = collectMessagesToSendLater;
        this.sendUnsentMessages(YES, completion);
    }
    didReceiveMessageForKey(key, message) {
        const sendResponseFunction = function (responseMessage, completion) {
            this._sendMessageForKey(CBSocketClient.responseMessageKey, responseMessage, message, NO, undefined, NO, nil, completion);
        }.bind(this);
        sendResponseFunction.sendIntermediateResponse = function (updateMessage, completion) {
            this._sendMessageForKey(CBSocketClient.responseMessageKey, updateMessage, message, YES, undefined, NO, nil, completion);
        }.bind(this);
        const sendUserBoundResponseFunction = function (responseMessage, completion) {
            this._sendMessageForKey(CBSocketClient.responseMessageKey, responseMessage, message, NO, undefined, YES, nil, completion);
        }.bind(this);
        sendUserBoundResponseFunction.sendIntermediateResponse = function (updateMessage, completion) {
            this._sendMessageForKey(CBSocketClient.responseMessageKey, updateMessage, message, YES, undefined, YES, nil, completion);
        }.bind(this);
        if (IS_SOCKET_ERROR(message.messageData)) {
            console.log("CBSocketClient did receive error message.");
            console.log(message.messageData);
        }
        this._callbackHolder.socketDidReceiveMessageForKey(key, message, sendResponseFunction);
    }
    addTargetForMessagesForKeys(keys, handlerFunction) {
        keys.forEach(function (key, index, array) {
            this.addTargetForMessagesForKey(key, handlerFunction);
        }.bind(this));
    }
    addTargetForMessagesForKey(key, handlerFunction) {
        this._callbackHolder.registerHandler(key, handlerFunction);
        if (IS_NOT(this._subscribedKeys[key])) {
            this._socket.on(key, function (message) {
                this.didReceiveMessageForKey(key, message);
            }.bind(this));
            this._subscribedKeys[key] = true;
        }
    }
    addTargetForOneMessageForKey(key, handlerFunction) {
        this._callbackHolder.registerOnetimeHandler(key, handlerFunction);
        if (IS_NOT(this._subscribedKeys[key])) {
            this._socket.on(key, function (message) {
                this.didReceiveMessageForKey(key, message);
            }.bind(this));
            this._subscribedKeys[key] = true;
        }
    }
}
CBSocketClient.responseMessageKey = "CBSocketResponseMessage";
CBSocketClient.multipleMessageKey = "CBSocketMultipleMessage";
CBSocketClient.disconnectionMessage = {
    _isCBSocketErrorMessage: YES,
    messageData: "Server disconnected"
};
CBSocketClient.completionPolicy = {
    "all": "all",
    "allDifferent": "allDifferent",
    "first": "first",
    "last": "last",
    "firstAndLast": "firstAndLast",
    "firstAndLastIfDifferent": "firstAndLastIfDifferent",
    "directOnly": "directOnly",
    "firstOnly": "firstOnly",
    "storedOrFirst": "storedOrFirst"
};
const SocketClient = new Proxy({ "name": "SocketClient" }, {
    get(target, key) {
        const result = (messageData, completionPolicy, isUserBound, intermediateResultFunction) => CBCore.sharedInstance.socketClient.resultForMessageForKey(key, messageData, completionPolicy, isUserBound, intermediateResultFunction);
        return result;
    }
});
/// <reference path="./Custom components/CBSocketClient.ts" />
class CBCore extends UIObject {
    constructor() {
        super();
        this._socketClient = new CBSocketClient(this);
        this._functionsToCallForEachSocketClient = [];
        this._class = CBCore;
        this.superclass = UIObject;
    }
    static get sharedInstance() {
        if (!CBCore._sharedInstance) {
            CBCore._sharedInstance = new CBCore();
        }
        return CBCore._sharedInstance;
    }
    broadcastMessageInRootViewTree(message) {
        UICore.main.rootViewController.view.broadcastEventInSubtree(message);
    }
    get socketClient() {
        return this._socketClient;
    }
    set languageKey(languageKey) {
        if (IS_NOT(languageKey)) {
            localStorage.removeItem("CBLanguageKey");
        }
        localStorage.setItem("CBLanguageKey", JSON.stringify(languageKey));
        this.didSetLanguageKey();
    }
    get languageKey() {
        const result = FIRST(localStorage.getItem("CBLanguageKey"), LanguageService.defaultLanguageKey).replace("\"", "").replace("\"", "");
        return result;
    }
    didSetLanguageKey() {
        UIRoute.currentRoute.routeWithComponent("settings", { "language": this.languageKey }, YES).applyByReplacingCurrentRouteInHistory();
    }
    readTextFileFromURL(fileURL) {
        var result = new Promise((resolve, reject) => {
            var rawFile = new XMLHttpRequest();
            rawFile.open("GET", fileURL, false);
            rawFile.onreadystatechange = function () {
                if (rawFile.readyState === 4) {
                    if (rawFile.status === 200 || rawFile.status == 0) {
                        var allText = rawFile.responseText;
                        resolve(allText);
                        return;
                    }
                    reject(rawFile.response);
                }
            };
            rawFile.send(null);
        });
        return result;
    }
    reloadSocketConnection() {
        // @ts-ignore
        this.socketClient.socket.disconnect();
        this._socketClient = new CBSocketClient(this);
        this._socketClient._messagesToBeSent = [];
        this._functionsToCallForEachSocketClient.forEach(function (functionToCall, index, array) {
            functionToCall();
        });
    }
    callFunctionForEachSocketClient(functionToCall) {
        this._functionsToCallForEachSocketClient.push(functionToCall);
        functionToCall();
    }
}
/// <reference path="./UIView.ts" />
/// <reference path="./UIObject.ts" />
class UIViewController extends UIObject {
    constructor(view) {
        super();
        this.view = view;
        this.loadIntrospectionVariables();
        this._UIViewController_constructorArguments = { "view": view };
        this._initInstanceVariables();
        this.loadSubviews();
        this.updateViewConstraints();
        this.updateViewStyles();
        this._layoutViewSubviews();
    }
    loadIntrospectionVariables() {
        this._class = UIViewController;
        this.superclass = UIObject;
    }
    _initInstanceVariables() {
        this.view = this._UIViewController_constructorArguments.view;
        this.view._viewControllerLayoutFunction = this.layoutViewsManually.bind(this);
        this.view._didLayoutSubviewsDelegateFunction = this.viewDidLayoutSubviews.bind(this);
        this.view._didReceiveBroadcastEventDelegateFunction = this.viewDidReceiveBroadcastEvent.bind(this);
        this.childViewControllers = [];
        this.parentViewController = nil;
    }
    handleRouteRecursively(route) {
        this.handleRoute(route);
        this.childViewControllers.forEach(function (controller, index, array) {
            if (!route.isHandled) {
                controller.handleRouteRecursively(route);
            }
        });
    }
    handleRoute(route) {
        return __awaiter(this, void 0, void 0, function* () {
        });
    }
    loadSubviews() {
    }
    viewWillAppear() {
        return __awaiter(this, void 0, void 0, function* () {
        });
    }
    viewDidAppear() {
        return __awaiter(this, void 0, void 0, function* () {
        });
    }
    viewWillDisappear() {
        return __awaiter(this, void 0, void 0, function* () {
        });
    }
    viewDidDisappear() {
        return __awaiter(this, void 0, void 0, function* () {
        });
    }
    updateViewConstraints() {
    }
    updateViewStyles() {
    }
    layoutViewsManually() {
    }
    _layoutViewSubviews() {
        this.view.layoutSubviews();
        this.viewDidLayoutSubviews();
    }
    viewDidLayoutSubviews() {
        // this.childViewControllers.forEach(function (controller, index, controllers) {
        //     controller._layoutViewSubviews();
        // })
    }
    viewDidReceiveBroadcastEvent(event) {
    }
    hasChildViewController(viewController) {
        // This is for performance reasons
        if (!IS(viewController)) {
            return NO;
        }
        for (var i = 0; i < this.childViewControllers.length; i++) {
            const childViewController = this.childViewControllers[i];
            if (childViewController == viewController) {
                return YES;
            }
        }
        return NO;
    }
    addChildViewController(viewController) {
        if (!this.hasChildViewController(viewController)) {
            viewController.willMoveToParentViewController(this);
            this.childViewControllers.push(viewController);
            //this.view.addSubview(viewController.view);
            //viewController.didMoveToParentViewController(this);
        }
    }
    removeFromParentViewController() {
        const index = this.parentViewController.childViewControllers.indexOf(this);
        if (index > -1) {
            this.parentViewController.childViewControllers.splice(index, 1);
            //this.view.removeFromSuperview();
            this.parentViewController = nil;
        }
    }
    willMoveToParentViewController(parentViewController) {
    }
    didMoveToParentViewController(parentViewController) {
        this.parentViewController = parentViewController;
    }
    removeChildViewController(controller) {
        controller.viewWillDisappear();
        if (IS(controller.parentViewController)) {
            controller.removeFromParentViewController();
        }
        if (IS(controller.view)) {
            controller.view.removeFromSuperview();
        }
        controller.viewDidDisappear();
    }
    addChildViewControllerInContainer(controller, containerView) {
        controller.viewWillAppear();
        this.addChildViewController(controller);
        containerView.addSubview(controller.view);
        controller.didMoveToParentViewController(this);
        controller.viewDidAppear();
    }
    addChildViewControllerInDialogView(controller, dialogView) {
        controller.viewWillAppear();
        this.addChildViewController(controller);
        dialogView.view = controller.view;
        var originalDismissFunction = dialogView.dismiss.bind(dialogView);
        dialogView.dismiss = animated => {
            originalDismissFunction(animated);
            this.removeChildViewController(controller);
        };
        controller.didMoveToParentViewController(this);
        controller.viewDidAppear();
    }
}
///<reference path="UIView.ts"/>
class UIDialogView extends UIView {
    constructor(elementID, viewHTMLElement) {
        super(elementID, viewHTMLElement);
        this._view = nil;
        this.animationDuration = 0.25;
        this._zIndex = 100;
        this.isVisible = NO;
        this.dismissesOnTapOutside = YES;
        this.addTargetForControlEvent(UIView.controlEvent.PointerTap, function (sender, event) {
            this.didDetectTapOutside(sender, event);
        }.bind(this));
        this.backgroundColor = UIColor.colorWithRGBA(0, 10, 25).colorWithAlpha(0.75); //CBColor.primaryContentColor.colorWithAlpha(0.75)
        this.zIndex = this._zIndex;
    }
    didDetectTapOutside(sender, event) {
        if (event.target == this.viewHTMLElement && this.dismissesOnTapOutside) {
            this.dismiss(this._appearedAnimated);
        }
    }
    set zIndex(zIndex) {
        this._zIndex = zIndex;
        this.style.zIndex = "" + zIndex;
    }
    get zIndex() {
        return this._zIndex;
    }
    set view(view) {
        this._view.removeFromSuperview();
        this._view = view;
        this.addSubview(view);
    }
    get view() {
        return this._view;
    }
    willAppear(animated = NO) {
        if (animated) {
            this.style.opacity = "0";
        }
        this.style.height = "";
        this._frame = null;
    }
    animateAppearing() {
        this.style.opacity = "1";
    }
    animateDisappearing() {
        this.style.opacity = "0";
    }
    showInView(containerView, animated) {
        animated = (animated && !IS_FIREFOX);
        this._appearedAnimated = animated;
        this.willAppear(animated);
        containerView.addSubview(this);
        if (animated) {
            this.layoutSubviews();
            UIView.animateViewOrViewsWithDurationDelayAndFunction(this, this.animationDuration, 0, undefined, function () {
                this.animateAppearing();
            }.bind(this), nil);
        }
        else {
            this.setNeedsLayout();
        }
        this.isVisible = YES;
    }
    showInRootView(animated) {
        this.showInView(UICore.main.rootViewController.view, animated);
    }
    dismiss(animated) {
        animated = (animated && !IS_FIREFOX);
        if (animated == undefined) {
            animated = this._appearedAnimated;
        }
        if (animated) {
            UIView.animateViewOrViewsWithDurationDelayAndFunction(this, this.animationDuration, 0, undefined, function () {
                this.animateDisappearing();
            }.bind(this), function () {
                if (this.isVisible == NO) {
                    this.removeFromSuperview();
                }
            }.bind(this));
        }
        else {
            this.removeFromSuperview();
        }
        this.isVisible = NO;
    }
    didReceiveBroadcastEvent(event) {
        super.didReceiveBroadcastEvent(event);
        if (event.name == UICore.broadcastEventName.WindowDidResize) {
            this.setNeedsLayout();
        }
    }
    layoutSubviews() {
        if (!IS(this.view)) {
            return;
        }
        //this.frame = this.superview.bounds;
        this.setPosition(0, 0, 0, 0, 0, "100%");
        this.setPosition(0, 0, 0, 0, UIView.pageHeight, "100%");
        const bounds = this.bounds;
        const margin = 20;
        //this.view.centerInContainer();
        this.view.style.position = "relative";
        // this.view.style.maxHeight = "" + (bounds.height - margin * 2).integerValue + "px";
        // this.view.style.maxWidth = "" + (bounds.width - margin * 2).integerValue + "px";
        // var viewIntrinsicRectangle = this.view.intrinsicContentSize();
        // this.view.frame = new UIRectangle((bounds.width - viewIntrinsicRectangle.width)*0.5,  )
        super.layoutSubviews();
    }
}
/// <reference path="./UIView.ts" />
class UIBaseButton extends UIView {
    constructor(elementID, elementType, initViewData) {
        super(elementID, nil, elementType, initViewData);
        this._selected = NO;
        this._highlighted = NO;
        this._isToggleable = NO;
        this._class = UIButton;
        this.superclass = UIView;
        this.initViewStateControl();
    }
    initViewStateControl() {
        this.class.superclass = UIView;
        // Instance variables
        this._isPointerInside = NO;
        const setHovered = function () {
            this.hovered = YES;
        }.bind(this);
        this.addTargetForControlEvent(UIView.controlEvent.PointerHover, setHovered);
        const setNotHovered = function () {
            this.hovered = NO;
        }.bind(this);
        this.addTargetForControlEvents([
            UIView.controlEvent.PointerLeave, UIView.controlEvent.PointerCancel, UIView.controlEvent.MultipleTouches
        ], setNotHovered);
        var highlightingTime;
        const setHighlighted = function () {
            this.highlighted = YES;
            highlightingTime = Date.now();
        }.bind(this);
        this.addTargetForControlEvent(UIView.controlEvent.PointerDown, setHighlighted);
        this.addTargetForControlEvent(UIView.controlEvent.PointerEnter, setHighlighted);
        const setNotHighlighted = function () {
            this.highlighted = NO;
        }.bind(this);
        const setNotHighlightedWithMinimumDuration = function () {
            const minimumDurationInMilliseconds = 50;
            const elapsedTime = Date.now() - highlightingTime;
            if (minimumDurationInMilliseconds < elapsedTime) {
                this.highlighted = NO;
            }
            else {
                setTimeout(function () {
                    this.highlighted = NO;
                }.bind(this), minimumDurationInMilliseconds - elapsedTime);
            }
        }.bind(this);
        this.addTargetForControlEvents([
            UIView.controlEvent.PointerLeave, UIView.controlEvent.PointerCancel, UIView.controlEvent.MultipleTouches
        ], setNotHighlighted);
        this.addTargetForControlEvent(UIView.controlEvent.PointerUp, setNotHighlightedWithMinimumDuration);
        // Handle enter key press
        this.addTargetForControlEvent(UIView.controlEvent.EnterDown, function () {
            setHighlighted();
            setNotHighlightedWithMinimumDuration();
        });
        this.addTargetForControlEvent(UIView.controlEvent.Focus, function (sender, event) {
            this.focused = YES;
        }.bind(this));
        this.addTargetForControlEvent(UIView.controlEvent.Blur, function (sender, event) {
            this.focused = NO;
        }.bind(this));
        this.updateContentForCurrentState();
        this.pausesPointerEvents = YES;
        this.tabIndex = 1;
        this.style.cursor = "pointer";
        //this.style.outline = "none";
        this.nativeSelectionEnabled = NO;
        this.addTargetForControlEvents([
            UIView.controlEvent.EnterDown, UIView.controlEvent.PointerUpInside
        ], function (sender, event) {
            if (this.isToggleable) {
                this.toggleSelectedState();
            }
        }.bind(this));
    }
    set hovered(hovered) {
        this._hovered = hovered;
        this.updateContentForCurrentState();
    }
    get hovered() {
        return this._hovered;
    }
    set highlighted(highlighted) {
        this._highlighted = highlighted;
        this.updateContentForCurrentState();
    }
    get highlighted() {
        return this._highlighted;
    }
    set focused(focused) {
        this._focused = focused;
        if (focused) {
            this.focus();
        }
        else {
            this.blur();
        }
        this.updateContentForCurrentState();
    }
    get focused() {
        return this._focused;
    }
    set selected(selected) {
        this._selected = selected;
        this.updateContentForCurrentState();
    }
    get selected() {
        return this._selected;
    }
    updateContentForCurrentState() {
        var updateFunction = this.updateContentForNormalState;
        if (this.selected && this.highlighted) {
            updateFunction = this.updateContentForSelectedAndHighlightedState;
        }
        else if (this.selected) {
            updateFunction = this.updateContentForSelectedState;
        }
        else if (this.focused) {
            updateFunction = this.updateContentForFocusedState;
        }
        else if (this.highlighted) {
            updateFunction = this.updateContentForHighlightedState;
        }
        else if (this.hovered) {
            updateFunction = this.updateContentForHoveredState;
        }
        if (!IS(updateFunction)) {
            this.backgroundColor = UIColor.nilColor;
        }
        else {
            updateFunction.call(this);
        }
    }
    updateContentForNormalState() {
    }
    updateContentForHoveredState() {
        this.updateContentForNormalState();
    }
    updateContentForFocusedState() {
        this.updateContentForHoveredState();
    }
    updateContentForHighlightedState() {
    }
    updateContentForSelectedState() {
    }
    updateContentForSelectedAndHighlightedState() {
        this.updateContentForSelectedState();
    }
    set enabled(enabled) {
        super.enabled = enabled;
        this.updateContentForCurrentEnabledState();
    }
    get enabled() {
        return super.enabled;
    }
    updateContentForCurrentEnabledState() {
        if (this.enabled) {
            this.alpha = 1;
        }
        else {
            this.alpha = 0.5;
        }
        this.userInteractionEnabled = this.enabled;
    }
    addStyleClass(styleClassName) {
        super.addStyleClass(styleClassName);
        if (this.styleClassName != styleClassName) {
            this.updateContentForCurrentState.call(this);
        }
    }
    didReceiveBroadcastEvent(event) {
        super.didReceiveBroadcastEvent(event);
        if (event.name == UIView.broadcastEventName.PageDidScroll || event.name ==
            UIView.broadcastEventName.AddedToViewTree) {
            this.hovered = NO;
            this.highlighted = NO;
        }
    }
    toggleSelectedState() {
        this.selected = !this.selected;
    }
    set isToggleable(isToggleable) {
        this._isToggleable = isToggleable;
    }
    get isToggleable() {
        return this._isToggleable;
    }
    layoutSubviews() {
        super.layoutSubviews();
        const bounds = this.bounds;
    }
    sendControlEventForKey(eventKey, nativeEvent) {
        if (eventKey == UIView.controlEvent.PointerUpInside && !this.highlighted) {
            // Do not send the event in this case
            //super.sendControlEventForKey(eventKey, nativeEvent);
            const asd = 1;
        }
        else {
            super.sendControlEventForKey(eventKey, nativeEvent);
        }
    }
    static getEventCoordinatesInDocument(touchOrMouseEvent) {
        // http://www.quirksmode.org/js/events_properties.html
        var posx = 0;
        var posy = 0;
        var e = touchOrMouseEvent;
        if (!e) {
            e = window.event;
        }
        if (e.pageX || e.pageY) {
            posx = e.pageX;
            posy = e.pageY;
        }
        else if (e.clientX || e.clientY) {
            posx = e.clientX + document.body.scrollLeft
                + document.documentElement.scrollLeft;
            posy = e.clientY + document.body.scrollTop
                + document.documentElement.scrollTop;
        }
        // posx and posy contain the mouse position relative to the document
        const coordinates = { "x": posx, "y": posy };
        return coordinates;
    }
    static getElementPositionInDocument(el) {
        //https://www.kirupa.com/html5/getting_mouse_click_position.htm
        var xPosition = 0;
        var yPosition = 0;
        while (el) {
            if (el.tagName == "BODY") {
                // Coordinates in document are coordinates in body, therefore subtracting the scroll position of the body is not needed
                //      // deal with browser quirks with body/window/document and page scroll
                //      var xScrollPos = el.scrollLeft || document.documentElement.scrollLeft;
                //      var yScrollPos = el.scrollTop || document.documentElement.scrollTop;
                //
                //      xPosition += (el.offsetLeft - xScrollPos + el.clientLeft);
                //      yPosition += (el.offsetTop - yScrollPos + el.clientTop);
            }
            else {
                xPosition += (el.offsetLeft - el.scrollLeft + el.clientLeft);
                yPosition += (el.offsetTop - el.scrollTop + el.clientTop);
            }
            el = el.offsetParent;
        }
        return {
            x: xPosition,
            y: yPosition
        };
    }
    static convertCoordinatesFromDocumentToElement(x, y, element) {
        const elementPositionInDocument = this.getElementPositionInDocument(element);
        const coordinatesInElement = { "x": x - elementPositionInDocument.x, "y": y - elementPositionInDocument.y };
        return coordinatesInElement;
    }
    static getEventCoordinatesInElement(touchOrMouseEvent, element) {
        const coordinatesInDocument = this.getEventCoordinatesInDocument(touchOrMouseEvent);
        const coordinatesInElement = this.convertCoordinatesFromDocumentToElement(coordinatesInDocument.x, coordinatesInDocument.y, element);
        return coordinatesInElement;
    }
}
/// <reference path="./UIBaseButton.ts" />
class UILink extends UIBaseButton {
    constructor(elementID, initViewData = nil) {
        super(elementID, "a", initViewData);
        this._class = UILink;
        this.superclass = UIBaseButton;
        this.stopsPointerEventPropagation = NO;
        this.pausesPointerEvents = NO;
    }
    initView(elementID, viewHTMLElement, initViewData) {
        super.initView(elementID, viewHTMLElement, initViewData);
        this.class.superclass = UIBaseButton;
        // Instance variables
        //this.style.position = "relative"
        viewHTMLElement.onclick = this.blur.bind(this);
    }
    get colors() {
        return this._colors;
    }
    set colors(value) {
        this._colors = value;
    }
    get viewHTMLElement() {
        return super.viewHTMLElement;
    }
    set text(text) {
        this.viewHTMLElement.textContent = text;
    }
    get text() {
        return this.viewHTMLElement.textContent;
    }
    set target(target) {
        this.viewHTMLElement.setAttribute("href", target);
    }
    get target() {
        const result = this.viewHTMLElement.getAttribute("href");
        return result;
    }
    set targetRouteForCurrentState(targetRouteForCurrentState) {
        this._targetRouteForCurrentState = targetRouteForCurrentState;
        this.updateTarget();
    }
    get targetRouteForCurrentState() {
        return this._targetRouteForCurrentState;
    }
    _targetRouteForCurrentState() {
        const result = UIRoute.currentRoute.routeByRemovingComponentsOtherThanOnesNamed(["settings"]);
        return result;
    }
    didReceiveBroadcastEvent(event) {
        super.didReceiveBroadcastEvent(event);
        if (event.name == UICore.broadcastEventName.RouteDidChange) {
            this.updateTarget();
        }
    }
    wasAddedToViewTree() {
        super.wasAddedToViewTree();
        this.updateTarget();
    }
    updateTarget() {
        const route = this.targetRouteForCurrentState();
        if (route instanceof UIRoute) {
            this.target = route.linkRepresentation;
            return;
        }
        this.target = route;
    }
    layoutSubviews() {
        super.layoutSubviews();
        const bounds = this.bounds;
    }
}
/// <reference path="./UILink.ts" />
class UILinkButton extends UILink {
    //link = new UILink(this.elementID + "Link");
    constructor(elementID, elementType, titleType) {
        super(elementID, { "elementType": elementType, "titleType": titleType });
        this._class = UILinkButton;
        this.superclass = UILink;
        this.button.addTargetForControlEvents([
            UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
        ], function (sender, event) {
            const asd = 1;
            window.location = this.target;
        }.bind(this));
        // this.link.hidden = YES;
        // this.addSubview(this.link);
    }
    initView(elementID, viewHTMLElement, initViewData) {
        super.initView(elementID, viewHTMLElement, initViewData);
        this.class.superclass = UILink;
        // Instance variables
        this.button = new UIButton(this.elementID + "Button", initViewData.elementType, initViewData.titleType);
        this.addSubview(this.button);
        this.style.position = "absolute";
    }
    get titleLabel() {
        return this.button.titleLabel;
    }
    get imageView() {
        return this.button.imageView;
    }
    set colors(colors) {
        this.button.colors = colors;
    }
    get colors() {
        return this.button.colors;
    }
    get viewHTMLElement() {
        return super.viewHTMLElement;
    }
    set target(target) {
        this.viewHTMLElement.setAttribute("href", target);
    }
    get target() {
        const result = this.viewHTMLElement.getAttribute("href");
        return result;
    }
    layoutSubviews() {
        super.layoutSubviews();
        const bounds = this.bounds;
        this.button.frame = bounds;
        this.button.layoutSubviews();
    }
}
/// <reference path="./UICore/UIViewController.ts" />
/// <reference path="./UICore/UIDialogView.ts" />
/// <reference path="./UICore/UILinkButton.ts" />
class RootViewController extends UIViewController {
    constructor(view) {
        // Calling super
        super(view);
        // Here are some suggested conventions that are used in UICore
        // Instance variables, it is good to initialize to nil or empty function, not leave as undefined to avoid
        // if blocks
        // this._firstView = nil;
        // this._secondView = nil;
        // this._testView = nil;
        // this._button = nil;
        // The nil object avoids unneccessary crashes by allowing you to call any function or access any variable on it, returning nil
        // Define properties with get and set functions so they can be accessed and set like variables
        // Name variables that should be private, like property variables, with a _ sign, this also holds for private functions
        // Avoid accessing variables and functions named with _ from outside as this creates strong coupling and hinders stability
        // Code for further setup if necessary
    }
    loadIntrospectionVariables() {
        super.loadIntrospectionVariables();
        this.superclass = UIViewController;
    }
    loadSubviews() {
        this.view.backgroundColor = UIColor.colorWithRGBA(225, 225, 225);
        document.body.style.backgroundColor = this.view.backgroundColor.stringValue;
        this.topBarView = new CellView("TopBarView").copyValuesRecursivelyFromObject({
            titleLabel: {
                text: "GETELEC - General Tool for Electron Emission Calculations",
                textAlignment: UITextView.textAlignment.center
            },
            colors: {
                background: { normal: new UIColor("#282c11") },
                titleLabel: { normal: UIColor.whiteColor }
            }
        });
        this.view.addSubview(this.topBarView);
        // Bottom bar
        this.bottomBarView = new BottomBarView("BottomBarView").copyValuesRecursivelyFromObject({
            style: { overflow: "hidden" }
        });
        this.view.addSubview(this.bottomBarView);
        this.backToMainButton = new CBFlatButton();
        this.backToMainButton.titleLabel.text = "&#8592;";
        this.backToMainButton.colors = {
            titleLabel: {
                normal: UIColor.whiteColor,
                highlighted: UIColor.whiteColor.colorWithAlpha(0.75),
                selected: UIColor.whiteColor.colorWithAlpha(0.5)
            },
            background: {
                normal: UIColor.transparentColor,
                highlighted: UIColor.transparentColor,
                selected: UIColor.transparentColor
            }
        };
        this.backToMainButton.calculateAndSetViewFrame = function () {
            this.setPosition(0, nil, 0, 0, nil, 50);
        };
        this.topBarView.addSubview(this.backToMainButton);
        this.backToMainButton.addControlEventTarget.EnterDown.PointerUpInside = (sender, event) => UIRoute.currentRoute.routeByRemovingComponentsOtherThanOnesNamed(["settings"]).apply();
    }
    handleRoute(route) {
        const _super = Object.create(null, {
            handleRoute: { get: () => super.handleRoute }
        });
        return __awaiter(this, void 0, void 0, function* () {
            _super.handleRoute.call(this, route);
            this.backToMainButton.hidden = NO;
            if (IS(route.componentWithName(GETELECViewController.routeComponentName))) {
                // Show article view
                if (!IS(this.useGETELECViewController)) {
                    this.useGETELECViewController = new GETELECViewController(new UIView("GETELECView"));
                }
                this.contentViewController = this.useGETELECViewController;
            }
            else if (IS(route.componentWithName(GETELECTableViewController.routeComponentName))) {
                // Show article view
                if (!IS(this.useGETELECViewController)) {
                    this.useGETELECTableViewController = new GETELECTableViewController(new UIView("GETELECTableView"));
                }
                this.contentViewController = this.useGETELECTableViewController;
            }
            else {
                // Main view controller
                if (!IS(this.mainViewController)) {
                    this.mainViewController = new UIViewController(new UIView("MainView"));
                    this.mainViewController.view.backgroundColor = UIColor.whiteColor;
                    this.mainViewController.view.viewHTMLElement.innerHTML = yield CBCore.sharedInstance.readTextFileFromURL("/frontpage.html");
                    //     LanguageService.stringForKey(
                    //     "frontPageHTML",
                    //     LanguageService.currentLanguageKey,
                    //     "This is the front page."
                    // )
                }
                this.contentViewController = this.mainViewController;
                this.backToMainButton.hidden = YES;
            }
        });
    }
    set contentViewController(controller) {
        if (this.contentViewController == controller) {
            return;
        }
        if (this.contentViewController) {
            this.removeChildViewController(this.contentViewController);
        }
        this._contentViewController = controller;
        this.addChildViewControllerInContainer(controller, this.view);
        this._layoutViewSubviews();
        this.view.setNeedsLayout();
    }
    get contentViewController() {
        return this._contentViewController || nil;
    }
    static get paddingLength() {
        return 20;
    }
    get paddingLength() {
        return this.class.paddingLength;
    }
    layoutViewsManually() {
        super.layoutViewsManually();
        // View bounds
        var bounds = this.view.bounds;
        this.topBarView.frame = bounds.rectangleWithHeight(50);
        this.contentViewController.view.frame = this.topBarView.frame.rectangleForNextRow(this.paddingLength, FIRST_OR_NIL(wrapInNil(this.contentViewController).intrinsicViewContentHeight([bounds.width, 1250].min()), wrapInNil(this.contentViewController).view.intrinsicContentHeight([bounds.width, 1250].min()))).rectangleWithWidth([bounds.width, 1250].min(), 0.5);
        this.contentViewController.view.style.boxShadow = "0 3px 6px 0 rgba(0, 0, 0, 0.1)";
        this.contentViewController.view.setMargins(0, 0, this.paddingLength, 0);
        //this.view.setPaddings(0, 0, this.paddingLength, 0);
        var bottomBarHeight = Math.max(100, this.bottomBarView.intrinsicContentHeight(this.contentViewController.view.frame.width));
        this.bottomBarView.frame = this.contentViewController.view.frame.rectangleWithY([
            this.contentViewController.view.frame.max.y + this.paddingLength * 2,
            UIView.pageHeight - bottomBarHeight
        ].max()).rectangleWithHeight(bottomBarHeight).rectangleWithWidth(this.view.bounds.width).rectangleWithX(0);
    }
}
// @ts-ignore
class UIRoute extends Array {
    constructor(hash) {
        super();
        this._isHandled = NO;
        this.completedComponents = [];
        if (!hash || !hash.startsWith) {
            return;
        }
        if (hash.startsWith("#")) {
            hash = hash.slice(1);
        }
        hash = decodeURIComponent(hash);
        const components = hash.split("]");
        components.forEach(function (component, index, array) {
            const componentName = component.split("[")[0];
            const parameters = {};
            if (!componentName) {
                return;
            }
            const parametersString = component.split("[")[1] || "";
            const parameterPairStrings = parametersString.split(",") || [];
            parameterPairStrings.forEach(function (pairString, index, array) {
                const keyAndValueArray = pairString.split(":");
                const key = decodeURIComponent(keyAndValueArray[0]);
                const value = decodeURIComponent(keyAndValueArray[1]);
                if (key) {
                    parameters[key] = value;
                }
            });
            this.push({
                name: componentName,
                parameters: parameters
            });
        }, this);
    }
    static get currentRoute() {
        return new UIRoute(window.location.hash);
    }
    apply() {
        window.location.hash = this.stringRepresentation;
    }
    applyByReplacingCurrentRouteInHistory() {
        window.location.replace(this.linkRepresentation);
    }
    copy() {
        var result = new UIRoute();
        result = Object.assign(result, this);
        return result;
    }
    childRoute() {
        var result = this.copy();
        result.completedComponents.forEach(function (component, index, completedComponents) {
            var indexInResult = result.indexOf(component);
            if (indexInResult > -1) {
                result.splice(indexInResult, 1);
            }
        });
        result.completedComponents = [];
        result.parentRoute = this;
        return result;
    }
    routeByRemovingComponentsOtherThanOnesNamed(componentNames) {
        const result = this.copy();
        const indexesToRemove = [];
        result.forEach(function (component, index, array) {
            if (!componentNames.contains(component.name)) {
                indexesToRemove.push(index);
            }
        });
        indexesToRemove.forEach(function (indexToRemove, index, array) {
            result.removeElementAtIndex(indexToRemove);
        });
        return result;
    }
    routeByRemovingComponentNamed(componentName) {
        const result = this.copy();
        const componentIndex = result.findIndex(function (component, index) {
            return (component.name == componentName);
        });
        if (componentIndex != -1) {
            result.splice(componentIndex, 1);
        }
        return result;
    }
    routeByRemovingParameterInComponent(componentName, parameterName, removeComponentIfEmpty = NO) {
        var result = this.copy();
        var parameters = result.componentWithName(componentName).parameters;
        if (IS_NOT(parameters)) {
            parameters = {};
        }
        delete parameters[parameterName];
        result = result.routeWithComponent(componentName, parameters);
        if (removeComponentIfEmpty && Object.keys(parameters).length == 0) {
            result = result.routeByRemovingComponentNamed(componentName);
        }
        return result;
    }
    routeBySettingParameterInComponent(componentName, parameterName, valueToSet) {
        var result = this.copy();
        if (IS_NIL(valueToSet) || IS_NIL(parameterName)) {
            return result;
        }
        var parameters = result.componentWithName(componentName).parameters;
        if (IS_NOT(parameters)) {
            parameters = {};
        }
        parameters[parameterName] = valueToSet;
        result = result.routeWithComponent(componentName, parameters);
        return result;
    }
    routeWithViewControllerComponent(viewController, parameters, extendParameters = NO) {
        return this.routeWithComponent(viewController.routeComponentName, parameters, extendParameters);
    }
    routeWithComponent(name, parameters, extendParameters = NO) {
        const result = this.copy();
        var component = result.componentWithName(name);
        if (IS_NOT(component)) {
            component = {
                name: name,
                parameters: {}
            };
            result.push(component);
        }
        if (IS_NOT(parameters)) {
            parameters = {};
        }
        if (extendParameters) {
            component.parameters = Object.assign(component.parameters, parameters);
        }
        else {
            component.parameters = parameters;
        }
        return result;
    }
    navigateBySettingComponent(name, parameters, extendParameters = NO) {
        this.routeWithComponent(name, parameters, extendParameters).apply();
    }
    componentWithViewController(viewController) {
        return this.componentWithName(viewController.routeComponentName);
    }
    componentWithName(name) {
        var result = nil;
        this.forEach(function (component, index, self) {
            if (component.name == name) {
                result = component;
            }
        });
        return result;
    }
    didcompleteComponent(component) {
        const self = this;
        const index = self.indexOf(component, 0);
        if (index > -1) {
            self.completedComponents.push(self.splice(index, 1)[0]);
            //self.completedComponents.push(component);
        }
    }
    set isHandled(isHandled) {
        this._isHandled = isHandled;
    }
    get isHandled() {
        return (this._isHandled || (this.length == 0 && this.completedComponents.length != 0));
    }
    get linkRepresentation() {
        return "#" + this.stringRepresentation;
    }
    get stringRepresentation() {
        var result = "";
        this.completedComponents.forEach(function (component, index, self) {
            result = result + component.name;
            const parameters = component.parameters;
            result = result + "[";
            Object.keys(parameters).forEach(function (key, index, keys) {
                if (index) {
                    result = result + ",";
                }
                result = result + encodeURIComponent(key) + ":" + encodeURIComponent(parameters[key]);
            });
            result = result + "]";
        });
        this.forEach(function (component, index, self) {
            result = result + component.name;
            const parameters = component.parameters;
            result = result + "[";
            Object.keys(parameters).forEach(function (key, index, keys) {
                if (index) {
                    result = result + ",";
                }
                result = result + encodeURIComponent(key) + ":" + encodeURIComponent(parameters[key]);
            });
            result = result + "]";
        });
        return result;
    }
}
/// <reference path="UIView.ts" />
/// <reference path="UIViewController.ts" />
/// <reference path="../RootViewController.ts" />
/// <reference path="UIRoute.ts" />
class UICore extends UIObject {
    constructor(rootDivElementID, rootViewControllerClass) {
        super();
        this.rootViewController = nil;
        this._class = UICore;
        this.superclass = UIObject;
        UICore.RootViewControllerClass = rootViewControllerClass;
        UICore.main = this;
        const rootViewElement = document.getElementById(rootDivElementID);
        const rootView = new UIView(rootDivElementID, rootViewElement);
        rootView.pausesPointerEvents = NO; //YES;
        if (UICore.RootViewControllerClass) {
            if (!(UICore.RootViewControllerClass.prototype instanceof UIViewController) ||
                UICore.RootViewControllerClass === UIViewController) {
                console.log("Error, UICore.RootViewControllerClass must be a or a subclass of UIViewController, falling back to UIViewController.");
                UICore.RootViewControllerClass = UIViewController;
            }
            this.rootViewController = new UICore.RootViewControllerClass(rootView);
        }
        else {
            this.rootViewController = new UIViewController(rootView);
        }
        this.rootViewController.viewWillAppear();
        this.rootViewController.viewDidAppear();
        this.rootViewController.view.addTargetForControlEvent(UIView.controlEvent.PointerUpInside, function (sender, event) {
            document.activeElement.blur();
        });
        const windowDidResize = function () {
            // Doing layout two times to prevent page scrollbars from confusing the layout
            this.rootViewController._layoutViewSubviews();
            UIView.layoutViewsIfNeeded();
            this.rootViewController._layoutViewSubviews();
            //UIView.layoutViewsIfNeeded()
            this.rootViewController.view.broadcastEventInSubtree({
                name: UICore.broadcastEventName.WindowDidResize,
                parameters: nil
            });
        };
        window.addEventListener("resize", windowDidResize.bind(this));
        const didScroll = function () {
            //code
            this.rootViewController.view.broadcastEventInSubtree({
                name: UIView.broadcastEventName.PageDidScroll,
                parameters: nil
            });
        }.bind(this);
        window.addEventListener("scroll", didScroll, false);
        const hashDidChange = function () {
            //code
            this.rootViewController.handleRouteRecursively(UIRoute.currentRoute);
            this.rootViewController.view.broadcastEventInSubtree({
                name: UICore.broadcastEventName.RouteDidChange,
                parameters: nil
            });
        }.bind(this);
        window.addEventListener("hashchange", hashDidChange.bind(this), false);
        hashDidChange();
    }
    static loadClass(className) {
        if (window[className]) {
            return;
        }
        document.writeln("<script type='text/javascript' src='dist/UICore/" + className + ".js'></script>");
    }
}
UICore.languageService = nil;
UICore.broadcastEventName = {
    "RouteDidChange": "RouteDidChange",
    "WindowDidResize": "WindowDidResize"
};
UICore.RootViewControllerClass = nil;
const IS_FIREFOX = navigator.userAgent.toLowerCase().indexOf("firefox") > -1;
const IS_SAFARI = /^((?!chrome|android).)*safari/i.test(navigator.userAgent);
Array.prototype.indexOf || (Array.prototype.indexOf = function (d, e) {
    var a;
    if (null == this) {
        throw new TypeError("\"this\" is null or not defined");
    }
    const c = Object(this), b = c.length >>> 0;
    if (0 === b) {
        return -1;
    }
    a = +e || 0;
    Infinity === Math.abs(a) && (a = 0);
    if (a >= b) {
        return -1;
    }
    for (a = Math.max(0 <= a ? a : b - Math.abs(a), 0); a < b;) {
        if (a in c && c[a] === d) {
            return a;
        }
        a++;
    }
    return -1;
});
/// <reference path="UIView.ts" />
/// <reference path="UICore.ts" />
class UITextView extends UIView {
    constructor(elementID, textViewType = UITextView.type.paragraph, viewHTMLElement = null) {
        super(elementID, viewHTMLElement, textViewType);
        this._textColor = UITextView.defaultTextColor;
        this._isSingleLine = YES;
        this.textPrefix = "";
        this.textSuffix = "";
        this._notificationAmount = 0;
        this._minFontSize = nil;
        this._maxFontSize = nil;
        this._automaticFontSizeSelection = NO;
        this.changesOften = NO;
        this._class = UITextView;
        this.superclass = UIView;
        this.text = "";
        this.style.overflow = "hidden";
        this.style.textOverflow = "ellipsis";
        this.isSingleLine = YES;
        this.textColor = this.textColor;
        this.userInteractionEnabled = YES;
        if (textViewType == UITextView.type.textArea) {
            this.pausesPointerEvents = YES;
            this.addTargetForControlEvent(UIView.controlEvent.PointerUpInside, function (sender, event) {
                sender.focus();
            });
        }
    }
    static _determinePXAndPTRatios() {
        const o = document.createElement("div");
        o.style.width = "1000pt";
        document.body.appendChild(o);
        UITextView._ptToPx = o.clientWidth / 1000;
        document.body.removeChild(o);
        UITextView._pxToPt = 1 / UITextView._ptToPx;
    }
    get textAlignment() {
        const result = this.style.textAlign;
        return result;
    }
    set textAlignment(textAlignment) {
        this._textAlignment = textAlignment;
        this.style.textAlign = textAlignment;
    }
    get textColor() {
        const result = this._textColor;
        return result;
    }
    set textColor(color) {
        this._textColor = color || UITextView.defaultTextColor;
        this.style.color = this._textColor.stringValue;
    }
    get isSingleLine() {
        return this._isSingleLine;
    }
    set isSingleLine(isSingleLine) {
        this._isSingleLine = isSingleLine;
        if (isSingleLine) {
            this.style.whiteSpace = "pre";
            return;
        }
        this.style.whiteSpace = "pre-wrap";
    }
    get notificationAmount() {
        return this._notificationAmount;
    }
    set notificationAmount(notificationAmount) {
        if (this._notificationAmount == notificationAmount) {
            return;
        }
        this._notificationAmount = notificationAmount;
        this.text = this.text;
        this.setNeedsLayoutUpToRootView();
        this.notificationAmountDidChange(notificationAmount);
    }
    notificationAmountDidChange(notificationAmount) {
    }
    get text() {
        return (this._text || this.viewHTMLElement.innerHTML);
    }
    set text(text) {
        this._text = text;
        var notificationText = "";
        if (this.notificationAmount) {
            notificationText = "<span style=\"color: " + UITextView.notificationTextColor.stringValue + ";\">" +
                (" (" + this.notificationAmount + ")").bold() + "</span>";
        }
        if (this.viewHTMLElement.innerHTML != this.textPrefix + text + this.textSuffix + notificationText) {
            this.viewHTMLElement.innerHTML = this.textPrefix + FIRST(text, "") + this.textSuffix + notificationText;
        }
        this.setNeedsLayout();
    }
    set innerHTML(innerHTML) {
        this.text = innerHTML;
    }
    get innerHTML() {
        return this.viewHTMLElement.innerHTML;
    }
    setText(key, defaultString, parameters) {
        this.setInnerHTML(key, defaultString, parameters);
    }
    get fontSize() {
        const style = window.getComputedStyle(this.viewHTMLElement, null).fontSize;
        const result = (parseFloat(style) * UITextView._pxToPt);
        return result;
    }
    set fontSize(fontSize) {
        this.style.fontSize = "" + fontSize + "pt";
    }
    useAutomaticFontSize(minFontSize = nil, maxFontSize = nil) {
        this._automaticFontSizeSelection = YES;
        this._minFontSize = minFontSize;
        this._maxFontSize = maxFontSize;
        this.setNeedsLayout();
    }
    static automaticallyCalculatedFontSize(bounds, currentSize, currentFontSize, minFontSize, maxFontSize) {
        minFontSize = FIRST(minFontSize, 1);
        maxFontSize = FIRST(maxFontSize, 100000000000);
        const heightMultiplier = bounds.height / (currentSize.height + 1);
        const widthMultiplier = bounds.width / (currentSize.width + 1);
        var multiplier = heightMultiplier;
        if (heightMultiplier > widthMultiplier) {
            multiplier = widthMultiplier;
        }
        const maxFittingFontSize = currentFontSize * multiplier;
        if (maxFittingFontSize > maxFontSize) {
            return maxFontSize;
        }
        if (minFontSize > maxFittingFontSize) {
            return minFontSize;
        }
        return maxFittingFontSize;
    }
    didReceiveBroadcastEvent(event) {
        super.didReceiveBroadcastEvent(event);
    }
    willMoveToSuperview(superview) {
        super.willMoveToSuperview(superview);
    }
    layoutSubviews() {
        super.layoutSubviews();
        if (this._automaticFontSizeSelection) {
            this.fontSize = UITextView.automaticallyCalculatedFontSize(new UIRectangle(0, 0, 1 *
                this.viewHTMLElement.offsetHeight, 1 *
                this.viewHTMLElement.offsetWidth), this.intrinsicContentSize(), this.fontSize, this._minFontSize, this._maxFontSize);
        }
    }
    intrinsicContentHeight(constrainingWidth = 0) {
        if (this.changesOften) {
            return super.intrinsicContentHeight(constrainingWidth);
        }
        const keyPath = (this.viewHTMLElement.innerHTML + "_csf_" + this.computedStyle.font).replace(new RegExp("\\.", "g"), "_") + "." +
            ("" + constrainingWidth).replace(new RegExp("\\.", "g"), "_");
        var result = UITextView._intrinsicHeightCache.valueForKeyPath(keyPath);
        if (IS_LIKE_NULL(result)) {
            result = super.intrinsicContentHeight(constrainingWidth);
            UITextView._intrinsicHeightCache.setValueForKeyPath(keyPath, result);
        }
        return result;
    }
    intrinsicContentWidth(constrainingHeight = 0) {
        if (this.changesOften) {
            return super.intrinsicContentWidth(constrainingHeight);
        }
        const keyPath = (this.viewHTMLElement.innerHTML + "_csf_" + this.computedStyle.font).replace(new RegExp("\\.", "g"), "_") + "." +
            ("" + constrainingHeight).replace(new RegExp("\\.", "g"), "_");
        var result = UITextView._intrinsicWidthCache.valueForKeyPath(keyPath);
        if (IS_LIKE_NULL(result)) {
            result = super.intrinsicContentWidth(constrainingHeight);
            UITextView._intrinsicWidthCache.setValueForKeyPath(keyPath, result);
        }
        return result;
    }
    intrinsicContentSize() {
        // This works but is slow
        const result = this.intrinsicContentSizeWithConstraints(nil, nil);
        return result;
    }
}
UITextView.defaultTextColor = UIColor.blackColor;
UITextView.notificationTextColor = UIColor.redColor;
UITextView._intrinsicHeightCache = new UIObject();
UITextView._intrinsicWidthCache = new UIObject();
UITextView.type = {
    "paragraph": "p",
    "header1": "h1",
    "header2": "h2",
    "header3": "h3",
    "header4": "h4",
    "header5": "h5",
    "header6": "h6",
    "textArea": "textarea",
    "textField": "input",
    "span": "span",
    "label": "label"
};
UITextView.textAlignment = {
    "left": "left",
    "center": "center",
    "right": "right",
    "justify": "justify"
};
UITextView._determinePXAndPTRatios();
// /**
//  * Uses canvas.measureText to compute and return the width of the given text of given font in pixels.
//  * 
//  * @param {String} text The text to be rendered.
//  * @param {String} font The css font descriptor that text is to be rendered with (e.g. "bold 14px verdana").
//  * 
//  * @see https://stackoverflow.com/questions/118241/calculate-text-width-with-javascript/21015393#21015393
//  */
// function getTextMetrics(text, font) {
//     // re-use canvas object for better performance
//     var canvas = getTextMetrics.canvas || (getTextMetrics.canvas = document.createElement("canvas"));
//     var context = canvas.getContext("2d");
//     context.font = font;
//     var metrics = context.measureText(text);
//     return metrics;
// }
/// <reference path="./UICore/UIViewController.ts" />
///<reference path="UICore/UITextView.ts"/>
class GETELECTableViewController extends UIViewController {
    constructor(view) {
        // Calling super
        super(view);
        // Code for further setup if necessary
    }
    loadIntrospectionVariables() {
        super.loadIntrospectionVariables();
        this.superclass = UIViewController;
    }
    loadSubviews() {
        this.view.backgroundColor = UIColor.whiteColor;
        this.titleLabel = new UITextView(this.view.elementID + "TitleLabel", UITextView.type.header3);
        this.titleLabel.textAlignment = UITextView.textAlignment.left;
        this.titleLabel.nativeSelectionEnabled = NO;
        this.titleLabel.isSingleLine = NO;
        this.view.addSubview(this.titleLabel);
        this.inputTextArea = new UITextArea(this.view.elementID + "InputTextArea");
        this.inputTextArea.placeholderText = "Input your data here.";
        this.inputTextArea.changesOften = YES;
        this.view.addSubview(this.inputTextArea);
        this.loadDataButton = new CBButton();
        this.loadDataButton.titleLabel.text = "Load data";
        this.view.addSubview(this.loadDataButton);
        this.loadDataButton.enabled = NO;
        this.resultsLabel = new UIView();
        this.resultsLabel.innerHTML = "Results";
        this.resultsLabel.hidden = YES;
        this.view.addSubview(this.resultsLabel);
        this.chartView = new UIView("ResultDataChart", nil, "canvas");
        this.chartView.hidden = YES;
        this.view.addSubview(this.chartView);
        this.inputTextArea.addControlEventTarget.TextChange = (sender, event) => {
            //console.log(sender);
            this.view.setNeedsLayoutUpToRootView();
            this.loadDataButton.enabled = IS(this.inputTextArea.text);
        };
        // @ts-ignore
        const ctx = this.chartView.viewHTMLElement.getContext("2d");
        // @ts-ignore
        this.chart = new Chart(ctx, {
            // The type of chart we want to create
            //type: "line",
            type: "scatter",
            // The data for our dataset
            data: {
                //labels: ['January', 'February', 'March', 'April', 'May', 'June'],
                datasets: [
                    {
                        label: "Input Data",
                        fill: false,
                        lineTension: 0,
                        //backgroundColor: UIColor.transparentColor,
                        borderColor: "rgb(125, 115, 252)" //,
                        //data: [{ x: 0, y: 0 }, { x: 2, y: 10 }, { x: 5, y: 5 }, { x: 7, y: 2 }, { x: 8, y: 20 }, {
                        // x: 11, y: 45 }]
                    },
                    {
                        label: "Fitted theoretical line",
                        fill: false,
                        lineTension: 0,
                        backgroundColor: "rgb(255, 99, 132)",
                        borderColor: "rgb(255, 99, 132)" //,
                        //data: [{ x: 0, y: 0 }, { x: 2, y: 10 }, { x: 5, y: 5 }, { x: 7, y: 2 }, { x: 8, y: 20 }, {
                        // x: 11, y: 45 }]
                    }
                ]
            },
            // Configuration options go here
            options: {
                fill: false,
                lineTension: 0,
                scales: {
                    x: {
                        type: "linear",
                        position: "bottom"
                    },
                    y: {
                        type: "logarithmic",
                        position: "bottom"
                    },
                    yAxes: [{
                            type: 'logarithmic',
                            ticks: {
                                min: 1,
                                max: 1000000,
                                callback: function (value, index, values) {
                                    if (value === 1000000)
                                        return "1M";
                                    if (value === 100000)
                                        return "100K";
                                    if (value === 10000)
                                        return "10K";
                                    if (value === 1000)
                                        return "1K";
                                    if (value === 100)
                                        return "100";
                                    if (value === 10)
                                        return "10";
                                    if (value === 1)
                                        return "1";
                                    return null;
                                }
                            }
                        }]
                },
                tooltips: {
                    mode: "interpolate",
                    intersect: true
                },
                plugins: {
                    crosshair: {
                        line: {
                            color: "#F66",
                            width: 1 // crosshair line width
                        },
                        sync: {
                            enabled: true,
                            group: 1,
                            suppressTooltips: false // suppress tooltips when showing a synced tracer
                        },
                        zoom: {
                            enabled: NO,
                            zoomboxBackgroundColor: "rgba(66,133,244,0.2)",
                            zoomboxBorderColor: "#48F",
                            zoomButtonText: "Reset Zoom",
                            zoomButtonClass: "reset-zoom" // reset zoom button class
                        },
                        callbacks: {
                            beforeZoom: function (start, end) {
                                return true;
                            },
                            afterZoom: function (start, end) {
                            }
                        }
                    }
                }
                // responsive: false,
                // maintainAspectRatio: false,
                // showScale: true
            }
        });
        this.loadDataButton.addControlEventTarget.EnterDown.PointerUpInside = (sender, event) => __awaiter(this, void 0, void 0, function* () {
            //console.log(sender)
            try {
                CBDialogViewShower.showActionIndicatorDialog("Loading.");
                const stringSocketClientResult = yield SocketClient.PerformFitFun({
                    inputData: this.inputTextArea.text
                });
                if (IS(stringSocketClientResult.errorResult)) {
                    console.log(stringSocketClientResult.errorResult);
                    CBDialogViewShower.alert("Failed to analyze data. " + JSON.stringify(stringSocketClientResult.errorResult));
                    return;
                }
                const result = JSON.parse(stringSocketClientResult.result);
                //console.log(result)
                this.resultsLabel.innerHTML = FIRST(result.resultHTMLString, "")
                    + IF(IS_NOT(result.resultHTMLString))(RETURNER(`Radius: ${result.Radius} beta: ${result.beta} sigma_Aeff: ${result.sigma_Aeff}`))();
                const pointPoints = this.pointObjectsFromValues(result.xplot_mrk, result.yplot_mrk);
                const linePoints = this.pointObjectsFromValues(result.xplot_line, result.yplot_line);
                this.chart.data.datasets[0].data = pointPoints;
                this.chart.data.datasets[0].showLine = NO;
                this.chart.data.datasets[1].data = linePoints;
                this.chart.data.datasets[1].showLine = YES;
                this.chart.data.datasets[1].pointRadius = 0;
                this.chart.data.datasets[1].interpolate = YES;
                this.chart.options.scales.xAxes[0].scaleLabel.display = true;
                this.chart.options.scales.xAxes[0].scaleLabel.labelString = "" + FIRST(result.xAxisUnit, "");
                this.chart.options.scales.yAxes[0].scaleLabel.display = true;
                this.chart.options.scales.yAxes[0].scaleLabel.labelString = "" + FIRST(result.yAxisUnit, "");
                this.resultsLabel.hidden = NO;
                this.chartView.hidden = NO;
                this.chart.update();
                this.view.setNeedsLayoutUpToRootView();
            }
            catch (exception) {
                console.log(exception);
                CBDialogViewShower.alert("Failed to display results. " + JSON.stringify(exception));
            }
            CBDialogViewShower.hideActionIndicatorDialog();
        });
    }
    pointObjectsFromValues(xValues, yValues) {
        var resultPoints = [];
        for (var i = 0; i < xValues.length; i++) {
            var x = xValues[i];
            var y = yValues[i];
            resultPoints.push({ x: x, y: y });
        }
        return resultPoints;
    }
    handleRoute(route) {
        const _super = Object.create(null, {
            handleRoute: { get: () => super.handleRoute }
        });
        return __awaiter(this, void 0, void 0, function* () {
            _super.handleRoute.call(this, route);
            const inquiryComponent = route.componentWithViewController(GETELECTableViewController);
            this.titleLabel.text = "Enter your data and press the load button. TABLE VIEW";
            route.didcompleteComponent(inquiryComponent);
        });
    }
    layoutViewsManually() {
        super.layoutViewsManually();
        const padding = RootViewController.paddingLength;
        const labelHeight = padding * 1.25;
        // View bounds
        var bounds = this.view.bounds;
        this.view.setPaddings(0, 0, padding, 0);
        this.titleLabel.frame = bounds.rectangleWithInsets(padding, padding, 0, padding)
            .rectangleWithHeight(this.titleLabel.intrinsicContentHeight(bounds.width));
        this.inputTextArea.frame = this.titleLabel.frame.rectangleForNextRow(padding * 0.25, this.inputTextArea.intrinsicContentHeight(this.titleLabel.frame.width) + 5);
        this.loadDataButton.frame = this.inputTextArea.frame.rectangleForNextRow(padding, labelHeight * 2);
        this.loadDataButton.setMargins(0, 0, padding, 0);
        this.resultsLabel.frame = this.loadDataButton.frame.rectangleForNextRow(padding);
        this.chartView.frame = this.resultsLabel.frame.rectangleForNextRow(padding, this.loadDataButton.frame.width * 0.5);
        //.rectangleWithInsets(padding * 5 * 0, padding * 5, 0, 0)
        //.rectangleByAddingX(-padding)
        //this.chartView.setPadding(padding);
        this.chartView.setMaxSizes(this.chartView.frame.height, this.chartView.frame.width);
    }
    intrinsicViewContentHeight(constrainingWidth = 0) {
        const padding = RootViewController.paddingLength;
        const labelHeight = padding * 1.5;
        var result = padding + this.titleLabel.intrinsicContentHeight(constrainingWidth) + padding * 0.25 +
            labelHeight * 2 + padding * 0.5 + padding +
            this.inputTextArea.intrinsicContentHeight(constrainingWidth);
        if (IS_NOT(this.chartView.hidden)) {
            result = result + padding + (constrainingWidth - padding * 2) * 0.5 + padding + labelHeight * 2;
        }
        return result;
    }
}
GETELECTableViewController.routeComponentName = "fit_funasdasd";
GETELECTableViewController.ParameterIdentifierName = { "ID": "ID" };
// @ts-ignore
Chart.Interaction.modes["interpolate"] = function (chart, e, options) {
    // This function has a separate license
    // MIT License
    //
    // Copyright (c) 2018 Abel Heinsbroek
    //
    // Permission is hereby granted, free of charge, to any person obtaining a copy
    // of this software and associated documentation files (the "Software"), to deal
    // in the Software without restriction, including without limitation the rights
    // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    // copies of the Software, and to permit persons to whom the Software is
    // furnished to do so, subject to the following conditions:
    //
    //     The above copyright notice and this permission notice shall be included in all
    // copies or substantial portions of the Software.
    //
    //     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    //     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    //     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    // SOFTWARE.
    var items = [];
    for (var datasetIndex = 0; datasetIndex < chart.data.datasets.length; datasetIndex++) {
        // check for interpolate setting
        if (!chart.data.datasets[datasetIndex].interpolate) {
            continue;
        }
        var meta = chart.getDatasetMeta(datasetIndex);
        // do not interpolate hidden charts
        if (meta.hidden) {
            continue;
        }
        var xScale = chart.scales[meta.xAxisID];
        var yScale = chart.scales[meta.yAxisID];
        var xValue = xScale.getValueForPixel(e.x);
        var data = chart.data.datasets[datasetIndex].data;
        var index = data.findIndex(function (o) {
            return o.x >= xValue;
        });
        if (index === -1) {
            continue;
        }
        // linear interpolate value
        var prev = data[index - 1];
        var next = data[index];
        if (prev && next) {
            var slope = (next.y - prev.y) / (next.x - prev.x);
            var interpolatedValue = prev.y + (xValue - prev.x) * slope;
        }
        if (chart.data.datasets[datasetIndex].steppedLine && prev) {
            interpolatedValue = prev.y;
        }
        if (isNaN(interpolatedValue)) {
            continue;
        }
        var yPosition = yScale.getPixelForValue(interpolatedValue);
        // do not interpolate values outside of the axis limits
        if (isNaN(yPosition)) {
            continue;
        }
        // create a 'fake' event point
        var fakePoint = {
            value: interpolatedValue,
            xValue: xValue,
            tooltipPosition: function () {
                return this._model;
            },
            hasValue: function () {
                return true;
            },
            _model: {
                x: e.x,
                y: yPosition
            },
            _datasetIndex: datasetIndex,
            _index: items.length,
            _xScale: {
                getLabelForIndex: function (indx) {
                    let xValue = items[indx].xValue;
                    xValue = xValue.toPrecision(5);
                    return xValue;
                }
            },
            _yScale: {
                getLabelForIndex: function (indx) {
                    let value = items[indx].value.toPrecision(5);
                    return value;
                }
            },
            _chart: chart
        };
        items.push(fakePoint);
    }
    // add other, not interpolated, items
    // @ts-ignore
    var xItems = Chart.Interaction.modes.x(chart, e, options);
    xItems.forEach((item, index, array) => {
        if (!chart.data.datasets[item._datasetIndex].interpolate) {
            items.push(item);
        }
    });
    return items;
};
/// <reference path="./UICore/UIViewController.ts" />
///<reference path="UICore/UITextView.ts"/>
class GETELECViewController extends UIViewController {
    constructor(view) {
        // Calling super
        super(view);
        // Code for further setup if necessary
    }
    loadIntrospectionVariables() {
        super.loadIntrospectionVariables();
        this.superclass = UIViewController;
    }
    loadSubviews() {
        this.view.backgroundColor = UIColor.whiteColor;
        this.titleLabel = new UITextView(this.view.elementID + "TitleLabel", UITextView.type.header3);
        this.titleLabel.textAlignment = UITextView.textAlignment.left;
        this.titleLabel.nativeSelectionEnabled = NO;
        this.titleLabel.isSingleLine = NO;
        this.view.addSubview(this.titleLabel);
        this.inputTextArea = new UITextArea(this.view.elementID + "InputTextArea");
        this.inputTextArea.placeholderText = "Input your data here.";
        this.inputTextArea.changesOften = YES;
        this.view.addSubview(this.inputTextArea);
        this.loadDataButton = new CBButton();
        this.loadDataButton.titleLabel.text = "Load data";
        this.view.addSubview(this.loadDataButton);
        this.loadDataButton.enabled = NO;
        this.resultsLabel = new UIView();
        this.resultsLabel.innerHTML = "Results";
        this.resultsLabel.hidden = YES;
        this.view.addSubview(this.resultsLabel);
        this.chartView = new UIView("ResultDataChart", nil, "canvas");
        this.chartView.hidden = YES;
        this.view.addSubview(this.chartView);
        this.inputTextArea.addControlEventTarget.TextChange = (sender, event) => {
            //console.log(sender);
            this.view.setNeedsLayoutUpToRootView();
            this.loadDataButton.enabled = IS(this.inputTextArea.text);
        };
        // @ts-ignore
        const ctx = this.chartView.viewHTMLElement.getContext("2d");
        // @ts-ignore
        this.chart = new Chart(ctx, {
            // The type of chart we want to create
            //type: "line",
            type: "scatter",
            // The data for our dataset
            data: {
                //labels: ['January', 'February', 'March', 'April', 'May', 'June'],
                datasets: [
                    {
                        label: "Input Data",
                        fill: false,
                        lineTension: 0,
                        //backgroundColor: UIColor.transparentColor,
                        borderColor: "rgb(125, 115, 252)" //,
                        //data: [{ x: 0, y: 0 }, { x: 2, y: 10 }, { x: 5, y: 5 }, { x: 7, y: 2 }, { x: 8, y: 20 }, {
                        // x: 11, y: 45 }]
                    },
                    {
                        label: "Fitted line",
                        fill: false,
                        lineTension: 0,
                        backgroundColor: "rgb(255, 99, 132)",
                        borderColor: "rgb(255, 99, 132)" //,
                        //data: [{ x: 0, y: 0 }, { x: 2, y: 10 }, { x: 5, y: 5 }, { x: 7, y: 2 }, { x: 8, y: 20 }, {
                        // x: 11, y: 45 }]
                    }
                ]
            },
            // Configuration options go here
            options: {
                fill: false,
                lineTension: 0,
                scales: {
                    x: {
                        type: "linear",
                        position: "bottom"
                    },
                    y: {
                        type: "logarithmic",
                        position: "bottom"
                    },
                    yAxes: [
                        {
                            type: "logarithmic",
                            ticks: {
                                min: 0,
                                max: 1000000,
                                callback: function (value, index, values) {
                                    if (value === 1000000) {
                                        return "1M";
                                    }
                                    if (value === 100000) {
                                        return "100K";
                                    }
                                    if (value === 10000) {
                                        return "10K";
                                    }
                                    if (value === 1000) {
                                        return "1K";
                                    }
                                    if (value === 100) {
                                        return "100";
                                    }
                                    if (value === 10) {
                                        return "10";
                                    }
                                    if (value === 1) {
                                        return "1";
                                    }
                                    if (value === 0) {
                                        return "0";
                                    }
                                    return null;
                                }
                            }
                        }
                    ]
                },
                tooltips: {
                    mode: "interpolate",
                    intersect: true
                },
                plugins: {
                    crosshair: {
                        line: {
                            color: "#F66",
                            width: 1 // crosshair line width
                        },
                        sync: {
                            enabled: true,
                            group: 1,
                            suppressTooltips: false // suppress tooltips when showing a synced tracer
                        },
                        zoom: {
                            enabled: NO,
                            zoomboxBackgroundColor: "rgba(66,133,244,0.2)",
                            zoomboxBorderColor: "#48F",
                            zoomButtonText: "Reset Zoom",
                            zoomButtonClass: "reset-zoom" // reset zoom button class
                        },
                        callbacks: {
                            beforeZoom: function (start, end) {
                                return true;
                            },
                            afterZoom: function (start, end) {
                            }
                        }
                    }
                }
                // responsive: false,
                // maintainAspectRatio: false,
                // showScale: true
            }
        });
        this.loadDataButton.addControlEventTarget.EnterDown.PointerUpInside = (sender, event) => __awaiter(this, void 0, void 0, function* () {
            //console.log(sender)
            try {
                CBDialogViewShower.showActionIndicatorDialog("Loading.");
                const stringSocketClientResult = yield SocketClient.PerformFitFun({
                    inputData: this.inputTextArea.text
                });
                if (IS(stringSocketClientResult.errorResult)) {
                    console.log(stringSocketClientResult.errorResult);
                    CBDialogViewShower.alert("Failed to analyze data. " + JSON.stringify(stringSocketClientResult.errorResult));
                    CBDialogViewShower.hideActionIndicatorDialog();
                    return;
                }
                const result = JSON.parse(stringSocketClientResult.result);
                //console.log(result)
                this.resultsLabel.innerHTML = FIRST(result.resultHTMLString, "") +
                    IF(IS_NOT(result.resultHTMLString))(RETURNER(`Radius: ${result.Radius} &beta: ${result.beta} &sigma Aeff: ${result.sigma_Aeff}`))
                        .ELSE(RETURNER(""));
                const pointPoints = this.pointObjectsFromValues(result.xplot_mrk, result.yplot_mrk);
                const linePoints = this.pointObjectsFromValues(result.xplot_line, result.yplot_line);
                this.chart.data.datasets[0].data = pointPoints;
                this.chart.data.datasets[0].showLine = NO;
                this.chart.data.datasets[1].data = linePoints;
                this.chart.data.datasets[1].showLine = YES;
                this.chart.data.datasets[1].pointRadius = 0;
                this.chart.data.datasets[1].interpolate = YES;
                this.chart.options.scales.xAxes[0].scaleLabel.display = true;
                this.chart.options.scales.xAxes[0].scaleLabel.labelString = "" + FIRST(result.xAxisUnit, "");
                this.chart.options.scales.yAxes[0].scaleLabel.display = true;
                this.chart.options.scales.yAxes[0].scaleLabel.labelString = "" + FIRST(result.yAxisUnit, "");
                this.resultsLabel.hidden = NO;
                this.chartView.hidden = NO;
                this.chart.update();
                this.view.setNeedsLayoutUpToRootView();
            }
            catch (exception) {
                console.log(exception);
                CBDialogViewShower.alert("Failed to display results. " + JSON.stringify(exception));
            }
            CBDialogViewShower.hideActionIndicatorDialog();
        });
    }
    pointObjectsFromValues(xValues, yValues) {
        var resultPoints = [];
        for (var i = 0; i < xValues.length; i++) {
            var x = xValues[i];
            var y = yValues[i];
            resultPoints.push({ x: x, y: y });
        }
        return resultPoints;
    }
    handleRoute(route) {
        const _super = Object.create(null, {
            handleRoute: { get: () => super.handleRoute }
        });
        return __awaiter(this, void 0, void 0, function* () {
            _super.handleRoute.call(this, route);
            const inquiryComponent = route.componentWithViewController(GETELECViewController);
            this.titleLabel.text = "Enter your data and press the load button.";
            route.didcompleteComponent(inquiryComponent);
        });
    }
    layoutViewsManually() {
        super.layoutViewsManually();
        const padding = RootViewController.paddingLength;
        const labelHeight = padding * 1.25;
        // View bounds
        var bounds = this.view.bounds;
        this.view.setPaddings(0, 0, padding, 0);
        this.titleLabel.frame = bounds.rectangleWithInsets(padding, padding, 0, padding)
            .rectangleWithHeight(this.titleLabel.intrinsicContentHeight(bounds.width - padding * 2));
        this.inputTextArea.frame = this.titleLabel.frame.rectangleForNextRow(padding * 0.25, this.inputTextArea.intrinsicContentHeight(this.titleLabel.frame.width) + 5);
        this.loadDataButton.frame = this.inputTextArea.frame.rectangleForNextRow(padding, labelHeight * 2);
        this.loadDataButton.setMargins(0, 0, padding, 0);
        this.resultsLabel.frame = this.loadDataButton.frame.rectangleForNextRow(padding);
        this.chartView.frame = this.resultsLabel.frame.rectangleForNextRow(padding, this.loadDataButton.frame.width * 0.5);
        //.rectangleWithInsets(padding * 5 * 0, padding * 5, 0, 0)
        //.rectangleByAddingX(-padding)
        //this.chartView.setPadding(padding);
        this.chartView.setMaxSizes(this.chartView.frame.height, this.chartView.frame.width);
    }
    intrinsicViewContentHeight(constrainingWidth = 0) {
        const padding = RootViewController.paddingLength;
        const labelHeight = padding * 1.5;
        var result = padding + this.titleLabel.intrinsicContentHeight(constrainingWidth) + padding * 0.25 +
            labelHeight * 2 + padding * 0.5 + padding +
            this.inputTextArea.intrinsicContentHeight(constrainingWidth);
        if (IS_NOT(this.chartView.hidden)) {
            result = result + padding + (constrainingWidth - padding * 2) * 0.5 + padding + labelHeight * 2;
        }
        return result;
    }
}
GETELECViewController.routeComponentName = "fit_fun";
GETELECViewController.ParameterIdentifierName = { "ID": "ID" };
// @ts-ignore
Chart.Interaction.modes["interpolate"] = function (chart, e, options) {
    // This function has a separate license
    // MIT License
    //
    // Copyright (c) 2018 Abel Heinsbroek
    //
    // Permission is hereby granted, free of charge, to any person obtaining a copy
    // of this software and associated documentation files (the "Software"), to deal
    // in the Software without restriction, including without limitation the rights
    // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    // copies of the Software, and to permit persons to whom the Software is
    // furnished to do so, subject to the following conditions:
    //
    //     The above copyright notice and this permission notice shall be included in all
    // copies or substantial portions of the Software.
    //
    //     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    //     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    //     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    // SOFTWARE.
    var items = [];
    for (var datasetIndex = 0; datasetIndex < chart.data.datasets.length; datasetIndex++) {
        // check for interpolate setting
        if (!chart.data.datasets[datasetIndex].interpolate) {
            continue;
        }
        var meta = chart.getDatasetMeta(datasetIndex);
        // do not interpolate hidden charts
        if (meta.hidden) {
            continue;
        }
        var xScale = chart.scales[meta.xAxisID];
        var yScale = chart.scales[meta.yAxisID];
        var xValue = xScale.getValueForPixel(e.x);
        var data = chart.data.datasets[datasetIndex].data;
        var index = data.findIndex(function (o) {
            return o.x >= xValue;
        });
        if (index === -1) {
            continue;
        }
        // linear interpolate value
        var prev = data[index - 1];
        var next = data[index];
        if (prev && next) {
            var slope = (next.y - prev.y) / (next.x - prev.x);
            var interpolatedValue = prev.y + (xValue - prev.x) * slope;
        }
        if (chart.data.datasets[datasetIndex].steppedLine && prev) {
            interpolatedValue = prev.y;
        }
        if (isNaN(interpolatedValue)) {
            continue;
        }
        var yPosition = yScale.getPixelForValue(interpolatedValue);
        // do not interpolate values outside of the axis limits
        if (isNaN(yPosition)) {
            continue;
        }
        // create a 'fake' event point
        var fakePoint = {
            value: interpolatedValue,
            xValue: xValue,
            tooltipPosition: function () {
                return this._model;
            },
            hasValue: function () {
                return true;
            },
            _model: {
                x: e.x,
                y: yPosition
            },
            _datasetIndex: datasetIndex,
            _index: items.length,
            _xScale: {
                getLabelForIndex: function (indx) {
                    let xValue = items[indx].xValue;
                    xValue = xValue.toPrecision(5);
                    return xValue;
                }
            },
            _yScale: {
                getLabelForIndex: function (indx) {
                    let value = items[indx].value.toPrecision(5);
                    return value;
                }
            },
            _chart: chart
        };
        items.push(fakePoint);
    }
    // add other, not interpolated, items
    // @ts-ignore
    var xItems = Chart.Interaction.modes.x(chart, e, options);
    xItems.forEach((item, index, array) => {
        if (!chart.data.datasets[item._datasetIndex].interpolate) {
            items.push(item);
        }
    });
    return items;
};
/// <reference path="./UIBaseButton.ts" />
class UIButton extends UIBaseButton {
    constructor(elementID, elementType, titleType = UITextView.type.span) {
        super(elementID, elementType, { "titleType": titleType });
        this.usesAutomaticTitleFontSize = NO;
        this.minAutomaticFontSize = nil;
        this.maxAutomaticFontSize = 25;
        this._class = UIButton;
        this.superclass = UIBaseButton;
    }
    initView(elementID, viewHTMLElement, initViewData) {
        this.class.superclass = UIBaseButton;
        // Instance variables
        this.colors = {
            titleLabel: {
                normal: UIColor.whiteColor,
                highlighted: UIColor.whiteColor,
                selected: UIColor.whiteColor
            },
            background: {
                normal: UIColor.blueColor,
                highlighted: UIColor.greenColor,
                selected: UIColor.redColor
            }
        };
        this._imageView = new UIImageView(elementID + "ImageView");
        this._imageView.hidden = YES;
        this.addSubview(this.imageView);
        this.imageView.fillMode = UIImageView.fillMode.aspectFitIfLarger;
        if (IS_NOT_NIL(initViewData.titleType)) {
            this._titleLabel = new UITextView(elementID + "TitleLabel", initViewData.titleType);
            this.titleLabel.style.whiteSpace = "nowrap";
            this.addSubview(this.titleLabel);
            this.titleLabel.userInteractionEnabled = NO;
        }
        this.contentPadding = 10;
        this.imageView.userInteractionEnabled = NO;
        this.titleLabel.textAlignment = UITextView.textAlignment.center;
        this.titleLabel.nativeSelectionEnabled = NO;
    }
    get contentPadding() {
        return this._contentPadding.integerValue;
    }
    set contentPadding(contentPadding) {
        this._contentPadding = contentPadding;
        this.setNeedsLayout();
    }
    set hovered(hovered) {
        this._hovered = hovered;
        this.updateContentForCurrentState();
    }
    get hovered() {
        return this._hovered;
    }
    set highlighted(highlighted) {
        this._highlighted = highlighted;
        this.updateContentForCurrentState();
    }
    get highlighted() {
        return this._highlighted;
    }
    set focused(focused) {
        this._focused = focused;
        if (focused) {
            this.focus();
        }
        else {
            this.blur();
        }
        this.updateContentForCurrentState();
    }
    get focused() {
        return this._focused;
    }
    set selected(selected) {
        this._selected = selected;
        this.updateContentForCurrentState();
    }
    get selected() {
        return this._selected;
    }
    updateContentForCurrentState() {
        var updateFunction = this.updateContentForNormalState;
        if (this.selected && this.highlighted) {
            updateFunction = this.updateContentForSelectedAndHighlightedState;
        }
        else if (this.selected) {
            updateFunction = this.updateContentForSelectedState;
        }
        else if (this.focused) {
            updateFunction = this.updateContentForFocusedState;
        }
        else if (this.highlighted) {
            updateFunction = this.updateContentForHighlightedState;
        }
        else if (this.hovered) {
            updateFunction = this.updateContentForHoveredState;
        }
        if (!IS(updateFunction)) {
            this.titleLabel.textColor = UIColor.nilColor;
            this.backgroundColor = UIColor.nilColor;
        }
        else {
            updateFunction.call(this);
        }
        this.updateContentForCurrentEnabledState();
    }
    updateContentForNormalState() {
        this.backgroundColor = this.colors.background.normal;
        this.titleLabel.textColor = this.colors.titleLabel.normal;
    }
    updateContentForHoveredState() {
        this.updateContentForNormalState();
        if (this.colors.background.hovered) {
            this.backgroundColor = this.colors.background.hovered;
        }
        if (this.colors.titleLabel.hovered) {
            this.titleLabel.textColor = this.colors.titleLabel.hovered;
        }
    }
    updateContentForFocusedState() {
        this.updateContentForHoveredState();
        if (this.colors.background.focused) {
            this.backgroundColor = this.colors.background.focused;
        }
        if (this.colors.titleLabel.focused) {
            this.titleLabel.textColor = this.colors.titleLabel.focused;
        }
    }
    updateContentForHighlightedState() {
        this.backgroundColor = this.colors.background.highlighted;
        this.titleLabel.textColor = this.colors.titleLabel.highlighted;
    }
    updateContentForSelectedState() {
        this.backgroundColor = this.colors.background.selected;
        this.titleLabel.textColor = this.colors.titleLabel.selected;
    }
    updateContentForSelectedAndHighlightedState() {
        this.updateContentForSelectedState();
        if (this.colors.background.selectedAndHighlighted) {
            this.backgroundColor = this.colors.background.selectedAndHighlighted;
        }
        if (this.colors.titleLabel.selectedAndHighlighted) {
            this.titleLabel.textColor = this.colors.titleLabel.selectedAndHighlighted;
        }
    }
    set enabled(enabled) {
        super.enabled = enabled;
        this.updateContentForCurrentState();
    }
    get enabled() {
        return super.enabled;
    }
    updateContentForCurrentEnabledState() {
        if (this.enabled) {
            this.alpha = 1;
        }
        else {
            this.alpha = 0.5;
        }
        this.userInteractionEnabled = this.enabled;
    }
    addStyleClass(styleClassName) {
        super.addStyleClass(styleClassName);
        if (this.styleClassName != styleClassName) {
            this.updateContentForCurrentState.call(this);
        }
    }
    get titleLabel() {
        return this._titleLabel;
    }
    get imageView() {
        return this._imageView;
    }
    layoutSubviews() {
        super.layoutSubviews();
        var bounds = this.bounds;
        this.hoverText = this.titleLabel.text;
        // Image only if text is not present
        if (IS_NOT(this.imageView.hidden) && !IS(this.titleLabel.text)) {
            this.imageView.frame = bounds;
        }
        // Text only if image is not present
        if (IS(this.imageView.hidden) && IS(this.titleLabel.text)) {
            var titleElement = this.titleLabel.viewHTMLElement;
            this.titleLabel.style.left = this.contentPadding;
            this.titleLabel.style.right = this.contentPadding;
            // this.titleLabel.style.marginLeft = ""
            // this.titleLabel.style.right = this.contentPadding
            this.titleLabel.style.top = "50%";
            this.titleLabel.style.transform = "translateY(-50%)";
            this.titleLabel.frame = new UIRectangle(nil, nil, nil, nil);
            if (this.usesAutomaticTitleFontSize) {
                var hidden = this.titleLabel.hidden;
                this.titleLabel.hidden = YES;
                this.titleLabel.fontSize = 15;
                this.titleLabel.fontSize = UITextView.automaticallyCalculatedFontSize(new UIRectangle(0, 0, this.bounds.height, 1 *
                    this.titleLabel.viewHTMLElement.offsetWidth), this.titleLabel.intrinsicContentSize(), this.titleLabel.fontSize, this.minAutomaticFontSize, this.maxAutomaticFontSize);
                this.titleLabel.hidden = hidden;
            }
        }
        // Image and text both present
        if (IS_NOT(this.imageView.hidden) && IS(this.titleLabel.text)) {
            const imageShareOfWidth = 0.25;
            bounds = bounds.rectangleWithInset(this.contentPadding);
            const imageFrame = bounds.copy();
            imageFrame.width = bounds.height - this.contentPadding * 0.5;
            this.imageView.frame = imageFrame;
            var titleElement = this.titleLabel.viewHTMLElement;
            this.titleLabel.style.left = imageFrame.max.x + this.contentPadding;
            this.titleLabel.style.right = this.contentPadding;
            this.titleLabel.style.top = "50%";
            this.titleLabel.style.transform = "translateY(-50%)";
            if (this.usesAutomaticTitleFontSize) {
                var hidden = this.titleLabel.hidden;
                this.titleLabel.hidden = YES;
                this.titleLabel.fontSize = 15;
                this.titleLabel.fontSize = UITextView.automaticallyCalculatedFontSize(new UIRectangle(0, 0, this.bounds.height, 1 *
                    this.titleLabel.viewHTMLElement.offsetWidth), this.titleLabel.intrinsicContentSize(), this.titleLabel.fontSize, this.minAutomaticFontSize, this.maxAutomaticFontSize);
                this.titleLabel.hidden = hidden;
            }
        }
        this.applyClassesAndStyles();
    }
    initViewStyleSelectors() {
        this.initStyleSelector("." + this.styleClassName, "background-color: lightblue;");
        // var selectorWithoutImage = "." + this.styleClassName + " ." + this.imageView.styleClassName + " + ." + this.titleLabel.styleClassName;
        // this.initStyleSelector(
        //     selectorWithoutImage,
        //     "left: " + this.contentPadding + ";" +
        //     "right: " + this.contentPadding + ";" +
        //     "top: 50%;" +
        //     "transform: translateY(-50%);");
    }
}
/// <reference path="../UICore/UIButton.ts" />
class CBButton extends UIButton {
    constructor(elementID, elementType) {
        super(elementID, elementType);
        this._class = CBButton;
        this.superclass = UIButton;
    }
    initView(elementID, viewHTMLElement, initViewData) {
        super.initView(elementID, viewHTMLElement, initViewData);
        //this.style.borderRadius = "2px"
        this.style.outline = "none";
        this.colors.titleLabel.normal = UIColor.whiteColor;
        this.setBackgroundColorsWithNormalColor(UIColor.blueColor);
        this.colors.titleLabel.selected = UIColor.blueColor;
    }
    setBackgroundColorsWithNormalColor(normalBackgroundColor) {
        this.colors.background.normal = normalBackgroundColor;
        this.colors.background.hovered = UIColor.colorWithRGBA(40, 168, 183); // normalBackgroundColor.colorByMultiplyingRGB(0.85)
        this.colors.background.focused = normalBackgroundColor; // normalBackgroundColor.colorByMultiplyingRGB(0.5)
        this.colors.background.highlighted = UIColor.colorWithRGBA(48, 196, 212); // normalBackgroundColor.colorByMultiplyingRGB(0.7)
        this.colors.background.selected = UIColor.whiteColor;
        this.updateContentForCurrentState();
    }
    updateContentForNormalState() {
        super.updateContentForNormalState();
        this.setBorder(0, 0);
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.24)"
    }
    updateContentForHoveredState() {
        super.updateContentForHoveredState();
        this.setBorder(0, 0);
        //this.titleLabel.textColor = UIColor.whiteColor.colorByMultiplyingRGB(0.85);
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.18)"
    }
    updateContentForFocusedState() {
        super.updateContentForFocusedState();
        this.setBorder(0, 1, UIColor.blueColor);
        //this.titleLabel.textColor = UIColor.whiteColor.colorByMultiplyingRGB(0.85);
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.18)"
    }
    updateContentForHighlightedState() {
        super.updateContentForHighlightedState();
        this.setBorder(0, 0);
        //this.titleLabel.textColor = UIColor.whiteColor.colorByMultiplyingRGB(0.7);
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.12)"
    }
    updateContentForCurrentEnabledState() {
        super.updateContentForCurrentEnabledState();
        if (IS_NOT(this.enabled)) {
            this.titleLabel.textColor = new UIColor("#adadad");
            this.backgroundColor = new UIColor("#e5e5e5");
            this.alpha = 1;
        }
    }
}
///<reference path="../UICore/UIButton.ts"/>
class CBCheckbox extends UIView {
    constructor(elementID) {
        super(elementID, nil, "label");
        this._class = CBCheckbox;
        this.superclass = UIView;
    }
    initView(elementID, viewHTMLElement) {
        super.initView(elementID, viewHTMLElement);
        this.checkbox = new UIView(elementID + "Checkbox", nil, "input");
        this.checkbox.viewHTMLElement.setAttribute("type", "checkbox");
        this.addSubview(this.checkbox);
        this.titleLabel = new UITextView(elementID + "TitleLabel", UITextView.type.span);
        this.addSubview(this.titleLabel);
        this.checkbox.userInteractionEnabled = NO;
        this.titleLabel.userInteractionEnabled = NO;
        this.titleLabel.style.overflow = "visible";
        this.titleLabel.style.lineHeight = "1.5";
        this.style.cursor = "pointer";
        this.viewHTMLElement.onchange = (event) => {
            this.sendControlEventForKey(CBCheckbox.controlEvent.SelectionChange, event);
        };
    }
    get isStatic() {
        const result = this.checkbox.viewHTMLElement.classList.contains("staticCheckbox");
        return result;
    }
    set isStatic(isStatic) {
        if (isStatic) {
            this.checkbox.viewHTMLElement.classList.add("staticCheckbox");
        }
        else {
            this.checkbox.viewHTMLElement.classList.remove("staticCheckbox");
        }
    }
    get selected() {
        return this.checkbox.viewHTMLElement.checked;
    }
    set selected(selected) {
        this.checkbox.viewHTMLElement.checked = selected;
    }
    layoutSubviews() {
        super.layoutSubviews();
        const titleWidth = this.bounds.width - 35;
        this.titleLabel.setMaxSizes(nil, titleWidth);
        this.titleLabel.textPrefix = "<span style=\"position: absolute; overflow: hidden; left: 0; top: 0;text-overflow: ellipsis; white-space: pre; padding-left: 35px; width: " +
            titleWidth + "px;\">";
        this.titleLabel.textSuffix = "</span>";
        this.hoverText = this.titleLabel.text.replace(this.titleLabel.textPrefix, "").replace(this.titleLabel.textSuffix, "");
    }
}
CBCheckbox.controlEvent = Object.assign({}, UIView.controlEvent, {
    "SelectionChange": "SelectionChange"
});
/// <reference path="../UICore/UIColor.ts" />
class CBColor extends UIColor {
    constructor(stringValue) {
        super(stringValue);
        this._class = CBColor;
        this.superclass = UIColor;
    }
    static get primaryTintColor() {
        return new CBColor("rgba(52, 127, 230, 1)");
    }
    static get secondaryTintColor() {
        return new CBColor("rgba(0, 196, 212, 1)");
    }
    static get primaryContentColor() {
        return new CBColor("rgb(35, 35, 35)");
    }
    static get greenTintColor() {
        return new CBColor("rgba(51, 188, 125, 1)");
    }
    static get redTintColor() {
        return new CBColor("rgba(236, 88, 111, 1)");
    }
}
/// <reference path="./CBButton.ts" />
class CBFlatButton extends CBButton {
    constructor(elementID, elementType) {
        super(elementID, elementType);
        this._class = CBFlatButton;
        this.superclass = CBButton;
    }
    initView(elementID, viewHTMLElement, initViewData) {
        super.initView(elementID, viewHTMLElement, initViewData);
        this.colors = {
            titleLabel: {
                normal: UIColor.blueColor,
                highlighted: UIColor.blueColor,
                selected: UIColor.blueColor
            },
            background: {
                normal: UIColor.transparentColor,
                hovered: new UIColor("#F8F8F8"),
                highlighted: new UIColor("#ebebeb"),
                selected: new UIColor("#ebebeb")
            }
        };
    }
    set titleLabelColor(titleLabelColor) {
        this.colors.titleLabel.normal = titleLabelColor;
        this.colors.titleLabel.highlighted = titleLabelColor;
        this.colors.titleLabel.selected = titleLabelColor;
        this.updateContentForCurrentState();
    }
    get titleLabelColor() {
        const result = this.colors.titleLabel.normal;
        return result;
    }
    updateContentForNormalState() {
        UIButton.prototype.updateContentForNormalState.call(this);
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.24)"
    }
    updateContentForHoveredState() {
        UIButton.prototype.updateContentForHoveredState.call(this);
        //this.titleLabel.textColor = UIColor.whiteColor.colorByMultiplyingRGB(0.85);
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.18)";
    }
    updateContentForFocusedState() {
        UIButton.prototype.updateContentForFocusedState.call(this);
        //this.titleLabel.textColor = UIColor.whiteColor.colorByMultiplyingRGB(0.85);
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.18)"
    }
    updateContentForHighlightedState() {
        UIButton.prototype.updateContentForHighlightedState.call(this);
        //this.titleLabel.textColor = UIColor.whiteColor.colorByMultiplyingRGB(0.7);
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.12)"
    }
}
CBFlatButton.colors = {
    titleLabel: {
        normal: UIColor.blueColor,
        highlighted: UIColor.blueColor,
        selected: UIColor.blueColor
    },
    background: {
        normal: UIColor.transparentColor,
        hovered: new UIColor("#F8F8F8"),
        highlighted: new UIColor("#ebebeb"),
        selected: new UIColor("#ebebeb")
    }
};
/// <reference path="../UICore/UIView.ts" />
class RowView extends UIView {
    constructor(elementID, cells = [], cellWidths = []) {
        super(elementID);
        this.padding = 0;
        this._rowHeight = 50;
        this._class = RowView;
        this.superclass = UIView;
        this._cells = cells;
        this._cellWeights = cellWidths;
    }
    get cells() {
        return this._cells;
    }
    set cells(cells) {
        const previousCells = this.cells;
        // The cells are added to this._cells in this.addCell()
        previousCells.copy().forEach(function (cell, index, array) {
            if (!cells.contains(cell)) {
                cell.removeFromSuperview();
                this._cells.removeElement(cell);
            }
        }.bind(this));
        cells.copy().forEach(function (cell, index, array) {
            if (!IS(cell.superview)) {
                this.addCell(cell, 1, index);
            }
        }, this);
        this._previousLayoutBounds = nil;
        this.setNeedsLayout();
    }
    removeCellAtIndex(index) {
        const remainingCells = this.cells.copy();
        remainingCells.removeElementAtIndex(index);
        this.cellWeights.removeElementAtIndex(index);
        this.cells = remainingCells;
    }
    removeFirstCell() {
        this.removeCellAtIndex(0);
    }
    removeLastCell() {
        this.removeCellAtIndex(this.cells.length - 1);
    }
    addCell(cell, weight = 1, index = this.cells.length) {
        if (this.cells.contains(cell) && IS_NOT_NIL(cell)) {
            return;
        }
        this.cells.insertElementAtIndex(index, cell);
        this.cellWeights.insertElementAtIndex(index, weight);
        this.addSubview(cell);
        this.setNeedsLayout();
    }
    get cellWeights() {
        return this._cellWeights;
    }
    set cellWeights(widths) {
        this._cellWeights = widths;
        this._previousLayoutBounds = nil;
        this.setNeedsLayout();
    }
    get cellWidths() {
        return this._cellWidths;
    }
    set cellWidths(widths) {
        this._cellWidths = widths;
        this._previousLayoutBounds = nil;
        this.setNeedsLayout();
    }
    get rowHeight() {
        var result = IF(this._rowHeight)(() => this._rowHeight)
            .ELSE(() => this.cells.map((value, index, array) => value.intrinsicContentHeight(this.bounds.width)).max());
        return result; //this.bounds.height
    }
    layoutSubviews() {
        const bounds = this.bounds;
        if (bounds.isEqualTo(this._previousLayoutBounds)) {
            return;
        }
        super.layoutSubviews();
        this._previousLayoutBounds = bounds;
        bounds.distributeViewsAlongWidth(this._cells, this._cellWeights, this.padding, this._cellWidths);
        this.cells.forEach(function (cell, index, array) {
            cell.frame = cell.frame.rectangleWithHeight(this.rowHeight);
        }.bind(this));
    }
}
/// <reference path="./CBFlatButton.ts" />
/// <reference path="./RowView.ts" />
/// <reference path="../UICore/UICoreExtensions.ts" />
class CBDialogView extends UIView {
    constructor(elementID) {
        super(elementID);
        this.titleLabel = nil;
        this.titleRow = nil;
        this.placeholderLabel = nil;
        this._view = nil;
        this.questionLabel = nil;
        this.yesButton = nil;
        this.noButton = nil;
        this.cancelButton = nil;
        this.noButtonDismissesDialog = YES;
        this._class = CBDialogView;
        this.superclass = UIView;
        this.updateContent();
    }
    initView(elementID, viewHTMLElement) {
        super.initView(elementID, viewHTMLElement);
        // this.initPlaceholderLabelIfNeeded()
        // this.initTitleLabelIfNeeded()
        // this.initTitleRow()
        // this.titleRow.cellWeights = [1, 0.75, 1.2, 1]
        // var cells = this.titleRow.cells as CellView[];
        // var cellTitles = ["TIME:", "COMPANY:", "COMMENTS:", "PRICE:"]
        // cells.forEach(function (cell, index, array) {
        //     cell.titleLabel.text = cellTitles[index];
        //     cell.titleLabel.useSairaFont()
        // })
        // this.initQuestionLabelIfNeeded()
        // this.initYesNoButtonsIfNeeded()
        this.style.borderRadius = "5px";
    }
    set view(view) {
        this.view.removeFromSuperview();
        this._view = view;
        this.addSubview(view);
    }
    get view() {
        return this._view;
    }
    initCancelButtonIfNeeded() {
        if (IS_NOT(this.cancelButton)) {
            this.cancelButton = new CBFlatButton(this.elementID + "CancelButton");
            this.cancelButton.titleLabel.text = "OK";
            this.cancelButton.titleLabelColor = UIColor.greenColor;
            this.cancelButton.titleLabel.style.fontWeight = "600";
            this.addSubview(this.cancelButton);
            this.cancelButton.addTargetForControlEvents([
                UIButton.controlEvent.PointerUpInside, UIButton.controlEvent.EnterDown
            ], function (sender, event) {
                this.dialogView.dismiss();
            }.bind(this));
        }
    }
    initYesNoButtonsIfNeeded() {
        if (IS_NOT(this.yesButton)) {
            this.yesButton = new CBFlatButton(this.elementID + "YesButton");
            this.yesButton.titleLabel.setText("cBDDialogViewYES", "YES");
            this.yesButton.titleLabelColor = UIColor.greenColor;
            this.yesButton.titleLabel.style.fontWeight = "600";
            this.addSubview(this.yesButton);
            this.noButton = new CBFlatButton(this.elementID + "NoButton");
            this.noButton.titleLabel.setText("cBDDialogViewNO", "NO");
            this.noButton.titleLabelColor = UIColor.redColor;
            this.noButton.titleLabel.style.fontWeight = "600";
            this.addSubview(this.noButton);
            // this.yesButton.addTargetForControlEvents([UIButton.controlEvent.PointerUpInside, UIButton.controlEvent.EnterDown], function (this: CBDialogView, sender: UIButton, event: Event) {
            //     this.dialogView.dismiss()
            // }.bind(this))
            this.noButton.addTargetForControlEvents([
                UIButton.controlEvent.PointerUpInside, UIButton.controlEvent.EnterDown
            ], function (sender, event) {
                if (this.noButtonDismissesDialog) {
                    this.dialogView.dismiss();
                }
            }.bind(this));
        }
    }
    initQuestionLabelIfNeeded() {
        if (IS_NOT(this.questionLabel)) {
            this.questionLabel = new UITextView(this.elementID + "QuestionLabel", UITextView.type.header3);
            this.questionLabel.text = LanguageService.stringForKey("cBDDialogViewAreYouSure", LanguageService.currentLanguageKey, "Are you sure?");
            this.questionLabel.textAlignment = UITextView.textAlignment.center;
            this.addSubview(this.questionLabel);
        }
    }
    initTitleRow() {
        if (IS_NOT(this.titleRow)) {
            this.titleRow = new RowView(this.elementID + "TitleRow");
            this.addSubview(this.titleRow);
        }
    }
    initTitleLabelIfNeeded() {
        if (IS_NOT(this.titleLabel)) {
            this.titleLabel = new UITextView(this.elementID + "TitleLabel", UITextView.type.header4);
            this.titleLabel.text = "Accept offer";
            this.titleLabel.style.fontSize = "24";
            this.addSubview(this.titleLabel);
        }
    }
    initPlaceholderLabelIfNeeded() {
        if (IS_NOT(this.placeholderLabel)) {
            this.placeholderLabel = new UITextView(this.elementID + "PlaceholderLabel", UITextView.type.header4);
            this.placeholderLabel.text = "No offers have been made yet.";
            this.addSubview(this.placeholderLabel);
        }
    }
    updateContent() {
        function stringFromValue(value) {
            if (IS(value)) {
                return value;
            }
            return "-";
        }
        this.layoutSubviews();
    }
    layoutSubviews() {
        super.layoutSubviews();
        const padding = RootViewController.paddingLength;
        const labelHeight = padding * 0.75;
        const bounds = this.bounds.rectangleWithInsets(padding * 0.5, padding * 0.5, 0, 0);
        const topObject = { frame: bounds.rectangleWithHeight(0) };
        this.titleLabel.frame = bounds.rectangleWithHeight(this.titleLabel.intrinsicContentHeight(bounds.width)).rectangleWithY(bounds.y +
            padding * 0.5);
        this.titleLabel.style.marginLeft = "" + (padding * 0.5) + "px";
        this.titleRow.frame = FIRST_OR_NIL(this.titleLabel, topObject).frame.rectangleForNextRow(0, this.titleRow.intrinsicContentHeight(bounds.width));
        this.view.frame = FIRST_OR_NIL(this.titleRow, this.titleLabel, topObject).frame.rectangleForNextRow(padding, this.view.intrinsicContentHeight(bounds.width)).rectangleWithWidth(bounds.width).rectangleWithX(bounds.x);
        this.questionLabel.frame = FIRST_OR_NIL(this.view, this.titleRow, this.titleLabel, topObject).frame.rectangleForNextRow(padding, this.questionLabel.intrinsicContentHeight(bounds.width)).rectangleWithWidth(bounds.width).rectangleWithX(bounds.x);
        const buttons = [this.yesButton, this.noButton];
        const buttonsFrame = FIRST_OR_NIL(this.questionLabel, this.view, this.titleRow, this.titleLabel, topObject).frame.rectangleForNextRow(padding, labelHeight *
            2)
            .rectangleWithWidth([buttons.everyElement.titleLabel.intrinsicContentWidth().max(), 120].max() * 2 + padding * 3, 0.5);
        buttonsFrame.distributeViewsAlongWidth(buttons);
        this.cancelButton.frame = buttonsFrame.rectangleWithWidth(150, 0.5);
        this.placeholderLabel.frame = bounds.rectangleWithInsets(padding * 0.5, padding * 0.5, 0, 0);
        if (IS(this.view)) {
            this.placeholderLabel.hidden = YES;
            // this.titleLabel.hidden = NO
            // this.titleRow.hidden = NO
        }
        else {
            this.placeholderLabel.hidden = NO;
            // this.titleLabel.hidden = YES
            // this.titleRow.hidden = YES
        }
    }
    intrinsicContentHeight(constrainingWidth = 0) {
        const padding = RootViewController.paddingLength;
        const labelHeight = padding * 0.75;
        function addPaddingIfNeeded(value) {
            return IF(value)(function () {
                return value + padding;
            })();
        }
        var result = padding;
        result = result + addPaddingIfNeeded(this.titleLabel.intrinsicContentHeight(constrainingWidth));
        result = result + addPaddingIfNeeded(this.titleRow.intrinsicContentHeight(constrainingWidth));
        result = result + addPaddingIfNeeded(this.view.intrinsicContentHeight(constrainingWidth) + padding);
        result = result + addPaddingIfNeeded(this.questionLabel.intrinsicContentHeight(constrainingWidth));
        result = result + IF(FIRST_OR_NIL(this.yesButton, this.cancelButton))(RETURNER(labelHeight))();
        return result;
    }
}
/// <reference path="../UICore/UICore.ts" />
/// <reference path="../UICore/UIRoute.ts" />
class LanguageService {
    static useStoredLanguageValues(values = {}) {
        LanguageService.languages = JSON.parse(JSON.stringify(LanguageService.languageValues)).objectByCopyingValuesRecursivelyFromObject(values);
        LanguageService.broadcastLanguageChangeEvent();
    }
    static broadcastLanguageChangeEvent() {
        UICore.main.rootViewController.view.broadcastEventInSubtree({
            name: UIView.broadcastEventName.LanguageChanged,
            parameters: {}
        });
    }
    static get defaultLanguageKey() {
        return "en";
    }
    static get currentLanguageKey() {
        if (!LanguageService._currentLanguageKey) {
            LanguageService.updateCurrentLanguageKey();
        }
        return LanguageService._currentLanguageKey;
    }
    static updateCurrentLanguageKey(route = UIRoute.currentRoute) {
        var result = route.componentWithName("settings").parameters.language;
        if (IS_NOT(result)) {
            result = LanguageService.defaultLanguageKey;
        }
        LanguageService._currentLanguageKey = result;
    }
    get currentLanguageKey() {
        const result = LanguageService.currentLanguageKey;
        return result;
    }
    static stringForKey(key, languageKey, defaultString, parameters) {
        var result;
        if (IS(key) && LanguageService.languages[languageKey] &&
            IS_DEFINED(LanguageService.languages[languageKey][key])) {
            result = LanguageService.languages[languageKey][key];
        }
        else {
            result = defaultString;
        }
        if (IS(parameters)) {
            const parameterKeys = Object.keys(parameters);
            parameterKeys.forEach(function (key, index, array) {
                const keyString = "%" + key + "%";
                const parameter = parameters[key];
                var parameterString;
                if (parameter instanceof Object) {
                    parameterString = UICore.languageService.stringForCurrentLanguage(parameter);
                }
                else {
                    parameterString = parameter;
                }
                result = result.replace(new RegExp(keyString, "g"), parameterString);
            });
        }
        return result;
    }
    stringForKey(key, languageKey, defaultString, parameters) {
        return LanguageService.stringForKey(key, languageKey, defaultString, parameters);
    }
    static localizedTextObjectForKey(key, defaultString = key, parameters) {
        const result = {};
        LanguageService.languages.forEach(function (languageObject, languageKey) {
            result[languageKey] = LanguageService.stringForKey(key, languageKey, defaultString, parameters);
        });
        return result;
    }
    localizedTextObjectForKey(key, defaultString, parameters) {
        const result = LanguageService.localizedTextObjectForKey(key, defaultString, parameters);
        return result;
    }
    static localizedTextObjectForText(text) {
        if (IS_NOT(text)) {
            return nil;
        }
        const result = {
            [LanguageService.defaultLanguageKey]: text
        };
        return result;
    }
    localizedTextObjectForText(text) {
        const result = LanguageService.localizedTextObjectForText(text);
        return result;
    }
    static stringForCurrentLanguage(localizedTextObject) {
        if (!LanguageService || !localizedTextObject) {
            const asd = 1;
        }
        localizedTextObject = FIRST_OR_NIL(localizedTextObject);
        var result = localizedTextObject[LanguageService.currentLanguageKey];
        if (IS_NOT(result)) {
            result = localizedTextObject[LanguageService.defaultLanguageKey];
        }
        if (IS_NOT(result)) {
            result = localizedTextObject["en"];
        }
        if (IS_NOT(result)) {
            result = "";
        }
        return result;
    }
    stringForCurrentLanguage(localizedTextObject) {
        return LanguageService.stringForCurrentLanguage(localizedTextObject);
    }
}
LanguageService.languageValues = {
    "en": {
        "languageName": "English",
        "languageNameShort": "ENG",
        "topBarTitle": "Carrum",
        "frontPageHTML": "Asdasd",
        "bottomBarContent": "PUT BOTTOM BAR STUFF IN HERE"
    },
    "est": {
        "languageName": "Eesti keel",
        "languageNameShort": "EST",
        "topBarTitle": "Carrum"
    }
};
LanguageService.languages = JSON.parse(JSON.stringify(LanguageService.languageValues));
UICore.languageService = LanguageService;
/// <reference path="../UICore/UICoreExtensions.ts" />
/// <reference path="../UICore/UICore.ts" />
/// <reference path="../UICore/UIColor.ts" />
/// <reference path="./CBDialogView.ts" />
/// <reference path="../Custom components/LanguageService.ts" />
class CBDialogViewShower extends UIObject {
    constructor(elementID) {
        super();
        this._class = CBDialogViewShower;
        this.superclass = UIObject;
        this.dialogView = new UIDialogView();
        this.dialogView.view = new CBDialogView(elementID);
        this.dialogView.view.backgroundColor = UIColor.whiteColor;
        const dialogLayoutFunction = this.dialogView.layoutSubviews.bind(this.dialogView);
        this.dialogView.layoutSubviews = function () {
            dialogLayoutFunction();
            this.dialogView.view.frame = new UIRectangle(0, 0, window.innerHeight, window.innerWidth).rectangleWithHeight(this.dialogView.view.intrinsicContentHeight(this.getDialogWidth()), 0.5).rectangleWithWidth(this.getDialogWidth(), 0.5);
            this.dialogView.frame = UICore.main.rootViewController.view.bounds;
        }.bind(this);
        this.dialogView.view.yesButton.addTargetForControlEvents([
            UIButton.controlEvent.PointerUpInside, UIButton.controlEvent.EnterDown
        ], function (sender, event) {
            this.yesButtonWasPressed();
        }.bind(this));
        this.dialogView.view.noButton.addTargetForControlEvents([
            UIButton.controlEvent.PointerUpInside, UIButton.controlEvent.EnterDown
        ], function (sender, event) {
            this.noButtonWasPressed();
        }.bind(this));
    }
    getDialogWidth() {
        const padding = RootViewController.paddingLength;
        const labelHeight = padding * 0.75;
        var result = 250;
        const width = this.dialogView.view.titleLabel.intrinsicContentWidth() + padding * 2;
        result = Math.max(result, this.dialogView.view.view.intrinsicContentWidth(this.dialogView.view.view.viewHTMLElement.naturalHeight || 1000000000));
        result = Math.max(result, width);
        result = Math.min(result, 1000);
        const dialogMaxWidth = (this.dialogView.superview ||
            { "bounds": new UIRectangle(0, 0, 0, result) }).bounds.width;
        result = Math.min(result, dialogMaxWidth);
        return result;
    }
    yesButtonWasPressed() {
    }
    noButtonWasPressed() {
    }
    cancelButtonWasPressed() {
    }
    showQuestionDialogInRootView(titleTextObject, questionTextObject) {
        this.dialogView.view.initTitleLabelIfNeeded();
        this.dialogView.view.titleLabel.localizedTextObject = titleTextObject;
        this.dialogView.view.initQuestionLabelIfNeeded();
        if (IS(questionTextObject)) {
            this.dialogView.view.questionLabel.localizedTextObject = questionTextObject;
        }
        this.dialogView.view.initYesNoButtonsIfNeeded();
        this.dialogView.view.noButton.addTargetForControlEvents([
            UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
        ], function (sender, event) {
            this.noButtonWasPressed();
        }.bind(this));
        this.dialogView.view.yesButton.addTargetForControlEvents([
            UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
        ], function (sender, event) {
            this.yesButtonWasPressed();
        }.bind(this));
        this.dialogView.showInRootView(YES);
    }
    showMessageDialogInRootView(titleTextObject) {
        this.dialogView.view.initTitleLabelIfNeeded();
        this.dialogView.view.titleLabel.localizedTextObject = titleTextObject;
        this.dialogView.view.initCancelButtonIfNeeded();
        this.dialogView.view.cancelButton.addTargetForControlEvents([
            UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
        ], function (sender, event) {
            this.cancelButtonWasPressed();
        }.bind(this));
        this.dialogView.showInRootView(YES);
        this.dialogView.view.cancelButton.focus();
    }
    showDialogInRootView(view) {
        this.dialogView.view.view = view;
        this.dialogView.view.initCancelButtonIfNeeded();
        this.dialogView.view.cancelButton.addTargetForControlEvents([
            UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
        ], function (sender, event) {
            this.cancelButtonWasPressed();
        }.bind(this));
        this.dialogView.showInRootView(YES);
        this.dialogView.view.cancelButton.focus();
    }
    showImageDialogInRootView(imageURL, deleteImageCallback) {
        var loadingLabel = new UITextView();
        loadingLabel.text = "Loading image.";
        loadingLabel.textAlignment = UITextView.textAlignment.center;
        this.dialogView.view.view = loadingLabel;
        const imageView = new UIImageView();
        imageView.imageSource = imageURL;
        imageView.viewHTMLElement.onload = (event) => {
            this.dialogView.view.view = imageView;
            imageView.setNeedsLayoutUpToRootView();
        };
        imageView.fillMode = UIImageView.fillMode.aspectFitIfLarger;
        if (IS(deleteImageCallback)) {
            this.dialogView.view.initYesNoButtonsIfNeeded();
            this.dialogView.view.yesButton.titleLabel.text = "Close";
            this.dialogView.view.noButton.titleLabel.text = "Delete";
            this.dialogView.view.noButtonDismissesDialog = NO;
            this.dialogView.view.noButton.addTargetForControlEvents([
                UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
            ], function (sender, event) {
                const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(function () {
                    //this.dialogView.dismiss()
                }.bind(this));
                var textObject = LanguageService.localizedTextObjectForText("Delete this image.");
                dialogShower.showQuestionDialogInRootView(textObject);
                dialogShower.yesButtonWasPressed = function () {
                    deleteImageCallback();
                    dialogShower.dialogView.dismiss();
                };
            }.bind(this));
            this.dialogView.view.yesButton.addTargetForControlEvents([
                UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
            ], function (sender, event) {
                this.dialogView.dismiss();
            }.bind(this));
        }
        else {
            this.dialogView.view.initCancelButtonIfNeeded();
            this.dialogView.view.cancelButton.titleLabel.text = "Close";
        }
        this.dialogView.showInRootView(YES);
        this.dialogView.view.cancelButton.focus();
    }
    showActionIndicatorDialogInRootView(message) {
        const actionIndicator = new UIActionIndicator();
        this.dialogView.zIndex = 150;
        this.dialogView.view.view = actionIndicator;
        actionIndicator.style.minHeight = "100px";
        this.dialogView.view.initQuestionLabelIfNeeded();
        this.dialogView.view.questionLabel.text = message;
        actionIndicator.start();
        this.dialogView.view.backgroundColor = CBColor.primaryContentColor.colorWithAlpha(0.5);
        this.dialogView.view.questionLabel.textColor = UIColor.whiteColor;
        this.dialogView.dismissesOnTapOutside = NO;
        CBDialogViewShower.currentActionIndicatorDialogViewShower = this;
        this.dialogView.showInRootView(NO);
        this.dialogView.view.cancelButton.focus();
    }
    static showNextDialog() {
        (CBDialogViewShower.nextShowDialogFunctions.firstElement || nil)();
        CBDialogViewShower.nextShowDialogFunctions.removeElementAtIndex(0);
    }
    static alert(text, dismissCallback = nil) {
        const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback);
        const textObject = LanguageService.localizedTextObjectForText(text);
        const showDialogFunction = dialogShower.showMessageDialogInRootView.bind(dialogShower, textObject);
        CBDialogViewShower._showDialogWithFunction(showDialogFunction, dialogShower);
    }
    static localizedAlert(textObject, dismissCallback = nil) {
        const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback);
        const showDialogFunction = dialogShower.showMessageDialogInRootView.bind(dialogShower, textObject);
        CBDialogViewShower._showDialogWithFunction(showDialogFunction, dialogShower);
    }
    static showQuestionDialog(questionText, dismissCallback = nil) {
        const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback);
        const textObject = LanguageService.localizedTextObjectForText(questionText);
        const showDialogFunction = dialogShower.showQuestionDialogInRootView.bind(dialogShower, textObject);
        CBDialogViewShower._showDialogWithFunction(showDialogFunction, dialogShower);
        return dialogShower;
    }
    static showImageDialog(imageURL, deleteImageCallback = nil, dismissCallback = nil) {
        var dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback);
        var showDialogFunction = dialogShower.showImageDialogInRootView.bind(dialogShower, imageURL, deleteImageCallback);
        CBDialogViewShower._showDialogWithFunction(showDialogFunction, dialogShower);
        return dialogShower;
    }
    static showDialog(view, dismissCallback = nil) {
        const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback);
        const showDialogFunction = dialogShower.showDialogInRootView.bind(dialogShower, view);
        CBDialogViewShower._showDialogWithFunction(showDialogFunction, dialogShower);
        return dialogShower;
    }
    static showActionIndicatorDialog(message, dismissCallback = nil) {
        if (IS(CBDialogViewShower.currentActionIndicatorDialogViewShower)) {
            CBDialogViewShower.currentActionIndicatorDialogViewShower.dialogView.view.questionLabel.text = message;
            CBDialogViewShower.currentActionIndicatorDialogViewShower.dialogView.view.setNeedsLayoutUpToRootView();
            return;
        }
        const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback);
        dialogShower.showActionIndicatorDialogInRootView(message);
        return dialogShower;
    }
    static hideActionIndicatorDialog() {
        CBDialogViewShower.currentActionIndicatorDialogViewShower.dialogView.dismiss();
        CBDialogViewShower.currentActionIndicatorDialogViewShower = nil;
    }
    static _dialogShowerWithDismissCallback(dismissCallback) {
        const dialogShower = new CBDialogViewShower();
        const dismissFunction = dialogShower.dialogView.dismiss.bind(dialogShower.dialogView);
        dialogShower.dialogView.dismiss = function () {
            dismissFunction();
            dismissCallback();
            CBDialogViewShower.currentDialogViewShower = null;
            CBDialogViewShower.showNextDialog();
        };
        return dialogShower;
    }
    static _showDialogWithFunction(showDialogFunction, dialogShower) {
        if (IS(CBDialogViewShower.currentDialogViewShower)) {
            CBDialogViewShower.nextShowDialogFunctions.push(showDialogFunction);
        }
        else {
            CBDialogViewShower.currentDialogViewShower = dialogShower;
            showDialogFunction();
        }
    }
}
CBDialogViewShower.nextShowDialogFunctions = [];
CBDialogViewShower.currentActionIndicatorDialogViewShower = nil;
/// <reference path="../UICore/UILinkButton.ts" />
class CBLinkButton extends UILinkButton {
    constructor(elementID, elementType) {
        super(elementID, elementType);
        this._class = CBLinkButton;
        this.superclass = UILinkButton;
    }
    initView(elementID, viewHTMLElement, initViewData) {
        super.initView(elementID, viewHTMLElement, initViewData);
        this.button.removeFromSuperview();
        this.button = new CBButton(this.elementID + "Button", initViewData.elementType);
        this.addSubview(this.button);
    }
}
/// <reference path="../UICore/UIView.ts" />
class CellView extends UIButton {
    //titleLabel: UITextView;
    constructor(elementID, titleLabelType = UITextView.type.span) {
        super(elementID, undefined, titleLabelType);
        this._isAButton = NO;
        this.leftInset = 0;
        this.rightInset = 0;
        this._class = CellView;
        this.superclass = UIView;
        // this.titleLabel = new UITextView(this.elementID + "TitleLabel", titleLabelType);
        // this.addSubview(this.titleLabel);
        this.updateForCurrentIsAButtonState();
    }
    set isAButton(isAButton) {
        this._isAButton = isAButton;
        this.updateForCurrentIsAButtonState();
    }
    get isAButton() {
        return this._isAButton;
    }
    updateForCurrentIsAButtonState() {
        if (this._isAButton) {
            this.style.cursor = "pointer";
            this.titleLabel.userInteractionEnabled = NO;
            this.titleLabel.nativeSelectionEnabled = NO;
            this.titleLabel.textAlignment = UITextView.textAlignment.center;
            this.nativeSelectionEnabled = NO;
            this.style.outline = "";
            this.colors = {
                titleLabel: {
                    normal: UIColor.blueColor,
                    highlighted: UIColor.blueColor,
                    selected: UIColor.blueColor
                },
                background: {
                    normal: UIColor.transparentColor,
                    hovered: new UIColor("#F8F8F8"),
                    highlighted: new UIColor("#ebebeb"),
                    selected: new UIColor("#ebebeb")
                }
            };
        }
        else {
            this.style.cursor = "";
            this.titleLabel.userInteractionEnabled = YES;
            this.titleLabel.nativeSelectionEnabled = YES;
            this.titleLabel.textAlignment = UITextView.textAlignment.left;
            this.nativeSelectionEnabled = YES;
            this.style.outline = "none";
            this.colors = {
                titleLabel: {
                    normal: UIColor.blackColor,
                    highlighted: UIColor.blackColor,
                    selected: UIColor.blackColor
                },
                background: {
                    normal: UIColor.transparentColor,
                    highlighted: UIColor.transparentColor,
                    selected: UIColor.transparentColor
                }
            };
        }
        this.updateContentForCurrentState();
    }
    initRightImageViewIfNeeded() {
        if (this._rightImageView) {
            return;
        }
        this._rightImageView = new UIImageView(this.elementID + "RightImageView");
        this._rightImageView.userInteractionEnabled = NO;
    }
    set rightImageSource(imageSource) {
        if (IS(imageSource)) {
            this.initRightImageViewIfNeeded();
            this._rightImageView.imageSource = imageSource;
            this.addSubview(this._rightImageView);
        }
        else {
            this._rightImageView.removeFromSuperview();
        }
    }
    get rightImageSource() {
        var result = nil;
        if (this._rightImageView) {
            result = this._rightImageView.imageSource;
        }
        return result;
    }
    layoutSubviews() {
        super.layoutSubviews();
        const padding = RootViewController.paddingLength;
        const labelHeight = padding;
        const bounds = this.bounds;
        this.titleLabel.centerYInContainer();
        this.titleLabel.style.left = "" + (padding * 0.5 + this.leftInset).integerValue + "px";
        this.titleLabel.style.right = "" + (padding * 0.5 + this.rightInset).integerValue + "px";
        this.titleLabel.style.maxHeight = "100%";
        this.titleLabel.style.overflow = "hidden";
        //this.titleLabel.style.whiteSpace = "nowrap";
        if (this._rightImageView && this._rightImageView.superview == this) {
            // var imageHeight = bounds.height - padding;
            // this._rightImageView.frame = new UIRectangle(bounds.width - imageHeight - padding * 0.5, padding * 0.5, imageHeight, imageHeight);
            this._rightImageView.frame = bounds.rectangleWithInsets(this.leftInset, padding * 0.5 +
                this.rightInset, 0, 0).rectangleWithWidth(24, 1).rectangleWithHeight(24, 0.5);
            this.titleLabel.style.right = "" +
                (padding * 0.5 + this.rightInset + this._rightImageView.frame.width).integerValue + "px";
        }
    }
}
///<reference path="../UICore/UIView.ts"/>
class SearchTextField extends UIView {
    constructor(elementID) {
        super(elementID);
        this._class = SearchTextField;
        this.superclass = UIView;
    }
    initView(elementID, viewHTMLElement) {
        super.initView(elementID, viewHTMLElement);
        viewHTMLElement.classList.add("input", "input--search");
        viewHTMLElement.innerHTML = "\
            <button type=\"button\" class=\"input__button\">\
                <i class=\"material-icons\">search</i>\
            </button>\
            <input type=\"search\" class=\"input__field\" placeholder=\"Search\">";
        this._textField = new UIView(nil, this.textFieldElement);
        this._searchButton = new UIView(nil, this.searchButtonElement);
        this._searchButton.addTargetForControlEvents([
            UIView.controlEvent.EnterDown, UIView.controlEvent.PointerUpInside
        ], function (sender, event) {
            this.performSearch();
        }.bind(this));
        this._textField.addTargetForControlEvent(UIView.controlEvent.EnterDown, function (sender, event) {
            this.performSearch();
        }.bind(this));
        this._textField.viewHTMLElement.oninput = function (event) {
            this._textField.sendControlEventForKey(UITextField.controlEvent.TextChange, event);
        }.bind(this);
    }
    get searchButtonElement() {
        return this.viewHTMLElement.querySelector("button");
    }
    get textFieldElement() {
        return this.viewHTMLElement.querySelector("input");
    }
    set placeholderText(placeholderText) {
        this.textFieldElement.setAttribute("placeholder", placeholderText);
    }
    get placeholderText() {
        return this.textFieldElement.getAttribute("placeholder");
    }
    get text() {
        return this.textFieldElement.value;
    }
    set text(text) {
        this.textFieldElement.value = text;
    }
    focus() {
        this.textFieldElement.focus();
    }
    blur() {
        this.textFieldElement.blur();
    }
    performSearch() {
    }
    layoutSubviews() {
        super.layoutSubviews();
    }
}
///<reference path="UIObject.ts"/>
class UIKeyValueStringFilter extends UIObject {
    constructor(useSeparateWebWorker = NO) {
        super();
        this._isThreadClosed = NO;
        this._webWorker = UIKeyValueStringFilter._sharedWebWorker;
        this._class = UIKeyValueStringFilter;
        this.superclass = UIObject;
        if (useSeparateWebWorker) {
            this._webWorker = new Worker("compiledScripts/UIKeyValueStringFilterWebWorker.js");
        }
        UIKeyValueStringFilter._instanceNumber = UIKeyValueStringFilter._instanceNumber + 1;
        this._instanceNumber = UIKeyValueStringFilter._instanceNumber;
    }
    get instanceIdentifier() {
        return this._instanceNumber;
    }
    filterData(filteringString, data, excludedData, dataKeyPath, identifier, completion) {
        if (this._isThreadClosed) {
            return;
        }
        const startTime = Date.now();
        const instanceIdentifier = this.instanceIdentifier;
        this._webWorker.onmessage = function (message) {
            if (message.data.instanceIdentifier == instanceIdentifier) {
                console.log("Filtering took " + (Date.now() - startTime) + " ms.");
                completion(message.data.filteredData, message.data.filteredIndexes, message.data.identifier);
            }
        };
        try {
            this._webWorker.postMessage({
                "filteringString": filteringString,
                "data": data,
                "excludedData": excludedData,
                "dataKeyPath": dataKeyPath,
                "identifier": identifier,
                "instanceIdentifier": instanceIdentifier
            });
        }
        catch (exception) {
            completion([], [], identifier);
        }
    }
    closeThread() {
        this._isThreadClosed = YES;
        if (this._webWorker != UIKeyValueStringFilter._sharedWebWorker) {
            this._webWorker.terminate();
        }
    }
}
UIKeyValueStringFilter._sharedWebWorker = new Worker("compiledScripts/UIKeyValueStringFilterWebWorker.js");
UIKeyValueStringFilter._instanceNumber = -1;
/// <reference path="./UIView.ts" />
class UINativeScrollView extends UIView {
    constructor(elementID, viewHTMLElement) {
        super(elementID, viewHTMLElement);
        this.animationDuration = 0;
        this._class = UINativeScrollView;
        this.superclass = UIView;
        this.style.cssText = this.style.cssText + "-webkit-overflow-scrolling: touch;";
        this.style.overflow = "auto";
        // this.scrollsX = YES;
        // this.scrollsY = YES;
        this.viewHTMLElement.addEventListener("scroll", function (event) {
            this.didScrollToPosition(new UIPoint(this.viewHTMLElement.scrollLeft, this.viewHTMLElement.scrollTop));
            this.broadcastEventInSubtree({
                name: UIView.broadcastEventName.PageDidScroll,
                parameters: nil
            });
        }.bind(this));
    }
    didScrollToPosition(offsetPosition) {
    }
    get scrollsX() {
        const result = (this.style.overflowX == "scroll");
        return result;
    }
    set scrollsX(scrolls) {
        if (scrolls) {
            this.style.overflowX = "scroll";
        }
        else {
            this.style.overflowX = "hidden";
        }
    }
    get scrollsY() {
        const result = (this.style.overflowY == "scroll");
        return result;
    }
    set scrollsY(scrolls) {
        if (scrolls) {
            this.style.overflowY = "scroll";
        }
        else {
            this.style.overflowY = "hidden";
        }
    }
    get contentOffset() {
        const result = new UIPoint(this.viewHTMLElement.scrollLeft, this.viewHTMLElement.scrollTop);
        return result;
    }
    set contentOffset(offsetPoint) {
        if (this.animationDuration) {
            this.scrollXTo(this.viewHTMLElement, offsetPoint.x, this.animationDuration);
            this.scrollYTo(this.viewHTMLElement, offsetPoint.y, this.animationDuration);
            return;
        }
        this.viewHTMLElement.scrollLeft = offsetPoint.x;
        this.viewHTMLElement.scrollTop = offsetPoint.y;
    }
    scrollToBottom() {
        this.contentOffset = new UIPoint(this.contentOffset.x, this.scrollSize.height - this.frame.height);
    }
    scrollToTop() {
        this.contentOffset = new UIPoint(this.contentOffset.x, 0);
    }
    get isScrolledToBottom() {
        return this.contentOffset.isEqualTo(new UIPoint(this.contentOffset.x, this.scrollSize.height -
            this.frame.height));
    }
    get isScrolledToTop() {
        return this.contentOffset.isEqualTo(new UIPoint(this.contentOffset.x, 0));
    }
    scrollYTo(element, to, duration) {
        duration = duration * 1000;
        const start = element.scrollTop;
        const change = to - start;
        const increment = 10;
        const animateScroll = function (elapsedTime) {
            elapsedTime += increment;
            const position = this.easeInOut(elapsedTime, start, change, duration);
            element.scrollTop = position;
            if (elapsedTime < duration) {
                setTimeout(function () {
                    animateScroll(elapsedTime);
                }, increment);
            }
        }.bind(this);
        animateScroll(0);
    }
    scrollXTo(element, to, duration) {
        duration = duration * 1000;
        const start = element.scrollTop;
        const change = to - start;
        const increment = 10;
        const animateScroll = function (elapsedTime) {
            elapsedTime += increment;
            const position = this.easeInOut(elapsedTime, start, change, duration);
            element.scrollLeft = position;
            if (elapsedTime < duration) {
                setTimeout(function () {
                    animateScroll(elapsedTime);
                }, increment);
            }
        }.bind(this);
        animateScroll(0);
    }
    easeInOut(currentTime, start, change, duration) {
        currentTime /= duration / 2;
        if (currentTime < 1) {
            return change / 2 * currentTime * currentTime + start;
        }
        currentTime -= 1;
        return -change / 2 * (currentTime * (currentTime - 2) - 1) + start;
    }
}
/// <reference path="./UINativeScrollView.ts" />
class UITableView extends UINativeScrollView {
    constructor(elementID) {
        super(elementID);
        this.allRowsHaveEqualHeight = NO;
        this._visibleRows = [];
        this._firstLayoutVisibleRows = [];
        this._rowPositions = [];
        this._highestValidRowPositionIndex = 0;
        this._reusableViews = {};
        this._removedReusableViews = {};
        this._rowIDIndex = 0;
        this.reloadsOnLanguageChange = YES;
        this.sidePadding = 0;
        this._persistedData = [];
        this._needsDrawingOfVisibleRowsBeforeLayout = NO;
        this._isDrawVisibleRowsScheduled = NO;
        this.animationDuration = 0.25;
        this._class = UITableView;
        this.superclass = UINativeScrollView;
        this.scrollsX = NO;
    }
    initView(elementID, viewHTMLElement) {
        super.initView(elementID, viewHTMLElement);
        this._fullHeightView = new UIView();
        this._fullHeightView.hidden = YES;
        this._fullHeightView.userInteractionEnabled = NO;
        this.addSubview(this._fullHeightView);
    }
    loadData() {
        this._persistedData = [];
        this._calculatePositionsUntilIndex(this.numberOfRows() - 1);
        this._needsDrawingOfVisibleRowsBeforeLayout = YES;
        this.setNeedsLayout();
        // this.forEachViewInSubtree(function(view) {
        //     view.setNeedsLayout();
        // })
    }
    reloadData() {
        this._removeVisibleRows();
        this._removeAllReusableRows();
        this._rowPositions = [];
        this._highestValidRowPositionIndex = 0;
        this.loadData();
    }
    highlightChanges(previousData, newData) {
        previousData = previousData.map(function (dataPoint, index, array) {
            return JSON.stringify(dataPoint);
        });
        newData = newData.map(function (dataPoint, index, array) {
            return JSON.stringify(dataPoint);
        });
        const newIndexes = [];
        newData.forEach(function (value, index, array) {
            if (!previousData.contains(value)) {
                newIndexes.push(index);
            }
        });
        newIndexes.forEach(function (index) {
            if (this.isRowWithIndexVisible(index)) {
                this.highlightRowAsNew(this.viewForRowWithIndex(index));
            }
        }.bind(this));
    }
    highlightRowAsNew(row) {
    }
    invalidateSizeOfRowWithIndex(index, animateChange = NO) {
        if (this._rowPositions[index]) {
            this._rowPositions[index].isValid = NO;
        }
        this._highestValidRowPositionIndex = Math.min(this._highestValidRowPositionIndex, index - 1);
        // if (index == 0) {
        //     this._highestValidRowPositionIndex = 0;
        //     this._rowPositions = [];
        // }
        this._needsDrawingOfVisibleRowsBeforeLayout = YES;
        this._shouldAnimateNextLayout = animateChange;
    }
    _calculateAllPositions() {
        this._calculatePositionsUntilIndex(this.numberOfRows() - 1);
    }
    _calculatePositionsUntilIndex(maxIndex) {
        var validPositionObject = this._rowPositions[this._highestValidRowPositionIndex];
        if (!IS(validPositionObject)) {
            validPositionObject = {
                bottomY: 0,
                topY: 0,
                isValid: YES
            };
        }
        var previousBottomY = validPositionObject.bottomY;
        if (!this._rowPositions.length) {
            this._highestValidRowPositionIndex = -1;
        }
        for (var i = this._highestValidRowPositionIndex + 1; i <= maxIndex; i++) {
            var height;
            const rowPositionObject = this._rowPositions[i];
            if (IS((rowPositionObject || nil).isValid)) {
                height = rowPositionObject.bottomY - rowPositionObject.topY;
            }
            else {
                height = this.heightForRowWithIndex(i);
            }
            const positionObject = {
                bottomY: previousBottomY + height,
                topY: previousBottomY,
                isValid: YES
            };
            if (i < this._rowPositions.length) {
                this._rowPositions[i] = positionObject;
            }
            else {
                this._rowPositions.push(positionObject);
            }
            this._highestValidRowPositionIndex = i;
            previousBottomY = previousBottomY + height;
        }
    }
    indexesForVisibleRows(paddingRatio = 0.5) {
        const firstVisibleY = this.contentOffset.y - this.bounds.height * paddingRatio;
        const lastVisibleY = firstVisibleY + this.bounds.height * (1 + paddingRatio);
        const numberOfRows = this.numberOfRows();
        if (this.allRowsHaveEqualHeight) {
            const rowHeight = this.heightForRowWithIndex(0);
            var firstIndex = firstVisibleY / rowHeight;
            var lastIndex = lastVisibleY / rowHeight;
            firstIndex = Math.trunc(firstIndex);
            lastIndex = Math.trunc(lastIndex) + 1;
            firstIndex = Math.max(firstIndex, 0);
            lastIndex = Math.min(lastIndex, numberOfRows - 1);
            var result = [];
            for (var i = firstIndex; i < lastIndex + 1; i++) {
                result.push(i);
            }
            return result;
        }
        var accumulatedHeight = 0;
        var result = [];
        this._calculateAllPositions();
        const rowPositions = this._rowPositions;
        for (var i = 0; i < numberOfRows; i++) {
            const height = rowPositions[i].bottomY - rowPositions[i].topY; // this.heightForRowWithIndex(i)
            accumulatedHeight = accumulatedHeight + height;
            if (accumulatedHeight >= firstVisibleY) {
                result.push(i);
            }
            if (accumulatedHeight >= lastVisibleY) {
                break;
            }
        }
        return result;
    }
    _removeVisibleRows() {
        const visibleRows = [];
        this._visibleRows.forEach(function (row, index, array) {
            this._persistedData[row._UITableViewRowIndex] = this.persistenceDataItemForRowWithIndex(row._UITableViewRowIndex, row);
            row.removeFromSuperview();
            this._removedReusableViews[row._UITableViewReusabilityIdentifier].push(row);
        }, this);
        this._visibleRows = visibleRows;
    }
    _removeAllReusableRows() {
        this._reusableViews.forEach(function (rows) {
            rows.forEach(function (row, index, array) {
                this._persistedData[row._UITableViewRowIndex] = this.persistenceDataItemForRowWithIndex(row._UITableViewRowIndex, row);
                row.removeFromSuperview();
                this._markReusableViewAsUnused(row);
            }.bind(this));
        }.bind(this));
    }
    _markReusableViewAsUnused(row) {
        if (!this._removedReusableViews[row._UITableViewReusabilityIdentifier].contains(row)) {
            this._removedReusableViews[row._UITableViewReusabilityIdentifier].push(row);
        }
    }
    _drawVisibleRows() {
        if (!this.isMemberOfViewTree) {
            return;
        }
        const visibleIndexes = this.indexesForVisibleRows();
        const minIndex = visibleIndexes[0];
        const maxIndex = visibleIndexes[visibleIndexes.length - 1];
        const removedViews = [];
        const visibleRows = [];
        this._visibleRows.forEach(function (row, index, array) {
            if (row._UITableViewRowIndex < minIndex || row._UITableViewRowIndex > maxIndex) {
                //row.removeFromSuperview();
                this._persistedData[row._UITableViewRowIndex] = this.persistenceDataItemForRowWithIndex(row._UITableViewRowIndex, row);
                this._removedReusableViews[row._UITableViewReusabilityIdentifier].push(row);
                removedViews.push(row);
            }
            else {
                visibleRows.push(row);
            }
        }, this);
        this._visibleRows = visibleRows;
        visibleIndexes.forEach(function (rowIndex, index, array) {
            if (this.isRowWithIndexVisible(rowIndex)) {
                return;
            }
            const view = this.viewForRowWithIndex(rowIndex);
            //view._UITableViewRowIndex = rowIndex;
            this._firstLayoutVisibleRows.push(view);
            this._visibleRows.push(view);
            this.addSubview(view);
        }, this);
        for (var i = 0; i < removedViews.length; i++) {
            var view = removedViews[i];
            if (this._visibleRows.indexOf(view) == -1) {
                //this._persistedData[view._UITableViewRowIndex] = this.persistenceDataItemForRowWithIndex(view._UITableViewRowIndex, view);
                view.removeFromSuperview();
                //this._removedReusableViews[view._UITableViewReusabilityIdentifier].push(view);
            }
        }
        //this.setNeedsLayout();
    }
    visibleRowWithIndex(rowIndex) {
        for (var i = 0; i < this._visibleRows.length; i++) {
            const row = this._visibleRows[i];
            if (row._UITableViewRowIndex == rowIndex) {
                return row;
            }
        }
        return nil;
    }
    isRowWithIndexVisible(rowIndex) {
        return IS(this.visibleRowWithIndex(rowIndex));
    }
    reusableViewForIdentifier(identifier, rowIndex) {
        if (!this._removedReusableViews[identifier]) {
            this._removedReusableViews[identifier] = [];
        }
        if (this._removedReusableViews[identifier] && this._removedReusableViews[identifier].length) {
            const view = this._removedReusableViews[identifier].pop();
            view._UITableViewRowIndex = rowIndex;
            Object.assign(view, this._persistedData[rowIndex] || this.defaultRowPersistenceDataItem());
            return view;
        }
        if (!this._reusableViews[identifier]) {
            this._reusableViews[identifier] = [];
        }
        const newView = this.newReusableViewForIdentifier(identifier, this._rowIDIndex);
        this._rowIDIndex = this._rowIDIndex + 1;
        if (this._rowIDIndex > 40) {
            const asd = 1;
        }
        newView._UITableViewReusabilityIdentifier = identifier;
        newView._UITableViewRowIndex = rowIndex;
        Object.assign(newView, this._persistedData[rowIndex] || this.defaultRowPersistenceDataItem());
        this._reusableViews[identifier].push(newView);
        return newView;
    }
    // Functions that should be overridden to draw the correct content START
    newReusableViewForIdentifier(identifier, rowIDIndex) {
        const view = new UIButton(this.elementID + "Row" + rowIDIndex);
        view.stopsPointerEventPropagation = NO;
        view.pausesPointerEvents = NO;
        return view;
    }
    heightForRowWithIndex(index) {
        return 50;
    }
    numberOfRows() {
        return 10000;
    }
    defaultRowPersistenceDataItem() {
    }
    persistenceDataItemForRowWithIndex(rowIndex, row) {
    }
    viewForRowWithIndex(rowIndex) {
        const row = this.reusableViewForIdentifier("Row", rowIndex);
        row.titleLabel.text = "Row " + rowIndex;
        return row;
    }
    // Functions that should be overridden to draw the correct content END
    // Functions that trigger redrawing of the content
    didScrollToPosition(offsetPosition) {
        super.didScrollToPosition(offsetPosition);
        this.forEachViewInSubtree(function (view) {
            view._isPointerValid = NO;
        });
        if (!this._isDrawVisibleRowsScheduled) {
            this._isDrawVisibleRowsScheduled = YES;
            UIView.runFunctionBeforeNextFrame(function () {
                this._calculateAllPositions();
                this._drawVisibleRows();
                this.setNeedsLayout();
                this._isDrawVisibleRowsScheduled = NO;
            }.bind(this));
        }
    }
    wasAddedToViewTree() {
        this.loadData();
    }
    setFrame(rectangle, zIndex, performUncheckedLayout) {
        const frame = this.frame;
        super.setFrame(rectangle, zIndex, performUncheckedLayout);
        if (frame.isEqualTo(rectangle) && !performUncheckedLayout) {
            return;
        }
        this._needsDrawingOfVisibleRowsBeforeLayout = YES;
    }
    didReceiveBroadcastEvent(event) {
        super.didReceiveBroadcastEvent(event);
        if (event.name == UIView.broadcastEventName.LanguageChanged && this.reloadsOnLanguageChange) {
            this.reloadData();
        }
    }
    _layoutAllRows(positions = this._rowPositions) {
        const bounds = this.bounds;
        this._visibleRows.forEach(function (row, index, array) {
            const frame = bounds.copy();
            const positionObject = positions[row._UITableViewRowIndex];
            frame.min.y = positionObject.topY;
            frame.max.y = positionObject.bottomY;
            row.frame = frame;
            row.style.width = "" + (bounds.width - this.sidePadding * 2).integerValue + "px";
            row.style.left = "" + this.sidePadding.integerValue + "px";
        }, this);
        this._fullHeightView.frame = bounds.rectangleWithHeight((positions.lastElement ||
            nil).bottomY).rectangleWithWidth(bounds.width * 0.5);
        this._firstLayoutVisibleRows = [];
    }
    _animateLayoutAllRows() {
        UIView.animateViewOrViewsWithDurationDelayAndFunction(this._visibleRows, this.animationDuration, 0, undefined, function () {
            this._layoutAllRows();
        }.bind(this), function () {
            // this._calculateAllPositions()
            // this._layoutAllRows()
        }.bind(this));
    }
    layoutSubviews() {
        const previousPositions = JSON.parse(JSON.stringify(this._rowPositions));
        const previousVisibleRowsLength = this._visibleRows.length;
        if (this._needsDrawingOfVisibleRowsBeforeLayout) {
            //this._calculateAllPositions()
            this._drawVisibleRows();
            this._needsDrawingOfVisibleRowsBeforeLayout = NO;
        }
        super.layoutSubviews();
        if (!this.numberOfRows() || !this.isMemberOfViewTree) {
            return;
        }
        if (this._shouldAnimateNextLayout) {
            // Need to do layout with the previous positions
            this._layoutAllRows(previousPositions);
            if (previousVisibleRowsLength < this._visibleRows.length) {
                UIView.runFunctionBeforeNextFrame(function () {
                    this._animateLayoutAllRows();
                }.bind(this));
            }
            else {
                this._animateLayoutAllRows();
            }
            this._shouldAnimateNextLayout = NO;
        }
        else {
            // if (this._needsDrawingOfVisibleRowsBeforeLayout) {
            //     this._drawVisibleRows();
            //     this._needsDrawingOfVisibleRowsBeforeLayout = NO;
            // }
            this._calculateAllPositions();
            this._layoutAllRows();
        }
    }
    intrinsicContentHeight(constrainingWidth = 0) {
        var result = 10000000000000000000000000;
        if (this._rowPositions.length) {
            result = this._rowPositions[this._rowPositions.length - 1].bottomY;
        }
        return result;
    }
}
///<reference path="UIView.ts"/>
class UIImageView extends UIView {
    constructor(elementID, viewHTMLElement = null) {
        super(elementID, viewHTMLElement, "img");
        this._hiddenWhenEmpty = NO;
        this._class = UIImageView;
        this.superclass = UIView;
        //this.actionIndicator = new UIActionIndicator(elementID + "ActionIndicator");
    }
    get viewHTMLElement() {
        return super.viewHTMLElement;
    }
    static dataURL(url, callback) {
        const xhr = new XMLHttpRequest();
        xhr.open("get", url);
        xhr.responseType = "blob";
        xhr.onload = function () {
            const fr = new FileReader();
            fr.onload = function () {
                callback(this.result);
            };
            fr.readAsDataURL(xhr.response); // async call
        };
        xhr.send();
    }
    static dataURLWithMaxSize(URLString, maxSize, completion) {
        const imageView = new UIImageView();
        imageView.imageSource = URLString;
        imageView.viewHTMLElement.onload = function () {
            const originalSize = imageView.intrinsicContentSize();
            var multiplier = maxSize / Math.max(originalSize.height, originalSize.width);
            multiplier = Math.min(1, multiplier);
            const result = imageView.getDataURL((originalSize.height * multiplier).integerValue, (originalSize.width *
                multiplier).integerValue);
            completion(result);
        };
    }
    static dataURLWithSizes(URLString, height, width, completion) {
        const imageView = new UIImageView();
        imageView.imageSource = URLString;
        imageView.viewHTMLElement.onload = function () {
            const result = imageView.getDataURL(height, width);
            completion(result);
        };
    }
    getDataURL(height, width) {
        const img = this.viewHTMLElement;
        // Create an empty canvas element
        const canvas = document.createElement("canvas");
        canvas.width = width;
        canvas.height = height;
        // Copy the image contents to the canvas
        const ctx = canvas.getContext("2d");
        ctx.drawImage(img, 0, 0, width, height);
        // Get the data-URL formatted image
        // Firefox supports PNG and JPEG. You could check img.src to
        // guess the original format, but be aware the using "image/jpg"
        // will re-encode the image.
        const dataURL = canvas.toDataURL("image/png");
        return dataURL;
        //return dataURL.replace(/^data:image\/(png|jpg);base64,/, "");
    }
    get imageSource() {
        return this.viewHTMLElement.src;
    }
    set imageSource(sourceString) {
        if (IS_NOT(sourceString)) {
            sourceString = "";
        }
        this.viewHTMLElement.src = sourceString;
        if (this.hiddenWhenEmpty) {
            this.hidden = IS_NOT(this.imageSource);
        }
        if (!sourceString || !sourceString.length) {
            //this.actionIndicator.stop();
            this.hidden = YES;
            return;
        }
        else {
            this.hidden = NO;
        }
        // this.superview.addSubview(this.actionIndicator);
        // this.actionIndicator.frame = this.frame;
        // this.actionIndicator.start();
        // this.actionIndicator.backgroundColor = UIColor.redColor
        // @ts-ignore
        this.viewHTMLElement.onload = function (event) {
            this.superview.setNeedsLayout();
            //this.actionIndicator.removeFromSuperview();
        }.bind(this);
    }
    setImageSource(key, defaultString) {
        const languageName = UICore.languageService.currentLanguageKey;
        this.imageSource = UICore.languageService.stringForKey(key, languageName, defaultString, nil);
    }
    didReceiveBroadcastEvent(event) {
        super.didReceiveBroadcastEvent(event);
        if (event.name == UIView.broadcastEventName.LanguageChanged || event.name ==
            UIView.broadcastEventName.AddedToViewTree) {
            this._setImageSourceFromKeyIfPossible();
        }
    }
    willMoveToSuperview(superview) {
        super.willMoveToSuperview(superview);
        this._setImageSourceFromKeyIfPossible();
    }
    _setImageSourceFromKeyIfPossible() {
        if (this._sourceKey && this._defaultSource) {
            this.setImageSource(this._sourceKey, this._defaultSource);
        }
    }
    get fillMode() {
        return this._fillMode;
    }
    set fillMode(fillMode) {
        this._fillMode = fillMode;
        this.style.objectFit = fillMode;
    }
    get hiddenWhenEmpty() {
        return this._hiddenWhenEmpty;
    }
    set hiddenWhenEmpty(hiddenWhenEmpty) {
        this._hiddenWhenEmpty = hiddenWhenEmpty;
        if (hiddenWhenEmpty) {
            this.hidden = IS_NOT(this.imageSource);
        }
    }
    didMoveToSuperview(superview) {
        super.didMoveToSuperview(superview);
    }
    layoutSubviews() {
        super.layoutSubviews();
    }
    intrinsicContentSize() {
        const result = new UIRectangle(0, 0, this.viewHTMLElement.naturalHeight, this.viewHTMLElement.naturalWidth);
        return result;
    }
    intrinsicContentSizeWithConstraints(constrainingHeight = 0, constrainingWidth = 0) {
        const heightRatio = constrainingHeight / this.viewHTMLElement.naturalHeight;
        const widthRatio = constrainingWidth / this.viewHTMLElement.naturalWidth;
        const multiplier = Math.max(heightRatio, widthRatio);
        const result = new UIRectangle(0, 0, this.viewHTMLElement.naturalHeight *
            multiplier, this.viewHTMLElement.naturalWidth * multiplier);
        return result;
    }
}
UIImageView.fillMode = {
    "stretchToFill": "fill",
    "aspectFit": "contain",
    "aspectFill": "cover",
    "center": "none",
    "aspectFitIfLarger": "scale-down"
};
/// <reference path="UITextView.ts" />
// @ts-ignore
class UITextField extends UITextView {
    constructor(elementID, viewHTMLElement = null, type = UITextView.type.textField) {
        super(elementID, type, viewHTMLElement);
        this._class = UITextField;
        this.superclass = UITextView;
        this.viewHTMLElement.setAttribute("type", "text");
        this.backgroundColor = UIColor.whiteColor;
        this.addTargetForControlEvent(UIView.controlEvent.PointerUpInside, function (sender, event) {
            sender.focus();
        });
        this.viewHTMLElement.oninput = (event) => {
            this.sendControlEventForKey(UITextField.controlEvent.TextChange, event);
        };
        this.style.webkitUserSelect = "text";
        this.nativeSelectionEnabled = YES;
        this.pausesPointerEvents = NO;
    }
    get addControlEventTarget() {
        return super.addControlEventTarget;
    }
    get viewHTMLElement() {
        return this._viewHTMLElement;
    }
    set text(text) {
        this.viewHTMLElement.value = text;
    }
    get text() {
        return this.viewHTMLElement.value;
    }
    set placeholderText(text) {
        this.viewHTMLElement.placeholder = text;
    }
    get placeholderText() {
        return this.viewHTMLElement.placeholder;
    }
    setPlaceholderText(key, defaultString) {
        this._placeholderTextKey = key;
        this._defaultPlaceholderText = defaultString;
        const languageName = UICore.languageService.currentLanguageKey;
        this.placeholderText = UICore.languageService.stringForKey(key, languageName, defaultString, nil);
    }
    didReceiveBroadcastEvent(event) {
        super.didReceiveBroadcastEvent(event);
        if (event.name == UIView.broadcastEventName.LanguageChanged || event.name ==
            UIView.broadcastEventName.AddedToViewTree) {
            this._setPlaceholderFromKeyIfPossible();
        }
    }
    willMoveToSuperview(superview) {
        super.willMoveToSuperview(superview);
        this._setPlaceholderFromKeyIfPossible();
    }
    _setPlaceholderFromKeyIfPossible() {
        if (this._placeholderTextKey && this._defaultPlaceholderText) {
            this.setPlaceholderText(this._placeholderTextKey, this._defaultPlaceholderText);
        }
    }
    get isSecure() {
        const result = (this.viewHTMLElement.type == "password");
        return result;
    }
    set isSecure(secure) {
        var type = "text";
        if (secure) {
            type = "password";
        }
        this.viewHTMLElement.type = type;
    }
}
UITextField.controlEvent = Object.assign({}, UIView.controlEvent, {
    "TextChange": "TextChange"
});
///<reference path="../UICore/UIButton.ts"/>
class SearchableDropdownRow extends UIButton {
    constructor(elementID) {
        super(elementID);
        this.type = SearchableDropdownRow.type.selectableItem;
        this._class = SearchableDropdownRow;
        this.superclass = UIButton;
    }
    initView(elementID, viewHTMLElement, initViewData) {
        super.initView(elementID, viewHTMLElement, initViewData);
        this._checkbox = new CBCheckbox(elementID + "Checkbox");
        this._checkbox.userInteractionEnabled = NO;
        this.addSubview(this._checkbox);
        this.style.outline = "none";
    }
    set titleText(titleText) {
        this.titleLabel.text = titleText;
        this._checkbox.titleLabel.text = titleText;
    }
    get titleText() {
        return this.titleLabel.text;
    }
    set selected(selected) {
        super.selected = selected;
        this._checkbox.selected = selected;
    }
    get selected() {
        return super.selected;
    }
    set focused(focused) {
        this._focused = focused;
        this.updateContentForCurrentState();
    }
    get focused() {
        return this._focused;
    }
    updateContentForNormalState() {
        if (this.type == SearchableDropdownRow.type.sectionTitle) {
            this.backgroundColor = UIColor.transparentColor;
            this.titleLabel.textColor = CBColor.primaryContentColor.colorWithAlpha(0.5);
            this.style.borderTop = "1px solid rgba(0, 0, 0, 0.3)";
            this.titleLabel.style.marginLeft = "";
            this.textSuffix = "";
            if (this._checkbox) {
                this._checkbox.hidden = YES;
            }
            this.titleLabel.hidden = NO;
        }
        else if (this.type == SearchableDropdownRow.type.selectedItem) {
            this.backgroundColor = UIColor.transparentColor;
            this.titleLabel.textColor = CBColor.primaryTintColor;
            this.style.borderTop = "";
            //this.titleLabel.style.marginLeft = "20px"
            this.textSuffix = "";
            if (this._checkbox) {
                this._checkbox.hidden = NO;
                this._checkbox.titleLabel.textColor = CBColor.primaryContentColor;
            }
            this.titleLabel.hidden = YES;
        }
        else if (this.type == SearchableDropdownRow.type.customItem) {
            this.backgroundColor = UIColor.transparentColor;
            this.titleLabel.textColor = CBColor.primaryTintColor;
            this.style.borderTop = "1px solid rgba(0, 0, 0, 0.3)";
            this.titleLabel.style.marginLeft = "";
            this.textSuffix = LanguageService.stringForKey("searchableDropdownCustomItem", LanguageService.currentLanguageKey, "-Custom item");
            if (this._checkbox) {
                this._checkbox.hidden = YES;
            }
            this.titleLabel.hidden = NO;
        }
        else {
            this.backgroundColor = UIColor.transparentColor;
            this.titleLabel.textColor = CBColor.primaryTintColor;
            this.style.borderTop = "";
            this.titleLabel.style.marginLeft = "20px";
            this.textSuffix = "";
            if (this._checkbox) {
                this._checkbox.hidden = YES;
            }
            this.titleLabel.hidden = NO;
        }
        this.userInteractionEnabled = YES;
    }
    get textSuffix() {
        const result = this.titleLabel.textSuffix;
        return result;
    }
    set textSuffix(textSuffix) {
        this.titleLabel.textSuffix = textSuffix;
        this._checkbox.titleLabel.textSuffix = textSuffix;
    }
    updateContentForHoveredState() {
        this.backgroundColor = CBColor.primaryContentColor.colorWithAlpha(0.05);
    }
    updateContentForHighlightedState() {
        this.backgroundColor = CBColor.primaryContentColor.colorWithAlpha(0.08);
    }
    updateContentForFocusedState() {
        this.backgroundColor = CBColor.primaryContentColor.colorWithAlpha(0.08);
    }
    updateContentForSelectedState() {
        this.updateContentForNormalState();
        if (this.type == SearchableDropdownRow.type.selectableItem || this.type ==
            SearchableDropdownRow.type.customItem) {
            this._checkbox.hidden = NO;
            //this._checkbox.titleLabel.hidden = YES;
            this._checkbox.titleLabel.textColor = CBColor.primaryTintColor;
            this._titleLabel.hidden = YES;
        }
    }
    wasRemovedFromViewTree() {
        super.wasRemovedFromViewTree();
        this.highlighted = NO;
        this.hovered = NO;
    }
    layoutSubviews() {
        super.layoutSubviews();
        if (this._checkbox) {
            this._checkbox.style.top = "25%";
            this._checkbox.style.height = "30px";
            this._checkbox.style.left = "" + this.contentPadding + "px";
            this._checkbox.style.right = "" + this.contentPadding + "px";
        }
    }
}
SearchableDropdownRow.type = {
    "sectionTitle": "SectionTitle",
    "selectedItem": "SelectedItem",
    "selectableItem": "SelectableItem",
    "customItem": "CustomItem"
};
/// <reference path="../UICore/UIButton.ts" />
/// <reference path="../UICore/UIKeyValueStringFilter.ts" />
/// <reference path="../UICore/UITableView.ts" />
/// <reference path="../UICore/UIImageView.ts" />
/// <reference path="../UICore/UITextView.ts" />
/// <reference path="../UICore/UITextField.ts" />
/// <reference path="../UICore/UIDialogView.ts" />
/// <reference path="./SearchTextField.ts" />
/// <reference path="./CBCheckbox.ts" />
/// <reference path="./CBColor.ts" />
/// <reference path="./SearchableDropdownRow.ts" />
/// <reference path="../RootViewController.ts" />
/// <reference path="./LanguageService.ts" />
class SearchableDropdown extends UIButton {
    constructor(elementID) {
        super(elementID);
        this._data = [];
        this._filteredData = [];
        this._excludedData = [];
        this.tintColor = CBColor.primaryTintColor;
        this.selectedIndices = [];
        this._selectedData = [];
        this._drawingData = [];
        this._isDrawingDataValid = NO;
        this.isSingleSelection = NO;
        this.showsSelectedSectionInMultipleSelectionMode = NO;
        this.allowsCustomItem = NO;
        this.keepFocusedRowVisible = YES;
        this._class = SearchableDropdown;
        this.superclass = UIButton;
    }
    initView(elementID, viewHTMLElement, initViewData) {
        super.initView(elementID, viewHTMLElement, initViewData);
        //this.style.borderRadius = "2px";
        this._titleLabel.text = "Current value";
        this._titleLabel.textAlignment = UITextView.textAlignment.left;
        this.overflowLabel = new UITextView(elementID + "OverflowLabel");
        this.overflowLabel.textColor = CBColor.primaryContentColor;
        this.overflowLabel.textAlignment = UITextView.textAlignment.right;
        this.addSubview(this.overflowLabel);
        this._rightImageView = new UIImageView(this.elementID + "RightImageView");
        this._rightImageView.imageSource = "images/baseline-arrow_drop_down-24px.svg";
        this._rightImageView.userInteractionEnabled = NO;
        this.addSubview(this._rightImageView);
        this.setNeedsLayout();
        this._containerView = new UIView(elementID + "ContainerView");
        this._containerView.style.boxShadow = "0 9px 13px 0 rgba(0,0,0,0.26)";
        this._containerView.style.borderRadius = "2px";
        this._searchTextField = new SearchTextField(elementID + "SearchTextField");
        this._searchTextField.placeholderText = "Search";
        this._containerView.addSubview(this._searchTextField);
        this._searchTextField._textField.addTargetForControlEvent(UITextField.controlEvent.TextChange, function (sender, event) {
            this.updateFilteredData(this._searchTextField.text);
        }.bind(this));
        this._searchTextField._textField.addTargetForControlEvent(UIView.controlEvent.EscDown, function (sender, event) {
            if (IS(this._searchTextField.text)) {
                this._searchTextField.text = "";
                this.updateFilteredData("");
            }
            else {
                this._dialogView.dismiss(YES);
            }
        }.bind(this));
        this._searchTextField._textField.addTargetForControlEvent(UIView.controlEvent.DownArrowDown, function (sender, event) {
            if (this.focusedRowIndex < (this.drawingData.length - 1)) {
                this.focusedRowIndex = this.focusedRowIndex + 1;
            }
        }.bind(this));
        this._searchTextField._textField.addTargetForControlEvent(UIView.controlEvent.UpArrowDown, function (sender, event) {
            if (this.focusedRowIndex > 0) {
                this.focusedRowIndex = this.focusedRowIndex - 1;
            }
        }.bind(this));
        this._searchTextField._textField.addTargetForControlEvent(UIView.controlEvent.EnterDown, function (sender, event) {
            const isTouchDevice = "ontouchstart" in document.documentElement;
            if (isTouchDevice) {
                this._searchTextField.blur();
                return;
            }
            const datapoint = this.drawingData[this.focusedRowIndex];
            const alreadySelected = this.selectedDataContains(datapoint);
            if (alreadySelected) {
                this.selectedData.removeElement(datapoint);
            }
            else if (this.isSingleSelection) {
                this.selectedIndices = [this.focusedRowIndex];
                this.selectedData = [datapoint];
                // view.selected = YES;
                this.selectionDidChange(this.selectedData);
                this.updateContentForCurrentSelection();
                this._dialogView.dismiss();
            }
            else {
                this.selectedData.push(datapoint);
            }
        }.bind(this));
        this._tableView = new UITableView(elementID + "TableView");
        this._containerView.addSubview(this._tableView);
        this._tableView.backgroundColor = UIColor.whiteColor;
        this._dialogView = new UIDialogView(elementID + "DialogView");
        this._dialogView.view = this._containerView;
        this._dialogView.backgroundColor = UIColor.transparentColor;
        this.addTargetForControlEvents([
            UIView.controlEvent.PointerUpInside, UIView.controlEvent.EnterDown
        ], function (sender, event) {
            if (this._dialogView.isVisible) {
                this.closeDropdown();
            }
            else {
                this.openDropdown();
            }
        }.bind(this));
        this._dialogView.addTargetForControlEvent(UIView.controlEvent.PointerDown, function (sender, event) {
            if (sender.viewHTMLElement == event.target) {
                sender.dismiss();
            }
        });
        const dialogLayoutFunction = this._dialogView.layoutSubviews;
        this._dialogView.layoutSubviews = function () {
            //dialogLayoutFunction.call(this._dialogView);
            this._dialogView.frame = this.rootView.bounds; //this.superview.rectangleInView(this.superview.bounds, this.rootView);
            const padding = RootViewController.paddingLength;
            const labelHeight = padding;
            const searchTextFieldHeight = this.bounds.height;
            this._containerView.frame = this.superview.rectangleInView(this.frame, this.rootView)
                .rectangleWithHeight(this.expandedContainerViewHeight);
            this._searchTextField.frame = this._containerView.bounds.rectangleWithHeight(searchTextFieldHeight).rectangleWithInsets(0, 16, 0, 0);
            this._tableView.frame = this._containerView.bounds.rectangleWithInsets(0, 0, 0, searchTextFieldHeight);
            //this._tableView.style.maxHeight = "" + this._tableView.intrinsicContentHeight() + "px";
        }.bind(this);
        // var animateDialogAppearing = this._dialogView.animateAppearing;
        // this._dialogView.animateAppearing = function (this: SearchableDropdown) {
        //     animateDialogAppearing.call(this._dialogView);
        //     this.style.zIndex = "101";
        //     if (this.highlightLabel) {
        //         this.highlightLabel.textColor = UIColor.whiteColor;
        //         this.highlightLabel.style.zIndex = "101";
        //         this.highlightLabel.userInteractionEnabled = NO;
        //     }
        //     this._searchTextField.focus();
        // }.bind(this)
        // var animateDialogDisappearing = this._dialogView.animateDisappearing;
        // this._dialogView.animateDisappearing = function (this: SearchableDropdown) {
        //     animateDialogDisappearing.call(this._dialogView);
        //     this.style.zIndex = "0";
        //     if (this.highlightLabel) {
        //         this.highlightLabel.textColor = CBColor.primaryContentColor;
        //         this.highlightLabel.style.zIndex = "0";
        //         this.highlightLabel.userInteractionEnabled = YES;
        //     }
        //     this._searchTextField.blur();
        // }.bind(this)
        this._tableView.numberOfRows = function () {
            var result = this.drawingData.length;
            if (IS_NOT(this.isSingleSelection) && this.showsSelectedSectionInMultipleSelectionMode) {
                result = result + this.selectedData.length;
            }
            return result;
        }.bind(this);
        const newReusableViewForIdentifierFunction = this._tableView.newReusableViewForIdentifier.bind(this._tableView);
        this._tableView.newReusableViewForIdentifier = function (identifier, rowIndex) {
            //if (identifier == "SelectedItemRow") {
            const view = new SearchableDropdownRow(elementID + identifier + rowIndex);
            view.stopsPointerEventPropagation = NO;
            view.pausesPointerEvents = NO;
            return view;
            // }
            // return newReusableViewForIdentifierFunction(identifier, rowIndex);
        }.bind(this);
        const viewForSelectedItemRow = function (index) {
            const view = this._tableView.reusableViewForIdentifier("SelectedItemRow", index);
            view.titleLabel.text = LanguageService.stringForCurrentLanguage(this.selectedData[index].title);
            view.selected = YES;
            return view;
        }.bind(this);
        this._tableView.viewForRowWithIndex = function (index) {
            const view = this._tableView.reusableViewForIdentifier("SubjectView", index);
            view.style.borderBottomColor = "";
            view.style.borderBottomStyle = "";
            view.style.borderBottomWidth = "";
            const rowWasHovered = function (sender, event) {
                this.focusedRowIndex = index;
            }.bind(this);
            view.addTargetForControlEvent(UIButton.controlEvent.PointerHover, rowWasHovered);
            view.removeTargetForControlEvent(UIButton.controlEvent.PointerHover, view._SearchableDropdownRowWasHoveredFunction);
            view._SearchableDropdownRowWasHoveredFunction = rowWasHovered;
            view.focused = (this.focusedRowIndex == index);
            if (!this.isSingleSelection && this.showsSelectedSectionInMultipleSelectionMode) {
                if (index < this.selectedData.length) {
                    view.type = SearchableDropdownRow.type.selectedItem;
                    view.titleText = LanguageService.stringForCurrentLanguage(this.selectedData[index].title);
                    view.selected = YES;
                    view.updateContentForCurrentState();
                    if (index == this.selectedData.length - 1) {
                        view.style.borderBottomColor = UIColor.colorWithRGBA(100, 100, 100).stringValue;
                        view.style.borderBottomStyle = "solid";
                        view.style.borderBottomWidth = "1px";
                    }
                    var viewWasTapped = function (sender, event) {
                        this.selectedIndices.removeElementAtIndex(index);
                        const selectedItem = this.selectedData[index];
                        this.selectedData.removeElement(selectedItem);
                        view.selected = NO;
                        this.selectionDidChange(this.selectedData);
                        //this.performFunctionWithDelay(0.15, function (this: SearchableDropdown<T>) {
                        this.updateContentForCurrentSelection();
                        //}.bind(this))
                        // // Unselecting main body item if possible
                        // var visibleRows = this._tableView._visibleRows as SearchableDropdownRow[]
                        // var visibleIndices = this._tableView.indexesForVisibleRows()
                        // visibleRows.forEach(function (this: SearchableDropdown<T>, visibleRow: SearchableDropdownRow, visibleRowIndex: number, visibleRowsArray: SearchableDropdownRow[]) {
                        //     var rowDataIndex = visibleIndices[visibleRowIndex]
                        //     visibleRow.selected = this.selectedData.contains(this.drawingData[rowDataIndex])
                        // }.bind(this))
                        this._searchTextField.focus();
                        if (view.viewWasTappedFunction) {
                            view.removeTargetForControlEvents([
                                UIView.controlEvent.EnterDown, UIView.controlEvent.PointerTap
                            ], view.viewWasTappedFunction);
                        }
                    }.bind(this);
                    // Removing previous event target if possible
                    if (view.viewWasTappedFunction) {
                        view.removeTargetForControlEvents([
                            UIView.controlEvent.EnterDown, UIView.controlEvent.PointerTap
                        ], view.viewWasTappedFunction);
                    }
                    // Adding event target
                    view.addTargetForControlEvents([
                        UIView.controlEvent.EnterDown, UIView.controlEvent.PointerTap
                    ], viewWasTapped);
                    view.viewWasTappedFunction = viewWasTapped;
                    return view;
                }
                index = index - this.selectedData.length;
            }
            // Datapoint
            const datapoint = this.drawingData[index];
            if (IS_NOT(datapoint)) {
                return;
            }
            // Setting different style for section title rows
            if (datapoint.isADropdownDataSection) {
                view.type = SearchableDropdownRow.type.sectionTitle;
                view.userInteractionEnabled = NO;
            }
            else {
                view.type = SearchableDropdownRow.type.selectableItem;
                view.userInteractionEnabled = YES;
            }
            if (datapoint._id == (this._customItem || nil)._id) {
                view.type = SearchableDropdownRow.type.customItem;
            }
            view.updateContentForNormalState();
            view.updateContentForCurrentState();
            // Setting row title
            view.titleText = LanguageService.stringForCurrentLanguage(datapoint.title);
            view.titleLabel.textAlignment = UITextView.textAlignment.left;
            // Highlighting of row
            //view.highlighted = (this.highlightedRowIndex == index);
            // Selecting of rows
            view.selected = this.selectedRowIdentifiers.contains(datapoint._id); //this.selectedIndices.contains(index); //this.selectedDataContains(datapoint);
            var viewWasTapped = function (sender, event) {
                if (view.selected) {
                    this.selectedIndices.removeElement(index);
                    this.selectedData.removeElement(datapoint);
                }
                else {
                    if (this.isSingleSelection) {
                        this.selectedIndices = [index];
                        this.selectedData = [datapoint];
                        // view.selected = YES;
                        this.selectionDidChange(this.selectedData);
                        this.updateContentForCurrentSelection();
                        this._dialogView.dismiss();
                        return;
                    }
                    else {
                        this.selectedIndices.push(index);
                        this.selectedData.push(datapoint);
                    }
                }
                //view.selected = !view.selected;
                const selectedData = this.selectedData;
                if (!view.selected) {
                    view.selected = YES;
                    this.performFunctionWithDelay(0.25, function () {
                        this.selectionDidChange(selectedData);
                        this.updateContentForCurrentSelection();
                        if (this.showsSelectedSectionInMultipleSelectionMode) {
                            this._tableView.contentOffset = this._tableView.contentOffset.pointByAddingY(view.frame.height);
                        }
                    }.bind(this));
                }
                else {
                    view._checkbox.selected = NO;
                    this.selectionDidChange(selectedData);
                    this.performFunctionWithDelay(0.15, function () {
                        view.selected = NO;
                        this.updateContentForCurrentSelection();
                        if (this.showsSelectedSectionInMultipleSelectionMode) {
                            this._tableView.contentOffset = this._tableView.contentOffset.pointByAddingY(-view.frame.height);
                        }
                    }.bind(this));
                }
                this._searchTextField.focus();
            }.bind(this);
            // Removing previous event target if possible
            if (view.viewWasTappedFunction) {
                view.removeTargetForControlEvents([
                    UIView.controlEvent.EnterDown, UIView.controlEvent.PointerUpInside
                ], view.viewWasTappedFunction);
            }
            // Adding event target
            view.addTargetForControlEvents([
                UIView.controlEvent.EnterDown, UIView.controlEvent.PointerUpInside
            ], viewWasTapped);
            view.viewWasTappedFunction = viewWasTapped;
            return view;
        }.bind(this);
        this._keyValueStringFilter = new UIKeyValueStringFilter();
    }
    openDropdown() {
        this._dialogView.showInView(this.rootView, YES);
        this._searchTextField.focus();
    }
    closeDropdown() {
        this._dialogView.dismiss(YES);
    }
    boundsDidChange() {
        super.boundsDidChange();
        this.setNeedsLayout();
    }
    get focusedRowIndex() {
        return this._focusedRowIndex;
    }
    set focusedRowIndex(focusedRowIndex) {
        const previousFocusedRowIndex = this.focusedRowIndex;
        this._focusedRowIndex = focusedRowIndex;
        if (previousFocusedRowIndex != focusedRowIndex) {
            this._tableView.visibleRowWithIndex(previousFocusedRowIndex).focused = NO;
            const focusedRow = this._tableView.visibleRowWithIndex(this.focusedRowIndex);
            focusedRow.focused = YES;
            if (!this.keepFocusedRowVisible) {
                return;
            }
            var contentOffset = this._tableView.contentOffset;
            if (focusedRow.frame.y < contentOffset.y) {
                contentOffset.y = focusedRow.frame.y;
            }
            if (focusedRow.frame.max.y > (contentOffset.y + this._tableView.bounds.height)) {
                contentOffset = contentOffset.pointByAddingY(-(contentOffset.y + this._tableView.bounds.height -
                    focusedRow.frame.max.y));
            }
            const animationDuration = this._tableView.animationDuration;
            this._tableView.animationDuration = 0;
            this._tableView.contentOffset = contentOffset;
            this._tableView.animationDuration = animationDuration;
        }
    }
    set expandedContainerViewHeight(expandedContainerViewHeight) {
        this._expandedContainerViewHeight = expandedContainerViewHeight;
        this._dialogView.setNeedsLayout();
    }
    get expandedContainerViewHeight() {
        if (IS(this._expandedContainerViewHeight)) {
            return this._expandedContainerViewHeight;
        }
        const padding = RootViewController.paddingLength;
        const labelHeight = padding;
        const result = this.superview.bounds.height - this.frame.max.y - padding;
        return result;
    }
    selectedDataContains(datapoint) {
        for (var i = 0; i < this.selectedData.length; i++) {
            const value = this.selectedData[i];
            if (value._id == datapoint._id) {
                return YES;
            }
        }
        return NO;
    }
    updateContentForNormalState() {
        this.style.borderBottom = "1px solid rgba(0,0,0,0.12)";
        this.titleLabel.textColor = CBColor.primaryContentColor;
        this.backgroundColor = UIColor.transparentColor;
        this.style.borderBottomColor = CBColor.primaryContentColor.colorWithAlpha(0.12).stringValue;
    }
    updateContentForHighlightedState() {
        this.style.borderBottomWidth = "2px";
        this.style.borderBottomColor = this.tintColor.stringValue;
    }
    selectionDidChange(selectedRows) {
        // Update title to indicate the selection
        this.updateTitleWithSelection(selectedRows);
        this.sendControlEventForKey(SearchableDropdown.controlEvent.SelectionDidChange, nil);
    }
    updateContentForCurrentSelection() {
        // Update selected rows in the top of the table view
        this._tableView.reloadData();
        this.setNeedsLayout();
    }
    get placeholderText() {
        if (IS_UNDEFINED(this._placeholderText)) {
            this._placeholderText = "Not selected";
        }
        return this._placeholderText;
    }
    set placeholderText(placeholderText) {
        this._placeholderText = placeholderText;
        this.updateTitleWithSelection(this.selectedData);
    }
    setPlaceholderText(key, defaultString, parameters) {
        this.placeholderLocalizedTextObject = LanguageService.localizedTextObjectForKey(key, defaultString, parameters);
        //this.titleLabel.setText(key, defaultString, parameters);
    }
    get placeholderLocalizedTextObject() {
        if (IS_UNDEFINED(this._placeholderLocalizedTextObject)) {
            this._placeholderLocalizedTextObject = LanguageService.localizedTextObjectForText("Not selected");
        }
        return this._placeholderLocalizedTextObject;
    }
    set placeholderLocalizedTextObject(placeholderText) {
        this._placeholderLocalizedTextObject = placeholderText;
        this.updateTitleWithSelection(this.selectedData);
    }
    updateTitleWithSelection(selectedRows) {
        this.titleLabel.localizedTextObject = this.placeholderLocalizedTextObject;
        if (selectedRows && selectedRows.length) {
            const maxWidth = this.titleLabel.bounds.width;
            this.titleLabel.localizedTextObject = nil;
            this.titleLabel.text = "";
            var stopLooping = NO;
            selectedRows.forEach(function (selectedDatapoint, index, array) {
                if (stopLooping) {
                    return;
                }
                var selectedString = LanguageService.stringForCurrentLanguage(selectedDatapoint.title);
                if (index) {
                    selectedString = ", " + selectedString;
                }
                const previousText = this.titleLabel.text;
                this.titleLabel.text = this.titleLabel.text + selectedString;
                this.overflowLabel.text = "+" + (array.length - index - 1);
                if (index == array.length - 1) {
                    this.overflowLabel.text = "";
                }
                if (index && (this.bounds.width -
                    (this.overflowLabel.intrinsicContentWidth() + this.titleLabel.intrinsicContentWidth()) - 20) < 0) {
                    this.titleLabel.text = previousText;
                    this.overflowLabel.text = "+" + (array.length - index - 2 * 0);
                    stopLooping = YES;
                }
            }, this);
        }
    }
    updateFilteredData(filteringString) {
        this._filteredData = [];
        this.data.forEach(function (sectionOrRow, index, array) {
            if (LanguageService.stringForCurrentLanguage(sectionOrRow.title).toLowerCase().contains(filteringString.toLowerCase())) {
                this.filteredData.push(sectionOrRow);
            }
            else if (sectionOrRow.isADropdownDataSection) {
                this._keyValueStringFilter.filterData(filteringString, sectionOrRow.rowsData, this._excludedData, "title." + LanguageService.currentLanguageKey, sectionOrRow, function (filteredData, filteredIndexes, sectionFromThread) {
                    if (filteredData.length) {
                        this.filteredData.push({
                            _id: sectionFromThread._id,
                            title: sectionFromThread.title,
                            rowsData: filteredData,
                            isADropdownDataSection: sectionFromThread.isADropdownDataSection,
                            isADropdownDataRow: sectionFromThread.isADropdownDataRow,
                            attachedObject: sectionFromThread.attachedObject,
                            itemCode: sectionFromThread.itemCode,
                            dropdownCode: sectionFromThread.dropdownCode
                        });
                        if (this.allowsCustomItem && this._searchTextField.text && this._customItem) {
                            this.filteredData.removeElement(this._customItem);
                            this.filteredData.push(this._customItem);
                        }
                        this._isDrawingDataValid = NO;
                        this._tableView.reloadData();
                    }
                }.bind(this));
            }
        }.bind(this));
        if (this.allowsCustomItem && this._searchTextField.text) {
            this.filteredData.removeElement(this._customItem);
            this.initCustomItemWithTitle(this._searchTextField.text);
            this.filteredData.push(this._customItem);
        }
        if (this.filteredData.length) {
            this.focusedRowIndex = 0;
        }
        else {
            this.focusedRowIndex = null;
        }
        this._isDrawingDataValid = NO;
        this._tableView.reloadData();
    }
    initCustomItemWithTitle(title) {
        if (IS_NOT(title)) {
            this._customItem = undefined;
        }
        else {
            this._customItem = {
                _id: "" + MAKE_ID(),
                title: LanguageService.localizedTextObjectForText(title),
                rowsData: [],
                isADropdownDataSection: NO,
                isADropdownDataRow: YES,
                attachedObject: undefined,
                itemCode: "custom_item",
                dropdownCode: null
            };
        }
    }
    selectItemOrCustomItemWithTitle(title) {
        if (IS_NOT(title)) {
            this._customItem = undefined;
        }
        else {
            var item = this.drawingData.find(function (dataItem, index, array) {
                return (LanguageService.stringForCurrentLanguage(dataItem.title) == title);
            });
            if (this.allowsCustomItem && IS_NOT(item)) {
                this._searchTextField.text = title;
                this.updateFilteredData(title);
                item = this._customItem;
            }
            if (IS_NOT(this.isSingleSelection)) {
                if (IS_NOT(this.selectedDataContains(item))) {
                    const selectedItemCodes = this.selectedItemCodes.copy();
                    selectedItemCodes.push(item.itemCode);
                    this.selectedItemCodes = selectedItemCodes;
                }
                return;
            }
            this.selectedItemCodes = [item.itemCode];
        }
    }
    set data(data) {
        this._data = data;
        this.updateFilteredData(this._searchTextField.text);
    }
    get data() {
        return this._data;
    }
    set filteredData(data) {
        this._filteredData = data;
        this._isDrawingDataValid = NO;
    }
    get filteredData() {
        return this._filteredData;
    }
    get drawingData() {
        if (this._isDrawingDataValid) {
            return this._drawingData;
        }
        const result = [];
        this._filteredData.forEach(function (section) {
            result.push({
                _id: section._id,
                title: section.title,
                // @ts-ignore
                rowsData: null,
                isADropdownDataSection: section.isADropdownDataSection,
                isADropdownDataRow: NO,
                attachedObject: section.attachedObject,
                itemCode: section.itemCode,
                dropdownCode: section.dropdownCode
            });
            if (section.rowsData) {
                section.rowsData.forEach(function (row) {
                    result.push(row);
                }.bind(this));
            }
        }.bind(this));
        this._drawingData = result;
        this._isDrawingDataValid = YES;
        return result;
    }
    get selectedData() {
        return this._selectedData;
    }
    set selectedData(selectedData) {
        this._selectedData = selectedData;
    }
    clearSelection() {
        this.selectedData = [];
        this.selectedIndices = [];
        this.updateTitleWithSelection(this.selectedData);
        this.updateContentForCurrentSelection();
        this.selectionDidChange(this.selectedData);
    }
    get selectedItemCodes() {
        return this.selectedData.map(function (item) {
            return item.itemCode;
        });
    }
    set selectedItemCodes(selectedItemCodes) {
        const selectedData = [];
        const selectedIndices = [];
        this._drawingData.forEach(function (item, index, array) {
            if (selectedItemCodes.contains(item.itemCode)) {
                selectedData.push(item);
                selectedIndices.push(index);
            }
        });
        this.selectedData = selectedData;
        this.selectedIndices = selectedIndices;
        this.updateTitleWithSelection(this.selectedData);
        this.updateContentForCurrentSelection();
        this.selectionDidChange(this.selectedData);
    }
    get selectedRowIdentifiers() {
        const result = [];
        this.selectedData.forEach(function (selectedDatapoint) {
            result.push(selectedDatapoint._id);
        }.bind(this));
        return result;
    }
    wasAddedToViewTree() {
        super.wasAddedToViewTree();
        this.setNeedsLayout();
    }
    layoutSubviews() {
        super.layoutSubviews();
        const bounds = this.bounds;
        const padding = RootViewController.paddingLength;
        const labelHeight = padding;
        this.updateTitleWithSelection(this.selectedData);
        if (this._rightImageView) {
            // var imageHeight = bounds.height - padding;
            // this._rightImageView.frame = new UIRectangle(bounds.width - imageHeight - padding * 0.5, padding * 0.5, imageHeight, imageHeight);
            this._rightImageView.frame = bounds.rectangleWithInsets(0, padding *
                0.5, 0, 0).rectangleWithWidth(24, 1).rectangleWithHeight(24, 0.5);
        }
        if (this.overflowLabel) {
            //this.overflowLabel.frame = bounds.rectangleWithInsets(padding, padding, 0, 0).rectangleWithInsets(0, this._rightImageView.bounds.width || 0, 0, 0).rectangleWithWidth(36, 1);
            //this.titleLabel.style.maxWidth = "" + (bounds.width - (bounds.width - this.titleLabel.viewHTMLElement.offsetLeft - (this.overflowLabel.frame.x || 0)))
            this.overflowLabel.style.width = "36px";
            this.overflowLabel.style.right = "32px";
            this.overflowLabel.centerYInContainer();
            this.titleLabel.style.marginRight = "60px";
            //this.overflowLabel.centerYInContainer();
            //this.titleLabel.frame = this.titleLabel.frame.rectangleWithWidth(bounds.width - (bounds.width - this.titleLabel.frame.x - (this.overflowLabel.frame.x || 0)))
        }
    }
}
SearchableDropdown.controlEvent = Object.assign({}, UIView.controlEvent, {
    "SelectionDidChange": "SelectionDidChange"
});
///<reference path="UIView.ts"/>
class UIActionIndicator extends UIView {
    constructor(elementID) {
        super(elementID);
        this._size = 50;
        this._class = UIActionIndicator;
        this.superclass = UIView;
        this.indicatorView = new UIView(this.elementID + "IndicatorView");
        this.addSubview(this.indicatorView);
        this.indicatorView.viewHTMLElement.classList.add("LukeHaasLoader");
        this.hidden = YES;
    }
    set size(size) {
        this._size = size;
        this.setNeedsLayoutUpToRootView();
    }
    get size() {
        return this._size;
    }
    set hidden(hidden) {
        super.hidden = hidden;
        if (hidden) {
            this.indicatorView.removeFromSuperview();
        }
        else {
            this.addSubview(this.indicatorView);
        }
    }
    start() {
        this.hidden = NO;
    }
    stop() {
        this.hidden = YES;
    }
    layoutSubviews() {
        super.layoutSubviews();
        const bounds = this.bounds;
        //this.indicatorView.centerInContainer();
        this.indicatorView.style.height = "" + this._size.integerValue + "px";
        this.indicatorView.style.width = "" + this._size.integerValue + "px";
        const minSize = Math.min(this.bounds.height, this.bounds.width);
        this.indicatorView.style.maxHeight = "" + minSize.integerValue + "px";
        this.indicatorView.style.maxWidth = "" + minSize.integerValue + "px";
        const size = Math.min(this._size, minSize);
        this.indicatorView.style.left = "" + ((bounds.width - size) * 0.5 - 11).integerValue + "px";
        this.indicatorView.style.top = "" + ((bounds.height - size) * 0.5 - 11).integerValue + "px";
    }
}
///<reference path="UIView.ts"/>
class UIDateTimeInput extends UIView {
    constructor(elementID, type = UIDateTimeInput.type.DateTime) {
        super(elementID, nil, "input");
        this._class = UIDateTimeInput;
        this.superclass = UIView;
        this.viewHTMLElement.setAttribute("type", type);
        this.viewHTMLElement.onchange = (event) => {
            this.sendControlEventForKey(UIDateTimeInput.controlEvent.ValueChange, event);
        };
    }
    get addControlEventTarget() {
        return super.addControlEventTarget;
    }
    get date() {
        const result = new Date(this.viewHTMLElement.value);
        return result;
    }
}
UIDateTimeInput.controlEvent = Object.assign({}, UIView.controlEvent, {
    "ValueChange": "ValueChange"
});
UIDateTimeInput.type = {
    "Date": "date",
    "Time": "time",
    "DateTime": "datetime"
};
UIDateTimeInput.format = {
    "European": "DD-MM-YYYY",
    "ISOComputer": "YYYY-MM-DD",
    "American": "MM/DD/YYYY"
};
///<reference path="UIObject.ts"/>
class UIKeyValueStringSorter extends UIObject {
    constructor(useSeparateWebWorker = NO) {
        super();
        this._isThreadClosed = NO;
        this._webWorker = UIKeyValueStringSorter._sharedWebWorker;
        this._class = UIKeyValueStringSorter;
        this.superclass = UIObject;
        if (useSeparateWebWorker) {
            this._webWorker = new Worker("compiledScripts/UIKeyValueStringSorterWebWorker.js");
        }
        UIKeyValueStringSorter._instanceNumber = UIKeyValueStringSorter._instanceNumber + 1;
        this._instanceNumber = UIKeyValueStringSorter._instanceNumber;
    }
    get instanceIdentifier() {
        return this._instanceNumber;
    }
    sortData(data, sortingInstructions, identifier, completion) {
        if (this._isThreadClosed) {
            return;
        }
        const startTime = Date.now();
        const instanceIdentifier = this.instanceIdentifier;
        this._webWorker.onmessage = function (message) {
            if (message.data.instanceIdentifier == instanceIdentifier) {
                console.log("Sorting " + data.length + " items took " + (Date.now() - startTime) + " ms.");
                completion(message.data.sortedData, message.data.sortedIndexes, message.data.identifier);
            }
        };
        try {
            this._webWorker.postMessage({
                "data": data,
                "sortingInstructions": sortingInstructions,
                "identifier": identifier,
                "instanceIdentifier": instanceIdentifier
            });
        }
        catch (exception) {
            completion([], [], identifier);
        }
    }
    sortedData(data, sortingInstructions, identifier = MAKE_ID()) {
        const result = new Promise((resolve, reject) => {
            this.sortData(data, sortingInstructions, identifier, (sortedData, sortedIndexes, sortedIdentifier) => {
                if (sortedIdentifier == identifier) {
                    resolve({
                        sortedData: sortedData,
                        sortedIndexes: sortedIndexes,
                        identifier: sortedIdentifier
                    });
                }
            });
        });
        return result;
    }
    closeThread() {
        this._isThreadClosed = YES;
        if (this._webWorker != UIKeyValueStringSorter._sharedWebWorker) {
            this._webWorker.terminate();
        }
    }
}
UIKeyValueStringSorter._sharedWebWorker = new Worker("compiledScripts/UIKeyValueStringSorterWebWorker.js");
UIKeyValueStringSorter._instanceNumber = -1;
UIKeyValueStringSorter.dataType = {
    "string": "string"
};
UIKeyValueStringSorter.direction = {
    "descending": "descending",
    "ascending": "ascending"
};
///<reference path="UIObject.ts"/>
class UILayoutGrid extends UIObject {
    constructor(frame) {
        super();
        this._subframes = [];
        this._class = UILayoutGrid;
        this.superclass = UIObject;
        this._frame = frame;
    }
    splitXInto(numberOfFrames) {
        if (this._subframes.length == 0) {
            for (var i = 0; i < numberOfFrames; i++) {
                const asd = 1;
            }
        }
    }
}
/// <reference path="./UIView.ts" />
class UIScrollView extends UIView {
    constructor(elementID, viewHTMLElement) {
        super(elementID, viewHTMLElement);
        this._contentOffset = new UIPoint(0, 0);
        this._contentScale = 1;
        this._scrollEnabled = YES;
        this._class = UIScrollView;
        this.superclass = UIView;
        this.containerView = new UIView(elementID + "ContainerView");
        super.addSubview(this.containerView);
        this.style.overflow = "hidden";
        this.pausesPointerEvents = NO; //YES;
        this.addTargetForControlEvent(UIView.controlEvent.PointerDown, function () {
            this._pointerDown = YES;
        }.bind(this));
        this.addTargetForControlEvent(UIView.controlEvent.PointerUp, function () {
            this._pointerDown = NO;
            this._previousClientPoint = null;
            scrollStopped();
        }.bind(this));
        function scrollStopped() {
            // Handle paging if needed
        }
        this.addTargetForControlEvent(UIView.controlEvent.PointerMove, function (sender, event) {
            if (!(this._pointerDown && this._scrollEnabled && this._enabled)) {
                return;
            }
            const currentClientPoint = new UIPoint(nil, nil);
            if (window.MouseEvent && event instanceof MouseEvent) {
                currentClientPoint.x = event.clientX;
                currentClientPoint.y = event.clientY;
            }
            if (window.TouchEvent && event instanceof TouchEvent) {
                const touchEvent = event;
                if (touchEvent.touches.length != 1) {
                    this._pointerDown = NO;
                    this._previousClientPoint = null;
                    scrollStopped();
                    return;
                }
                currentClientPoint.x = touchEvent.touches[0].clientX;
                currentClientPoint.y = touchEvent.touches[0].clientY;
            }
            if (!this._previousClientPoint) {
                this._previousClientPoint = currentClientPoint;
                return;
            }
            const changePoint = currentClientPoint.copy().subtract(this._previousClientPoint);
            if (this.containerView.bounds.width <= this.bounds.width) {
                changePoint.x = 0;
            }
            if (0 < this.contentOffset.x + changePoint.x) {
                changePoint.x = -this.contentOffset.x;
            }
            if (this.contentOffset.x + changePoint.x < -this.bounds.width) {
                changePoint.x = -this.bounds.width - this.contentOffset.x;
            }
            if (this.containerView.bounds.height <= this.bounds.height) {
                changePoint.y = 0;
            }
            if (0 < this.contentOffset.y + changePoint.y) {
                changePoint.y = -this.contentOffset.y;
            }
            if (this.contentOffset.y + changePoint.y < -this.bounds.height) {
                changePoint.y = -this.bounds.height - this.contentOffset.y;
            }
            this.contentOffset = this.contentOffset.add(changePoint);
            this._previousClientPoint = currentClientPoint;
        }.bind(this));
    }
    invalidateIntrinsicContentFrame() {
        this._intrinsicContentFrame = nil;
    }
    get contentOffset() {
        return this._contentOffset;
    }
    set contentOffset(offset) {
        this._contentOffset = offset;
        this.setNeedsLayout();
    }
    layoutSubviews() {
        super.layoutSubviews();
        // var intrinsicContentFrame = this._intrinsicContentFrame;
        // if (!IS(intrinsicContentFrame)) {
        //     intrinsicContentFrame = this.containerView.intrinsicContentSizeWithConstraints();   
        // }
        // intrinsicContentFrame.offsetByPoint(this.contentOffset);
        // intrinsicContentFrame.height = this.containerView.viewHTMLElement.scrollHeight;
        // intrinsicContentFrame.width = this.containerView.viewHTMLElement.scrollWidth;
        // this.containerView.frame = intrinsicContentFrame;
        this.containerView.frame = this.containerView.bounds.offsetByPoint(this.contentOffset);
    }
    // get _subviews() {
    //     return super.subviews;
    // }
    // set _subviews(subviews: UIView[]) {
    //     super.subviews = subviews;
    // }
    // get subviews() {
    //     return this.containerView.subviews;
    // }
    // set subviews(subviews: UIView[]) {
    //     this.containerView.subviews = subviews;
    //     this.invalidateIntrinsicContentFrame();
    // }
    hasSubview(view) {
        return this.containerView.hasSubview(view);
    }
    addSubview(view) {
        this.containerView.addSubview(view);
        this.invalidateIntrinsicContentFrame();
    }
}
///<reference path="UIView.ts"/>
class UISlideScrollerView extends UIView {
    constructor(elementID, viewHTMLElement) {
        super(elementID, viewHTMLElement);
        this._targetIndex = 0;
        this._isAnimating = NO;
        this._isAnimationOngoing = NO;
        this._animationTimer = nil;
        this._slideViews = [];
        this.wrapAround = YES;
        this.animationDuration = 0.35;
        this.animationDelay = 2;
        this._currentPageIndex = 0;
        this._class = UIScrollView;
        this.superclass = UIView;
        this._scrollView = new UIScrollView(elementID + "ScrollView");
        this.addSubview(this._scrollView);
        this._scrollView._scrollEnabled = NO;
        this._scrollView.addTargetForControlEvent(UIView.controlEvent.PointerMove, function (sender, event) {
            if (event instanceof MouseEvent) {
                this._animationTimer.invalidate();
            }
        }.bind(this));
        this._scrollView.addTargetForControlEvent(UIView.controlEvent.PointerLeave, function () {
            if (this._isAnimating && event instanceof MouseEvent) {
                this.startAnimating();
            }
        }.bind(this));
        // Touch events
        this._scrollView.addTargetForControlEvent(UIView.controlEvent.PointerDown, function (sender, event) {
            if (event instanceof TouchEvent) {
                this._animationTimer.invalidate();
            }
        }.bind(this));
        this._scrollView.addTargetForControlEvents([
            UIView.controlEvent.PointerUp, UIView.controlEvent.PointerCancel
        ], function (sender, event) {
            if (event instanceof TouchEvent && this._isAnimating) {
                this.startAnimating();
            }
        }.bind(this));
        // Page indicator
        this.pageIndicatorsView = new UIView(elementID + "PageIndicatorsView");
        this.addSubview(this.pageIndicatorsView);
    }
    buttonForPageIndicatorWithIndex(index) {
        const result = new UIButton(this.viewHTMLElement.id + "PageIndicatorButton" + index);
        result.addTargetForControlEvents([
            UIView.controlEvent.PointerUpInside, UIView.controlEvent.EnterUp
        ], function (sender, event) {
            this.scrollToPageWithIndex(index, YES);
            if (this._isAnimating) {
                this.startAnimating();
            }
        }.bind(this));
        result.addTargetForControlEvent(UIView.controlEvent.PointerMove, function () {
            this._animationTimer.invalidate();
        }.bind(this));
        result.updateContentForNormalState = function () {
            result.backgroundColor = UIColor.blueColor;
            result.titleLabel.textColor = UIColor.whiteColor;
        };
        result.frame = new UIRectangle(nil, nil, 20, 50);
        // result.style.height = "20px";
        // result.style.width = "50px";
        result.style.display = "table-cell";
        result.style.position = "relative";
        // var resultContent = new UIView(result.viewHTMLElement.id + "Content");
        // resultContent.backgroundColor = UIColor.whiteColor;
        // resultContent.centerYInContainer();
        // resultContent.style.height = "10px";
        // resultContent.style.height = "100%";
        // resultContent.style.borderRadius = "5px";
        // result.addSubview(resultContent);
        return result;
    }
    addSlideView(view) {
        this.slideViews.push(view);
        this.updateSlideViews();
    }
    set slideViews(views) {
        this._slideViews = views;
        this.updateSlideViews();
    }
    get slideViews() {
        return this._slideViews;
    }
    get currentPageIndex() {
        const result = this._currentPageIndex;
        return result;
    }
    set currentPageIndex(index) {
        this._currentPageIndex = index;
        this._slideViews[index].willAppear();
        //this._scrollView.contentOffset.x = -this._slideViews[index].frame.min.x; //-this.bounds.width * index;
        //this._scrollView.contentOffset.x = Math.round(this._scrollView.contentOffset.x);
        this._scrollView.contentOffset = this._scrollView.contentOffset.pointWithX(-this._slideViews[index].frame.min.x);
        this.pageIndicatorsView.subviews.forEach(function (button, index, array) {
            button.selected = NO;
        });
        this.pageIndicatorsView.subviews[index].selected = YES;
    }
    scrollToPreviousPage(animated) {
        if (this.slideViews.length == 0) {
            return;
        }
        var targetIndex = this.currentPageIndex;
        if (this.wrapAround) {
            targetIndex = (this.currentPageIndex - 1) % (this.slideViews.length);
        }
        else if (this.currentPageIndex - 1 < this.slideViews.length) {
            targetIndex = this.currentPageIndex - 1;
        }
        else {
            return;
        }
        this.scrollToPageWithIndex(targetIndex, animated);
    }
    scrollToNextPage(animated) {
        if (this.slideViews.length == 0) {
            return;
        }
        var targetIndex = this.currentPageIndex;
        if (this.wrapAround) {
            targetIndex = (this.currentPageIndex + 1) % (this.slideViews.length);
        }
        else if (this.currentPageIndex + 1 < this.slideViews.length) {
            targetIndex = this.currentPageIndex + 1;
        }
        else {
            return;
        }
        this.scrollToPageWithIndex(targetIndex, animated);
    }
    scrollToPageWithIndex(targetIndex, animated = YES) {
        this._targetIndex = targetIndex;
        // this._slideViews[this.currentPageIndex]._shouldLayout = NO;
        // this._slideViews[this._targetIndex]._shouldLayout = YES;
        //this._slideViews[this._targetIndex].hidden = NO;
        this.willScrollToPageWithIndex(targetIndex);
        this._isAnimationOngoing = YES;
        //var previousView = this._slideViews[this.currentPageIndex];
        if (animated) {
            UIView.animateViewOrViewsWithDurationDelayAndFunction(this._scrollView.containerView, this.animationDuration, 0, undefined, function () {
                this.currentPageIndex = targetIndex;
            }.bind(this), function () {
                this.didScrollToPageWithIndex(targetIndex);
                this._isAnimationOngoing = NO;
                //previousView.hidden = YES;
            }.bind(this));
        }
        else {
            this.currentPageIndex = targetIndex;
            this.didScrollToPageWithIndex(targetIndex);
            //previousView.hidden = YES;
        }
    }
    willScrollToPageWithIndex(index) {
        const targetView = this.slideViews[index];
        if (IS(targetView) && targetView.willAppear && targetView.willAppear instanceof Function) {
            targetView.willAppear();
        }
    }
    didScrollToPageWithIndex(index) {
    }
    startAnimating() {
        this._isAnimating = YES;
        this._animationTimer.invalidate();
        this._animationTimer = new UITimer(this.animationDelay + this.animationDuration, YES, function () {
            this.scrollToNextPage(YES);
        }.bind(this));
    }
    stopAnimating() {
        this._isAnimating = NO;
        this._animationTimer.invalidate();
    }
    updateSlideViews() {
        this._scrollView.containerView.subviews.slice().forEach(function (subview, index, array) {
            subview.removeFromSuperview();
        });
        this.pageIndicatorsView.subviews.slice().forEach(function (subview, index, array) {
            subview.removeFromSuperview();
        });
        this._slideViews.forEach(function (view, index, array) {
            this._scrollView.addSubview(view);
            this.pageIndicatorsView.addSubview(this.buttonForPageIndicatorWithIndex(index));
        }.bind(this));
    }
    didReceiveBroadcastEvent(event) {
        super.didReceiveBroadcastEvent(event);
        if (event.name == UICore.broadcastEventName.WindowDidResize) {
            this.currentPageIndex = this.currentPageIndex;
        }
    }
    set frame(frame) {
        super.frame = frame;
        this.currentPageIndex = this.currentPageIndex;
    }
    get frame() {
        return super.frame;
    }
    layoutSubviews() {
        super.layoutSubviews();
        if (this.bounds.isEqualTo(this._previousLayoutBounds)) {
            return;
        }
        const bounds = this.bounds;
        this._previousLayoutBounds = bounds;
        this._scrollView.frame = bounds;
        this._scrollView.containerView.frame = bounds.rectangleWithWidth(bounds.width *
            this.slideViews.length).performFunctionWithSelf(function (self) {
            self.offsetByPoint(this._scrollView.contentOffset);
            return self;
        }.bind(this));
        this._slideViews.forEach(function (view, index, array) {
            view.frame = bounds.rectangleWithX((this.bounds.width + 1) * index);
        }.bind(this));
        this.layoutPageIndicators();
    }
    layoutPageIndicators() {
        this.pageIndicatorsView.centerXInContainer();
        this.pageIndicatorsView.style.bottom = "20px";
        this.pageIndicatorsView.style.height = "20px";
        this.pageIndicatorsView.style.display = "table-row";
    }
    removeFromSuperview() {
        super.removeFromSuperview();
        this.stopAnimating();
    }
}
///<reference path="UIObject.ts"/>
class UIStringFilter extends UIObject {
    constructor(useSeparateWebWorker = NO) {
        super();
        this._isThreadClosed = NO;
        this._webWorker = UIStringFilter._sharedWebWorker;
        this._class = UIStringFilter;
        this.superclass = UIObject;
        if (useSeparateWebWorker) {
            this._webWorker = new Worker("compiledScripts/UIStringFilterWebWorker.js");
        }
        UIStringFilter._instanceNumber = UIStringFilter._instanceNumber + 1;
        this._instanceNumber = UIStringFilter._instanceNumber;
    }
    get instanceIdentifier() {
        return this._instanceNumber;
    }
    filterData(filteringString, data, excludedData, identifier, completion) {
        if (this._isThreadClosed) {
            return;
        }
        //var startTime = Date.now();
        const instanceIdentifier = this.instanceIdentifier;
        this._webWorker.onmessage = function (message) {
            if (message.data.instanceIdentifier == instanceIdentifier) {
                //console.log("Filtering took " + (Date.now() - startTime) + " ms");
                completion(message.data.filteredData, message.data.filteredIndexes, message.data.identifier);
            }
        };
        this._webWorker.postMessage({
            "filteringString": filteringString,
            "data": data,
            "excludedData": excludedData,
            "identifier": identifier,
            "instanceIdentifier": instanceIdentifier
        });
    }
    filteredData(filteringString, data, excludedData = [], identifier = MAKE_ID()) {
        const result = new Promise((resolve, reject) => {
            this.filterData(filteringString, data, excludedData, identifier, (filteredData, filteredIndexes, filteredIdentifier) => {
                if (filteredIdentifier == identifier) {
                    resolve({
                        filteredData: filteredData,
                        filteredIndexes: filteredIndexes,
                        identifier: filteredIdentifier
                    });
                }
            });
        });
        return result;
    }
    closeThread() {
        this._isThreadClosed = YES;
        if (this._webWorker != UIStringFilter._sharedWebWorker) {
            this._webWorker.terminate();
        }
    }
}
UIStringFilter._sharedWebWorker = new Worker("compiledScripts/UIStringFilterWebWorker.js");
UIStringFilter._instanceNumber = -1;
/// <reference path="UITextField.ts" />
class UITextArea extends UITextField {
    constructor(elementID, viewHTMLElement = null) {
        super(elementID, viewHTMLElement, UITextView.type.textArea);
        this._class = UITextArea;
        this.superclass = UITextField;
        this.viewHTMLElement.removeAttribute("type");
        this.style.overflow = "auto";
        this.style.webkitUserSelect = "text";
        this.pausesPointerEvents = NO;
    }
    get addControlEventTarget() {
        return super.addControlEventTarget;
    }
    // @ts-ignore
    get viewHTMLElement() {
        // @ts-ignore
        return super.viewHTMLElement;
    }
}
///<reference path="UIObject.ts"/>
class UITimer extends UIObject {
    constructor(interval, repeats, target) {
        super();
        this.interval = interval;
        this.repeats = repeats;
        this.target = target;
        this.isValid = YES;
        this.superclass = UIObject;
        this.schedule();
    }
    schedule() {
        const callback = function () {
            if (this.repeats == NO) {
                this.invalidate();
            }
            this.target();
        }.bind(this);
        this._intervalID = window.setInterval(callback, this.interval * 1000);
    }
    reschedule() {
        this.invalidate();
        this.schedule();
    }
    fire() {
        if (this.repeats == NO) {
            this.invalidate();
        }
        else {
            this.reschedule();
        }
        this.target();
    }
    invalidate() {
        if (this.isValid) {
            clearInterval(this._intervalID);
            this.isValid = NO;
        }
    }
}

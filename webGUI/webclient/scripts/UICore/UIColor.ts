///<reference path="UIObject.ts"/>





interface UIColorDescriptor {
    
    red: number;
    green: number;
    blue: number;
    alpha?: number;
    
}





class UIColor extends UIObject {
    
    
    constructor(public stringValue: string) {
        
        super()
        
        this._class = UIColor
        this.superclass = UIObject
        
        
        
    }
    
    
    
    toString() {
        return this.stringValue
    }
    
    static get redColor() {
        return new UIColor("red")
    }
    
    static get blueColor() {
        return new UIColor("blue")
    }
    
    static get greenColor() {
        return new UIColor("green")
    }
    
    static get yellowColor() {
        return new UIColor("yellow")
    }
    
    static get blackColor() {
        return new UIColor("black")
    }
    
    static get brownColor() {
        return new UIColor("brown")
    }
    
    static get whiteColor() {
        return new UIColor("white")
    }
    
    static get greyColor() {
        return new UIColor("grey")
    }
    
    static get lightGreyColor() {
        return new UIColor("lightgrey")
    }
    
    static get transparentColor() {
        return new UIColor("transparent")
    }
    
    static get undefinedColor() {
        return new UIColor("")
    }
    
    static get nilColor() {
        return new UIColor("")
    }
    
    
    
    static nameToHex(name: string) {
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
        }[name.toLowerCase()]
    }
    
    
    
    static hexToDescriptor(c: string): UIColorDescriptor {
        if (c[0] === "#") {
            c = c.substr(1)
        }
        const r = parseInt(c.slice(0, 2), 16)
        const g = parseInt(c.slice(2, 4), 16)
        const b = parseInt(c.slice(4, 6), 16)
        const a = parseInt(c.slice(6, 8), 16)
        
        const result = { "red": r, "green": g, "blue": b, "alpha": a }
        
        return result
        
        //return 'rgb(' + r + ',' + g + ',' + b + ')';
        
    }
    
    static rgbToDescriptor(colorString: string) {
        
        
        if (colorString.startsWith("rgba(")) {
            
            colorString = colorString.slice(5, colorString.length - 1)
            
        }
        
        if (colorString.startsWith("rgb(")) {
            
            colorString = colorString.slice(4, colorString.length - 1) + ", 0"
            
        }
        
        
        const components = colorString.split(",")
        
        
        
        const result = {
            "red": Number(components[0]),
            "green": Number(components[1]),
            "blue": Number(components[2]),
            "alpha": Number(components[3])
        }
        
        
        return result
        
        
    }
    
    
    
    get colorDescriptor(): UIColorDescriptor {
        
        var descriptor
        
        const colorHEXFromName = UIColor.nameToHex(this.stringValue)
        
        if (this.stringValue.startsWith("rgb")) {
            
            descriptor = UIColor.rgbToDescriptor(this.stringValue)
            
        }
        else if (colorHEXFromName) {
            
            descriptor = UIColor.hexToDescriptor(colorHEXFromName)
            
        }
        else {
            
            descriptor = UIColor.hexToDescriptor(this.stringValue)
            
        }
        
        return descriptor
        
    }
    
    
    
    colorWithRed(red: number) {
        
        
        const descriptor = this.colorDescriptor
        
        const result = new UIColor("rgba(" + red + "," + descriptor.green + "," + descriptor.blue + "," +
            descriptor.alpha + ")")
        
        return result
        
    }
    
    colorWithGreen(green: number) {
        
        
        const descriptor = this.colorDescriptor
        
        const result = new UIColor("rgba(" + descriptor.red + "," + green + "," + descriptor.blue + "," +
            descriptor.alpha + ")")
        
        return result
        
    }
    
    colorWithBlue(blue: number) {
        
        
        const descriptor = this.colorDescriptor
        
        const result = new UIColor("rgba(" + descriptor.red + "," + descriptor.green + "," + blue + "," +
            descriptor.alpha + ")")
        
        return result
        
    }
    
    
    colorWithAlpha(alpha: number) {
        
        
        const descriptor = this.colorDescriptor
        
        const result = new UIColor("rgba(" + descriptor.red + "," + descriptor.green + "," + descriptor.blue + "," +
            alpha + ")")
        
        return result
        
    }
    
    
    
    static colorWithRGBA(red: number, green: number, blue: number, alpha: number = 1) {
        
        
        const result = new UIColor("rgba(" + red + "," + green + "," + blue + "," + alpha + ")")
        
        return result
        
        
    }
    
    static colorWithDescriptor(descriptor: UIColorDescriptor) {
        
        
        const result = new UIColor("rgba(" + descriptor.red.toFixed(0) + "," + descriptor.green.toFixed(0) + "," +
            descriptor.blue.toFixed(0) + "," + this.defaultAlphaToOne(descriptor.alpha) + ")")
        
        return result
        
    }
    
    
    private static defaultAlphaToOne(value = 1) {
        if (value != value) {
            value = 1
        }
        return value
    }
    
    
    
    colorByMultiplyingRGB(multiplier: number) {
        
        const descriptor = this.colorDescriptor
        
        descriptor.red = descriptor.red * multiplier
        descriptor.green = descriptor.green * multiplier
        descriptor.blue = descriptor.blue * multiplier
        
        const result = UIColor.colorWithDescriptor(descriptor)
        
        return result
        
    }
    
    
    
    
    
}















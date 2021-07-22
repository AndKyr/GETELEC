/// <reference path="../UICore/UIColor.ts" />





class CBColor extends UIColor {
    
    
    
    
    
    constructor(stringValue: string) {
        
        super(stringValue)
        
        
        this._class = CBColor
        this.superclass = UIColor
        
        
        
    }
    
    
    
    
    
    static get primaryTintColor() {
        
        return new CBColor("rgba(52, 127, 230, 1)")
        
    }
    
    static get secondaryTintColor() {
        
        return new CBColor("rgba(0, 196, 212, 1)")
        
    }
    
    
    static get primaryContentColor() {
        
        return new CBColor("rgb(35, 35, 35)")
        
    }
    
    static get greenTintColor() {
        
        return new CBColor("rgba(51, 188, 125, 1)");
        
    }
    
    static get redTintColor() {
        
        return new CBColor("rgba(236, 88, 111, 1)");
        
    }
    
    
    
    
    
}






















/// <reference path="./UIObject.ts" />




class UIPoint extends UIObject {
    
    constructor(public x: number, public y: number) {
        
        super()
        
        this._class = UIPoint
        this.superclass = UIObject
        
    }
    
    
    
    
    
    copy() {
        return new UIPoint(this.x, this.y)
    }
    
    
    isEqualTo(point: UIPoint) {
    
        const result = (this.x == point.x && this.y == point.y)
    
        return result
        
    }
    
    scale(zoom: number) {
        const x = this.x
        const y = this.y
        this.x = x * zoom
        this.y = y * zoom
        return this
    }
    
    
    
    add(v: UIPoint) {
        this.x = this.x + v.x
        this.y = this.y + v.y
        return this
    }
    
    subtract(v: UIPoint) {
        this.x = this.x - v.x
        this.y = this.y - v.y
        return this
    }
    
    to(b: UIPoint) {
        const a = this
        const ab = b.copy().add(a.copy().scale(-1))
        return ab
    }
    
    pointWithX(x: number) {
        const result = this.copy()
        result.x = x
        return result
    }
    
    pointWithY(y: number) {
        const result = this.copy()
        result.y = y
        return result
    }
    
    pointByAddingX(x: number) {
        return this.pointWithX(this.x + x)
    }
    
    pointByAddingY(y: number) {
        return this.pointWithY(this.y + y)
    }
    
    
    get length() {
        var result = this.x * this.x + this.y * this.y
        result = Math.sqrt(result)
        return result
    }
    
    
    didChange(b) {
        
        // Callback to be set by delegate
        
    }
    
    
    
    
    
}























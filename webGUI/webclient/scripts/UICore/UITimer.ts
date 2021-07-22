///<reference path="UIObject.ts"/>









class UITimer extends UIObject {
    
    
    
    
    
    _intervalID: number
    
    isValid: boolean = YES
    
    
    constructor(public interval: number, public repeats: boolean, public target: Function) {
        
        super()
        
        this.superclass = UIObject
        
        
        this.schedule()
        
        
        
    }
    
    
    
    schedule() {
    
        const callback = function () {
            if (this.repeats == NO) {
                this.invalidate()
            }
            this.target()
        }.bind(this)
    
        this._intervalID = window.setInterval(callback, this.interval * 1000)
        
    }
    
    
    reschedule() {
        
        this.invalidate()
        this.schedule()
        
    }
    
    
    
    fire() {
        if (this.repeats == NO) {
            this.invalidate()
        }
        else {
            this.reschedule()
        }
        this.target()
    }
    
    invalidate() {
        
        if (this.isValid) {
            
            clearInterval(this._intervalID)
            
            this.isValid = NO
            
        }
        
        
        
    }
    
    
    
    
    
}



























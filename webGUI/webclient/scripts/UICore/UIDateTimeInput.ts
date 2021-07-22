///<reference path="UIView.ts"/>





class UIDateTimeInput extends UIView {
    
    
    
    
    
    constructor(elementID: string, type: string = UIDateTimeInput.type.DateTime) {
        
        super(elementID, nil, "input")
        
        this._class = UIDateTimeInput
        this.superclass = UIView
        
        
        this.viewHTMLElement.setAttribute("type", type)
        
        this.viewHTMLElement.onchange = (event) => {
            this.sendControlEventForKey(UIDateTimeInput.controlEvent.ValueChange, event)
        }
        
        
        
        
        
    }
    
    
    
    
    
    static controlEvent = Object.assign({}, UIView.controlEvent, {
        
        "ValueChange": "ValueChange"
        
    })
    
    get addControlEventTarget(): UIViewAddControlEventTargetObject<typeof UIDateTimeInput.controlEvent> {
        
        return super.addControlEventTarget as any
        
    }
    
    
    
    
    
    static type = {
        
        "Date": "date",
        "Time": "time",
        "DateTime": "datetime"
        
        
    }
    
    
    static format = {
        
        "European": "DD-MM-YYYY",
        "ISOComputer": "YYYY-MM-DD",
        "American": "MM/DD/YYYY"
        
    }
    
    
    
    
    
    get date() {
        
        const result = new Date((this.viewHTMLElement as HTMLInputElement).value)
        
        return result
        
    }
    
    
    
    
    
}

























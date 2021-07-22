/// <reference path="UITextField.ts" />




class UITextArea extends UITextField {
    
    
    
    
    
    constructor(elementID, viewHTMLElement = null) {
        
        super(elementID, viewHTMLElement, UITextView.type.textArea)
        
        this._class = UITextArea
        this.superclass = UITextField
        
        this.viewHTMLElement.removeAttribute("type")
        
        this.style.overflow = "auto"
        
        this.style.webkitUserSelect = "text"
        
        this.pausesPointerEvents = NO
        
    }
    
    
    get addControlEventTarget(): UIViewAddControlEventTargetObject<typeof UITextArea.controlEvent> {
        
        return super.addControlEventTarget as any;
        
    }
    
    // @ts-ignore
    get viewHTMLElement(): HTMLTextAreaElement {
        
        // @ts-ignore
        return super.viewHTMLElement
        
    }
    
    
    
    
    
}








































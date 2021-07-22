///<reference path="../UICore/UIView.ts"/>





class SearchTextField extends UIView {
    
    _textField: UIView
    _searchButton: UIView
    
    
    
    
    
    constructor(elementID: string) {
        
        super(elementID)
        
        
        this._class = SearchTextField
        this.superclass = UIView
        
        
        
        
    }
    
    
    
    initView(elementID: string, viewHTMLElement: HTMLElement) {
        
        super.initView(elementID, viewHTMLElement)
        
        
        viewHTMLElement.classList.add("input", "input--search")
        
        
        viewHTMLElement.innerHTML = "\
            <button type=\"button\" class=\"input__button\">\
                <i class=\"material-icons\">search</i>\
            </button>\
            <input type=\"search\" class=\"input__field\" placeholder=\"Search\">"
        
        
        
        
        this._textField = new UIView(nil, this.textFieldElement)
        this._searchButton = new UIView(nil, this.searchButtonElement)
        
        
        
        this._searchButton.addTargetForControlEvents([
            UIView.controlEvent.EnterDown, UIView.controlEvent.PointerUpInside
        ], function (this: SearchTextField, sender: UIView, event: Event) {
            
            this.performSearch()
            
        }.bind(this))
        
        this._textField.addTargetForControlEvent(
            UIView.controlEvent.EnterDown,
            function (this: SearchTextField, sender: UIView, event: Event) {
                
                this.performSearch()
                
            }.bind(this)
        )
        
        
        
        this._textField.viewHTMLElement.oninput = function (this: SearchTextField, event) {
            
            
            this._textField.sendControlEventForKey(UITextField.controlEvent.TextChange, event)
            
            
        }.bind(this)
        
        
        
        
        
    }
    
    
    
    
    
    get searchButtonElement() {
        
        return this.viewHTMLElement.querySelector("button") as HTMLButtonElement
        
        
    }
    
    get textFieldElement() {
        
        return this.viewHTMLElement.querySelector("input") as HTMLInputElement
        
        
    }
    
    
    
    set placeholderText(placeholderText: string) {
        
        this.textFieldElement.setAttribute("placeholder", placeholderText)
        
    }
    
    get placeholderText() {
        
        
        return this.textFieldElement.getAttribute("placeholder")
        
        
    }
    
    
    
    
    
    get text() {
        
        return this.textFieldElement.value
        
    }
    
    set text(text: string) {
        
        this.textFieldElement.value = text
        
    }
    
    
    
    focus() {
        
        this.textFieldElement.focus()
        
    }
    
    blur() {
        
        this.textFieldElement.blur()
        
    }
    
    
    
    
    
    performSearch(): any {
    
    
    
    }
    
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
        
        
        
        
        
    }
    
    
    
    
    
}
























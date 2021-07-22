///<reference path="../UICore/UIButton.ts"/>


class CBCheckbox extends UIView {
    
    
    checkbox: UIView
    titleLabel: UITextView
    
    
    
    
    
    constructor(elementID: string) {
        
        super(elementID, nil, "label")
        
        
        this._class = CBCheckbox
        this.superclass = UIView
        
        
        
        
    }
    
    
    
    
    
    initView(elementID: string, viewHTMLElement: HTMLElement) {
        
        super.initView(elementID, viewHTMLElement)
        
        
        this.checkbox = new UIView(elementID + "Checkbox", nil, "input")
        this.checkbox.viewHTMLElement.setAttribute("type", "checkbox")
        
        this.addSubview(this.checkbox)
        
        this.titleLabel = new UITextView(elementID + "TitleLabel", UITextView.type.span)
        
        this.addSubview(this.titleLabel)
        
        
        this.checkbox.userInteractionEnabled = NO
        this.titleLabel.userInteractionEnabled = NO
        
        this.titleLabel.style.overflow = "visible"
        this.titleLabel.style.lineHeight = "1.5"
        
        
        this.style.cursor = "pointer"
        
        this.viewHTMLElement.onchange = (event) => {
            this.sendControlEventForKey(CBCheckbox.controlEvent.SelectionChange, event)
        }
        
        
    }
    
    
    
    
    
    static controlEvent = Object.assign({}, UIView.controlEvent, {
        
        
        "SelectionChange": "SelectionChange"
        
        
    })
    
    
    
    
    
    get isStatic() {
        
        const result = this.checkbox.viewHTMLElement.classList.contains("staticCheckbox")
        
        return result
        
    }
    
    set isStatic(isStatic: boolean) {
        
        if (isStatic) {
            
            this.checkbox.viewHTMLElement.classList.add("staticCheckbox")
            
        }
        else {
            
            this.checkbox.viewHTMLElement.classList.remove("staticCheckbox")
            
        }
        
        
        
    }
    
    
    
    get selected() {
        
        return (this.checkbox.viewHTMLElement as HTMLInputElement).checked
        
    }
    
    set selected(selected: boolean) {
        
        (this.checkbox.viewHTMLElement as HTMLInputElement).checked = selected
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
        
        const titleWidth = this.bounds.width - 35
        this.titleLabel.setMaxSizes(nil, titleWidth)
        
        this.titleLabel.textPrefix = "<span style=\"position: absolute; overflow: hidden; left: 0; top: 0;text-overflow: ellipsis; white-space: pre; padding-left: 35px; width: " +
            titleWidth + "px;\">"
        
        this.titleLabel.textSuffix = "</span>"
        
        this.hoverText = this.titleLabel.text.replace(
            this.titleLabel.textPrefix,
            ""
        ).replace(this.titleLabel.textSuffix, "")
        
        
        
        
        
    }
    
    
    
    
    
}

































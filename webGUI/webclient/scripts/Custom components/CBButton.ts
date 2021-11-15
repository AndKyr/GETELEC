/// <reference path="../UICore/UIButton.ts" />





class CBButton extends UIButton {
    
    
    
    
    
    constructor(elementID?: string, elementType?: string) {
        
        super(elementID, elementType)
        
        this._class = CBButton
        this.superclass = UIButton
        
        
        
        
    }
    
    
    
    
    
    initView(elementID: string, viewHTMLElement: HTMLElement, initViewData) {
        
        super.initView(elementID, viewHTMLElement, initViewData)
        
        
        //this.style.borderRadius = "2px"
        
        this.style.outline = "none"
        
        
        this.colors.titleLabel.normal = UIColor.whiteColor
        this.setBackgroundColorsWithNormalColor(UIColor.colorWithRGBA(0, 79, 32))
        
        this.colors.titleLabel.selected = UIColor.blackColor
        
        
    }
    
    
    
    
    
    setBackgroundColorsWithNormalColor(normalBackgroundColor: UIColor) {
        
        this.colors.background.normal = normalBackgroundColor
        
        this.colors.background.hovered = UIColor.colorWithRGBA(87, 173, 122) // normalBackgroundColor.colorByMultiplyingRGB(0.85)
        
        this.colors.background.focused = normalBackgroundColor // normalBackgroundColor.colorByMultiplyingRGB(0.5)
        
        this.colors.background.highlighted = UIColor.colorWithRGBA(48, 196, 212) // normalBackgroundColor.colorByMultiplyingRGB(0.7)
        
        this.colors.background.selected = UIColor.whiteColor
        
        this.updateContentForCurrentState()
        
    }
    
    
    
    updateContentForNormalState() {
        
        super.updateContentForNormalState()
        
        this.setBorder(0, 0)
        
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.24)"
        
    }
    
    
    updateContentForHoveredState() {
        
        super.updateContentForHoveredState()
        
        this.setBorder(0, 0)
        
        //this.titleLabel.textColor = UIColor.whiteColor.colorByMultiplyingRGB(0.85);
        
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.18)"
        
    }
    
    updateContentForFocusedState() {
        
        super.updateContentForFocusedState()
        
        this.setBorder(0, 1, UIColor.blueColor)
        
        //this.titleLabel.textColor = UIColor.whiteColor.colorByMultiplyingRGB(0.85);
        
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.18)"
        
    }
    
    
    
    updateContentForHighlightedState() {
        
        super.updateContentForHighlightedState()
        
        this.setBorder(0, 0)
        
        //this.titleLabel.textColor = UIColor.whiteColor.colorByMultiplyingRGB(0.7);
        
        //this.style.boxShadow = "0 2px 2px 0 rgba(0,0,0,0.12)"
        
    }
    
    
    
    
    
    updateContentForCurrentEnabledState() {
        
        super.updateContentForCurrentEnabledState()
        
        if (IS_NOT(this.enabled)) {
            
            this.titleLabel.textColor = new UIColor("#adadad")
            
            this.backgroundColor = new UIColor("#e5e5e5")
            
            this.alpha = 1
            
        }
        
    }
    
    
    
    
    
}




























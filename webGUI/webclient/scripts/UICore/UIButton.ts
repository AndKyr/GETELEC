/// <reference path="./UIBaseButton.ts" />



interface UIButtonColorSpecifier {
    
    titleLabel: UIButtonElementColorSpecifier;
    background: UIButtonElementColorSpecifier;
    
}


interface UIButtonElementColorSpecifier {
    
    normal: UIColor;
    hovered?: UIColor;
    highlighted: UIColor;
    focused?: UIColor;
    selected: UIColor;
    selectedAndHighlighted?: UIColor;
    
}



class UIButton extends UIBaseButton {
    
    _contentPadding: any
    _titleLabel: UITextView
    _imageView: UIImageView
    
    usesAutomaticTitleFontSize = NO
    minAutomaticFontSize: number = nil
    maxAutomaticFontSize: number = 25
    
    colors: UIButtonColorSpecifier
    
    
    constructor(elementID: string, elementType?: string, titleType = UITextView.type.span) {
        
        super(elementID, elementType, { "titleType": titleType })
        
        this._class = UIButton
        this.superclass = UIBaseButton
        
        
    }
    
    
    
    
    
    initView(elementID, viewHTMLElement, initViewData: { titleType: string }) {
        
        this.class.superclass = UIBaseButton
        
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
            
        }
        
        
        this._imageView = new UIImageView(elementID + "ImageView")
        this._imageView.hidden = YES
        this.addSubview(this.imageView)
        
        this.imageView.fillMode = UIImageView.fillMode.aspectFitIfLarger
        
        
        if (IS_NOT_NIL(initViewData.titleType)) {
            
            this._titleLabel = new UITextView(elementID + "TitleLabel", initViewData.titleType)
            this.titleLabel.style.whiteSpace = "nowrap"
            this.addSubview(this.titleLabel)
            
            this.titleLabel.userInteractionEnabled = NO
            
        }
        
        this.contentPadding = 10
        
        this.imageView.userInteractionEnabled = NO
        this.titleLabel.textAlignment = UITextView.textAlignment.center
        this.titleLabel.nativeSelectionEnabled = NO
        
    }
    
    
    
    get contentPadding() {
        
        return this._contentPadding.integerValue
        
    }
    
    set contentPadding(contentPadding) {
        
        this._contentPadding = contentPadding
        
        this.setNeedsLayout()
        
    }
    
    
    
    
    public set hovered(hovered: boolean) {
        this._hovered = hovered
        this.updateContentForCurrentState()
    }
    
    public get hovered(): boolean {
        return this._hovered
    }
    
    public set highlighted(highlighted: boolean) {
        this._highlighted = highlighted
        this.updateContentForCurrentState()
    }
    
    public get highlighted(): boolean {
        return this._highlighted
    }
    
    public set focused(focused: boolean) {
        this._focused = focused
        if (focused) {
            this.focus()
        }
        else {
            this.blur()
        }
        this.updateContentForCurrentState()
    }
    
    public get focused(): boolean {
        return this._focused
    }
    
    public set selected(selected: boolean) {
        this._selected = selected
        this.updateContentForCurrentState()
    }
    
    public get selected(): boolean {
        return this._selected
    }
    
    
    
    
    
    updateContentForCurrentState() {
    
        var updateFunction: Function = this.updateContentForNormalState
        if (this.selected && this.highlighted) {
            updateFunction = this.updateContentForSelectedAndHighlightedState
        }
        else if (this.selected) {
            updateFunction = this.updateContentForSelectedState
        }
        else if (this.focused) {
            updateFunction = this.updateContentForFocusedState
        }
        else if (this.highlighted) {
            updateFunction = this.updateContentForHighlightedState
        }
        else if (this.hovered) {
            updateFunction = this.updateContentForHoveredState
        }
        
        if (!IS(updateFunction)) {
            this.titleLabel.textColor = UIColor.nilColor
            this.backgroundColor = UIColor.nilColor
        }
        else {
            updateFunction.call(this)
        }
        
        this.updateContentForCurrentEnabledState();
        
    }
    
    updateContentForNormalState() {
        
        this.backgroundColor = this.colors.background.normal
        this.titleLabel.textColor = this.colors.titleLabel.normal
        
    }
    
    updateContentForHoveredState() {
        
        this.updateContentForNormalState()
        
        if (this.colors.background.hovered) {
            this.backgroundColor = this.colors.background.hovered
        }
        
        if (this.colors.titleLabel.hovered) {
            this.titleLabel.textColor = this.colors.titleLabel.hovered
        }
        
    }
    
    updateContentForFocusedState() {
        
        this.updateContentForHoveredState()
        
        if (this.colors.background.focused) {
            this.backgroundColor = this.colors.background.focused
        }
        
        if (this.colors.titleLabel.focused) {
            this.titleLabel.textColor = this.colors.titleLabel.focused
        }
        
    }
    
    updateContentForHighlightedState() {
        
        this.backgroundColor = this.colors.background.highlighted
        this.titleLabel.textColor = this.colors.titleLabel.highlighted
        
    }
    
    updateContentForSelectedState() {
        
        this.backgroundColor = this.colors.background.selected
        this.titleLabel.textColor = this.colors.titleLabel.selected
        
    }
    
    updateContentForSelectedAndHighlightedState() {
        
        this.updateContentForSelectedState()
        
        if (this.colors.background.selectedAndHighlighted) {
            this.backgroundColor = this.colors.background.selectedAndHighlighted
        }
        
        if (this.colors.titleLabel.selectedAndHighlighted) {
            this.titleLabel.textColor = this.colors.titleLabel.selectedAndHighlighted
        }
        
    }
    
    
    set enabled(enabled: boolean) {
        
        super.enabled = enabled
        
        this.updateContentForCurrentState()
        
    }
    
    get enabled() {
        
        return super.enabled
        
    }
    
    updateContentForCurrentEnabledState() {
        
        if (this.enabled) {
            this.alpha = 1
        }
        else {
            this.alpha = 0.5
        }
        
        this.userInteractionEnabled = this.enabled
        
    }
    
    
    
    addStyleClass(styleClassName: string) {
        
        super.addStyleClass(styleClassName)
        
        if (this.styleClassName != styleClassName) {
            
            this.updateContentForCurrentState.call(this)
            
        }
        
    }
    
    
    get titleLabel() {
        
        return this._titleLabel
        
    }
    
    get imageView() {
        
        return this._imageView
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
    
        var bounds = this.bounds
    
        this.hoverText = this.titleLabel.text
        
        // Image only if text is not present
        if (IS_NOT(this.imageView.hidden) && !IS(this.titleLabel.text)) {
            
            
            this.imageView.frame = bounds
            
            
        }
        
        // Text only if image is not present
        if (IS(this.imageView.hidden) && IS(this.titleLabel.text)) {
            
            var titleElement = this.titleLabel.viewHTMLElement
            
            
            this.titleLabel.style.left = this.contentPadding
            this.titleLabel.style.right = this.contentPadding
            // this.titleLabel.style.marginLeft = ""
            // this.titleLabel.style.right = this.contentPadding
            this.titleLabel.style.top = "50%"
            this.titleLabel.style.transform = "translateY(-50%)"
            this.titleLabel.frame = new UIRectangle(nil, nil, nil, nil)
            
            if (this.usesAutomaticTitleFontSize) {
                
                var hidden = this.titleLabel.hidden
                
                this.titleLabel.hidden = YES
                
                this.titleLabel.fontSize = 15
                
                this.titleLabel.fontSize = UITextView.automaticallyCalculatedFontSize(
                    new UIRectangle(
                        0,
                        0,
                        this.bounds.height,
                        1 *
                        this.titleLabel.viewHTMLElement.offsetWidth
                    ),
                    this.titleLabel.intrinsicContentSize(),
                    this.titleLabel.fontSize,
                    this.minAutomaticFontSize,
                    this.maxAutomaticFontSize
                )
                
                this.titleLabel.hidden = hidden
                
            }
            
            
        }
        
        // Image and text both present
        if (IS_NOT(this.imageView.hidden) && IS(this.titleLabel.text)) {
    
            const imageShareOfWidth = 0.25
    
            bounds = bounds.rectangleWithInset(this.contentPadding)
    
            const imageFrame = bounds.copy()
            imageFrame.width = bounds.height - this.contentPadding * 0.5
            this.imageView.frame = imageFrame
            
            var titleElement = this.titleLabel.viewHTMLElement
            
            this.titleLabel.style.left = imageFrame.max.x + this.contentPadding
            this.titleLabel.style.right = this.contentPadding
            this.titleLabel.style.top = "50%"
            this.titleLabel.style.transform = "translateY(-50%)"
            
            if (this.usesAutomaticTitleFontSize) {
                
                var hidden = this.titleLabel.hidden
                
                this.titleLabel.hidden = YES
                
                this.titleLabel.fontSize = 15
                
                this.titleLabel.fontSize = UITextView.automaticallyCalculatedFontSize(
                    new UIRectangle(
                        0,
                        0,
                        this.bounds.height,
                        1 *
                        this.titleLabel.viewHTMLElement.offsetWidth
                    ),
                    this.titleLabel.intrinsicContentSize(),
                    this.titleLabel.fontSize,
                    this.minAutomaticFontSize,
                    this.maxAutomaticFontSize
                )
                
                this.titleLabel.hidden = hidden
                
            }
            
        }
        
        this.applyClassesAndStyles()
        
    }
    
    initViewStyleSelectors() {
        
        this.initStyleSelector("." + this.styleClassName, "background-color: lightblue;")
        
        // var selectorWithoutImage = "." + this.styleClassName + " ." + this.imageView.styleClassName + " + ." + this.titleLabel.styleClassName;
        
        // this.initStyleSelector(
        //     selectorWithoutImage,
        //     "left: " + this.contentPadding + ";" +
        //     "right: " + this.contentPadding + ";" +
        //     "top: 50%;" +
        //     "transform: translateY(-50%);");
        
    }
    
    
    
    
    
}






























































































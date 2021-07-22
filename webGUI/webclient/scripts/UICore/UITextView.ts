/// <reference path="UIView.ts" />
/// <reference path="UICore.ts" />





class UITextView extends UIView {
    
    
    _textColor: UIColor = UITextView.defaultTextColor
    _textAlignment: string
    
    _isSingleLine = YES
    
    textPrefix = ""
    textSuffix = ""
    
    _notificationAmount = 0
    
    _minFontSize: number = nil
    _maxFontSize: number = nil
    
    _automaticFontSizeSelection = NO
    
    changesOften = NO
    
    
    static defaultTextColor = UIColor.blackColor
    static notificationTextColor = UIColor.redColor
    
    static _intrinsicHeightCache: { [x: string]: { [x: string]: number; }; } & UIObject = new UIObject() as any
    static _intrinsicWidthCache: { [x: string]: { [x: string]: number; }; } & UIObject = new UIObject() as any
    
    
    static _ptToPx: number
    static _pxToPt: number
    _text: string
    
    
    constructor(elementID?: string, textViewType = UITextView.type.paragraph, viewHTMLElement = null) {
        
        super(elementID, viewHTMLElement, textViewType)
        
        this._class = UITextView
        this.superclass = UIView
        
        this.text = ""
        
        this.style.overflow = "hidden"
        this.style.textOverflow = "ellipsis"
        this.isSingleLine = YES
        
        this.textColor = this.textColor
        
        this.userInteractionEnabled = YES
        
        if (textViewType == UITextView.type.textArea) {
            
            this.pausesPointerEvents = YES
            
            this.addTargetForControlEvent(UIView.controlEvent.PointerUpInside, function (sender, event) {
                
                sender.focus()
                
            })
            
            
        }
        
        
    }
    
    
    
    
    
    static _determinePXAndPTRatios() {
    
        const o = document.createElement("div")
        o.style.width = "1000pt"
        document.body.appendChild(o)
        UITextView._ptToPx = o.clientWidth / 1000
        document.body.removeChild(o)
        UITextView._pxToPt = 1 / UITextView._ptToPx
        
    }
    
    
    
    
    
    static type = {
        
        "paragraph": "p",
        "header1": "h1",
        "header2": "h2",
        "header3": "h3",
        "header4": "h4",
        "header5": "h5",
        "header6": "h6",
        "textArea": "textarea",
        "textField": "input",
        "span": "span",
        "label": "label"
        
    }
    
    
    
    static textAlignment = {
        
        "left": "left",
        "center": "center",
        "right": "right",
        "justify": "justify"
        
    }
    
    get textAlignment() {
        const result = this.style.textAlign
        return result
    }
    
    set textAlignment(textAlignment: string) {
        this._textAlignment = textAlignment
        this.style.textAlign = textAlignment
    }
    
    
    
    get textColor() {
        const result = this._textColor
        return result
    }
    
    set textColor(color: UIColor) {
        
        this._textColor = color || UITextView.defaultTextColor
        this.style.color = this._textColor.stringValue
        
    }
    
    
    get isSingleLine() {
        
        return this._isSingleLine
        
    }
    
    set isSingleLine(isSingleLine: boolean) {
        
        this._isSingleLine = isSingleLine
        
        
        if (isSingleLine) {
            
            this.style.whiteSpace = "pre"
            
            return
            
        }
        
        this.style.whiteSpace = "pre-wrap"
        
    }
    
    
    get notificationAmount() {
        
        return this._notificationAmount
        
    }
    
    set notificationAmount(notificationAmount: number) {
        
        if (this._notificationAmount == notificationAmount) {
            
            return
            
        }
        
        this._notificationAmount = notificationAmount
        
        this.text = this.text
        
        this.setNeedsLayoutUpToRootView()
        
        this.notificationAmountDidChange(notificationAmount)
        
    }
    
    notificationAmountDidChange(notificationAmount: number) {
        
        
        
        
        
    }
    
    
    
    
    
    get text() {
        
        return (this._text || this.viewHTMLElement.innerHTML)
        
    }
    
    set text(text) {
        
        this._text = text
    
        var notificationText = ""
    
        if (this.notificationAmount) {
            
            notificationText = "<span style=\"color: " + UITextView.notificationTextColor.stringValue + ";\">" +
                (" (" + this.notificationAmount + ")").bold() + "</span>"
            
        }
        
        if (this.viewHTMLElement.innerHTML != this.textPrefix + text + this.textSuffix + notificationText) {
            
            this.viewHTMLElement.innerHTML = this.textPrefix + FIRST(text, "") + this.textSuffix + notificationText
            
        }
        
        this.setNeedsLayout();
        
    }
    
    set innerHTML(innerHTML: string) {
        
        this.text = innerHTML
        
    }
    
    get innerHTML() {
        
        return this.viewHTMLElement.innerHTML
        
    }
    
    
    
    setText(key: string, defaultString: string, parameters?: { [x: string]: string | UILocalizedTextObject }) {
        
        this.setInnerHTML(key, defaultString, parameters)
        
    }
    
    
    get fontSize() {
    
        const style = window.getComputedStyle(this.viewHTMLElement, null).fontSize
    
        const result = (parseFloat(style) * UITextView._pxToPt)
    
        return result
        
    }
    
    set fontSize(fontSize: number) {
        
        
        this.style.fontSize = "" + fontSize + "pt"
        
        
    }
    
    
    
    
    
    useAutomaticFontSize(minFontSize: number = nil, maxFontSize: number = nil) {
        
        
        this._automaticFontSizeSelection = YES
        
        
        this._minFontSize = minFontSize
        
        this._maxFontSize = maxFontSize
        
        this.setNeedsLayout()
        
        
    }
    
    
    
    
    
    static automaticallyCalculatedFontSize(
        bounds: UIRectangle,
        currentSize: UIRectangle,
        currentFontSize: number,
        minFontSize?: number,
        maxFontSize?: number
    ) {
        
        minFontSize = FIRST(minFontSize, 1)
        
        maxFontSize = FIRST(maxFontSize, 100000000000)
    
    
        const heightMultiplier = bounds.height / (currentSize.height + 1)
    
        const widthMultiplier = bounds.width / (currentSize.width + 1)
    
    
    
    
        var multiplier = heightMultiplier
    
        if (heightMultiplier > widthMultiplier) {
            
            multiplier = widthMultiplier
            
            
            
        }
    
    
    
    
        const maxFittingFontSize = currentFontSize * multiplier
    
    
        if (maxFittingFontSize > maxFontSize) {
            
            return maxFontSize
            
        }
        
        if (minFontSize > maxFittingFontSize) {
            
            return minFontSize
            
        }
        
        
        return maxFittingFontSize
        
        
        
    }
    
    
    
    
    
    didReceiveBroadcastEvent(event: UIViewBroadcastEvent) {
        
        super.didReceiveBroadcastEvent(event)
        
    }
    
    
    willMoveToSuperview(superview: UIView) {
        
        super.willMoveToSuperview(superview)
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
        
        
        if (this._automaticFontSizeSelection) {
            
            this.fontSize = UITextView.automaticallyCalculatedFontSize(
                new UIRectangle(0, 0, 1 *
                    this.viewHTMLElement.offsetHeight, 1 *
                    this.viewHTMLElement.offsetWidth),
                this.intrinsicContentSize(),
                this.fontSize,
                this._minFontSize,
                this._maxFontSize
            )
            
            
            
            
        }
        
        
        
        
        
    }
    
    
    
    
    
    intrinsicContentHeight(constrainingWidth = 0) {
        
        if (this.changesOften) {
            
            return super.intrinsicContentHeight(constrainingWidth)
            
        }
    
        const keyPath = (this.viewHTMLElement.innerHTML + "_csf_" + this.computedStyle.font).replace(new RegExp("\\.", "g"), "_") + "." +
            ("" + constrainingWidth).replace(new RegExp("\\.", "g"), "_")
    
        var result = UITextView._intrinsicHeightCache.valueForKeyPath(keyPath)
    
        if (IS_LIKE_NULL(result)) {
            
            result = super.intrinsicContentHeight(constrainingWidth)
            
            UITextView._intrinsicHeightCache.setValueForKeyPath(keyPath, result)
            
        }
        
        
        return result
        
    }
    
    intrinsicContentWidth(constrainingHeight = 0) {
        
        if (this.changesOften) {
            
            return super.intrinsicContentWidth(constrainingHeight)
            
        }
    
        const keyPath = (this.viewHTMLElement.innerHTML + "_csf_" + this.computedStyle.font).replace(new RegExp("\\.", "g"), "_") + "." +
            ("" + constrainingHeight).replace(new RegExp("\\.", "g"), "_")
    
        var result = UITextView._intrinsicWidthCache.valueForKeyPath(keyPath)
    
        if (IS_LIKE_NULL(result)) {
            
            result = super.intrinsicContentWidth(constrainingHeight)
            
            UITextView._intrinsicWidthCache.setValueForKeyPath(keyPath, result)
            
        }
        
        
        return result
        
    }
    
    
    
    
    intrinsicContentSize() {
        
        // This works but is slow
        const result = this.intrinsicContentSizeWithConstraints(nil, nil)
    
        return result
        
    }
    
    
    
    
    
}





UITextView._determinePXAndPTRatios()





// /**
//  * Uses canvas.measureText to compute and return the width of the given text of given font in pixels.
//  * 
//  * @param {String} text The text to be rendered.
//  * @param {String} font The css font descriptor that text is to be rendered with (e.g. "bold 14px verdana").
//  * 
//  * @see https://stackoverflow.com/questions/118241/calculate-text-width-with-javascript/21015393#21015393
//  */
// function getTextMetrics(text, font) {
//     // re-use canvas object for better performance
//     var canvas = getTextMetrics.canvas || (getTextMetrics.canvas = document.createElement("canvas"));
//     var context = canvas.getContext("2d");
//     context.font = font;
//     var metrics = context.measureText(text);
//     return metrics;
// }

































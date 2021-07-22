///<reference path="UIView.ts"/>









class UIImageView extends UIView {
    
    
    //actionIndicator: UIActionIndicator;
    _sourceKey: string
    _defaultSource: string
    
    _fillMode: any
    
    _hiddenWhenEmpty = NO
    
    
    constructor(elementID?, viewHTMLElement = null) {
        
        super(elementID, viewHTMLElement, "img")
        
        this._class = UIImageView
        this.superclass = UIView
        
        //this.actionIndicator = new UIActionIndicator(elementID + "ActionIndicator");
        
        
    }
    
    
    
    
    static fillMode = {
        
        "stretchToFill": "fill",
        "aspectFit": "contain",
        "aspectFill": "cover",
        "center": "none",
        "aspectFitIfLarger": "scale-down"
        
    }
    
    
    
    
    
    get viewHTMLElement() {
        
        return super.viewHTMLElement as HTMLImageElement
        
    }
    
    
    
    
    
    
    
    
    
    static dataURL(url, callback) {
        const xhr = new XMLHttpRequest()
        xhr.open("get", url)
        xhr.responseType = "blob"
        xhr.onload = function () {
            const fr = new FileReader()
    
            fr.onload = function () {
                callback(this.result)
            }
            
            fr.readAsDataURL(xhr.response) // async call
        }
        
        xhr.send()
    }
    
    
    
    static dataURLWithMaxSize(URLString: string, maxSize: number, completion: (resultURLString: string) => void) {
    
        const imageView = new UIImageView()
        imageView.imageSource = URLString
        
        imageView.viewHTMLElement.onload = function () {
    
            const originalSize = imageView.intrinsicContentSize()
    
    
            var multiplier = maxSize / Math.max(originalSize.height, originalSize.width)
    
            multiplier = Math.min(1, multiplier)
    
    
            const result = imageView.getDataURL((originalSize.height * multiplier).integerValue, (originalSize.width *
                multiplier).integerValue)
    
            completion(result)
            
        }
        
    }
    
    
    static dataURLWithSizes(
        URLString: string,
        height: number,
        width: number,
        completion: (resultURLString: string) => void
    ) {
    
        const imageView = new UIImageView()
        imageView.imageSource = URLString
        
        imageView.viewHTMLElement.onload = function () {
    
            const result = imageView.getDataURL(height, width)
            completion(result)
            
        }
        
    }
    
    
    getDataURL(height?: number, width?: number) {
    
        const img = this.viewHTMLElement
    
        // Create an empty canvas element
        const canvas = document.createElement("canvas")
        canvas.width = width
        canvas.height = height
        
        // Copy the image contents to the canvas
        const ctx = canvas.getContext("2d")
        ctx.drawImage(img, 0, 0, width, height)
        
        // Get the data-URL formatted image
        // Firefox supports PNG and JPEG. You could check img.src to
        // guess the original format, but be aware the using "image/jpg"
        // will re-encode the image.
        const dataURL = canvas.toDataURL("image/png")
    
        return dataURL
        
        //return dataURL.replace(/^data:image\/(png|jpg);base64,/, "");
        
    }
    
    
    
    
    get imageSource() {
        
        return this.viewHTMLElement.src
        
    }
    
    set imageSource(sourceString: string) {
        
        if (IS_NOT(sourceString)) {
            sourceString = ""
        }
        
        this.viewHTMLElement.src = sourceString
        
        if (this.hiddenWhenEmpty) {
            this.hidden = IS_NOT(this.imageSource)
        }
        
        if (!sourceString || !sourceString.length) {
            
            //this.actionIndicator.stop();
            this.hidden = YES
            
            return
            
        }
        else {
            
            this.hidden = NO
            
        }
        
        // this.superview.addSubview(this.actionIndicator);
        // this.actionIndicator.frame = this.frame;
        // this.actionIndicator.start();
        // this.actionIndicator.backgroundColor = UIColor.redColor
        
        // @ts-ignore
        this.viewHTMLElement.onload = function (this: UIImageView, event: Event) {
            
            this.superview.setNeedsLayout()
            
            //this.actionIndicator.removeFromSuperview();
            
        }.bind(this)
        
    }
    
    
    
    setImageSource(key: string, defaultString: string) {
    
        const languageName = UICore.languageService.currentLanguageKey
        this.imageSource = UICore.languageService.stringForKey(key, languageName, defaultString, nil)
        
    }
    
    
    
    didReceiveBroadcastEvent(event: UIViewBroadcastEvent) {
        
        super.didReceiveBroadcastEvent(event)
        
        if (event.name == UIView.broadcastEventName.LanguageChanged || event.name ==
            UIView.broadcastEventName.AddedToViewTree) {
            
            this._setImageSourceFromKeyIfPossible()
            
        }
        
    }
    
    
    willMoveToSuperview(superview: UIView) {
        
        super.willMoveToSuperview(superview)
        
        this._setImageSourceFromKeyIfPossible()
        
    }
    
    _setImageSourceFromKeyIfPossible() {
        
        if (this._sourceKey && this._defaultSource) {
            
            this.setImageSource(this._sourceKey, this._defaultSource)
            
        }
        
    }
    
    
    
    get fillMode() {
        
        
        return this._fillMode
        
    }
    
    set fillMode(fillMode) {
        
        this._fillMode = fillMode;
        
        (this.style as any).objectFit = fillMode
        
    }
    
    
    get hiddenWhenEmpty() {
        return this._hiddenWhenEmpty
    }
    
    set hiddenWhenEmpty(hiddenWhenEmpty: boolean) {
        this._hiddenWhenEmpty = hiddenWhenEmpty
        if (hiddenWhenEmpty) {
            this.hidden = IS_NOT(this.imageSource)
        }
    }
    
    
    
    
    
    didMoveToSuperview(superview: UIView) {
        
        super.didMoveToSuperview(superview)
        
        
        
        
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
        
        
        
    }
    
    
    
    
    intrinsicContentSize() {
    
    
        const result = new UIRectangle(0, 0, this.viewHTMLElement.naturalHeight, this.viewHTMLElement.naturalWidth)
    
        return result
        
        
    }
    
    intrinsicContentSizeWithConstraints(constrainingHeight = 0, constrainingWidth = 0) {
    
        const heightRatio = constrainingHeight / this.viewHTMLElement.naturalHeight
    
        const widthRatio = constrainingWidth / this.viewHTMLElement.naturalWidth
    
        const multiplier = Math.max(heightRatio, widthRatio)
    
        const result = new UIRectangle(0, 0, this.viewHTMLElement.naturalHeight *
            multiplier, this.viewHTMLElement.naturalWidth * multiplier)
    
        return result
        
        
    }
    
    
    
    
    
}















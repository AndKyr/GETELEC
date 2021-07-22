/// <reference path="../UICore/UIView.ts" />





class CellView extends UIButton {
    
    
    _isAButton = NO
    _rightImageView: UIImageView
    
    leftInset = 0
    rightInset = 0
    
    //titleLabel: UITextView;
    
    constructor(elementID?: string, titleLabelType: string = UITextView.type.span) {
        
        super(elementID, undefined, titleLabelType)
        
        this._class = CellView
        this.superclass = UIView
        
        // this.titleLabel = new UITextView(this.elementID + "TitleLabel", titleLabelType);
        // this.addSubview(this.titleLabel);
        
        this.updateForCurrentIsAButtonState()
        
        
        
        
        
    }
    
    
    
    
    
    set isAButton(isAButton: boolean) {
        
        this._isAButton = isAButton
        
        this.updateForCurrentIsAButtonState()
        
    }
    
    get isAButton() {
        
        return this._isAButton
        
    }
    
    
    updateForCurrentIsAButtonState() {
        
        
        if (this._isAButton) {
            
            this.style.cursor = "pointer"
            
            this.titleLabel.userInteractionEnabled = NO
            
            this.titleLabel.nativeSelectionEnabled = NO
            
            this.titleLabel.textAlignment = UITextView.textAlignment.center
            
            this.nativeSelectionEnabled = NO
            
            this.style.outline = ""
            
            this.colors = {
                
                titleLabel: {
                    
                    normal: UIColor.blueColor,
                    highlighted: UIColor.blueColor,
                    selected: UIColor.blueColor
                    
                },
                background: {
                    
                    normal: UIColor.transparentColor,
                    hovered: new UIColor("#F8F8F8"),
                    highlighted: new UIColor("#ebebeb"),
                    selected: new UIColor("#ebebeb")
                    
                }
                
            }
            
            
        }
        else {
            
            this.style.cursor = ""
            
            this.titleLabel.userInteractionEnabled = YES
            
            this.titleLabel.nativeSelectionEnabled = YES
            
            this.titleLabel.textAlignment = UITextView.textAlignment.left
            
            this.nativeSelectionEnabled = YES
            
            this.style.outline = "none"
            
            this.colors = {
                
                titleLabel: {
                    
                    normal: UIColor.blackColor,
                    highlighted: UIColor.blackColor,
                    selected: UIColor.blackColor
                    
                },
                background: {
                    
                    normal: UIColor.transparentColor,
                    highlighted: UIColor.transparentColor,
                    selected: UIColor.transparentColor
                    
                }
                
            }
            
            
        }
        
        this.updateContentForCurrentState()
        
        
        
    }
    
    
    
    
    
    initRightImageViewIfNeeded() {
        
        if (this._rightImageView) {
            
            return
            
        }
        
        
        this._rightImageView = new UIImageView(this.elementID + "RightImageView")
        
        this._rightImageView.userInteractionEnabled = NO
        
        
        
        
    }
    
    
    set rightImageSource(imageSource: string) {
        
        
        
        if (IS(imageSource)) {
            
            this.initRightImageViewIfNeeded()
            
            this._rightImageView.imageSource = imageSource
            
            this.addSubview(this._rightImageView)
            
        }
        else {
            
            this._rightImageView.removeFromSuperview()
            
        }
        
        
    }
    
    get rightImageSource() {
    
        var result = nil
    
        if (this._rightImageView) {
            
            result = this._rightImageView.imageSource
            
        }
        
        return result
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
    
        const padding = RootViewController.paddingLength
        const labelHeight = padding
    
        const bounds = this.bounds
    
    
        this.titleLabel.centerYInContainer()
        
        this.titleLabel.style.left = "" + (padding * 0.5 + this.leftInset).integerValue + "px"
        this.titleLabel.style.right = "" + (padding * 0.5 + this.rightInset).integerValue + "px"
        
        this.titleLabel.style.maxHeight = "100%"
        
        this.titleLabel.style.overflow = "hidden"
        
        //this.titleLabel.style.whiteSpace = "nowrap";
        
        
        if (this._rightImageView && this._rightImageView.superview == this) {
            
            // var imageHeight = bounds.height - padding;
            // this._rightImageView.frame = new UIRectangle(bounds.width - imageHeight - padding * 0.5, padding * 0.5, imageHeight, imageHeight);
            
            this._rightImageView.frame = bounds.rectangleWithInsets(this.leftInset, padding * 0.5 +
                this.rightInset, 0, 0).rectangleWithWidth(24, 1).rectangleWithHeight(24, 0.5)
            
            this.titleLabel.style.right = "" +
                (padding * 0.5 + this.rightInset + this._rightImageView.frame.width).integerValue + "px"
            
        }
        
        
        
        
        
    }
    
    
    
    
}




















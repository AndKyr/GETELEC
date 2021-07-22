///<reference path="UIView.ts"/>





class UIDialogView<ViewType extends UIView = UIView> extends UIView {
    
    
    
    _view: ViewType = nil
    
    _appearedAnimated: boolean
    
    animationDuration: number = 0.25
    
    _zIndex: number = 100
    
    isVisible: boolean = NO
    
    dismissesOnTapOutside = YES
    
    constructor(elementID?: string, viewHTMLElement?: HTMLElement) {
        
        
        super(elementID, viewHTMLElement)
        
        this.addTargetForControlEvent(
            UIView.controlEvent.PointerTap,
            function (this: UIDialogView, sender: UIView, event: Event) {
                
                this.didDetectTapOutside(sender, event)
                
            }.bind(this)
        )
        
        this.backgroundColor = UIColor.colorWithRGBA(0, 10, 25).colorWithAlpha(0.75) //CBColor.primaryContentColor.colorWithAlpha(0.75)
        
        this.zIndex = this._zIndex
        
    }
    
    
    
    didDetectTapOutside(sender: UIView, event: Event) {
        
        if (event.target == this.viewHTMLElement && this.dismissesOnTapOutside) {
            this.dismiss(this._appearedAnimated)
        }
        
    }
    
    
    set zIndex(zIndex: number) {
        
        this._zIndex = zIndex
        this.style.zIndex = "" + zIndex
        
    }
    
    get zIndex() {
        
        return this._zIndex
        
    }
    
    
    
    set view(view: ViewType) {
        
        this._view.removeFromSuperview()
        
        this._view = view
        
        this.addSubview(view)
        
    }
    
    
    get view() {
        
        return this._view
        
    }
    
    
    
    willAppear(animated: boolean = NO) {
        
        if (animated) {
            
            this.style.opacity = "0"
            
        }
        
        this.style.height = ""
        
        this._frame = null
        
    }
    
    
    animateAppearing() {
        
        this.style.opacity = "1"
        
    }
    
    animateDisappearing() {
        
        this.style.opacity = "0"
        
    }
    
    
    
    showInView(containerView: UIView, animated: boolean) {
        
        
        animated = (animated && !IS_FIREFOX)
        
        this._appearedAnimated = animated
        
        this.willAppear(animated)
        
        
        containerView.addSubview(this)
        
        if (animated) {
            
            this.layoutSubviews()
            
            UIView.animateViewOrViewsWithDurationDelayAndFunction(
                this,
                this.animationDuration,
                0,
                undefined,
                function () {
                    
                    
                    this.animateAppearing()
                    
                    
                }.bind(this),
                nil
            )
            
            
        }
        else {
            
            this.setNeedsLayout()
            
        }
        
        this.isVisible = YES
        
    }
    
    
    
    showInRootView(animated: boolean) {
        
        this.showInView(UICore.main.rootViewController.view, animated)
        
    }
    
    
    
    dismiss(animated?: boolean) {
        
        
        animated = (animated && !IS_FIREFOX)
        
        if (animated == undefined) {
            
            animated = this._appearedAnimated
            
        }
        
        if (animated) {
            
            UIView.animateViewOrViewsWithDurationDelayAndFunction(
                this,
                this.animationDuration,
                0,
                undefined,
                function () {
                    
                    this.animateDisappearing()
                    
                }.bind(this),
                function () {
                    
                    if (this.isVisible == NO) {
                        
                        this.removeFromSuperview()
                        
                    }
                    
                }.bind(this)
            )
            
        }
        else {
            
            this.removeFromSuperview()
            
        }
        
        this.isVisible = NO
        
    }
    
    
    
    
    
    didReceiveBroadcastEvent(event: UIViewBroadcastEvent) {
        
        super.didReceiveBroadcastEvent(event)
        
        if (event.name == UICore.broadcastEventName.WindowDidResize) {
            
            this.setNeedsLayout()
            
        }
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        
        if (!IS(this.view)) {
            
            return
            
        }
        
        //this.frame = this.superview.bounds;
        
        this.setPosition(0, 0, 0, 0, 0, "100%")
        this.setPosition(0, 0, 0, 0, UIView.pageHeight, "100%")
        
        const bounds = this.bounds
        
        const margin = 20
        
        //this.view.centerInContainer();
        
        this.view.style.position = "relative"
        
        // this.view.style.maxHeight = "" + (bounds.height - margin * 2).integerValue + "px";
        // this.view.style.maxWidth = "" + (bounds.width - margin * 2).integerValue + "px";
        
        
        
        
        // var viewIntrinsicRectangle = this.view.intrinsicContentSize();
        // this.view.frame = new UIRectangle((bounds.width - viewIntrinsicRectangle.width)*0.5,  )
        
        super.layoutSubviews()
        
    }
    
    
    
    
    
}
























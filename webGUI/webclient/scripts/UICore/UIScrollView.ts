/// <reference path="./UIView.ts" />





class UIScrollView extends UIView {
    
    
    _contentOffset: UIPoint = new UIPoint(0, 0)
    _contentScale: number = 1
    
    containerView: UIView
    
    _pointerDown: boolean
    
    _scrollEnabled: boolean = YES
    
    _previousClientPoint: UIPoint
    
    _intrinsicContentFrame: UIRectangle
    
    constructor(elementID: string, viewHTMLElement?: HTMLElement) {
        
        super(elementID, viewHTMLElement)
        
        this._class = UIScrollView
        this.superclass = UIView
        
        
        this.containerView = new UIView(elementID + "ContainerView")
        
        super.addSubview(this.containerView)
        
        
        this.style.overflow = "hidden"
        
        this.pausesPointerEvents = NO //YES;
        
        
        this.addTargetForControlEvent(UIView.controlEvent.PointerDown, function () {
            
            this._pointerDown = YES
            
        }.bind(this))
        
        this.addTargetForControlEvent(UIView.controlEvent.PointerUp, function () {
            
            this._pointerDown = NO
            
            this._previousClientPoint = null
            
            scrollStopped()
            
        }.bind(this))
        
        
        function scrollStopped() {
            
            // Handle paging if needed
            
        }
        
        
        this.addTargetForControlEvent(UIView.controlEvent.PointerMove, function (sender: UIScrollView, event: Event) {
            
            if (!(this._pointerDown && this._scrollEnabled && this._enabled)) {
                
                return
                
            }
    
            const currentClientPoint = new UIPoint(nil, nil)
    
            if ((window as any).MouseEvent && event instanceof MouseEvent) {
                
                currentClientPoint.x = (event as MouseEvent).clientX
                currentClientPoint.y = (event as MouseEvent).clientY
                
            }
            
            if ((window as any).TouchEvent && event instanceof TouchEvent) {
    
                const touchEvent: TouchEvent = event
    
                if (touchEvent.touches.length != 1) {
                    
                    this._pointerDown = NO
                    this._previousClientPoint = null
                    
                    scrollStopped()
                    
                    return
                    
                }
                
                currentClientPoint.x = touchEvent.touches[0].clientX
                currentClientPoint.y = touchEvent.touches[0].clientY
                
            }
            
            if (!this._previousClientPoint) {
                
                this._previousClientPoint = currentClientPoint
                
                return
                
            }
    
            const changePoint = currentClientPoint.copy().subtract(this._previousClientPoint)
    
    
            if (this.containerView.bounds.width <= this.bounds.width) {
                changePoint.x = 0
            }
            if (0 < this.contentOffset.x + changePoint.x) {
                changePoint.x = -this.contentOffset.x
            }
            if (this.contentOffset.x + changePoint.x < -this.bounds.width) {
                changePoint.x = -this.bounds.width - this.contentOffset.x
            }
            
            if (this.containerView.bounds.height <= this.bounds.height) {
                changePoint.y = 0
            }
            if (0 < this.contentOffset.y + changePoint.y) {
                changePoint.y = -this.contentOffset.y
            }
            if (this.contentOffset.y + changePoint.y < -this.bounds.height) {
                changePoint.y = -this.bounds.height - this.contentOffset.y
            }
            
            this.contentOffset = this.contentOffset.add(changePoint)
            
            this._previousClientPoint = currentClientPoint
            
        }.bind(this))
        
        
        
        
    }
    
    
    
    
    
    invalidateIntrinsicContentFrame() {
        
        this._intrinsicContentFrame = nil
        
    }
    
    
    
    
    
    get contentOffset() {
        
        return this._contentOffset
        
    }
    
    set contentOffset(offset: UIPoint) {
        
        this._contentOffset = offset
        this.setNeedsLayout()
        
    }
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
        
        // var intrinsicContentFrame = this._intrinsicContentFrame;
        // if (!IS(intrinsicContentFrame)) {
        //     intrinsicContentFrame = this.containerView.intrinsicContentSizeWithConstraints();   
        // }
        // intrinsicContentFrame.offsetByPoint(this.contentOffset);
        // intrinsicContentFrame.height = this.containerView.viewHTMLElement.scrollHeight;
        // intrinsicContentFrame.width = this.containerView.viewHTMLElement.scrollWidth;
        // this.containerView.frame = intrinsicContentFrame;
        
        this.containerView.frame = this.containerView.bounds.offsetByPoint(this.contentOffset)
        
        
        
        
    }
    
    
    
    // get _subviews() {
    //     return super.subviews;
    // }
    
    // set _subviews(subviews: UIView[]) {
    //     super.subviews = subviews;
    // }
    
    // get subviews() {
    //     return this.containerView.subviews;
    // }
    
    // set subviews(subviews: UIView[]) {
    
    //     this.containerView.subviews = subviews;
    
    //     this.invalidateIntrinsicContentFrame();
    
    
    // }
    
    
    hasSubview(view: UIView) {
        return this.containerView.hasSubview(view)
    }
    
    addSubview(view: UIView) {
        
        this.containerView.addSubview(view)
        
        this.invalidateIntrinsicContentFrame()
        
        
    }
    
    
    
    
    
}






















































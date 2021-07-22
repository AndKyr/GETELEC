///<reference path="UIView.ts"/>









class UISlideScrollerView extends UIView {
    
    
    
    _previousLayoutBounds: UIRectangle
    _targetIndex: number = 0
    pageIndicatorsView: UIView
    _isAnimating: boolean = NO
    _isAnimationOngoing: boolean = NO
    _animationTimer: UITimer = nil
    _scrollView: UIScrollView
    _slideViews: UIView[] = []
    
    wrapAround: boolean = YES
    
    animationDuration: number = 0.35
    animationDelay: number = 2
    
    _currentPageIndex: number = 0
    
    constructor(elementID: string, viewHTMLElement?: HTMLElement) {
        
        super(elementID, viewHTMLElement)
        
        this._class = UIScrollView
        this.superclass = UIView
        
        
        this._scrollView = new UIScrollView(elementID + "ScrollView")
        
        this.addSubview(this._scrollView)
        
        
        this._scrollView._scrollEnabled = NO
        
        this._scrollView.addTargetForControlEvent(
            UIView.controlEvent.PointerMove,
            function (sender: UIView, event: Event) {
                
                if (event instanceof MouseEvent) {
                    this._animationTimer.invalidate()
                }
                
            }.bind(this)
        )
        
        this._scrollView.addTargetForControlEvent(UIView.controlEvent.PointerLeave, function () {
            
            if (this._isAnimating && event instanceof MouseEvent) {
                this.startAnimating()
            }
            
        }.bind(this))
        
        // Touch events
        this._scrollView.addTargetForControlEvent(UIView.controlEvent.PointerDown, function (sender, event) {
            
            if (event instanceof TouchEvent) {
                this._animationTimer.invalidate()
            }
            
        }.bind(this))
        
        this._scrollView.addTargetForControlEvents([
            UIView.controlEvent.PointerUp, UIView.controlEvent.PointerCancel
        ], function (sender, event) {
            
            if (event instanceof TouchEvent && this._isAnimating) {
                
                this.startAnimating()
                
            }
            
        }.bind(this))
        
        
        // Page indicator
        
        this.pageIndicatorsView = new UIView(elementID + "PageIndicatorsView")
        this.addSubview(this.pageIndicatorsView)
        
        
        
        
        
    }
    
    
    
    
    
    buttonForPageIndicatorWithIndex(index: number): UIButton {
    
        const result = new UIButton(this.viewHTMLElement.id + "PageIndicatorButton" + index)
    
        result.addTargetForControlEvents([
            UIView.controlEvent.PointerUpInside, UIView.controlEvent.EnterUp
        ], function (sender, event) {
            
            this.scrollToPageWithIndex(index, YES)
            
            if (this._isAnimating) {
                
                this.startAnimating()
                
            }
            
            
            
        }.bind(this))
        
        result.addTargetForControlEvent(UIView.controlEvent.PointerMove, function () {
            
            this._animationTimer.invalidate()
            
        }.bind(this))
        
        
        result.updateContentForNormalState = function () {
            
            result.backgroundColor = UIColor.blueColor
            result.titleLabel.textColor = UIColor.whiteColor
            
        }
        
        
        result.frame = new UIRectangle(nil, nil, 20, 50)
        
        // result.style.height = "20px";
        // result.style.width = "50px";
        result.style.display = "table-cell"
        result.style.position = "relative"
        
        
        // var resultContent = new UIView(result.viewHTMLElement.id + "Content");
        // resultContent.backgroundColor = UIColor.whiteColor;
        // resultContent.centerYInContainer();
        // resultContent.style.height = "10px";
        // resultContent.style.height = "100%";
        // resultContent.style.borderRadius = "5px";
        
        // result.addSubview(resultContent);
        
        
        return result
        
    }
    
    
    
    addSlideView(view: UIView) {
        
        this.slideViews.push(view)
        
        this.updateSlideViews()
        
    }
    
    set slideViews(views: UIView[]) {
        
        this._slideViews = views
        
        this.updateSlideViews()
        
    }
    
    get slideViews() {
        return this._slideViews
    }
    
    
    
    get currentPageIndex() {
    
        const result = this._currentPageIndex
    
        return result
        
    }
    
    set currentPageIndex(index: number) {
        
        this._currentPageIndex = index
        
        this._slideViews[index].willAppear()
        
        //this._scrollView.contentOffset.x = -this._slideViews[index].frame.min.x; //-this.bounds.width * index;
        //this._scrollView.contentOffset.x = Math.round(this._scrollView.contentOffset.x);
        
        this._scrollView.contentOffset = this._scrollView.contentOffset.pointWithX(-this._slideViews[index].frame.min.x)
        
        
        this.pageIndicatorsView.subviews.forEach(function (button: UIButton, index, array) {
            
            button.selected = NO
            
        });
        
        
        
        (this.pageIndicatorsView.subviews[index] as UIButton).selected = YES
        
        
        
    }
    
    
    
    scrollToPreviousPage(animated: boolean) {
        
        if (this.slideViews.length == 0) {
            return
        }
    
        var targetIndex = this.currentPageIndex
    
        if (this.wrapAround) {
            targetIndex = (this.currentPageIndex - 1) % (this.slideViews.length)
        }
        else if (this.currentPageIndex - 1 < this.slideViews.length) {
            targetIndex = this.currentPageIndex - 1
        }
        else {
            return
        }
        
        this.scrollToPageWithIndex(targetIndex, animated)
        
    }
    
    scrollToNextPage(animated: boolean) {
        
        if (this.slideViews.length == 0) {
            return
        }
    
        var targetIndex = this.currentPageIndex
    
        if (this.wrapAround) {
            targetIndex = (this.currentPageIndex + 1) % (this.slideViews.length)
        }
        else if (this.currentPageIndex + 1 < this.slideViews.length) {
            targetIndex = this.currentPageIndex + 1
        }
        else {
            return
        }
        
        this.scrollToPageWithIndex(targetIndex, animated)
        
    }
    
    
    
    scrollToPageWithIndex(targetIndex: number, animated: boolean = YES) {
        
        this._targetIndex = targetIndex
        
        // this._slideViews[this.currentPageIndex]._shouldLayout = NO;
        // this._slideViews[this._targetIndex]._shouldLayout = YES;
        
        //this._slideViews[this._targetIndex].hidden = NO;
        
        this.willScrollToPageWithIndex(targetIndex)
        
        this._isAnimationOngoing = YES
        
        //var previousView = this._slideViews[this.currentPageIndex];
        
        if (animated) {
            
            
            
            UIView.animateViewOrViewsWithDurationDelayAndFunction(
                this._scrollView.containerView,
                this.animationDuration,
                0,
                undefined,
                function (this: UISlideScrollerView) {
                    
                    
                    this.currentPageIndex = targetIndex
                    
                    
                    
                }.bind(this),
                function (this: UISlideScrollerView) {
                    
                    this.didScrollToPageWithIndex(targetIndex)
                    
                    this._isAnimationOngoing = NO
                    
                    //previousView.hidden = YES;
                    
                }.bind(this)
            )
            
        }
        else {
            
            
            this.currentPageIndex = targetIndex
            this.didScrollToPageWithIndex(targetIndex)
            
            //previousView.hidden = YES;
            
        }
        
    }
    
    
    
    willScrollToPageWithIndex(index: number) {
    
        const targetView = this.slideViews[index]
    
        if (IS(targetView) && (targetView as any).willAppear && (targetView as any).willAppear instanceof Function) {
            
            (targetView as any).willAppear()
            
        }
        
    }
    
    didScrollToPageWithIndex(index: number) {
        
        
        
    }
    
    
    
    startAnimating() {
        
        this._isAnimating = YES
        
        this._animationTimer.invalidate()
        
        this._animationTimer = new UITimer(this.animationDelay + this.animationDuration, YES, function () {
            
            this.scrollToNextPage(YES)
            
        }.bind(this))
        
    }
    
    stopAnimating() {
        
        this._isAnimating = NO
        this._animationTimer.invalidate()
        
    }
    
    
    
    
    
    updateSlideViews() {
        
        this._scrollView.containerView.subviews.slice().forEach(function (subview, index, array) {
            
            subview.removeFromSuperview()
            
        })
        
        this.pageIndicatorsView.subviews.slice().forEach(function (subview, index, array) {
            
            subview.removeFromSuperview()
            
        })
        
        this._slideViews.forEach(function (view, index, array) {
            
            this._scrollView.addSubview(view)
            
            this.pageIndicatorsView.addSubview(this.buttonForPageIndicatorWithIndex(index))
            
        }.bind(this))
        
    }
    
    
    didReceiveBroadcastEvent(event: UIViewBroadcastEvent) {
        
        super.didReceiveBroadcastEvent(event)
        
        if (event.name == UICore.broadcastEventName.WindowDidResize) {
            
            this.currentPageIndex = this.currentPageIndex
            
            
        }
        
        
    }
    
    
    set frame(frame: UIRectangle) {
        
        super.frame = frame
        
        this.currentPageIndex = this.currentPageIndex
        
    }
    
    get frame() {
        
        return super.frame
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
        
        if (this.bounds.isEqualTo(this._previousLayoutBounds)) {
            return
        }
    
        const bounds = this.bounds
    
        this._previousLayoutBounds = bounds
        
        this._scrollView.frame = bounds
        
        
        
        this._scrollView.containerView.frame = bounds.rectangleWithWidth(bounds.width *
            this.slideViews.length).performFunctionWithSelf(function (this: UISlideScrollerView, self: UIRectangle) {
            
            self.offsetByPoint(this._scrollView.contentOffset)
            
            return self
            
        }.bind(this))
        
        this._slideViews.forEach(function (view, index, array) {
            
            view.frame = bounds.rectangleWithX((this.bounds.width + 1) * index)
            
        }.bind(this))
        
        
        this.layoutPageIndicators()
        
        
    }
    
    
    layoutPageIndicators() {
        
        this.pageIndicatorsView.centerXInContainer()
        this.pageIndicatorsView.style.bottom = "20px"
        this.pageIndicatorsView.style.height = "20px"
        this.pageIndicatorsView.style.display = "table-row"
        
    }
    
    
    
    
    removeFromSuperview() {
        
        
        super.removeFromSuperview()
        
        this.stopAnimating()
        
        
    }
    
    
    
    
    
}









































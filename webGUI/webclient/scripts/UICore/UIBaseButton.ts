/// <reference path="./UIView.ts" />





class UIBaseButton extends UIView {
    
    _selected: boolean = NO
    _highlighted: boolean = NO
    
    _isPointerInside: boolean
    
    
    _isToggleable: boolean = NO
    _hovered: boolean
    _focused: boolean
    
    
    
    
    
    constructor(elementID: string, elementType?: string, initViewData?: any) {
        
        super(elementID, nil, elementType, initViewData)
        
        this._class = UIButton
        this.superclass = UIView
        
        this.initViewStateControl()
        
    }
    
    
    
    
    
    initViewStateControl() {
        
        this.class.superclass = UIView
        
        // Instance variables
        
        
        this._isPointerInside = NO
    
    
        const setHovered = function () {
            this.hovered = YES
        }.bind(this)
        this.addTargetForControlEvent(UIView.controlEvent.PointerHover, setHovered)
    
        const setNotHovered = function () {
        
            this.hovered = NO
        
        }.bind(this)
    
        this.addTargetForControlEvents([
            UIView.controlEvent.PointerLeave, UIView.controlEvent.PointerCancel, UIView.controlEvent.MultipleTouches
        ], setNotHovered)
    
    
        var highlightingTime
        const setHighlighted = function () {
            this.highlighted = YES
            highlightingTime = Date.now()
        }.bind(this)
        this.addTargetForControlEvent(UIView.controlEvent.PointerDown, setHighlighted)
        this.addTargetForControlEvent(UIView.controlEvent.PointerEnter, setHighlighted)
    
        const setNotHighlighted = function () {
            this.highlighted = NO
        }.bind(this)
        const setNotHighlightedWithMinimumDuration = function () {
            const minimumDurationInMilliseconds = 50
            const elapsedTime = Date.now() - highlightingTime
            if (minimumDurationInMilliseconds < elapsedTime) {
                this.highlighted = NO
            }
            else {
                setTimeout(function () {
                    this.highlighted = NO
                }.bind(this), minimumDurationInMilliseconds - elapsedTime)
            }
        }.bind(this)
        this.addTargetForControlEvents([
            UIView.controlEvent.PointerLeave, UIView.controlEvent.PointerCancel, UIView.controlEvent.MultipleTouches
        ], setNotHighlighted)
        this.addTargetForControlEvent(UIView.controlEvent.PointerUp, setNotHighlightedWithMinimumDuration)
        
        // Handle enter key press
        this.addTargetForControlEvent(UIView.controlEvent.EnterDown, function () {
            
            setHighlighted()
            setNotHighlightedWithMinimumDuration()
            
        })
        
        
        this.addTargetForControlEvent(
            UIView.controlEvent.Focus,
            function (this: UIBaseButton, sender: UIView, event: Event) {
                
                this.focused = YES
                
            }.bind(this)
        )
        
        this.addTargetForControlEvent(
            UIView.controlEvent.Blur,
            function (this: UIBaseButton, sender: UIView, event: Event) {
                
                this.focused = NO
                
            }.bind(this)
        )
        
        
        this.updateContentForCurrentState()
        
        this.pausesPointerEvents = YES
        this.tabIndex = 1
        
        this.style.cursor = "pointer"
        
        //this.style.outline = "none";
        
        
        this.nativeSelectionEnabled = NO
        
        
        this.addTargetForControlEvents([
            UIView.controlEvent.EnterDown, UIView.controlEvent.PointerUpInside
        ], function (this: UIBaseButton, sender, event) {
            
            if (this.isToggleable) {
                
                this.toggleSelectedState()
                
            }
            
        }.bind(this))
        
        
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
            this.backgroundColor = UIColor.nilColor
        }
        else {
            updateFunction.call(this)
        }
        
    }
    
    updateContentForNormalState() {
        
        
        
    }
    
    updateContentForHoveredState() {
        
        this.updateContentForNormalState()
        
    }
    
    updateContentForFocusedState() {
        
        this.updateContentForHoveredState()
        
    }
    
    updateContentForHighlightedState() {
        
        
        
    }
    
    updateContentForSelectedState() {
        
        
        
    }
    
    updateContentForSelectedAndHighlightedState() {
        
        this.updateContentForSelectedState()
        
    }
    
    
    set enabled(enabled: boolean) {
        
        super.enabled = enabled
        
        this.updateContentForCurrentEnabledState()
        
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
    
    
    
    
    
    didReceiveBroadcastEvent(event: UIViewBroadcastEvent) {
        
        super.didReceiveBroadcastEvent(event)
        
        if (event.name == UIView.broadcastEventName.PageDidScroll || event.name ==
            UIView.broadcastEventName.AddedToViewTree) {
            
            this.hovered = NO
            
            this.highlighted = NO
            
            
        }
        
        
    }
    
    
    
    
    
    toggleSelectedState() {
        
        
        this.selected = !this.selected
        
        
    }
    
    set isToggleable(isToggleable: boolean) {
        
        this._isToggleable = isToggleable
        
    }
    
    get isToggleable() {
        
        return this._isToggleable
        
    }
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
    
        const bounds = this.bounds
    
    
    
    
    }
    
    
    
    
    
    sendControlEventForKey(eventKey: string, nativeEvent: Event) {
        
        if (eventKey == UIView.controlEvent.PointerUpInside && !this.highlighted) {
            
            // Do not send the event in this case
            //super.sendControlEventForKey(eventKey, nativeEvent);
    
            const asd = 1
    
        }
        else {
            
            super.sendControlEventForKey(eventKey, nativeEvent)
            
        }
        
    }
    
    
    
    
    
    static getEventCoordinatesInDocument(touchOrMouseEvent) {
        // http://www.quirksmode.org/js/events_properties.html
        var posx = 0
        var posy = 0
        var e = touchOrMouseEvent
        if (!e) {
            e = window.event
        }
        if (e.pageX || e.pageY) {
            posx = e.pageX
            posy = e.pageY
        }
        else if (e.clientX || e.clientY) {
            posx = e.clientX + document.body.scrollLeft
                + document.documentElement.scrollLeft
            posy = e.clientY + document.body.scrollTop
                + document.documentElement.scrollTop
        }
        // posx and posy contain the mouse position relative to the document
    
        const coordinates = { "x": posx, "y": posy }
    
        return coordinates
        
    }
    
    
    
    static getElementPositionInDocument(el) {
        //https://www.kirupa.com/html5/getting_mouse_click_position.htm
        var xPosition = 0
        var yPosition = 0
    
        while (el) {
            if (el.tagName == "BODY") {
                
                // Coordinates in document are coordinates in body, therefore subtracting the scroll position of the body is not needed
                
                //      // deal with browser quirks with body/window/document and page scroll
                //      var xScrollPos = el.scrollLeft || document.documentElement.scrollLeft;
                //      var yScrollPos = el.scrollTop || document.documentElement.scrollTop;
                //
                //      xPosition += (el.offsetLeft - xScrollPos + el.clientLeft);
                //      yPosition += (el.offsetTop - yScrollPos + el.clientTop);
            }
            else {
                xPosition += (el.offsetLeft - el.scrollLeft + el.clientLeft)
                yPosition += (el.offsetTop - el.scrollTop + el.clientTop)
            }
            
            el = el.offsetParent
        }
        return {
            x: xPosition,
            y: yPosition
        }
    }
    
    static convertCoordinatesFromDocumentToElement(x, y, element) {
        const elementPositionInDocument = this.getElementPositionInDocument(element)
        const coordinatesInElement = { "x": x - elementPositionInDocument.x, "y": y - elementPositionInDocument.y }
        return coordinatesInElement
    }
    
    static getEventCoordinatesInElement(touchOrMouseEvent, element) {
        const coordinatesInDocument = this.getEventCoordinatesInDocument(touchOrMouseEvent)
        const coordinatesInElement = this.convertCoordinatesFromDocumentToElement(
            coordinatesInDocument.x,
            coordinatesInDocument.y,
            element
        )
        return coordinatesInElement
    }
    
    
    
    
    
}






























































































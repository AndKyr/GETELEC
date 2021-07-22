/// <reference path="./UIBaseButton.ts" />



class UILink extends UIBaseButton {
    
    private _colors: UIButtonColorSpecifier
    
    
    constructor(elementID: string, initViewData = nil) {
        
        super(elementID, "a", initViewData)
        
        this._class = UILink
        this.superclass = UIBaseButton
        
        this.stopsPointerEventPropagation = NO
        
        this.pausesPointerEvents = NO
        
        
        
    }
    
    
    
    
    
    initView(elementID, viewHTMLElement, initViewData: { titleType: string }) {
        
        super.initView(elementID, viewHTMLElement, initViewData)
        
        this.class.superclass = UIBaseButton
        
        // Instance variables
        
        
        //this.style.position = "relative"
        
        
        viewHTMLElement.onclick = this.blur.bind(this)
        
        
        
        
    }
    
    
    
    
    
    get colors(): UIButtonColorSpecifier {
        return this._colors
    }
    
    set colors(value: UIButtonColorSpecifier) {
        this._colors = value
    }
    
    get viewHTMLElement() {
        
        return super.viewHTMLElement as HTMLLinkElement
        
    }
    
    set text(text: string) {
        
        this.viewHTMLElement.textContent = text
        
    }
    
    get text() {
        
        return this.viewHTMLElement.textContent
        
    }
    
    
    set target(target: string) {
        
        this.viewHTMLElement.setAttribute("href", target)
        
    }
    
    get target() {
    
        const result = this.viewHTMLElement.getAttribute("href")
    
        return result
        
    }
    
    
    
    
    
    set targetRouteForCurrentState(targetRouteForCurrentState: () => (UIRoute | string)) {
        
        this._targetRouteForCurrentState = targetRouteForCurrentState
        
        this.updateTarget()
        
    }
    
    get targetRouteForCurrentState() {
        
        return this._targetRouteForCurrentState
        
    }
    
    
    _targetRouteForCurrentState() {
    
        const result = UIRoute.currentRoute.routeByRemovingComponentsOtherThanOnesNamed(["settings"]) as (UIRoute | string)
    
        return result
        
    }
    
    
    
    
    
    didReceiveBroadcastEvent(event: UIViewBroadcastEvent) {
        
        super.didReceiveBroadcastEvent(event)
        
        if (event.name == UICore.broadcastEventName.RouteDidChange) {
            
            this.updateTarget()
            
        }
        
    }
    
    
    
    
    
    wasAddedToViewTree() {
        
        super.wasAddedToViewTree()
        
        this.updateTarget()
        
        
    }
    
    
    
    
    
    updateTarget() {
        
        const route = this.targetRouteForCurrentState()
        
        if (route instanceof UIRoute) {
            
            this.target = route.linkRepresentation
            
            return
            
        }
        
        this.target = route
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
    
        const bounds = this.bounds
    
    }
    
    
    
    
    
}






























































































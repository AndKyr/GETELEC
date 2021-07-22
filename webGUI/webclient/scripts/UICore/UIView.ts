///// <reference path="../autolayout/src/AutoLayout.js" />
/// <reference path="./UIObject.ts" />
/// <reference path="./UIColor.ts" />
/// <reference path="./UIRectangle.ts" />

//import Attribute from "../autolayout/src/Attribute";

//import AutoLayout from "../autolayout/src/AutoLayout";





declare module AutoLayout {
    
    
    class Constraint {
        
        [key: string]: any
        
    }
    
    
    class View {
        
        [key: string]: any
        
    }
    
    
    class VisualFormat {
        
        static parse(arg0: any, arg1: any): any;
        
        [key: string]: any
        
    }
    
    
    enum Attribute {
        
        LEFT, RIGHT, BOTTOM, TOP, CENTERX, CENTERY, WIDTH, HEIGHT, ZINDEX, VARIABLE, NOTANATTRIBUTE
        
    }
    
    
    enum Relation {
        
        EQU, LEQ, GEQ
        
    }
    
    
}





interface LooseObject {
    [key: string]: any
}


interface ControlEventTargetsObject {
    
    [key: string]: Function[];
    
}


interface UIViewBroadcastEvent {
    
    name: string;
    parameters: {
        [key: string]: string | string[];
    }
    
}



type UIViewAddControlEventTargetObject<T extends typeof UIView.controlEvent> = {
    
    [K in keyof T]: ((sender: UIView, event: Event) => void) & Partial<UIViewAddControlEventTargetObject<T>>
    
}





class UIView extends UIObject {
    
    _nativeSelectionEnabled: boolean = YES
    _shouldLayout: boolean
    _UITableViewRowIndex: number
    _UITableViewReusabilityIdentifier: any
    _UIViewIntrinsicTemporaryWidth: string
    _UIViewIntrinsicTemporaryHeight: string
    _enabled: boolean = YES
    _frame: any
    _backgroundColor: UIColor = UIColor.transparentColor
    
    _viewHTMLElement: HTMLElement & LooseObject
    
    _innerHTMLKey: string
    _defaultInnerHTML: string
    _parameters: { [x: string]: (string | UILocalizedTextObject) }
    
    _localizedTextObject: UILocalizedTextObject = nil
    
    _controlEventTargets: ControlEventTargetsObject = {} //{ "PointerDown": Function[]; "PointerMove": Function[]; "PointerLeave": Function[]; "PointerEnter": Function[]; "PointerUpInside": Function[]; "PointerUp": Function[]; "PointerHover": Function[]; };
    _frameTransform: string
    _viewControllerLayoutFunction: () => void = nil
    _didLayoutSubviewsDelegateFunction: () => void
    _didReceiveBroadcastEventDelegateFunction: (event: UIViewBroadcastEvent) => void
    _updateLayoutFunction: any
    // @ts-ignore
    _constraints: any[] //AutoLayout.Constraint[];
    superview: UIView
    subviews: UIView[]
    _styleClasses: any[]
    _isHidden: boolean = NO
    
    pausesPointerEvents: boolean = NO
    stopsPointerEventPropagation: boolean = YES
    _isPointerInside: boolean
    _isPointerValid: boolean
    _initialPointerPosition: UIPoint
    _hasPointerDragged: boolean
    _pointerDragThreshold = 2
    
    ignoresTouches: boolean = NO
    ignoresMouse: boolean = NO
    
    
    static _UIViewIndex: number = -1
    _UIViewIndex: number
    
    static _viewsToLayout: UIView[] = []
    
    forceIntrinsicSizeZero: boolean = NO
    _touchEventTime: number
    
    static _pageScale = 1
    
    constructor(
        elementID: string = ("UIView" +
            UIView.nextIndex),
        viewHTMLElement: HTMLElement & LooseObject = null,
        elementType: string = null,
        initViewData?: any
    ) {
        
        super()
        
        this._class = UIView
        this.superclass = UIObject
        
        // Instance variables
        
        UIView._UIViewIndex = UIView.nextIndex
        this._UIViewIndex = UIView._UIViewIndex
        
        this._styleClasses = []
        // Object.defineProperty(this, "styleClasses", { get: this.styleClasses, set: this.setStyleClasses });
        // Object.defineProperty(this, "styleClassName", { get: this.styleClassName });
        
        this._initViewHTMLElement(elementID, viewHTMLElement, elementType)
        
        this.subviews = []
        this.superview = nil
        
        // Object.defineProperty(this, "elementID", { get: this.elementID });
        
        // Object.defineProperty(this, "constraints", { get: this.constraints, set: this.setConstraints });
        this._constraints = []
        
        this._updateLayoutFunction = nil
        
        //Object.defineProperty(this, "backgroundColor", { get: this.backgroundColor, set: this.setBackgroundColor });
        //this.backgroundColor = "transparent";
        
        // Object.defineProperty(this, "alpha", { get: this.alpha, set: this.setAlpha });
        
        // Object.defineProperty(this, "frame", { get: this.frame, set: this.setFrame });
        // Object.defineProperty(this, "bounds", { get: this.bounds, set: this.setBounds });
        
        // Object.defineProperty(this, "userInteractionEnabled", { get: this.userInteractionEnabled, set: this.setUserInteractionEnabled });
        
        // this._controlEventTargets = {
        //     "PointerDown": [],
        //     "PointerMove": [],
        //     "PointerLeave": [],
        //     "PointerEnter": [],
        //     "PointerUpInside": [],
        //     "PointerUp": [],
        //     "PointerHover": []
        // }
        
        this._didLayoutSubviewsDelegateFunction = function () {
        }
        
        this._frameTransform = ""
        
        this.initView(this.viewHTMLElement.id, this.viewHTMLElement, initViewData)
        
        this._initViewCSSSelectorsIfNeeded()
        
        this._loadUIEvents()
        
        
        this.setNeedsLayout()
        
        
    }
    
    
    static get nextIndex() {
        
        return UIView._UIViewIndex + 1
        
    }
    
    static get pageHeight() {
        
        const body = document.body
        const html = document.documentElement
        
        const height = Math.max(
            body.scrollHeight,
            body.offsetHeight,
            html.clientHeight,
            html.scrollHeight,
            html.offsetHeight
        )
        
        return height
        
    }
    
    static get pageWidth() {
        
        const body = document.body
        const html = document.documentElement
        
        const width = Math.max(body.scrollWidth, body.offsetWidth, html.clientWidth, html.scrollWidth, html.offsetWidth)
        
        return width
        
    }
    
    
    
    
    
    initView(elementID: string, viewHTMLElement: HTMLElement, initViewData?: any) {
    
    
    
    }
    
    
    
    
    
    centerInContainer() {
        this.style.left = "50%"
        this.style.top = "50%"
        this.style.transform = "translateX(-50%) translateY(-50%)"
    }
    
    centerXInContainer() {
        this.style.left = "50%"
        this.style.transform = "translateX(-50%)"
    }
    
    centerYInContainer() {
        this.style.top = "50%"
        this.style.transform = "translateY(-50%)"
    }
    
    
    
    _initViewHTMLElement(elementID, viewHTMLElement, elementType = "div") {
        
        
        if (!IS(elementType)) {
            
            elementType = "div"
            
        }
        
        if (!IS(viewHTMLElement)) {
            
            this._viewHTMLElement = this.createElement(elementID, elementType)
            
            this.style.position = "absolute"
            this.style.margin = "0"
            
        }
        else {
            
            this._viewHTMLElement = viewHTMLElement
            
        }
        
        if (IS(elementID)) {
            
            
            this.viewHTMLElement.id = elementID
            
            
        }
        
        
        this.viewHTMLElement.obeyAutolayout = YES
        
        this.viewHTMLElement.UIView = this
        
        this.addStyleClass(this.styleClassName)
        
    }
    
    
    
    set nativeSelectionEnabled(selectable: boolean) {
        this._nativeSelectionEnabled = selectable
        if (!selectable) {
            this.style.cssText = this.style.cssText +
                " -webkit-touch-callout: none; -webkit-user-select: none; -khtml-user-select: none; -moz-user-select: none; -ms-user-select: none; user-select: none;"
        }
        else {
            this.style.cssText = this.style.cssText +
                " -webkit-touch-callout: text; -webkit-user-select: text; -khtml-user-select: text; -moz-user-select: text; -ms-user-select: text; user-select: text;"
        }
    }
    
    
    get nativeSelectionEnabled() {
        return this._nativeSelectionEnabled
    }
    
    
    
    get styleClassName() {
        
        const result = "UICore_UIView_" + this.class.name
        
        return result
        
    }
    
    _initViewCSSSelectorsIfNeeded() {
        
        if (!this.class._areViewCSSSelectorsInitialized) {
            
            this.initViewStyleSelectors()
            
            this.class._areViewCSSSelectorsInitialized = YES
            
        }
        
    }
    
    initViewStyleSelectors() {
        
        // Override this in a subclass
        
    }
    
    initStyleSelector(selector, style) {
        
        const styleRules = UIView.getStyleRules(selector)
        
        if (!styleRules) {
            
            UIView.createStyleSelector(selector, style)
            
        }
        
    }
    
    
    createElement(elementID, elementType) {
        var result = document.getElementById(elementID)
        if (!result) {
            result = document.createElement(elementType)
        }
        return result
    }
    
    public get viewHTMLElement() {
        return this._viewHTMLElement
    }
    
    public get elementID() {
        
        return this.viewHTMLElement.id
        
    }
    
    
    setInnerHTML(key: string, defaultString: string, parameters?: { [x: string]: string | UILocalizedTextObject }) {
        
        this._innerHTMLKey = key
        this._defaultInnerHTML = defaultString
        this._parameters = parameters
        
        const languageName = UICore.languageService.currentLanguageKey
        const result = UICore.languageService.stringForKey(key, languageName, defaultString, parameters)
        
        this.innerHTML = result
        
    }
    
    
    _setInnerHTMLFromKeyIfPossible() {
        
        if (this._innerHTMLKey && this._defaultInnerHTML) {
            
            this.setInnerHTML(this._innerHTMLKey, this._defaultInnerHTML, this._parameters)
            
        }
        
    }
    
    _setInnerHTMLFromLocalizedTextObjectIfPossible() {
        
        if (IS(this._localizedTextObject)) {
            
            this.innerHTML = UICore.languageService.stringForCurrentLanguage(this._localizedTextObject)
            
        }
        
    }
    
    
    get localizedTextObject() {
        
        return this._localizedTextObject
        
    }
    
    set localizedTextObject(localizedTextObject: UILocalizedTextObject) {
        
        this._localizedTextObject = localizedTextObject
        
        this._setInnerHTMLFromLocalizedTextObjectIfPossible()
        
    }
    
    
    get innerHTML() {
        return this.viewHTMLElement.innerHTML
    }
    
    
    set innerHTML(innerHTML) {
        
        if (this.innerHTML != innerHTML) {
            
            this.viewHTMLElement.innerHTML = FIRST(innerHTML, "")
            
        }
        
    }
    
    
    
    set hoverText(hoverText: string) {
        this.viewHTMLElement.setAttribute("title", hoverText)
    }
    
    get hoverText() {
        return this.viewHTMLElement.getAttribute("title")
    }
    
    
    get scrollSize() {
        
        const result = new UIRectangle(0, 0, this.viewHTMLElement.scrollHeight, this.viewHTMLElement.scrollWidth)
        
        return result
        
    }
    
    
    
    get dialogView(): UIDialogView {
        if (!IS(this.superview)) {
            return nil
        }
        if (!(this instanceof UIDialogView)) {
            return this.superview.dialogView
        }
        return this
    }
    
    
    get rootView(): UIView {
        if (IS(this.superview)) {
            return this.superview.rootView
        }
        return this
    }
    
    
    public set enabled(enabled: boolean) {
        this._enabled = enabled
        this.updateContentForCurrentEnabledState()
    }
    
    public get enabled(): boolean {
        return this._enabled
    }
    
    updateContentForCurrentEnabledState() {
        this.hidden = !this.enabled
        this.userInteractionEnabled = this.enabled
    }
    
    
    
    
    public get tabIndex(): number {
        
        return Number(this.viewHTMLElement.getAttribute("tabindex"))
        
    }
    
    
    public set tabIndex(index: number) {
        
        this.viewHTMLElement.setAttribute("tabindex", "" + index)
        
    }
    
    
    
    
    
    get styleClasses() {
        
        return this._styleClasses
        
    }
    
    set styleClasses(styleClasses) {
        
        this._styleClasses = styleClasses
        
    }
    
    hasStyleClass(styleClass) {
        
        // This is for performance reasons
        if (!IS(styleClass)) {
            return NO
        }
        
        const index = this.styleClasses.indexOf(styleClass)
        if (index > -1) {
            return YES
        }
        return NO
        
    }
    
    
    
    addStyleClass(styleClass: string) {
        
        if (!IS(styleClass)) {
            return
        }
        
        if (!this.hasStyleClass(styleClass)) {
            this._styleClasses.push(styleClass)
        }
        
    }
    
    removeStyleClass(styleClass: string) {
        
        // This is for performance reasons
        if (!IS(styleClass)) {
            return
        }
        
        const index = this.styleClasses.indexOf(styleClass)
        if (index > -1) {
            
            this.styleClasses.splice(index, 1)
            
        }
        
        
    }
    
    
    
    static findViewWithElementID(elementID: string): UIView {
        const viewHTMLElement = document.getElementById(elementID)
        if (IS_NOT(viewHTMLElement)) {
            return nil
        }
        // @ts-ignore
        const result = viewHTMLElement.UIView
        return result
    }
    
    
    
    
    static createStyleSelector(selector, style) {
        
        return
        
        // @ts-ignore
        if (!document.styleSheets) {
            return
        }
        if (document.getElementsByTagName("head").length == 0) {
            return
        }
        
        var styleSheet
        var mediaType
        
        if (document.styleSheets.length > 0) {
            for (var i = 0, l: any = document.styleSheets.length; i < l; i++) {
                if (document.styleSheets[i].disabled) {
                    continue
                }
                const media = document.styleSheets[i].media
                mediaType = typeof media
                
                if (mediaType === "string") {
                    if (media as any === "" || ((media as any).indexOf("screen") !== -1)) {
                        styleSheet = document.styleSheets[i]
                    }
                }
                else if (mediaType == "object") {
                    if (media.mediaText === "" || (media.mediaText.indexOf("screen") !== -1)) {
                        styleSheet = document.styleSheets[i]
                    }
                }
                
                if (typeof styleSheet !== "undefined") {
                    break
                }
            }
        }
        
        if (typeof styleSheet === "undefined") {
            const styleSheetElement = document.createElement("style")
            styleSheetElement.type = "text/css"
            document.getElementsByTagName("head")[0].appendChild(styleSheetElement)
            
            for (i = 0; i < document.styleSheets.length; i++) {
                if (document.styleSheets[i].disabled) {
                    continue
                }
                styleSheet = document.styleSheets[i]
            }
            
            mediaType = typeof styleSheet.media
        }
        
        if (mediaType === "string") {
            for (var i = 0, l = styleSheet.rules.length; i < l; i++) {
                if (styleSheet.rules[i].selectorText && styleSheet.rules[i].selectorText.toLowerCase() ==
                    selector.toLowerCase()) {
                    styleSheet.rules[i].style.cssText = style
                    return
                }
            }
            styleSheet.addRule(selector, style)
        }
        else if (mediaType === "object") {
            
            var styleSheetLength = 0
            
            try {
                
                styleSheetLength = (styleSheet.cssRules) ? styleSheet.cssRules.length : 0
                
            } catch (error) {
                
            }
            
            
            for (var i = 0; i < styleSheetLength; i++) {
                if (styleSheet.cssRules[i].selectorText && styleSheet.cssRules[i].selectorText.toLowerCase() ==
                    selector.toLowerCase()) {
                    styleSheet.cssRules[i].style.cssText = style
                    return
                }
            }
            styleSheet.insertRule(selector + "{" + style + "}", styleSheetLength)
        }
    }
    
    static getStyleRules(selector) {
        var selector = selector.toLowerCase()
        for (var i = 0; i < document.styleSheets.length; i++) {
            const styleSheet = document.styleSheets[i] as any
            var styleRules
            
            try {
                
                styleRules = styleSheet.cssRules ? styleSheet.cssRules : styleSheet.rules
                
            } catch (error) {
                
            }
            
            return styleRules
        }
    }
    
    
    
    get style() {
        return this.viewHTMLElement.style
    }
    
    get computedStyle() {
        return getComputedStyle(this.viewHTMLElement)
    }
    
    public get hidden(): boolean {
        return this._isHidden
    }
    
    
    public set hidden(v: boolean) {
        
        this._isHidden = v
        
        if (this._isHidden) {
            this.style.visibility = "hidden"
        }
        else {
            this.style.visibility = "visible"
        }
        
        
    }
    
    static set pageScale(scale: number) {
        
        UIView._pageScale = scale
        
        const zoom = scale
        const width = 100 / zoom
        const viewHTMLElement = UICore.main.rootViewController.view.viewHTMLElement
        viewHTMLElement.style.transformOrigin = "left top"
        viewHTMLElement.style.transform = "scale(" + zoom + ")"
        viewHTMLElement.style.width = width + "%"
        
    }
    
    static get pageScale() {
        
        return UIView._pageScale
        
    }
    
    
    
    
    
    calculateAndSetViewFrame() {
        
        // Use this method to calculate the frame for the view itself
        
        // This can be used when adding subviews to existing views like buttons
        
    }
    
    
    
    
    
    public get frame(): UIRectangle {
        
        // var result = new UIRectangle(1 * this.viewHTMLElement.offsetLeft, 1 * this.viewHTMLElement.offsetTop, 1 * this.viewHTMLElement.offsetHeight, 1 * this.viewHTMLElement.offsetWidth);
        
        // result.zIndex = 1 * this.style.zIndex;
        
        var result = this._frame
        
        if (!result) {
            
            result = new UIRectangle(1 * this.viewHTMLElement.offsetLeft, 1 * this.viewHTMLElement.offsetTop, 1 *
                this.viewHTMLElement.offsetHeight, 1 * this.viewHTMLElement.offsetWidth)
            result.zIndex = 0
            
        }
        
        return result.copy()
        
    }
    
    public set frame(rectangle: UIRectangle) {
        
        if (IS(rectangle)) {
            this.setFrame(rectangle)
        }
        
    }
    
    setFrame(rectangle, zIndex = 0, performUncheckedLayout = NO) {
        
        
        const frame = this._frame || new UIRectangle(nil, nil, nil, nil)
        
        if (zIndex != undefined) {
            rectangle.zIndex = zIndex
        }
        this._frame = rectangle
        
        // This is useless because frames are copied
        // frame.didChange = function () {
        //     // Do nothing
        // }
        // rectangle.didChange = function () {
        //     this.frame = rectangle;
        // }.bind(this);
        
        if (frame && frame.isEqualTo(rectangle) && !performUncheckedLayout) {
            return
        }
        
        
        UIView._setAbsoluteSizeAndPosition(
            this.viewHTMLElement,
            rectangle.topLeft.x,
            rectangle.topLeft.y,
            rectangle.width,
            rectangle.height,
            rectangle.zIndex
        )
        
        
        if (frame.height != rectangle.height || frame.width != rectangle.width || performUncheckedLayout) {
            
            this.setNeedsLayout()
            
            this.boundsDidChange()
            
            //this.layoutSubviews();
            
        }
        
        
        
    }
    
    
    
    get bounds() {
        
        var result: UIRectangle
        
        // if (IS_NOT(this._frame) && this.style.height == "" && this.style.width  == "" && this.style.left == "" && this.style.right == "" && this.style.bottom == "" && this.style.top == "") {
        
        //     result = new UIRectangle(0, 0, 0, 0)
        
        // }
        // else
        if (IS_NOT(this._frame)) {
            
            result = new UIRectangle(0, 0, 1 * this.viewHTMLElement.offsetHeight, 1 * this.viewHTMLElement.offsetWidth)
            
        }
        else {
            
            result = this.frame.copy()
            
            result.x = 0
            result.y = 0
            
        }
        
        return result
        
    }
    
    set bounds(rectangle) {
        
        const frame = this.frame
        
        this.frame = new UIRectangle(frame.topLeft.x, frame.topLeft.y, rectangle.height, rectangle.width)
        
    }
    
    
    boundsDidChange() {
        
        
        
    }
    
    
    setPosition(
        left: number | string = nil,
        right: number | string = nil,
        bottom: number | string = nil,
        top: number | string = nil,
        height: number | string = nil,
        width: number | string = nil
    ) {
        
        const previousBounds = this.bounds
        
        this.setStyleProperty("left", left)
        this.setStyleProperty("right", right)
        this.setStyleProperty("bottom", bottom)
        this.setStyleProperty("top", top)
        this.setStyleProperty("height", height)
        this.setStyleProperty("width", width)
        
        const bounds = this.bounds
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout()
            this.boundsDidChange()
        }
        
    }
    
    setSizes(height?: number | string, width?: number | string) {
        
        const previousBounds = this.bounds
        
        this.setStyleProperty("height", height)
        this.setStyleProperty("width", width)
        
        const bounds = this.bounds
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout()
            this.boundsDidChange()
        }
        
    }
    
    setMinSizes(height?: number | string, width?: number | string) {
        
        const previousBounds = this.bounds
        
        this.setStyleProperty("minHeight", height)
        this.setStyleProperty("minWidth", width)
        
        const bounds = this.bounds
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout()
            this.boundsDidChange()
        }
        
    }
    
    setMaxSizes(height?: number | string, width?: number | string) {
        
        const previousBounds = this.bounds
        
        this.setStyleProperty("maxHeight", height)
        this.setStyleProperty("maxWidth", width)
        
        const bounds = this.bounds
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout()
            this.boundsDidChange()
        }
        
    }
    
    setMargin(margin?: number | string) {
        
        const previousBounds = this.bounds
        
        this.setStyleProperty("margin", margin)
        
        const bounds = this.bounds
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout()
            this.boundsDidChange()
        }
        
    }
    
    setMargins(left?: number | string, right?: number | string, bottom?: number | string, top?: number | string) {
        
        const previousBounds = this.bounds
        
        this.setStyleProperty("marginLeft", left)
        this.setStyleProperty("marginRight", right)
        this.setStyleProperty("marginBottom", bottom)
        this.setStyleProperty("marginTop", top)
        
        const bounds = this.bounds
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout()
            this.boundsDidChange()
        }
        
    }
    
    setPadding(padding?: number | string) {
        
        const previousBounds = this.bounds
        
        this.setStyleProperty("padding", padding)
        
        const bounds = this.bounds
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout()
            this.boundsDidChange()
        }
        
    }
    
    setPaddings(left?: number | string, right?: number | string, bottom?: number | string, top?: number | string) {
        
        const previousBounds = this.bounds
        
        this.setStyleProperty("paddingLeft", left)
        this.setStyleProperty("paddingRight", right)
        this.setStyleProperty("paddingBottom", bottom)
        this.setStyleProperty("paddingTop", top)
        
        const bounds = this.bounds
        if (bounds.height != previousBounds.height || bounds.width != previousBounds.width) {
            this.setNeedsLayout()
            this.boundsDidChange()
        }
        
    }
    
    
    
    
    
    setBorder(
        radius: number | string = nil,
        width: number | string = 1,
        color: UIColor = UIColor.blackColor,
        style: string = "solid"
    ) {
        
        this.setStyleProperty("borderStyle", style)
        
        this.setStyleProperty("borderRadius", radius)
        
        this.setStyleProperty("borderColor", color.stringValue)
        
        this.setStyleProperty("borderWidth", width)
        
    }
    
    
    
    
    
    setStyleProperty(propertyName: string, value?: number | string) {
        
        
        try {
            
            if (IS_NIL(value)) {
                return
            }
            if (IS_DEFINED(value) && (value as Number).isANumber) {
                value = "" + (value as number).integerValue + "px"
            }
            this.style[propertyName] = value
            
        } catch (exception) {
            
            console.log(exception)
            
        }
        
        
    }
    
    
    
    get userInteractionEnabled() {
        
        const result = (this.style.pointerEvents != "none")
        
        return result
        
    }
    
    set userInteractionEnabled(userInteractionEnabled) {
        
        if (userInteractionEnabled) {
            
            this.style.pointerEvents = ""
            
        }
        else {
            
            this.style.pointerEvents = "none"
            
        }
        
    }
    
    
    
    get backgroundColor() {
        return this._backgroundColor
    }
    
    set backgroundColor(backgroundColor: UIColor) {
        
        this._backgroundColor = backgroundColor
        
        this.style.backgroundColor = backgroundColor.stringValue
        
    }
    
    
    
    get alpha() {
        return 1 * (this.style.opacity as any)
    }
    
    set alpha(alpha) {
        this.style.opacity = "" + alpha
    }
    
    
    
    
    
    static animateViewOrViewsWithDurationDelayAndFunction(
        viewOrViews: UIView | HTMLElement | UIView[] | HTMLElement[],
        duration: number,
        delay: number,
        timingStyle = "cubic-bezier(0.25,0.1,0.25,1)",
        transformFunction: Function,
        transitioncompletionFunction: Function
    ) {
        
        
        function callTransitioncompletionFunction() {
            
            (transitioncompletionFunction || nil)();
            
            (viewOrViews as UIView[]).forEach(function (view, index, array) {
                
                view.animationDidFinish()
                
            })
            
        }
        
        
        if (IS_FIREFOX) {
            
            // Firefox does not fire the transition completion event properly
            new UIObject().performFunctionWithDelay(delay + duration, callTransitioncompletionFunction)
            
            
            
        }
        
        
        if (!(viewOrViews instanceof Array)) {
            viewOrViews = [viewOrViews] as any
        }
        
        const transitionStyles = []
        const transitionDurations = []
        const transitionDelays = []
        const transitionTimings = []
        
        for (var i = 0; i < (viewOrViews as any).length; i++) {
            
            var view = viewOrViews[i]
            
            if (view.viewHTMLElement) {
                
                view = view.viewHTMLElement
                
            }
            
            view.addEventListener("transitionend", transitionDidFinish, true)
            
            transitionStyles.push(view.style.transition)
            transitionDurations.push(view.style.transitionDuration)
            transitionDelays.push(view.style.transitionDelay)
            transitionTimings.push(view.style.transitionTimingFunction)
            
            view.style.transition = "all"
            view.style.transitionDuration = "" + duration + "s"
            view.style.transitionDelay = "" + delay + "s"
            view.style.transitionTimingFunction = timingStyle
            
        }
        
        
        
        transformFunction()
        
        
        const transitionObject = {
            
            "finishImmediately": finishTransitionImmediately,
            "didFinish": transitionDidFinishManually,
            "views": viewOrViews,
            "registrationTime": Date.now()
            
        }
        
        function finishTransitionImmediately() {
            for (var i = 0; i < (viewOrViews as any).length; i++) {
                var view = viewOrViews[i]
                if (view.viewHTMLElement) {
                    view = view.viewHTMLElement
                }
                view.style.transition = "all"
                view.style.transitionDuration = "" + duration + "s"
                view.style.transitionDelay = "" + delay + "s"
                view.style.transition = transitionStyles[i]
                view.style.transitionDuration = transitionDurations[i]
                view.style.transitionDelay = transitionDelays[i]
                view.style.transitionTimingFunction = transitionTimings[i]
            }
        }
        
        function transitionDidFinish(event) {
            var view = event.srcElement
            if (!view) {
                return
            }
            if (view.viewHTMLElement) {
                view = view.viewHTMLElement
            }
            view.style.transition = transitionStyles[i]
            view.style.transitionDuration = transitionDurations[i]
            view.style.transitionDelay = transitionDelays[i]
            view.style.transitionTimingFunction = transitionTimings[i]
            
            callTransitioncompletionFunction()
            
            view.removeEventListener("transitionend", transitionDidFinish, true)
            
        }
        
        function transitionDidFinishManually() {
            for (var i = 0; i < (viewOrViews as any).length; i++) {
                
                var view = viewOrViews[i]
                
                if (view.viewHTMLElement) {
                    view = view.viewHTMLElement
                }
                
                view.style.transition = transitionStyles[i]
                view.style.transitionDuration = transitionDurations[i]
                view.style.transitionDelay = transitionDelays[i]
                view.style.transitionTimingFunction = transitionTimings[i]
                
                view.removeEventListener("transitionend", transitionDidFinish, true)
                
            }
            
            
            
        }
        
        return transitionObject
        
    }
    
    
    
    
    
    animationDidFinish() {
        
        
        
    }
    
    
    
    
    
    static _transformAttribute = (("transform" in document.documentElement.style) ? "transform" : undefined) ||
        (("-webkit-transform" in document.documentElement.style) ? "-webkit-transform" : "undefined") ||
        (("-moz-transform" in document.documentElement.style) ? "-moz-transform" : "undefined") ||
        (("-ms-transform" in document.documentElement.style) ? "-ms-transform" : "undefined") ||
        (("-o-transform" in document.documentElement.style) ? "-o-transform" : "undefined")
    
    static _setAbsoluteSizeAndPosition(element, left, top, width, height, zIndex = 0) {
        
        // if (!UIView._transformAttribute) {
        
        //     UIView._transformAttribute = (('transform' in document.documentElement.style) ? 'transform' : undefined);
        //     UIView._transformAttribute = UIView._transformAttribute || (('-webkit-transform' in document.documentElement.style) ? '-webkit-transform' : 'undefined');
        //     UIView._transformAttribute = UIView._transformAttribute || (('-moz-transform' in document.documentElement.style) ? '-moz-transform' : 'undefined');
        //     UIView._transformAttribute = UIView._transformAttribute || (('-ms-transform' in document.documentElement.style) ? '-ms-transform' : 'undefined');
        //     UIView._transformAttribute = UIView._transformAttribute || (('-o-transform' in document.documentElement.style) ? '-o-transform' : 'undefined');
        
        // }
        
        if (!IS(element) || !element.obeyAutolayout && !element.getAttribute("obeyAutolayout")) {
            return
        }
        
        if (element.id == "mainView") {
            
            
            var asd = 1
            
        }
        
        if (IS(height)) {
            height = height.integerValue + "px"
        }
        
        if (IS(width)) {
            width = width.integerValue + "px"
        }
        
        var str = element.style.cssText
        
        const frameTransform = UIView._transformAttribute + ": translate3d(" + (1 * left).integerValue + "px, " +
            (1 * top).integerValue + "px, " + zIndex.integerValue + "px)"
        
        if (element.UIView) {
            str = str + frameTransform + ";"
        }
        else {
            element.UIView._frameTransform = frameTransform
        }
        
        if (height == nil) {
            str = str + " height: unset;"
        }
        else {
            str = str + " height:" + height + ";"
        }
        
        if (width == nil) {
            str = str + " width: unset;"
        }
        else {
            str = str + " width:" + width + ";"
        }
        
        if (element.id == "mainView") {
            
            
            var asd = 1
            
        }
        
        element.style.cssText = element.style.cssText + str
        
    }
    
    
    
    static performAutoLayout(parentElement, visualFormatArray, constraintsArray) {
        
        
        const view = new AutoLayout.View()
        
        if (IS(visualFormatArray) && IS(visualFormatArray.length)) {
            view.addConstraints(AutoLayout.VisualFormat.parse(visualFormatArray, { extended: true }))
        }
        
        if (IS(constraintsArray) && IS(constraintsArray.length)) {
            view.addConstraints(constraintsArray)
        }
        
        const elements = {}
        for (var key in view.subViews) {
            
            if (!view.subViews.hasOwnProperty(key)) {
                continue
            }
            
            var element = nil
            
            try {
                
                element = parentElement.querySelector("#" + key)
                
            } catch (error) {
                
                //console.log("Error occurred " + error);
                
            }
            
            if (element && !element.obeyAutolayout && !element.getAttribute("obeyAutolayout")) {
            
            
            
            }
            else if (element) {
                
                element.className += element.className ? " abs" : "abs"
                elements[key] = element
                
            }
            
        }
        
        var parentUIView = nil
        
        if (parentElement.UIView) {
            
            parentUIView = parentElement.UIView
            
        }
        
        const updateLayout = function () {
            view.setSize(
                parentElement ? parentElement.clientWidth : window.innerWidth,
                parentElement ? parentElement.clientHeight : window.innerHeight
            )
            for (key in view.subViews) {
                
                if (!view.subViews.hasOwnProperty(key)) {
                    continue
                }
                
                const subView = view.subViews[key]
                
                if (elements[key]) {
                    UIView._setAbsoluteSizeAndPosition(
                        elements[key],
                        subView.left,
                        subView.top,
                        subView.width,
                        subView.height
                    )
                }
            }
            
            parentUIView.didLayoutSubviews()
            
        }
        
        updateLayout()
        return updateLayout
        
    }
    
    
    static runFunctionBeforeNextFrame(step: () => void) {
        
        if (IS_SAFARI) {
            
            // This creates a microtask
            Promise.resolve().then(step)
            
        }
        else {
            
            window.requestAnimationFrame(step)
            
        }
        
    }
    
    
    static scheduleLayoutViewsIfNeeded() {
        
        UIView.runFunctionBeforeNextFrame(UIView.layoutViewsIfNeeded)
        
    }
    
    
    static layoutViewsIfNeeded() {
        for (var i = 0; i < UIView._viewsToLayout.length; i++) {
            const view = UIView._viewsToLayout[i]
            view.layoutIfNeeded()
        }
        UIView._viewsToLayout = []
    }
    
    
    setNeedsLayout() {
        
        if (this._shouldLayout) {
            return
        }
        
        this._shouldLayout = YES
        
        // Register view for layout before next frame
        UIView._viewsToLayout.push(this)
        
        if (UIView._viewsToLayout.length == 1) {
            UIView.scheduleLayoutViewsIfNeeded()
        }
        
    }
    
    
    get needsLayout() {
        
        return this._shouldLayout
        
    }
    
    
    layoutIfNeeded() {
        
        if (!this._shouldLayout) {
            return
        }
        
        this._shouldLayout = NO
        
        try {
            
            this.layoutSubviews()
            
        } catch (exception) {
            
            console.log(exception)
            
        }
        
    }
    
    
    layoutSubviews() {
        
        
        this._shouldLayout = NO
        
        // Autolayout
        //window.removeEventListener('resize', this._updateLayoutFunction);
        if (this.constraints.length) {
            
            this._updateLayoutFunction = UIView.performAutoLayout(this.viewHTMLElement, null, this.constraints)
            
        }
        
        //this._updateLayoutFunction = this.layoutSubviews.bind(this);
        
        //window.addEventListener('resize', this._updateLayoutFunction);
        
        this._viewControllerLayoutFunction()
        
        this.applyClassesAndStyles()
        
        for (var i = 0; i < this.subviews.length; i++) {
            
            const subview = this.subviews[i]
            
            subview.calculateAndSetViewFrame()
            
            //subview.layoutSubviews();
            
        }
        
        this.didLayoutSubviews()
        
    }
    
    
    
    applyClassesAndStyles() {
        
        
        
        //var classesString = "";
        
        for (var i = 0; i < this.styleClasses.length; i++) {
            
            const styleClass = this.styleClasses[i]
            
            if (styleClass) {
                
                this.viewHTMLElement.classList.add(styleClass)
                
            }
            
            
            
            //classesString = classesString + " " + styleClass;
            
        }
        
        
        //this.viewHTMLElement.className = classesString;
        
        
        
    }
    
    didLayoutSubviews() {
        
        this._didLayoutSubviewsDelegateFunction()
        
    }
    
    get constraints() {
        return this._constraints
    }
    
    set constraints(constraints) {
        this._constraints = constraints
    }
    
    addConstraint(constraint) {
        
        this.constraints.push(constraint)
        
    }
    
    
    
    addConstraintsWithVisualFormat(visualFormatArray) {
        
        this.constraints = this.constraints.concat(AutoLayout.VisualFormat.parse(
            visualFormatArray,
            { extended: true }
        ))
        
    }
    
    static constraintWithView(view, attribute, relation, toView, toAttribute, multiplier, constant, priority) {
        
        var UIViewObject = nil
        var viewID = null
        if (view) {
            if (view.isKindOfClass && view.isKindOfClass(UIView)) {
                UIViewObject = view
                view = view.viewHTMLElement
            }
            viewID = view.id
        }
        
        var toUIViewObject = nil
        var toViewID = null
        if (toView) {
            if (toView.isKindOfClass && view.isKindOfClass(UIView)) {
                toUIViewObject = toView
                toView = toView.viewHTMLElement
            }
            toViewID = toView.id
        }
        
        const constraint = {
            
            view1: viewID,
            attr1: attribute,
            relation: relation,
            view2: toViewID,
            attr2: toAttribute,
            multiplier: multiplier,
            constant: constant,
            priority: priority
            
        }
        
        return constraint
        
    }
    
    static constraintAttribute = {
        
        "left": AutoLayout.Attribute.LEFT,
        "right": AutoLayout.Attribute.RIGHT,
        "bottom": AutoLayout.Attribute.BOTTOM,
        "top": AutoLayout.Attribute.TOP,
        "centerX": AutoLayout.Attribute.CENTERX,
        "centerY": AutoLayout.Attribute.CENTERY,
        "height": AutoLayout.Attribute.HEIGHT,
        "width": AutoLayout.Attribute.WIDTH,
        "zIndex": AutoLayout.Attribute.ZINDEX,
        // Not sure what these are for
        "constant": AutoLayout.Attribute.NOTANATTRIBUTE,
        "variable": AutoLayout.Attribute.VARIABLE
        
    }
    
    static constraintRelation = {
        
        "equal": AutoLayout.Relation.EQU,
        "lessThanOrEqual": AutoLayout.Relation.LEQ,
        "greaterThanOrEqual": AutoLayout.Relation.GEQ
        
    }
    
    
    
    subviewWithID(viewID) {
        
        
        var resultHTMLElement = nil
        
        try {
            
            resultHTMLElement = this.viewHTMLElement.querySelector("#" + viewID)
            
        } catch (error) {
        
        }
        
        if (resultHTMLElement && resultHTMLElement.UIView) {
            return resultHTMLElement.UIView
        }
        return nil
    }
    
    
    
    rectangleContainingSubviews() {
        
        const center = this.bounds.center
        
        var result = new UIRectangle(center.x, center.y, 0, 0)
        
        for (var i = 0; i < this.subviews.length; i++) {
            
            const subview = this.subviews[i]
            
            var frame = subview.frame
            
            const rectangleContainingSubviews = subview.rectangleContainingSubviews()
            
            frame = frame.concatenateWithRectangle(rectangleContainingSubviews)
            
            result = result.concatenateWithRectangle(frame)
            
        }
        
        return result
        
    }
    
    
    
    
    hasSubview(view: UIView) {
        
        // This is for performance reasons
        if (!IS(view)) {
            return NO
        }
        
        for (var i = 0; i < this.subviews.length; i++) {
            const subview = this.subviews[i]
            if (subview == view) {
                return YES
            }
        }
        return NO
    }
    
    get viewBelowThisView() {
        const result: UIView = (this.viewHTMLElement.previousElementSibling as any || {}).UIView
        return result
    }
    
    get viewAboveThisView() {
        const result: UIView = (this.viewHTMLElement.nextElementSibling as any || {}).UIView
        return result
    }
    
    addSubview(view: UIView, aboveView?: UIView) {
        
        if (!this.hasSubview(view) && IS(view)) {
            
            view.willMoveToSuperview(this)
            
            if (IS(aboveView)) {
                this.viewHTMLElement.insertBefore(view.viewHTMLElement, aboveView.viewHTMLElement.nextSibling)
                this.subviews.insertElementAtIndex(this.subviews.indexOf(aboveView), view)
            }
            else {
                this.viewHTMLElement.appendChild(view.viewHTMLElement)
                this.subviews.push(view)
            }
            view.didMoveToSuperview(this)
            
            if (this.superview && this.isMemberOfViewTree) {
                
                view.broadcastEventInSubtree({
                    
                    name: UIView.broadcastEventName.AddedToViewTree,
                    parameters: nil
                    
                })
                
            }
            
            this.setNeedsLayout()
            
        }
        
    }
    
    addSubviews(views: UIView[]) {
        views.forEach(function (this: UIView, view: UIView, index, array) {
            this.addSubview(view)
        }, this)
    }
    
    
    moveToBottomOfSuperview() {
        
        if (IS(this.superview)) {
            
            const bottomView = this.superview.subviews.firstElement
            
            if (bottomView == this) {
                
                return
                
            }
            
            this.superview.subviews.removeElement(this)
            
            this.superview.subviews.insertElementAtIndex(0, this)
            
            this.superview.viewHTMLElement.insertBefore(this.viewHTMLElement, bottomView.viewHTMLElement)
            
            
        }
        
        
        
    }
    
    moveToTopOfSuperview() {
        
        if (IS(this.superview)) {
            
            const topView = this.superview.subviews.lastElement
            
            if (topView == this) {
                
                return
                
            }
            
            this.superview.subviews.removeElement(this)
            
            this.superview.subviews.push(this)
            
            this.superview.viewHTMLElement.appendChild(this.viewHTMLElement)
            
            
        }
        
        
        
    }
    
    
    removeFromSuperview() {
        if (IS(this.superview)) {
            
            this.forEachViewInSubtree(function (view) {
                
                view.blur()
                
            })
            
            const index = this.superview.subviews.indexOf(this)
            if (index > -1) {
                this.superview.subviews.splice(index, 1)
                this.superview.viewHTMLElement.removeChild(this.viewHTMLElement)
                this.superview = nil
                
                this.broadcastEventInSubtree({
                    
                    name: UIView.broadcastEventName.RemovedFromViewTree,
                    parameters: nil
                    
                })
                
            }
        }
    }
    
    
    willAppear() {
        
        
        
    }
    
    
    
    willMoveToSuperview(superview: UIView) {
        
        this._setInnerHTMLFromKeyIfPossible()
        
        this._setInnerHTMLFromLocalizedTextObjectIfPossible()
        
    }
    
    didMoveToSuperview(superview: UIView) {
        
        this.superview = superview
        
    }
    
    wasAddedToViewTree() {
        
    }
    
    wasRemovedFromViewTree() {
        
    }
    
    get isMemberOfViewTree() {
        var element = this.viewHTMLElement
        for (var i = 0; element; i = i) {
            if (element.parentElement && element.parentElement == document.body) {
                return YES
            }
            element = element.parentElement
        }
        return NO
    }
    
    
    get allSuperviews() {
        const result = []
        var view: UIView = this
        for (var i = 0; IS(view); i = i) {
            result.push(view)
            view = view.superview
        }
        return result
    }
    
    
    setNeedsLayoutOnAllSuperviews() {
        
        this.allSuperviews.reverse().forEach(function (view: UIView, index, array) {
            
            view.setNeedsLayout()
            
        })
        
    }
    
    
    setNeedsLayoutUpToRootView() {
        
        this.setNeedsLayoutOnAllSuperviews()
        
        this.setNeedsLayout()
        
    }
    
    
    focus() {
        
        this.viewHTMLElement.focus()
        
    }
    
    
    blur() {
        
        this.viewHTMLElement.blur()
        
    }
    
    
    
    
    
    _loadUIEvents() {
        
        //this.viewHTMLElement = nil;
        
        const isTouchEventClassDefined: boolean = NO || (window as any).TouchEvent
        
        const pauseEvent = (event: Event, forced = NO) => {
            
            if (this.pausesPointerEvents || forced) {
                
                if (event.stopPropagation) {
                    event.stopPropagation()
                }
                if (event.preventDefault) {
                    event.preventDefault()
                }
                event.cancelBubble = true
                event.returnValue = false
                return false
                
            }
            
            if (event.stopPropagation && this.stopsPointerEventPropagation) {
                event.stopPropagation()
            }
            
        }
        
        const onMouseDown = (event) => {
            
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                ((this.ignoresMouse || (IS(this._touchEventTime) && (Date.now() - this._touchEventTime) > 500)) &&
                    event instanceof MouseEvent)) {
                return
            }
            
            this.sendControlEventForKey(UIView.controlEvent.PointerDown, event)
            
            this._isPointerInside = YES
            this._isPointerValid = YES
            this._initialPointerPosition = new UIPoint(event.clientX, event.clientY)
            if (isTouchEventClassDefined && event instanceof TouchEvent) {
                
                this._touchEventTime = Date.now()
                
                this._initialPointerPosition = new UIPoint(event.touches[0].clientX, event.touches[0].clientY)
                
                if (event.touches.length > 1) {
                    
                    onTouchCancel(event)
                    
                    return
                    
                }
                
                
            }
            else {
                
                this._touchEventTime = nil
                
                pauseEvent(event)
                
            }
            
            
            this._hasPointerDragged = NO
            
        }
        
        const onTouchStart = onMouseDown as any
        
        const onmouseup = (event) => {
            
            if (!this._isPointerValid) {
                return
            }
            
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return
            }
            
            if (this._isPointerInside) {
                
                onPointerUpInside(event)
                
                if (!this._hasPointerDragged) {
                    
                    this.sendControlEventForKey(UIView.controlEvent.PointerTap, event)
                    
                }
                
                
            }
            
            // This has to be sent after the more specific event so that UIButton can ignore it when not highlighted
            this.sendControlEventForKey(UIView.controlEvent.PointerUp, event)
            
            pauseEvent(event)
            
        }
        
        const onTouchEnd = onmouseup
        
        // function onMouseEnter(event) {
        
        //     this.sendControlEventForKey(UIView.controlEvent.PointerEnter, event);
        
        //     this._isPointerInside = YES;
        
        //     pauseEvent(event);
        
        // }
        
        const onmouseout = (event) => {
            
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return
            }
            
            this.sendControlEventForKey(UIView.controlEvent.PointerLeave, event)
            
            this._isPointerInside = NO
            
            pauseEvent(event)
            
        }
        
        const onTouchLeave = onmouseout
        
        var onTouchCancel = function (event) {
            
            if (!this._isPointerValid) {
                return
            }
            
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return
            }
            
            this._isPointerValid = NO
            
            this.sendControlEventForKey(UIView.controlEvent.PointerCancel, event)
            
        }.bind(this)
        
        const onmouseover = (event) => {
            
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return
            }
            
            this.sendControlEventForKey(UIView.controlEvent.PointerHover, event)
            
            this._isPointerInside = YES
            
            this._isPointerValid = YES
            
            pauseEvent(event)
            
        }
        
        const onMouseMove = (event) => {
            
            if (!this._isPointerValid) {
                return
            }
            
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return
            }
            
            if (IS_NOT(this._initialPointerPosition)) {
                
                this._initialPointerPosition = new UIPoint(event.clientX, event.clientY)
                
            }
            
            if (new UIPoint(event.clientX, event.clientY).to(this._initialPointerPosition).length >
                this._pointerDragThreshold) {
                
                this._hasPointerDragged = YES
                
            }
            
            
            this.sendControlEventForKey(UIView.controlEvent.PointerMove, event)
            
            pauseEvent(event)
            
        }
        
        const onTouchMove = function (event: TouchEvent) {
            
            if (!this._isPointerValid) {
                return
            }
            
            if ((this.ignoresTouches && isTouchEventClassDefined && event instanceof TouchEvent) ||
                (this.ignoresMouse && event instanceof MouseEvent)) {
                return
            }
            
            if (event.touches.length > 1) {
                
                onTouchZoom(event)
                
                return
                
            }
            
            const touch = event.touches[0]
            
            if (new UIPoint(touch.clientX, touch.clientY).to(this._initialPointerPosition).length >
                this._pointerDragThreshold) {
                
                this._hasPointerDragged = YES
                
            }
            
            
            if (this._isPointerInside && this.viewHTMLElement !=
                document.elementFromPoint(touch.clientX, touch.clientY)) {
                
                this._isPointerInside = NO
                
                this.sendControlEventForKey(UIView.controlEvent.PointerLeave, event)
                
            }
            
            
            this.sendControlEventForKey(UIView.controlEvent.PointerMove, event)
            
            //pauseEvent(event);
            
            
        }
        
        var onTouchZoom = function onTouchZoom(event: TouchEvent) {
            
            this.sendControlEventForKey(UIView.controlEvent.MultipleTouches, event)
            
        }.bind(this)
        
        
        var onPointerUpInside = (event) => {
            pauseEvent(event)
            
            this.sendControlEventForKey(UIView.controlEvent.PointerUpInside, event)
        }
        
        function eventKeyIsEnter(event) {
            if (event.keyCode !== 13) {
                return NO
            }
            return YES
        }
        
        function eventKeyIsTab(event) {
            if (event.keyCode !== 9) {
                return NO
            }
            return YES
        }
        
        function eventKeyIsEsc(event) {
            var result = false
            if ("key" in event) {
                result = (event.key == "Escape" || event.key == "Esc")
            }
            else {
                result = (event.keyCode == 27)
            }
            return result
        }
        
        function eventKeyIsLeft(event) {
            if (event.keyCode != "37") {
                return NO
            }
            return YES
        }
        
        function eventKeyIsRight(event) {
            if (event.keyCode != "39") {
                return NO
            }
            return YES
        }
        
        function eventKeyIsDown(event) {
            if (event.keyCode != "40") {
                return NO
            }
            return YES
        }
        
        function eventKeyIsUp(event) {
            if (event.keyCode != "38") {
                return NO
            }
            return YES
        }
        
        const onKeyDown = function (event) {
            
            if (eventKeyIsEnter(event)) {
                
                this.sendControlEventForKey(UIView.controlEvent.EnterDown, event)
                
            }
            
            if (eventKeyIsEsc(event)) {
                
                this.sendControlEventForKey(UIView.controlEvent.EscDown, event)
                
            }
            
            if (eventKeyIsTab(event) && this._controlEventTargets.TabDown && this._controlEventTargets.TabDown.length) {
                
                this.sendControlEventForKey(UIView.controlEvent.TabDown, event)
                
                pauseEvent(event, YES)
                
            }
            
            if (eventKeyIsLeft(event)) {
                
                this.sendControlEventForKey(UIView.controlEvent.LeftArrowDown, event)
                
            }
            
            if (eventKeyIsRight(event)) {
                
                this.sendControlEventForKey(UIView.controlEvent.RightArrowDown, event)
                
            }
            
            if (eventKeyIsDown(event)) {
                
                this.sendControlEventForKey(UIView.controlEvent.DownArrowDown, event)
                
            }
            
            if (eventKeyIsUp(event)) {
                
                this.sendControlEventForKey(UIView.controlEvent.UpArrowDown, event)
                
            }
            
        }.bind(this)
        
        const onKeyUp = function (event) {
            
            if (eventKeyIsEnter(event)) {
                
                this.sendControlEventForKey(UIView.controlEvent.EnterUp, event)
                
            }
            
        }.bind(this)
        
        
        const onfocus = function (event: Event) {
            
            this.sendControlEventForKey(UIView.controlEvent.Focus, event)
            
        }.bind(this)
        
        const onblur = function (event: Event) {
            
            this.sendControlEventForKey(UIView.controlEvent.Blur, event)
            
        }.bind(this)
        
        
        // Mouse and touch start events
        this._viewHTMLElement.onmousedown = onMouseDown.bind(this)
        this._viewHTMLElement.ontouchstart = onTouchStart.bind(this)
        // this.viewHTMLElement.addEventListener("mousedown", onMouseDown.bind(this), false)
        // this.viewHTMLElement.addEventListener('touchstart', onTouchStart.bind(this), false)
        // //this.viewHTMLElement.addEventListener("mouseenter", onMouseEnter.bind(this), false);
        
        // Mouse and touch move events
        this._viewHTMLElement.onmousemove = onMouseMove.bind(this)
        this._viewHTMLElement.ontouchmove = onTouchMove.bind(this)
        // this.viewHTMLElement.addEventListener("mousemove", onMouseMove.bind(this), false)
        // this.viewHTMLElement.addEventListener('touchmove', onTouchMove.bind(this), false)
        
        //this.viewHTMLElement.addEventListener("mousewheel", onmousewheel.bind(this), false)
        
        this._viewHTMLElement.onmouseover = onmouseover.bind(this)
        // this.viewHTMLElement.addEventListener("mouseover", onmouseover.bind(this), false)
        
        // Mouse and touch end events
        this._viewHTMLElement.onmouseup = onmouseup.bind(this)
        this._viewHTMLElement.ontouchend = onTouchEnd.bind(this)
        this._viewHTMLElement.ontouchcancel = onTouchCancel.bind(this)
        // this.viewHTMLElement.addEventListener("mouseup", onmouseup.bind(this), false)
        // this.viewHTMLElement.addEventListener('touchend', onTouchEnd.bind(this), false)
        // this.viewHTMLElement.addEventListener('touchcancel', onTouchCancel.bind(this), false)
        
        this._viewHTMLElement.onmouseout = onmouseout.bind(this)
        // this.viewHTMLElement.addEventListener("mouseout", onmouseout.bind(this), false)
        this._viewHTMLElement.addEventListener("touchleave", onTouchLeave.bind(this), false)
        
        // this.viewHTMLElement.onkeydown = onkeydown
        // this.viewHTMLElement.onkeyup = onkeyup
        this._viewHTMLElement.addEventListener("keydown", onKeyDown, false)
        this._viewHTMLElement.addEventListener("keyup", onKeyUp, false)
        
        // Focus events
        this._viewHTMLElement.onfocus = onfocus
        this._viewHTMLElement.onblur = onblur
        // this.viewHTMLElement.addEventListener("focus", onfocus, true)
        // this.viewHTMLElement.addEventListener("blur", onblur, true)
        
        
    }
    
    
    
    
    
    public static controlEvent = {
        
        "PointerDown": "PointerDown",
        "PointerMove": "PointerMove",
        "PointerLeave": "PointerLeave",
        "PointerEnter": "PointerEnter",
        "PointerUpInside": "PointerUpInside",
        "PointerTap": "PointerTap",
        "PointerUp": "PointerUp",
        "MultipleTouches": "PointerZoom",
        "PointerCancel": "PointerCancel",
        "PointerHover": "PointerHover",
        "EnterDown": "EnterDown",
        "EnterUp": "EnterUp",
        "EscDown": "EscDown",
        "TabDown": "TabDown",
        "LeftArrowDown": "LeftArrowDown",
        "RightArrowDown": "RightArrowDown",
        "DownArrowDown": "DownArrowDown",
        "UpArrowDown": "UpArrowDown",
        "Focus": "Focus",
        "Blur": "Blur"
        
    }
    
    controlEvent = UIView.controlEvent
    
    
    
    get addControlEventTarget() {
        
        const eventKeys = []
        
        
        const result: UIViewAddControlEventTargetObject<typeof UIView.controlEvent> = new Proxy(
            (this.constructor as any).controlEvent,
            {
                
                get: (target, key, receiver) => {
                    
                    eventKeys.push(key)
                    
                    return result
                    
                },
                set: (target, key, value, receiver) => {
                    
                    eventKeys.push(key)
                    
                    this.addTargetForControlEvents(eventKeys, value)
                    
                    return true
                    
                }
                
            }
        )
        
        return result
        
    }
    
    
    
    
    
    addTargetForControlEvents(eventKeys: string[], targetFunction: (sender: UIView, event: Event) => void) {
        
        eventKeys.forEach(function (this: UIView, key: string, index: number, array: string[]) {
            
            this.addTargetForControlEvent(key, targetFunction)
            
        }, this)
        
    }
    
    
    
    addTargetForControlEvent(eventKey: string, targetFunction: (sender: UIView, event: Event) => void) {
        
        var targets = this._controlEventTargets[eventKey]
        
        if (!targets) {
            // @ts-ignore
            targets = []
            this._controlEventTargets[eventKey] = targets
        }
        
        if (targets.indexOf(targetFunction) == -1) {
            targets.push(targetFunction)
        }
        
    }
    
    removeTargetForControlEvent(eventKey: string, targetFunction: (sender: UIView, event: Event) => void) {
        const targets = this._controlEventTargets[eventKey]
        if (!targets) {
            return
        }
        const index = targets.indexOf(targetFunction)
        if (index != -1) {
            targets.splice(index, 1)
        }
    }
    
    removeTargetForControlEvents(eventKeys: string[], targetFunction: (sender: UIView, event: Event) => void) {
        
        eventKeys.forEach(function (key, index, array) {
            
            this.removeTargetForControlEvent(key, targetFunction)
            
        }, this)
        
    }
    
    sendControlEventForKey(eventKey: string, nativeEvent: Event) {
        var targets = this._controlEventTargets[eventKey]
        if (!targets) {
            return
        }
        targets = targets.copy()
        for (var i = 0; i < targets.length; i++) {
            const target = targets[i]
            target(this, nativeEvent)
        }
    }
    
    
    
    
    
    static broadcastEventName = {
        
        "LanguageChanged": "LanguageChanged",
        "RemovedFromViewTree": "RemovedFromViewTree",
        "AddedToViewTree": "AddedToViewTree",
        "PageDidScroll": "PageDidScroll"
        
    }
    
    
    broadcastEventInSubtree(event: UIViewBroadcastEvent) {
        
        this.forEachViewInSubtree(function (view) {
            
            view.didReceiveBroadcastEvent(event)
            
            if (view._didReceiveBroadcastEventDelegateFunction) {
                
                view._didReceiveBroadcastEventDelegateFunction(event)
                
            }
            
        })
        
        
    }
    
    didReceiveBroadcastEvent(event: UIViewBroadcastEvent) {
        
        if (event.name == UIView.broadcastEventName.PageDidScroll) {
            
            this._isPointerValid = NO
            
        }
        
        if (event.name == UIView.broadcastEventName.AddedToViewTree) {
            
            this.wasAddedToViewTree()
            
        }
        
        if (event.name == UIView.broadcastEventName.RemovedFromViewTree) {
            
            this.wasRemovedFromViewTree()
            
        }
        
        if (event.name == UIView.broadcastEventName.LanguageChanged || event.name ==
            UIView.broadcastEventName.AddedToViewTree) {
            
            this._setInnerHTMLFromKeyIfPossible()
            
            this._setInnerHTMLFromLocalizedTextObjectIfPossible()
            
        }
        
        
        
    }
    
    
    
    
    
    forEachViewInSubtree(functionToCall: (view: UIView) => void) {
        
        functionToCall(this)
        
        this.subviews.forEach(function (subview, index, array) {
            
            subview.forEachViewInSubtree(functionToCall)
            
        })
        
    }
    
    
    
    
    
    rectangleInView(rectangle: UIRectangle, view: UIView) {
        if (!view.isMemberOfViewTree || !this.isMemberOfViewTree) {
            return nil
        }
        
        const viewClientRectangle = view.viewHTMLElement.getBoundingClientRect()
        const viewLocation = new UIPoint(viewClientRectangle.left, viewClientRectangle.top)
        
        const selfClientRectangle = this.viewHTMLElement.getBoundingClientRect()
        const selfLocation = new UIPoint(selfClientRectangle.left, selfClientRectangle.top)
        
        const offsetPoint = selfLocation.subtract(viewLocation)
        
        return rectangle.copy().offsetByPoint(offsetPoint)
    }
    
    rectangleFromView(rectangle: UIRectangle, view: UIView) {
        return view.rectangleInView(rectangle, this)
    }
    
    
    
    
    
    intrinsicContentSizeWithConstraints(constrainingHeight: number = 0, constrainingWidth: number = 0) {
        
        // This works but is slow
        
        const result = new UIRectangle(0, 0, 0, 0)
        if (this.rootView.forceIntrinsicSizeZero) {
            return result
        }
        
        var temporarilyInViewTree = NO
        var nodeAboveThisView: Node
        if (!this.isMemberOfViewTree) {
            document.body.appendChild(this.viewHTMLElement)
            temporarilyInViewTree = YES
            nodeAboveThisView = this.viewHTMLElement.nextSibling
        }
        
        const height = this.style.height
        const width = this.style.width
        
        this.style.height = "" + constrainingHeight
        this.style.width = "" + constrainingWidth
        
        
        const left = this.style.left
        const right = this.style.right
        const bottom = this.style.bottom
        const top = this.style.top
        
        this.style.left = ""
        this.style.right = ""
        this.style.bottom = ""
        this.style.top = ""
        
        
        const resultHeight = this.viewHTMLElement.scrollHeight
        
        
        const whiteSpace = this.style.whiteSpace
        this.style.whiteSpace = "nowrap"
        
        const resultWidth = this.viewHTMLElement.scrollWidth
        
        this.style.whiteSpace = whiteSpace
        
        
        
        this.style.height = height
        this.style.width = width
        
        this.style.left = left
        this.style.right = right
        this.style.bottom = bottom
        this.style.top = top
        
        if (temporarilyInViewTree) {
            document.body.removeChild(this.viewHTMLElement)
            if (this.superview) {
                if (nodeAboveThisView) {
                    this.superview.viewHTMLElement.insertBefore(this.viewHTMLElement, nodeAboveThisView)
                }
                else {
                    this.superview.viewHTMLElement.appendChild(this.viewHTMLElement)
                }
            }
        }
        
        result.height = resultHeight
        result.width = resultWidth
        
        
        return result
        
    }
    
    
    
    
    
    intrinsicContentWidth(constrainingHeight: number = 0): number {
        
        const result = this.intrinsicContentSizeWithConstraints(constrainingHeight).width
        
        return result
        
    }
    
    intrinsicContentHeight(constrainingWidth: number = 0): number {
        
        const result = this.intrinsicContentSizeWithConstraints(undefined, constrainingWidth).height
        
        return result
        
        
    }
    
    intrinsicContentSize(): UIRectangle {
        
        return nil
        
    }
    
    
    
    
    
}





































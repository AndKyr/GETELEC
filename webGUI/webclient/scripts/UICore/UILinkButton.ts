/// <reference path="./UILink.ts" />



class UILinkButton extends UILink {
    
    
    button: UIButton
    
    //link = new UILink(this.elementID + "Link");
    
    
    
    
    
    constructor(elementID: string, elementType?: string, titleType?: string) {
        
        super(elementID, { "elementType": elementType, "titleType": titleType })
        
        this._class = UILinkButton
        this.superclass = UILink
        
        
        this.button.addTargetForControlEvents([
            UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
        ], function (this: UILinkButton, sender: UIButton, event: Event) {
    
            const asd = 1
    
            window.location = this.target as any
            
            
        }.bind(this))
        
        // this.link.hidden = YES;
        
        // this.addSubview(this.link);
        
        
    }
    
    
    
    
    
    initView(elementID, viewHTMLElement, initViewData: { titleType: string, elementType: string }) {
        
        super.initView(elementID, viewHTMLElement, initViewData)
        
        this.class.superclass = UILink
        
        // Instance variables
        
        this.button = new UIButton(this.elementID + "Button", initViewData.elementType, initViewData.titleType)
        
        this.addSubview(this.button)
        
        this.style.position = "absolute"
        
        
        
        
        
    }
    
    
    
    
    
    get titleLabel() {
        
        return this.button.titleLabel
        
    }
    
    get imageView() {
        
        return this.button.imageView
        
    }
    
    
    set colors(colors: UIButtonColorSpecifier) {
        
        this.button.colors = colors
        
    }
    
    
    get colors() {
        
        return this.button.colors
        
    }
    
    
    
    
    
    get viewHTMLElement() {
        
        return super.viewHTMLElement as HTMLLinkElement
        
    }
    
    
    set target(target: string) {
        
        this.viewHTMLElement.setAttribute("href", target)
        
    }
    
    get target() {
    
        const result = this.viewHTMLElement.getAttribute("href")
    
        return result
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
    
        const bounds = this.bounds
    
        this.button.frame = bounds
        
        this.button.layoutSubviews()
        
        
    }
    
    
    
    
    
}






























































































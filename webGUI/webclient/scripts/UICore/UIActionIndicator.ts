///<reference path="UIView.ts"/>





class UIActionIndicator extends UIView {
    
    
    indicatorView: UIView
    
    _size: number = 50
    
    constructor(elementID?: string) {
        
        super(elementID)
        
        this._class = UIActionIndicator
        this.superclass = UIView
        
        this.indicatorView = new UIView(this.elementID + "IndicatorView")
        
        this.addSubview(this.indicatorView)
        
        this.indicatorView.viewHTMLElement.classList.add("LukeHaasLoader")
        
        this.hidden = YES
        
        
    }
    
    
    
    
    
    set size(size: number) {
        
        this._size = size
        
        this.setNeedsLayoutUpToRootView()
        
    }
    
    get size() {
        
        return this._size
        
    }
    
    
    
    set hidden(hidden: boolean) {
        
        super.hidden = hidden
        
        if (hidden) {
            
            this.indicatorView.removeFromSuperview()
            
        }
        else {
            
            this.addSubview(this.indicatorView)
            
        }
        
    }
    
    start() {
        
        this.hidden = NO
        
    }
    
    stop() {
        
        this.hidden = YES
        
    }
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
        
        const bounds = this.bounds
        
        //this.indicatorView.centerInContainer();
        
        this.indicatorView.style.height = "" + this._size.integerValue + "px"
        this.indicatorView.style.width = "" + this._size.integerValue + "px"
        
        const minSize = Math.min(this.bounds.height, this.bounds.width)
        
        this.indicatorView.style.maxHeight = "" + minSize.integerValue + "px"
        this.indicatorView.style.maxWidth = "" + minSize.integerValue + "px"
        
        const size = Math.min(this._size, minSize)
        
        this.indicatorView.style.left = "" + ((bounds.width - size) * 0.5 - 11).integerValue + "px"
        this.indicatorView.style.top = "" + ((bounds.height - size) * 0.5 - 11).integerValue + "px"
        
        
        
        
    }
    
    
    
    
    
}





























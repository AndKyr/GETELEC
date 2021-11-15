/// <reference path="./UICore/UIViewController.ts" />
/// <reference path="./UICore/UIDialogView.ts" />
/// <reference path="./UICore/UILinkButton.ts" />





class RootViewController extends UIViewController {
    
    private useGETELECViewController: GETELECViewController
    private useGETELECTableViewController: GETELECTableViewController
    private mainViewController: UIViewController
    
    private _contentViewController: UIViewController
    private topBarView: CellView
    
    private backToMainButton: CBFlatButton
    
    bottomBarView: any
    
    constructor(view) {
        
        // Calling super
        super(view)
        
        // Here are some suggested conventions that are used in UICore
        
        // Instance variables, it is good to initialize to nil or empty function, not leave as undefined to avoid
        // if blocks
        // this._firstView = nil;
        // this._secondView = nil;
        // this._testView = nil;
        // this._button = nil;
        
        // The nil object avoids unneccessary crashes by allowing you to call any function or access any variable on it, returning nil
        
        // Define properties with get and set functions so they can be accessed and set like variables
        
        // Name variables that should be private, like property variables, with a _ sign, this also holds for private functions
        // Avoid accessing variables and functions named with _ from outside as this creates strong coupling and hinders stability
        
        // Code for further setup if necessary
        
    }
    
    
    loadIntrospectionVariables() {
        
        super.loadIntrospectionVariables()
        this.superclass = UIViewController
        
    }
    
    
    loadSubviews() {
        
        this.view.backgroundColor = UIColor.colorWithRGBA(225, 225, 225)
        
        document.body.style.backgroundColor = this.view.backgroundColor.stringValue
        
        this.topBarView = new CellView("TopBarView").copyValuesRecursivelyFromObject({
            titleLabel: {
                text: "",
                textAlignment: UITextView.textAlignment.center
            },
            colors: {
                background: { normal: new UIColor("#004f20") },
                titleLabel: { normal: UIColor.whiteColor }
            }
        })
        this.view.addSubview(this.topBarView)
        
        // Bottom bar
        this.bottomBarView = new BottomBarView("BottomBarView").copyValuesRecursivelyFromObject({
            style: { overflow: "hidden" }
        })
        //this.view.addSubview(this.bottomBarView)
        
        this.backToMainButton = new CBFlatButton()
        this.backToMainButton.titleLabel.text = "&#8592;"
        this.backToMainButton.colors = {
            titleLabel: {
                normal: UIColor.whiteColor,
                highlighted: UIColor.whiteColor.colorWithAlpha(0.75),
                selected: UIColor.whiteColor.colorWithAlpha(0.5)
            },
            background: {
                normal: UIColor.transparentColor,
                highlighted: UIColor.transparentColor,
                selected: UIColor.transparentColor
            }
        }
        
        this.backToMainButton.calculateAndSetViewFrame = function (this: CBFlatButton) {
            this.setPosition(0, nil, 0, 0, nil, 50)
        }
        
        this.topBarView.addSubview(this.backToMainButton)
        
        this.backToMainButton.addControlEventTarget.EnterDown.PointerUpInside = (
            sender,
            event
        ) => UIRoute.currentRoute.routeByRemovingComponentsOtherThanOnesNamed(["settings"]).apply()
        
    }
    
    
    
    
    
    async handleRoute(route: UIRoute) {
        
        super.handleRoute(route)
        
        this.backToMainButton.hidden = NO
        
        if (IS(route.componentWithName(GETELECViewController.routeComponentName))) {
            
            // Show article view
            if (!IS(this.useGETELECViewController)) {
                
                this.useGETELECViewController = new GETELECViewController(new UIView("GETELECView"))
                
            }
            
            this.contentViewController = this.useGETELECViewController
            
        }
        else if (IS(route.componentWithName(GETELECTableViewController.routeComponentName))) {
            
            // Show article view
            if (!IS(this.useGETELECViewController)) {
                
                this.useGETELECTableViewController = new GETELECTableViewController(new UIView("GETELECTableView"))
                
            }
            
            this.contentViewController = this.useGETELECTableViewController
            
        }
        else {
            
            // Main view controller
            if (!IS(this.mainViewController)) {
                
                this.mainViewController = new UIViewController(new UIView("MainView"))
                
                this.mainViewController.view.backgroundColor = UIColor.theChosenGreyColor;
                
                this.mainViewController.view.viewHTMLElement.innerHTML = await CBCore.sharedInstance.readTextFileFromURL("/frontpage.html");
                    
                //     LanguageService.stringForKey(
                //     "frontPageHTML",
                //     LanguageService.currentLanguageKey,
                //     "This is the front page."
                // )
                
            }
            
            
            this.contentViewController = this.mainViewController
            this.backToMainButton.hidden = YES
            
        }
        
    }
    
    
    set contentViewController(controller: UIViewController) {
        
        if (this.contentViewController == controller) {
            return
        }
        
        if (this.contentViewController) {
            this.removeChildViewController(this.contentViewController)
        }
        
        this._contentViewController = controller
        this.addChildViewControllerInContainer(controller, this.view)
        this._layoutViewSubviews()
        
        this.view.setNeedsLayout()
        
    }
    
    get contentViewController(): UIViewController {
        return this._contentViewController || nil
    }
    
    
    static get paddingLength() {
        return 20
    }
    
    get paddingLength() {
        return this.class.paddingLength
    }
    
    
    layoutViewsManually() {
        
        super.layoutViewsManually()
        
        // View bounds
        var bounds = this.view.bounds
        
        this.topBarView.frame = bounds.rectangleWithHeight(50)
        
        
        this.contentViewController.view.frame = this.topBarView.frame.rectangleForNextRow(
            this.paddingLength,
            FIRST_OR_NIL(
                wrapInNil(this.contentViewController as any).intrinsicViewContentHeight([bounds.width, 1250].min()),
                wrapInNil(this.contentViewController).view.intrinsicContentHeight([bounds.width, 1250].min())
            )
        ).rectangleWithWidth([bounds.width, 1250].min(), 0.5)
        
        this.contentViewController.view.style.boxShadow = "0 3px 6px 0 rgba(0, 0, 0, 0.1)"
        
        this.contentViewController.view.setMargins(0, 0, this.paddingLength, 0)
        
        //this.view.setPaddings(0, 0, this.paddingLength, 0);
        
        var bottomBarHeight = Math.max(
            100,
            this.bottomBarView.intrinsicContentHeight(this.contentViewController.view.frame.width)
        )
        
        this.bottomBarView.frame = this.contentViewController.view.frame.rectangleWithY(
            [
                this.contentViewController.view.frame.max.y + this.paddingLength * 2,
                UIView.pageHeight - bottomBarHeight
            ].max()
        ).rectangleWithHeight(
            bottomBarHeight
        ).rectangleWithWidth(this.view.bounds.width).rectangleWithX(0)
        
    }
    
    
}








































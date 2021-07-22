/// <reference path="./UIView.ts" />
/// <reference path="./UIObject.ts" />





class UIViewController extends UIObject {
    
    
    parentViewController: UIViewController
    childViewControllers: UIViewController[]
    _UIViewController_constructorArguments: { "view": UIView; }
    static readonly routeComponentName: string;
    static readonly ParameterIdentifierName: any;
    
    
    
    constructor(public view: UIView) {
        
        super()
        
        
        this.loadIntrospectionVariables()
        this._UIViewController_constructorArguments = { "view": view }
        this._initInstanceVariables()
        
        this.loadSubviews()
        this.updateViewConstraints()
        this.updateViewStyles()
        this._layoutViewSubviews()
        
        
        
    }
    
    
    
    loadIntrospectionVariables() {
        
        this._class = UIViewController
        this.superclass = UIObject
        
        
    }
    
    _initInstanceVariables() {
        
        
        
        this.view = this._UIViewController_constructorArguments.view
        
        this.view._viewControllerLayoutFunction = this.layoutViewsManually.bind(this)
        this.view._didLayoutSubviewsDelegateFunction = this.viewDidLayoutSubviews.bind(this)
        this.view._didReceiveBroadcastEventDelegateFunction = this.viewDidReceiveBroadcastEvent.bind(this)
        
        this.childViewControllers = []
        this.parentViewController = nil
        
        
        
    }
    
    
    
    
    
    handleRouteRecursively(route: UIRoute) {
        
        this.handleRoute(route)
        
        this.childViewControllers.forEach(function (controller, index, array) {
            
            if (!route.isHandled) {
                controller.handleRouteRecursively(route)
            }
            
        })
        
    }
    
    async handleRoute(route: UIRoute) {
        
        
        
    }
    
    
    
    loadSubviews() {
    
    
    
    
    
    }
    
    async viewWillAppear() {
    
    
    
    }
    
    
    async viewDidAppear() {
    
    
    
    }
    
    
    async viewWillDisappear() {
    
    
    
    }
    
    async viewDidDisappear() {
    
    
    
    }
    
    
    updateViewConstraints() {
    
    
    
    
    
    }
    
    updateViewStyles() {
    
    
    
    
    
    }
    
    layoutViewsManually() {
        
        
        
    }
    
    _layoutViewSubviews() {
        
        this.view.layoutSubviews()
        
        this.viewDidLayoutSubviews()
        
    }
    
    
    
    viewDidLayoutSubviews() {
        
        // this.childViewControllers.forEach(function (controller, index, controllers) {
        
        //     controller._layoutViewSubviews();
        
        // })
        
        
        
    }
    
    
    
    viewDidReceiveBroadcastEvent(event: UIViewBroadcastEvent) {
        
        
        
    }
    
    
    
    
    hasChildViewController(viewController) {
        
        // This is for performance reasons
        if (!IS(viewController)) {
            return NO
        }
        
        for (var i = 0; i < this.childViewControllers.length; i++) {
    
            const childViewController = this.childViewControllers[i]
    
            if (childViewController == viewController) {
                return YES
            }
            
        }
        
        return NO
        
    }
    
    addChildViewController(viewController) {
        if (!this.hasChildViewController(viewController)) {
            viewController.willMoveToParentViewController(this)
            this.childViewControllers.push(viewController)
            //this.view.addSubview(viewController.view);
            //viewController.didMoveToParentViewController(this);
        }
    }
    
    
    removeFromParentViewController() {
        const index = this.parentViewController.childViewControllers.indexOf(this)
        if (index > -1) {
            this.parentViewController.childViewControllers.splice(index, 1)
            //this.view.removeFromSuperview();
            this.parentViewController = nil
        }
    }
    
    willMoveToParentViewController(parentViewController) {
    
    }
    
    
    didMoveToParentViewController(parentViewController) {
        
        this.parentViewController = parentViewController
        
    }
    
    removeChildViewController(controller: UIViewController) {
        
        controller.viewWillDisappear()
        if (IS(controller.parentViewController)) {
            controller.removeFromParentViewController()
        }
        if (IS(controller.view)) {
            controller.view.removeFromSuperview()
        }
        controller.viewDidDisappear()
        
    }
    
    addChildViewControllerInContainer(controller: UIViewController, containerView: UIView) {
        
        controller.viewWillAppear()
        this.addChildViewController(controller)
        containerView.addSubview(controller.view)
        controller.didMoveToParentViewController(this)
        controller.viewDidAppear()
        
    }
    
    addChildViewControllerInDialogView(controller: UIViewController, dialogView: UIDialogView) {
        
        controller.viewWillAppear()
        this.addChildViewController(controller)
        dialogView.view = controller.view
        
        var originalDismissFunction = dialogView.dismiss.bind(dialogView)
        
        dialogView.dismiss = animated => {
            
            originalDismissFunction(animated)
            
            this.removeChildViewController(controller)
            
        }
        
        controller.didMoveToParentViewController(this)
        controller.viewDidAppear()
        
    }
    
    
    
    
    
}



















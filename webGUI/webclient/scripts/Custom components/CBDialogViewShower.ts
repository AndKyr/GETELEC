/// <reference path="../UICore/UICoreExtensions.ts" />
/// <reference path="../UICore/UICore.ts" />
/// <reference path="../UICore/UIColor.ts" />
/// <reference path="./CBDialogView.ts" />
/// <reference path="../Custom components/LanguageService.ts" />





class CBDialogViewShower extends UIObject {
    
    dialogView: UIDialogView<CBDialogView>
    static currentDialogViewShower: CBDialogViewShower
    static nextShowDialogFunctions: Function[] = []
    static currentActionIndicatorDialogViewShower: CBDialogViewShower = nil
    
    
    constructor(elementID?: string) {
        
        super()
        
        this._class = CBDialogViewShower
        this.superclass = UIObject
        
        
        this.dialogView = new UIDialogView<CBDialogView>()
        
        this.dialogView.view = new CBDialogView(elementID)
        
        this.dialogView.view.backgroundColor = UIColor.whiteColor
        
        const dialogLayoutFunction = this.dialogView.layoutSubviews.bind(this.dialogView)
        
        this.dialogView.layoutSubviews = function (this: CBDialogViewShower) {
            
            dialogLayoutFunction()
            
            this.dialogView.view.frame = new UIRectangle(
                0,
                0,
                window.innerHeight,
                window.innerWidth
            ).rectangleWithHeight(
                this.dialogView.view.intrinsicContentHeight(this.getDialogWidth()),
                0.5
            ).rectangleWithWidth(this.getDialogWidth(), 0.5)
            
            this.dialogView.frame = UICore.main.rootViewController.view.bounds
            
            
        }.bind(this)
        
        
        this.dialogView.view.yesButton.addTargetForControlEvents([
            UIButton.controlEvent.PointerUpInside, UIButton.controlEvent.EnterDown
        ], function (this: CBDialogViewShower, sender: UIButton, event: Event) {
            
            this.yesButtonWasPressed()
            
            
        }.bind(this))
        
        this.dialogView.view.noButton.addTargetForControlEvents([
            UIButton.controlEvent.PointerUpInside, UIButton.controlEvent.EnterDown
        ], function (this: CBDialogViewShower, sender: UIButton, event: Event) {
            
            this.noButtonWasPressed()
            
            
        }.bind(this))
        
        
    }
    
    
    
    
    
    getDialogWidth() {
        
        const padding = RootViewController.paddingLength
        const labelHeight = padding * 0.75
        
        
        var result = 250
        
        const width = this.dialogView.view.titleLabel.intrinsicContentWidth() + padding * 2
        
        result = Math.max(result, this.dialogView.view.view.intrinsicContentWidth(this.dialogView.view.view.viewHTMLElement.naturalHeight || 1000000000))
        
        result = Math.max(result, width)
        
        result = Math.min(result, 1000)
        
        const dialogMaxWidth = (this.dialogView.superview ||
            { "bounds": new UIRectangle(0, 0, 0, result) }).bounds.width
        
        result = Math.min(result, dialogMaxWidth)
        
        return result
        
        
        
    }
    
    
    yesButtonWasPressed() {
        
        
        
    }
    
    noButtonWasPressed() {
        
        
        
    }
    
    
    cancelButtonWasPressed() {
        
        
        
    }
    
    
    
    
    
    showQuestionDialogInRootView(titleTextObject?: CBLocalizedTextObject, questionTextObject?: CBLocalizedTextObject) {
        
        
        this.dialogView.view.initTitleLabelIfNeeded()
        
        this.dialogView.view.titleLabel.localizedTextObject = titleTextObject
        
        
        
        this.dialogView.view.initQuestionLabelIfNeeded()
        
        if (IS(questionTextObject)) {
            
            this.dialogView.view.questionLabel.localizedTextObject = questionTextObject
            
        }
        
        
        this.dialogView.view.initYesNoButtonsIfNeeded()
        
        
        this.dialogView.view.noButton.addTargetForControlEvents([
            UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
        ], function (this: CBDialogViewShower, sender: UIButton, event: Event) {
            
            this.noButtonWasPressed()
            
            
        }.bind(this))
        
        this.dialogView.view.yesButton.addTargetForControlEvents([
            UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
        ], function (this: CBDialogViewShower, sender: UIButton, event: Event) {
            
            this.yesButtonWasPressed()
            
            
        }.bind(this))
        
        
        this.dialogView.showInRootView(YES)
        
        
    }
    
    showMessageDialogInRootView(titleTextObject?: CBLocalizedTextObject) {
        
        
        this.dialogView.view.initTitleLabelIfNeeded()
        
        this.dialogView.view.titleLabel.localizedTextObject = titleTextObject
        
        
        this.dialogView.view.initCancelButtonIfNeeded()
        
        
        this.dialogView.view.cancelButton.addTargetForControlEvents([
            UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
        ], function (this: CBDialogViewShower, sender: UIButton, event: Event) {
            
            this.cancelButtonWasPressed()
            
            
        }.bind(this))
        
        
        this.dialogView.showInRootView(YES)
        
        this.dialogView.view.cancelButton.focus()
        
        
    }
    
    showDialogInRootView(view: UIView) {
        
        
        this.dialogView.view.view = view
        
        this.dialogView.view.initCancelButtonIfNeeded()
        
        
        this.dialogView.view.cancelButton.addTargetForControlEvents([
            UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
        ], function (this: CBDialogViewShower, sender: UIButton, event: Event) {
            
            this.cancelButtonWasPressed()
            
            
        }.bind(this))
        
        
        this.dialogView.showInRootView(YES)
        
        this.dialogView.view.cancelButton.focus()
        
        
    }
    
    
    
    
    
    showImageDialogInRootView(imageURL: string, deleteImageCallback?: Function) {
        
        
        var loadingLabel = new UITextView()
        
        loadingLabel.text = "Loading image."
        
        loadingLabel.textAlignment = UITextView.textAlignment.center
        
        this.dialogView.view.view = loadingLabel
        
        
        const imageView = new UIImageView()
        
        imageView.imageSource = imageURL
        
        imageView.viewHTMLElement.onload = (event) => {
            this.dialogView.view.view = imageView
            
            imageView.setNeedsLayoutUpToRootView()
        }
        
        imageView.fillMode = UIImageView.fillMode.aspectFitIfLarger
        
        
        if (IS(deleteImageCallback)) {
            
            this.dialogView.view.initYesNoButtonsIfNeeded()
            
            this.dialogView.view.yesButton.titleLabel.text = "Close"
            
            this.dialogView.view.noButton.titleLabel.text = "Delete"
            
            this.dialogView.view.noButtonDismissesDialog = NO
            
            this.dialogView.view.noButton.addTargetForControlEvents([
                UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
            ], function (this: CBDialogViewShower, sender: UIButton, event: Event) {
                
                const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(function (this: CBDialogViewShower) {
                    
                    //this.dialogView.dismiss()
                    
                }.bind(this))
                
                var textObject = LanguageService.localizedTextObjectForText("Delete this image.")
                
                dialogShower.showQuestionDialogInRootView(textObject)
                
                dialogShower.yesButtonWasPressed = function () {
                    
                    deleteImageCallback()
                    
                    dialogShower.dialogView.dismiss()
                    
                }
                
                
            }.bind(this))
            
            
            this.dialogView.view.yesButton.addTargetForControlEvents([
                UIButton.controlEvent.EnterDown, UIButton.controlEvent.PointerUpInside
            ], function (this: CBDialogViewShower, sender: UIButton, event: Event) {
                
                
                this.dialogView.dismiss()
                
                
            }.bind(this))
            
        }
        else {
            
            this.dialogView.view.initCancelButtonIfNeeded()
            
            this.dialogView.view.cancelButton.titleLabel.text = "Close"
            
        }
        
        
        
        
        this.dialogView.showInRootView(YES)
        
        this.dialogView.view.cancelButton.focus()
        
        
    }
    
    showActionIndicatorDialogInRootView(message: string) {
        
        
        
        const actionIndicator = new UIActionIndicator()
        
        
        this.dialogView.zIndex = 150
        
        this.dialogView.view.view = actionIndicator
        
        actionIndicator.style.minHeight = "100px"
        
        this.dialogView.view.initQuestionLabelIfNeeded()
        
        this.dialogView.view.questionLabel.text = message
        
        
        actionIndicator.start()
        
        this.dialogView.view.backgroundColor = CBColor.primaryContentColor.colorWithAlpha(0.5)
        
        this.dialogView.view.questionLabel.textColor = UIColor.whiteColor
        
        
        
        this.dialogView.dismissesOnTapOutside = NO
        
        
        
        
        
        CBDialogViewShower.currentActionIndicatorDialogViewShower = this
        
        
        this.dialogView.showInRootView(NO)
        
        this.dialogView.view.cancelButton.focus()
        
        
    }
    
    
    
    
    
    static showNextDialog() {
        
        (CBDialogViewShower.nextShowDialogFunctions.firstElement || nil)()
        CBDialogViewShower.nextShowDialogFunctions.removeElementAtIndex(0)
        
    }
    
    
    
    static alert(text: string, dismissCallback: Function = nil) {
        
        const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback)
        const textObject = LanguageService.localizedTextObjectForText(text)
        const showDialogFunction = dialogShower.showMessageDialogInRootView.bind(dialogShower, textObject)
        CBDialogViewShower._showDialogWithFunction(showDialogFunction, dialogShower)
        
    }
    
    
    
    static localizedAlert(textObject: CBLocalizedTextObject, dismissCallback: Function = nil) {
        
        const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback)
        const showDialogFunction = dialogShower.showMessageDialogInRootView.bind(dialogShower, textObject)
        CBDialogViewShower._showDialogWithFunction(showDialogFunction, dialogShower)
        
    }
    
    
    static showQuestionDialog(questionText: string, dismissCallback: Function = nil) {
        
        const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback)
        const textObject = LanguageService.localizedTextObjectForText(questionText)
        const showDialogFunction = dialogShower.showQuestionDialogInRootView.bind(dialogShower, textObject)
        CBDialogViewShower._showDialogWithFunction(showDialogFunction, dialogShower)
        
        return dialogShower
        
    }
    
    
    static showImageDialog(imageURL: string, deleteImageCallback: Function = nil, dismissCallback: Function = nil) {
        
        var dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback)
        
        var showDialogFunction = dialogShower.showImageDialogInRootView.bind(
            dialogShower,
            imageURL,
            deleteImageCallback
        )
        CBDialogViewShower._showDialogWithFunction(showDialogFunction, dialogShower)
        
        return dialogShower
        
        
    }
    
    static showDialog(view: UIView, dismissCallback: Function = nil) {
        
        const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback)
        const showDialogFunction = dialogShower.showDialogInRootView.bind(dialogShower, view)
        CBDialogViewShower._showDialogWithFunction(showDialogFunction, dialogShower)
    
        return dialogShower;
    
    }
    
    static showActionIndicatorDialog(message: string, dismissCallback: Function = nil) {
        
        
        if (IS(CBDialogViewShower.currentActionIndicatorDialogViewShower)) {
            
            
            CBDialogViewShower.currentActionIndicatorDialogViewShower.dialogView.view.questionLabel.text = message
            
            CBDialogViewShower.currentActionIndicatorDialogViewShower.dialogView.view.setNeedsLayoutUpToRootView()
            
            return
            
        }
        
        
        const dialogShower = CBDialogViewShower._dialogShowerWithDismissCallback(dismissCallback)
        
        dialogShower.showActionIndicatorDialogInRootView(message)
        
        return dialogShower
        
        
    }
    
    
    static hideActionIndicatorDialog() {
        
        CBDialogViewShower.currentActionIndicatorDialogViewShower.dialogView.dismiss()
        
        CBDialogViewShower.currentActionIndicatorDialogViewShower = nil
        
    }
    
    
    
    private static _dialogShowerWithDismissCallback(dismissCallback: Function) {
        
        const dialogShower = new CBDialogViewShower()
        const dismissFunction = dialogShower.dialogView.dismiss.bind(dialogShower.dialogView)
        
        dialogShower.dialogView.dismiss = function () {
            dismissFunction()
            dismissCallback()
            CBDialogViewShower.currentDialogViewShower = null
            CBDialogViewShower.showNextDialog()
        }
        
        return dialogShower
        
    }
    
    
    private static _showDialogWithFunction(showDialogFunction: any, dialogShower: CBDialogViewShower) {
        
        if (IS(CBDialogViewShower.currentDialogViewShower)) {
            
            CBDialogViewShower.nextShowDialogFunctions.push(showDialogFunction)
            
        }
        else {
            
            CBDialogViewShower.currentDialogViewShower = dialogShower
            showDialogFunction()
            
        }
        
    }
    
    
    
    
}







































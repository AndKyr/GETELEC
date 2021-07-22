/// <reference path="./CBFlatButton.ts" />
/// <reference path="./RowView.ts" />
/// <reference path="../UICore/UICoreExtensions.ts" />



class CBDialogView extends UIView {
    
    
    
    
    
    titleLabel: UITextView = nil
    titleRow: RowView = nil
    placeholderLabel: UITextView = nil
    
    _view: UIView = nil
    
    questionLabel: UITextView = nil
    yesButton: CBFlatButton = nil
    noButton: CBFlatButton = nil
    cancelButton: CBFlatButton = nil
    
    
    noButtonDismissesDialog = YES
    
    
    constructor(elementID?: string) {
        
        super(elementID)
        
        
        this._class = CBDialogView
        this.superclass = UIView
        
        
        
        
        
        this.updateContent()
        
        
        
    }
    
    
    
    
    
    initView(elementID: string, viewHTMLElement: HTMLElement) {
        
        
        super.initView(elementID, viewHTMLElement)
        
        
        
        // this.initPlaceholderLabelIfNeeded()
        
        // this.initTitleLabelIfNeeded()
        
        // this.initTitleRow()
        
        
        // this.titleRow.cellWeights = [1, 0.75, 1.2, 1]
        
        // var cells = this.titleRow.cells as CellView[];
        
        // var cellTitles = ["TIME:", "COMPANY:", "COMMENTS:", "PRICE:"]
        
        // cells.forEach(function (cell, index, array) {
        
        //     cell.titleLabel.text = cellTitles[index];
        //     cell.titleLabel.useSairaFont()
        
        // })
        
        
        
        
        
        // this.initQuestionLabelIfNeeded()
        
        
        // this.initYesNoButtonsIfNeeded()
        
        
        
        
        
        this.style.borderRadius = "5px"
        
        
        
    }
    
    
    
    
    
    set view(view: UIView) {
        
        this.view.removeFromSuperview()
        
        this._view = view
        
        this.addSubview(view)
        
    }
    
    get view() {
        
        return this._view
        
    }
    
    
    
    
    
    initCancelButtonIfNeeded() {
        
        if (IS_NOT(this.cancelButton)) {
            
            this.cancelButton = new CBFlatButton(this.elementID + "CancelButton")
            this.cancelButton.titleLabel.text = "OK"
            this.cancelButton.titleLabelColor = UIColor.greenColor
            this.cancelButton.titleLabel.style.fontWeight = "600"
            this.addSubview(this.cancelButton)
            
            this.cancelButton.addTargetForControlEvents([
                UIButton.controlEvent.PointerUpInside, UIButton.controlEvent.EnterDown
            ], function (this: CBDialogView, sender: UIButton, event: Event) {
                this.dialogView.dismiss()
            }.bind(this))
            
        }
        
    }
    
    initYesNoButtonsIfNeeded() {
        
        if (IS_NOT(this.yesButton)) {
            
            this.yesButton = new CBFlatButton(this.elementID + "YesButton")
            this.yesButton.titleLabel.setText("cBDDialogViewYES", "YES")
            this.yesButton.titleLabelColor = UIColor.greenColor
            this.yesButton.titleLabel.style.fontWeight = "600"
            this.addSubview(this.yesButton)
            
            this.noButton = new CBFlatButton(this.elementID + "NoButton")
            this.noButton.titleLabel.setText("cBDDialogViewNO", "NO")
            this.noButton.titleLabelColor = UIColor.redColor
            this.noButton.titleLabel.style.fontWeight = "600"
            this.addSubview(this.noButton)
            
            // this.yesButton.addTargetForControlEvents([UIButton.controlEvent.PointerUpInside, UIButton.controlEvent.EnterDown], function (this: CBDialogView, sender: UIButton, event: Event) {
            //     this.dialogView.dismiss()
            // }.bind(this))
            
            this.noButton.addTargetForControlEvents([
                UIButton.controlEvent.PointerUpInside, UIButton.controlEvent.EnterDown
            ], function (this: CBDialogView, sender: UIButton, event: Event) {
                if (this.noButtonDismissesDialog) {
                    this.dialogView.dismiss()
                }
            }.bind(this))
            
        }
        
        
        
    }
    
    initQuestionLabelIfNeeded() {
        
        if (IS_NOT(this.questionLabel)) {
            
            this.questionLabel = new UITextView(this.elementID + "QuestionLabel", UITextView.type.header3)
            this.questionLabel.text = LanguageService.stringForKey(
                "cBDDialogViewAreYouSure",
                LanguageService.currentLanguageKey,
                "Are you sure?"
            )
            this.questionLabel.textAlignment = UITextView.textAlignment.center
            this.addSubview(this.questionLabel)
            
        }
        
    }
    
    initTitleRow() {
        
        if (IS_NOT(this.titleRow)) {
            
            this.titleRow = new RowView(this.elementID + "TitleRow")
            this.addSubview(this.titleRow)
            
        }
        
    }
    
    initTitleLabelIfNeeded() {
        
        if (IS_NOT(this.titleLabel)) {
            
            this.titleLabel = new UITextView(this.elementID + "TitleLabel", UITextView.type.header4)
            this.titleLabel.text = "Accept offer"
            this.titleLabel.style.fontSize = "24"
            
            this.addSubview(this.titleLabel)
            
        }
        
    }
    
    initPlaceholderLabelIfNeeded() {
        
        if (IS_NOT(this.placeholderLabel)) {
            
            this.placeholderLabel = new UITextView(this.elementID + "PlaceholderLabel", UITextView.type.header4)
            this.placeholderLabel.text = "No offers have been made yet."
            this.addSubview(this.placeholderLabel)
            
        }
        
    }
    
    updateContent() {
        
        function stringFromValue(value) {
            if (IS(value)) {
                return value
            }
            return "-"
        }
        
        
        
        this.layoutSubviews()
        
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
    
        const padding = RootViewController.paddingLength
        const labelHeight = padding * 0.75
    
        const bounds = this.bounds.rectangleWithInsets(padding * 0.5, padding * 0.5, 0, 0)
    
        const topObject = { frame: bounds.rectangleWithHeight(0) }
    
    
        this.titleLabel.frame = bounds.rectangleWithHeight(this.titleLabel.intrinsicContentHeight(bounds.width)).rectangleWithY(
            bounds.y +
            padding * 0.5)
        this.titleLabel.style.marginLeft = "" + (padding * 0.5) + "px"
        
        this.titleRow.frame = FIRST_OR_NIL(this.titleLabel, topObject).frame.rectangleForNextRow(
            0,
            this.titleRow.intrinsicContentHeight(bounds.width)
        )
        
        
        
        
        
        this.view.frame = FIRST_OR_NIL(this.titleRow, this.titleLabel, topObject).frame.rectangleForNextRow(
            padding,
            this.view.intrinsicContentHeight(bounds.width)
        ).rectangleWithWidth(bounds.width).rectangleWithX(bounds.x)
        
        
        
        this.questionLabel.frame = FIRST_OR_NIL(
            this.view,
            this.titleRow,
            this.titleLabel,
            topObject
        ).frame.rectangleForNextRow(
            padding,
            this.questionLabel.intrinsicContentHeight(bounds.width)
        ).rectangleWithWidth(bounds.width).rectangleWithX(bounds.x)
    
    
    
    
        const buttons = [this.yesButton, this.noButton]
    
        const buttonsFrame = FIRST_OR_NIL(
            this.questionLabel,
            this.view,
            this.titleRow,
            this.titleLabel,
            topObject
        ).frame.rectangleForNextRow(padding, labelHeight *
            2)
            .rectangleWithWidth([(buttons.everyElement.titleLabel.intrinsicContentWidth() as any as number[]).max(), 120].max() * 2 + padding * 3, 0.5)
        
        buttonsFrame.distributeViewsAlongWidth(buttons)
        
        this.cancelButton.frame = buttonsFrame.rectangleWithWidth(150, 0.5)
        
        
        
        
        this.placeholderLabel.frame = bounds.rectangleWithInsets(padding * 0.5, padding * 0.5, 0, 0)
        
        
        
        if (IS(this.view)) {
            
            this.placeholderLabel.hidden = YES
            
            // this.titleLabel.hidden = NO
            // this.titleRow.hidden = NO
            
        }
        else {
            
            
            this.placeholderLabel.hidden = NO
            
            // this.titleLabel.hidden = YES
            // this.titleRow.hidden = YES
            
            
            
        }
        
        
        
        
        
    }
    
    
    
    intrinsicContentHeight(constrainingWidth: number = 0) {
    
    
        const padding = RootViewController.paddingLength
        const labelHeight = padding * 0.75
    
        function addPaddingIfNeeded(value: number) {
            
            return IF<number>(value)(function () {
                return value + padding
            })()
            
        }
    
        var result = padding
    
        result = result + addPaddingIfNeeded(this.titleLabel.intrinsicContentHeight(constrainingWidth))
        
        result = result + addPaddingIfNeeded(this.titleRow.intrinsicContentHeight(constrainingWidth))
        
        result = result + addPaddingIfNeeded(this.view.intrinsicContentHeight(constrainingWidth) + padding)
        
        result = result + addPaddingIfNeeded(this.questionLabel.intrinsicContentHeight(constrainingWidth))
        
        result = result + IF<number>(FIRST_OR_NIL(this.yesButton, this.cancelButton))(RETURNER(labelHeight))()
        
        
        
        
        return result
        
    }
    
    
    
    
    
}



































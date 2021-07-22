/// <reference path="../UICore/UIButton.ts" />
/// <reference path="../UICore/UIKeyValueStringFilter.ts" />
/// <reference path="../UICore/UITableView.ts" />
/// <reference path="../UICore/UIImageView.ts" />
/// <reference path="../UICore/UITextView.ts" />
/// <reference path="../UICore/UITextField.ts" />
/// <reference path="../UICore/UIDialogView.ts" />
/// <reference path="./SearchTextField.ts" />
/// <reference path="./CBCheckbox.ts" />
/// <reference path="./CBColor.ts" />
/// <reference path="./SearchableDropdownRow.ts" />
/// <reference path="../RootViewController.ts" />
/// <reference path="./LanguageService.ts" />





interface CBLocalizedTextObject {
    
    [key: string]: string
    
}


interface CBDropdownData<T> {
    
    _id: string;
    name?: CBLocalizedTextObject;
    dropdownCode: string;
    data: CBDropdownDataItem<T>[];
    
}


interface CBDropdownDataItem<T> {
    
    _id: string;
    title: CBLocalizedTextObject;
    rowsData?: CBDropdownDataItem<T>[]
    isADropdownDataSection: boolean;
    isADropdownDataRow: boolean;
    
    attachedObject: T
    
    itemCode: string;
    dropdownCode: string;
    
}


interface CBSubjectDropdownItemAttachedObject {
    
    isAProduct?: boolean;
    isAService?: boolean;
    
}




class SearchableDropdown<T> extends UIButton {
    
    
    
    
    
    _keyValueStringFilter: UIKeyValueStringFilter
    _containerView: UIView
    _searchTextField: SearchTextField
    _dialogView: UIDialogView
    _tableView: UITableView
    _rightImageView: UIImageView
    
    highlightLabel: UITextView
    
    _data: CBDropdownDataItem<T>[] = []
    _filteredData: CBDropdownDataItem<T>[] = []
    
    _excludedData: string[] = []
    
    
    tintColor: UIColor = CBColor.primaryTintColor
    
    
    overflowLabel: UITextView
    selectedIndices: number[] = []
    _selectedData: CBDropdownDataItem<T>[] = []
    _drawingData: CBDropdownDataItem<T>[] = []
    
    _isDrawingDataValid = NO
    _placeholderText: string
    _expandedContainerViewHeight: number
    
    isSingleSelection = NO
    showsSelectedSectionInMultipleSelectionMode = NO
    
    allowsCustomItem = NO
    _customItem: CBDropdownDataItem<any>
    _focusedRowIndex: number
    
    keepFocusedRowVisible = YES
    _placeholderLocalizedTextObject: UILocalizedTextObject
    
    
    constructor(elementID: string) {
        
        super(elementID)
        
        
        
        this._class = SearchableDropdown
        this.superclass = UIButton
        
        
        
        
        
    }
    
    
    
    initView(elementID, viewHTMLElement, initViewData) {
        
        super.initView(elementID, viewHTMLElement, initViewData)
        
        
        //this.style.borderRadius = "2px";
        
        
        this._titleLabel.text = "Current value"
        this._titleLabel.textAlignment = UITextView.textAlignment.left
        
        this.overflowLabel = new UITextView(elementID + "OverflowLabel")
        this.overflowLabel.textColor = CBColor.primaryContentColor
        this.overflowLabel.textAlignment = UITextView.textAlignment.right
        this.addSubview(this.overflowLabel)
        
        
        this._rightImageView = new UIImageView(this.elementID + "RightImageView")
        this._rightImageView.imageSource = "images/baseline-arrow_drop_down-24px.svg"
        this._rightImageView.userInteractionEnabled = NO
        this.addSubview(this._rightImageView)
        
        
        this.setNeedsLayout()
        
        
        this._containerView = new UIView(elementID + "ContainerView")
        this._containerView.style.boxShadow = "0 9px 13px 0 rgba(0,0,0,0.26)"
        this._containerView.style.borderRadius = "2px"
        
        this._searchTextField = new SearchTextField(elementID + "SearchTextField")
        this._searchTextField.placeholderText = "Search"
        this._containerView.addSubview(this._searchTextField)
        
        this._searchTextField._textField.addTargetForControlEvent(
            UITextField.controlEvent.TextChange,
            function (this: SearchableDropdown<T>, sender, event) {
                
                this.updateFilteredData(this._searchTextField.text)
                
            }.bind(this)
        )
        
        this._searchTextField._textField.addTargetForControlEvent(
            UIView.controlEvent.EscDown,
            function (this: SearchableDropdown<T>, sender, event) {
                
                if (IS(this._searchTextField.text)) {
                    
                    this._searchTextField.text = ""
                    this.updateFilteredData("")
                    
                }
                else {
                    
                    this._dialogView.dismiss(YES)
                    
                }
                
                
            }.bind(this)
        )
        
        this._searchTextField._textField.addTargetForControlEvent(
            UIView.controlEvent.DownArrowDown,
            function (this: SearchableDropdown<T>, sender: UITextField, event: Event) {
                
                if (this.focusedRowIndex < (this.drawingData.length - 1)) {
                    
                    this.focusedRowIndex = this.focusedRowIndex + 1
                    
                }
                
            }.bind(this)
        )
        
        this._searchTextField._textField.addTargetForControlEvent(
            UIView.controlEvent.UpArrowDown,
            function (this: SearchableDropdown<T>, sender: UITextField, event: Event) {
                
                if (this.focusedRowIndex > 0) {
                    
                    this.focusedRowIndex = this.focusedRowIndex - 1
                    
                }
                
            }.bind(this)
        )
        
        this._searchTextField._textField.addTargetForControlEvent(
            UIView.controlEvent.EnterDown,
            function (this: SearchableDropdown<T>, sender, event) {
    
                const isTouchDevice = "ontouchstart" in document.documentElement
    
                if (isTouchDevice) {
                    
                    this._searchTextField.blur()
                    
                    return
                    
                }
    
    
                const datapoint = this.drawingData[this.focusedRowIndex]
    
    
                const alreadySelected = this.selectedDataContains(datapoint)
    
                if (alreadySelected) {
                    
                    this.selectedData.removeElement(datapoint)
                    
                }
                else if (this.isSingleSelection) {
                    
                    
                    this.selectedIndices = [this.focusedRowIndex]
                    
                    this.selectedData = [datapoint]
                    
                    // view.selected = YES;
                    
                    this.selectionDidChange(this.selectedData)
                    
                    this.updateContentForCurrentSelection()
                    
                    this._dialogView.dismiss()
                    
                    
                }
                else {
                    
                    this.selectedData.push(datapoint)
                    
                }
                
                
                
                
                
            }.bind(this)
        )
        
        this._tableView = new UITableView(elementID + "TableView")
        this._containerView.addSubview(this._tableView)
        this._tableView.backgroundColor = UIColor.whiteColor
        
        this._dialogView = new UIDialogView(elementID + "DialogView")
        this._dialogView.view = this._containerView
        this._dialogView.backgroundColor = UIColor.transparentColor
        
        this.addTargetForControlEvents([
            UIView.controlEvent.PointerUpInside, UIView.controlEvent.EnterDown
        ], function (this: SearchableDropdown<T>, sender, event) {
            if (this._dialogView.isVisible) {
                this.closeDropdown()
            }
            else {
                this.openDropdown()
            }
        }.bind(this))
        
        
        this._dialogView.addTargetForControlEvent(
            UIView.controlEvent.PointerDown,
            function (sender: UIDialogView, event: Event) {
                
                if (sender.viewHTMLElement == event.target) {
                    sender.dismiss()
                }
                
            }
        )
    
    
        const dialogLayoutFunction = this._dialogView.layoutSubviews
    
        this._dialogView.layoutSubviews = function (this: SearchableDropdown<T>) {
            
            //dialogLayoutFunction.call(this._dialogView);
            
            this._dialogView.frame = this.rootView.bounds //this.superview.rectangleInView(this.superview.bounds, this.rootView);
    
            const padding = RootViewController.paddingLength
            const labelHeight = padding
    
            const searchTextFieldHeight = this.bounds.height
    
            this._containerView.frame = this.superview.rectangleInView(this.frame, this.rootView)
                .rectangleWithHeight(this.expandedContainerViewHeight)
            
            this._searchTextField.frame = this._containerView.bounds.rectangleWithHeight(searchTextFieldHeight).rectangleWithInsets(
                0,
                16,
                0,
                0
            )
            
            this._tableView.frame = this._containerView.bounds.rectangleWithInsets(0, 0, 0, searchTextFieldHeight)
            
            //this._tableView.style.maxHeight = "" + this._tableView.intrinsicContentHeight() + "px";
            
            
        }.bind(this)
        
        
        // var animateDialogAppearing = this._dialogView.animateAppearing;
        // this._dialogView.animateAppearing = function (this: SearchableDropdown) {
        //     animateDialogAppearing.call(this._dialogView);
        //     this.style.zIndex = "101";
        //     if (this.highlightLabel) {
        //         this.highlightLabel.textColor = UIColor.whiteColor;
        //         this.highlightLabel.style.zIndex = "101";
        //         this.highlightLabel.userInteractionEnabled = NO;
        //     }
        
        //     this._searchTextField.focus();
        
        // }.bind(this)
        
        // var animateDialogDisappearing = this._dialogView.animateDisappearing;
        // this._dialogView.animateDisappearing = function (this: SearchableDropdown) {
        //     animateDialogDisappearing.call(this._dialogView);
        //     this.style.zIndex = "0";
        //     if (this.highlightLabel) {
        //         this.highlightLabel.textColor = CBColor.primaryContentColor;
        //         this.highlightLabel.style.zIndex = "0";
        //         this.highlightLabel.userInteractionEnabled = YES;
        //     }
        
        //     this._searchTextField.blur();
        
        // }.bind(this)
        
        
        
        
        
        this._tableView.numberOfRows = function (this: SearchableDropdown<T>) {
    
            var result = this.drawingData.length
    
            if (IS_NOT(this.isSingleSelection) && this.showsSelectedSectionInMultipleSelectionMode) {
                
                result = result + this.selectedData.length
                
            }
            
            
            return result
            
        }.bind(this)
    
    
        const newReusableViewForIdentifierFunction = this._tableView.newReusableViewForIdentifier.bind(this._tableView)
    
        this._tableView.newReusableViewForIdentifier = function (
            this: SearchableDropdown<T>,
            identifier: string,
            rowIndex: number
        ) {
            
            //if (identifier == "SelectedItemRow") {
    
            const view = new SearchableDropdownRow(elementID + identifier + rowIndex)
    
            view.stopsPointerEventPropagation = NO
            view.pausesPointerEvents = NO
            
            return view
            
            // }
            
            // return newReusableViewForIdentifierFunction(identifier, rowIndex);
            
        }.bind(this)
    
    
        const viewForSelectedItemRow = function (this: SearchableDropdown<T>, index: number) {
        
            const view = this._tableView.reusableViewForIdentifier("SelectedItemRow", index) as CBCheckbox
        
            view.titleLabel.text = LanguageService.stringForCurrentLanguage(this.selectedData[index].title)
        
            view.selected = YES
        
            return view
        
        }.bind(this)
    
        
    
    
    
        this._tableView.viewForRowWithIndex = function (this: SearchableDropdown<T>, index: number) {
    
    
    
    
            const view = this._tableView.reusableViewForIdentifier("SubjectView", index) as SearchableDropdownRow
    
    
            view.style.borderBottomColor = ""
            view.style.borderBottomStyle = ""
            view.style.borderBottomWidth = ""
    
    
    
    
            const rowWasHovered = function (this: SearchableDropdown<T>, sender: SearchableDropdownRow, event: Event) {
        
                this.focusedRowIndex = index
        
            }.bind(this)
    
            view.addTargetForControlEvent(UIButton.controlEvent.PointerHover, rowWasHovered)
            
            view.removeTargetForControlEvent(
                UIButton.controlEvent.PointerHover,
                view._SearchableDropdownRowWasHoveredFunction
            )
            
            view._SearchableDropdownRowWasHoveredFunction = rowWasHovered
            
            view.focused = (this.focusedRowIndex == index)
            
            
            
            if (!this.isSingleSelection && this.showsSelectedSectionInMultipleSelectionMode) {
                
                if (index < this.selectedData.length) {
                    
                    
                    view.type = SearchableDropdownRow.type.selectedItem
                    
                    view.titleText = LanguageService.stringForCurrentLanguage(this.selectedData[index].title)
                    
                    view.selected = YES
                    
                    view.updateContentForCurrentState()
                    
                    
                    if (index == this.selectedData.length - 1) {
                        
                        view.style.borderBottomColor = UIColor.colorWithRGBA(100, 100, 100).stringValue
                        
                        view.style.borderBottomStyle = "solid"
                        
                        view.style.borderBottomWidth = "1px"
                        
                    }
                    
                    var viewWasTapped = function (this: SearchableDropdown<T>, sender, event) {
                        
                        
                        this.selectedIndices.removeElementAtIndex(index)
    
                        const selectedItem = this.selectedData[index]
    
                        this.selectedData.removeElement(selectedItem)
                        
                        
                        view.selected = NO
                        
                        this.selectionDidChange(this.selectedData)
                        
                        //this.performFunctionWithDelay(0.15, function (this: SearchableDropdown<T>) {
                        
                        this.updateContentForCurrentSelection()
                        
                        //}.bind(this))
                        
                        // // Unselecting main body item if possible
                        // var visibleRows = this._tableView._visibleRows as SearchableDropdownRow[]
                        
                        // var visibleIndices = this._tableView.indexesForVisibleRows()
                        
                        // visibleRows.forEach(function (this: SearchableDropdown<T>, visibleRow: SearchableDropdownRow, visibleRowIndex: number, visibleRowsArray: SearchableDropdownRow[]) {
                        
                        //     var rowDataIndex = visibleIndices[visibleRowIndex]
                        
                        //     visibleRow.selected = this.selectedData.contains(this.drawingData[rowDataIndex])
                        
                        // }.bind(this))
                        
                        
                        
                        
                        this._searchTextField.focus()
                        
                        
                        if ((view as any).viewWasTappedFunction) {
                            view.removeTargetForControlEvents([
                                UIView.controlEvent.EnterDown, UIView.controlEvent.PointerTap
                            ], (view as any).viewWasTappedFunction)
                        }
                        
                        
                    }.bind(this)
                    
                    
                    // Removing previous event target if possible
                    if ((view as any).viewWasTappedFunction) {
                        view.removeTargetForControlEvents([
                            UIView.controlEvent.EnterDown, UIView.controlEvent.PointerTap
                        ], (view as any).viewWasTappedFunction)
                    }
                    
                    
                    // Adding event target
                    view.addTargetForControlEvents([
                        UIView.controlEvent.EnterDown, UIView.controlEvent.PointerTap
                    ], viewWasTapped);
                    (view as any).viewWasTappedFunction = viewWasTapped
                    
                    
                    return view
                    
                    
                    
                }
                
                index = index - this.selectedData.length
                
            }
            
            
            
            // Datapoint
            const datapoint = this.drawingData[index]
    
            if (IS_NOT(datapoint)) {
                return
            }
            
            // Setting different style for section title rows
            if (datapoint.isADropdownDataSection) {
                
                view.type = SearchableDropdownRow.type.sectionTitle
                
                view.userInteractionEnabled = NO
                
            }
            else {
                
                view.type = SearchableDropdownRow.type.selectableItem
                
                view.userInteractionEnabled = YES
                
            }
            
            if (datapoint._id == (this._customItem || nil)._id) {
                
                view.type = SearchableDropdownRow.type.customItem
                
            }
            
            view.updateContentForNormalState()
            
            view.updateContentForCurrentState()
            
            
            // Setting row title
            view.titleText = LanguageService.stringForCurrentLanguage(datapoint.title)
            view.titleLabel.textAlignment = UITextView.textAlignment.left
            
            
            // Highlighting of row
            //view.highlighted = (this.highlightedRowIndex == index);
            
            
            // Selecting of rows
            
            view.selected = this.selectedRowIdentifiers.contains(datapoint._id) //this.selectedIndices.contains(index); //this.selectedDataContains(datapoint);
            
            var viewWasTapped = function (this: SearchableDropdown<T>, sender, event) {
                
                if (view.selected) {
                    
                    this.selectedIndices.removeElement(index)
                    
                    this.selectedData.removeElement(datapoint)
                    
                }
                else {
                    
                    
                    if (this.isSingleSelection) {
                        
                        this.selectedIndices = [index]
                        
                        this.selectedData = [datapoint]
                        
                        // view.selected = YES;
                        
                        this.selectionDidChange(this.selectedData)
                        
                        this.updateContentForCurrentSelection()
                        
                        this._dialogView.dismiss()
                        
                        return
                        
                    }
                    else {
                        
                        this.selectedIndices.push(index)
                        
                        this.selectedData.push(datapoint)
                        
                        
                        
                    }
                    
                    
                    
                }
                
                //view.selected = !view.selected;
    
                const selectedData = this.selectedData
    
                if (!view.selected) {
                    
                    
                    view.selected = YES
                    
                    this.performFunctionWithDelay(0.25, function (this: SearchableDropdown<T>) {
                        
                        
                        this.selectionDidChange(selectedData)
                        
                        this.updateContentForCurrentSelection()
                        
                        if (this.showsSelectedSectionInMultipleSelectionMode) {
                            
                            this._tableView.contentOffset = this._tableView.contentOffset.pointByAddingY(view.frame.height)
                            
                        }
                        
                        
                    }.bind(this))
                    
                    
                }
                else {
                    
                    view._checkbox.selected = NO
                    
                    this.selectionDidChange(selectedData)
                    
                    this.performFunctionWithDelay(0.15, function (this: SearchableDropdown<T>) {
                        
                        
                        view.selected = NO
                        
                        this.updateContentForCurrentSelection()
                        
                        if (this.showsSelectedSectionInMultipleSelectionMode) {
                            
                            this._tableView.contentOffset = this._tableView.contentOffset.pointByAddingY(-view.frame.height)
                            
                        }
                        
                        
                        
                    }.bind(this))
                    
                    
                    
                }
                
                
                
                this._searchTextField.focus()
                
            }.bind(this)
            
            // Removing previous event target if possible
            if ((view as any).viewWasTappedFunction) {
                view.removeTargetForControlEvents([
                    UIView.controlEvent.EnterDown, UIView.controlEvent.PointerUpInside
                ], (view as any).viewWasTappedFunction)
            }
            
            // Adding event target
            view.addTargetForControlEvents([
                UIView.controlEvent.EnterDown, UIView.controlEvent.PointerUpInside
            ], viewWasTapped);
            (view as any).viewWasTappedFunction = viewWasTapped
            
            return view
            
        }.bind(this)
        
        this._keyValueStringFilter = new UIKeyValueStringFilter()
        
    }
    
    
    
    
    
    openDropdown() {
        this._dialogView.showInView(this.rootView, YES)
        this._searchTextField.focus()
    }
    
    closeDropdown() {
        this._dialogView.dismiss(YES)
    }
    
    boundsDidChange() {
        
        super.boundsDidChange()
        
        this.setNeedsLayout()
        
    }
    
    
    
    
    
    
    
    
    get focusedRowIndex() {
        
        return this._focusedRowIndex
        
    }
    
    set focusedRowIndex(focusedRowIndex: number) {
    
        const previousFocusedRowIndex = this.focusedRowIndex
    
        this._focusedRowIndex = focusedRowIndex
        
        
        if (previousFocusedRowIndex != focusedRowIndex) {
            
            (this._tableView.visibleRowWithIndex(previousFocusedRowIndex) as SearchableDropdownRow).focused = NO
    
            const focusedRow = this._tableView.visibleRowWithIndex(this.focusedRowIndex) as SearchableDropdownRow
    
            focusedRow.focused = YES
            
            if (!this.keepFocusedRowVisible) {
                
                return
                
            }
    
            var contentOffset = this._tableView.contentOffset
    
            if (focusedRow.frame.y < contentOffset.y) {
                
                contentOffset.y = focusedRow.frame.y
                
            }
            
            if (focusedRow.frame.max.y > (contentOffset.y + this._tableView.bounds.height)) {
                
                contentOffset = contentOffset.pointByAddingY(-(contentOffset.y + this._tableView.bounds.height -
                    focusedRow.frame.max.y))
                
            }
    
            const animationDuration = this._tableView.animationDuration
    
            this._tableView.animationDuration = 0
            
            this._tableView.contentOffset = contentOffset
            
            this._tableView.animationDuration = animationDuration
            
        }
        
        
    }
    
    
    set expandedContainerViewHeight(expandedContainerViewHeight: number) {
        
        
        this._expandedContainerViewHeight = expandedContainerViewHeight
        
        this._dialogView.setNeedsLayout()
        
        
    }
    
    
    get expandedContainerViewHeight() {
        
        if (IS(this._expandedContainerViewHeight)) {
            
            return this._expandedContainerViewHeight
            
        }
    
        const padding = RootViewController.paddingLength
        const labelHeight = padding
    
        const result = this.superview.bounds.height - this.frame.max.y - padding
    
    
        return result
        
    }
    
    
    
    
    
    selectedDataContains(datapoint: CBDropdownDataItem<T>) {
        for (var i = 0; i < this.selectedData.length; i++) {
            const value = this.selectedData[i]
            if (value._id == datapoint._id) {
                return YES
            }
        }
        return NO
    }
    
    
    
    
    
    updateContentForNormalState() {
        
        
        this.style.borderBottom = "1px solid rgba(0,0,0,0.12)"
        
        this.titleLabel.textColor = CBColor.primaryContentColor
        this.backgroundColor = UIColor.transparentColor
        
        this.style.borderBottomColor = CBColor.primaryContentColor.colorWithAlpha(0.12).stringValue
        
        
    }
    
    updateContentForHighlightedState() {
        
        this.style.borderBottomWidth = "2px"
        
        this.style.borderBottomColor = this.tintColor.stringValue
        
    }
    
    
    static controlEvent = Object.assign({}, UIView.controlEvent, {
        
        "SelectionDidChange": "SelectionDidChange"
        
    })
    
    
    
    selectionDidChange(selectedRows: CBDropdownDataItem<T>[]) {
        
        // Update title to indicate the selection
        this.updateTitleWithSelection(selectedRows)
        
        this.sendControlEventForKey(SearchableDropdown.controlEvent.SelectionDidChange, nil)
        
    }
    
    updateContentForCurrentSelection() {
        
        // Update selected rows in the top of the table view
        this._tableView.reloadData()
        
        this.setNeedsLayout()
        
    }
    
    
    get placeholderText() {
        
        if (IS_UNDEFINED(this._placeholderText)) {
            
            
            this._placeholderText = "Not selected"
            
            
        }
        
        return this._placeholderText
        
    }
    
    set placeholderText(placeholderText: string) {
        
        this._placeholderText = placeholderText
        
        this.updateTitleWithSelection(this.selectedData)
        
    }
    
    
    setPlaceholderText(key: string, defaultString: string, parameters?: { [x: string]: string }) {
        
        this.placeholderLocalizedTextObject = LanguageService.localizedTextObjectForKey(key, defaultString, parameters)
        
        //this.titleLabel.setText(key, defaultString, parameters);
        
    }
    
    get placeholderLocalizedTextObject() {
        
        if (IS_UNDEFINED(this._placeholderLocalizedTextObject)) {
            
            
            this._placeholderLocalizedTextObject = LanguageService.localizedTextObjectForText("Not selected")
            
            
        }
        
        return this._placeholderLocalizedTextObject
        
    }
    
    set placeholderLocalizedTextObject(placeholderText: UILocalizedTextObject) {
        
        this._placeholderLocalizedTextObject = placeholderText
        
        this.updateTitleWithSelection(this.selectedData)
        
    }
    
    
    updateTitleWithSelection(selectedRows: CBDropdownDataItem<T>[]) {
        
        
        this.titleLabel.localizedTextObject = this.placeholderLocalizedTextObject
        
        if (selectedRows && selectedRows.length) {
    
            const maxWidth = this.titleLabel.bounds.width
            this.titleLabel.localizedTextObject = nil
            this.titleLabel.text = ""
            var stopLooping = NO
    
            selectedRows.forEach(function (this: SearchableDropdown<T>, selectedDatapoint, index, array) {
                
                if (stopLooping) {
                    return
                }
    
                var selectedString = LanguageService.stringForCurrentLanguage(selectedDatapoint.title)
                if (index) {
                    selectedString = ", " + selectedString
                }
                const previousText = this.titleLabel.text
                this.titleLabel.text = this.titleLabel.text + selectedString
                
                this.overflowLabel.text = "+" + (array.length - index - 1)
                if (index == array.length - 1) {
                    this.overflowLabel.text = ""
                }
                
                if (index && (this.bounds.width -
                    (this.overflowLabel.intrinsicContentWidth() + this.titleLabel.intrinsicContentWidth()) - 20) < 0) {
                    
                    this.titleLabel.text = previousText
                    this.overflowLabel.text = "+" + (array.length - index - 2 * 0)
                    
                    stopLooping = YES
                }
                
            }, this)
            
        }
        
    }
    
    
    
    
    
    updateFilteredData(filteringString: string) {
        
        this._filteredData = []
        
        this.data.forEach(function (
            this: SearchableDropdown<T>,
            sectionOrRow: CBDropdownDataItem<T>,
            index: number,
            array: CBDropdownDataItem<T>[]
        ) {
            
            if (LanguageService.stringForCurrentLanguage(sectionOrRow.title).toLowerCase().contains(filteringString.toLowerCase())) {
                
                this.filteredData.push(sectionOrRow)
                
            }
            else if (sectionOrRow.isADropdownDataSection) {
                
                this._keyValueStringFilter.filterData(
                    filteringString,
                    sectionOrRow.rowsData,
                    this._excludedData,
                    "title." + LanguageService.currentLanguageKey,
                    sectionOrRow,
                    function (
                        this: SearchableDropdown<T>,
                        filteredData,
                        filteredIndexes,
                        sectionFromThread: CBDropdownDataItem<T>
                    ) {
                        
                        if (filteredData.length) {
                            
                            this.filteredData.push({
                                
                                _id: sectionFromThread._id,
                                title: sectionFromThread.title,
                                rowsData: filteredData,
                                isADropdownDataSection: sectionFromThread.isADropdownDataSection,
                                isADropdownDataRow: sectionFromThread.isADropdownDataRow,
                                
                                attachedObject: sectionFromThread.attachedObject,
                                
                                itemCode: sectionFromThread.itemCode,
                                dropdownCode: sectionFromThread.dropdownCode
                                
                            })
                            
                            if (this.allowsCustomItem && this._searchTextField.text && this._customItem) {
                                
                                this.filteredData.removeElement(this._customItem)
                                
                                this.filteredData.push(this._customItem)
                                
                            }
                            
                            this._isDrawingDataValid = NO
                            
                            this._tableView.reloadData()
                            
                        }
                        
                    }.bind(this)
                )
                
            }
            
        }.bind(this))
        
        
        if (this.allowsCustomItem && this._searchTextField.text) {
            
            this.filteredData.removeElement(this._customItem)
            
            this.initCustomItemWithTitle(this._searchTextField.text)
            
            this.filteredData.push(this._customItem)
            
        }
        
        
        if (this.filteredData.length) {
            
            this.focusedRowIndex = 0
            
        }
        else {
            
            this.focusedRowIndex = null
            
        }
        
        this._isDrawingDataValid = NO
        
        this._tableView.reloadData()
        
    }
    
    
    
    
    
    initCustomItemWithTitle(title: string) {
        
        if (IS_NOT(title)) {
            
            this._customItem = undefined
            
        }
        else {
            
            this._customItem = {
                
                _id: "" + MAKE_ID(),
                title: LanguageService.localizedTextObjectForText(title),
                rowsData: [],
                isADropdownDataSection: NO,
                isADropdownDataRow: YES,
                
                attachedObject: undefined,
                
                itemCode: "custom_item",
                dropdownCode: null
                
            }
            
        }
        
    }
    
    
    selectItemOrCustomItemWithTitle(title: string) {
        
        if (IS_NOT(title)) {
            
            this._customItem = undefined
            
            
            
        }
        else {
    
            var item = this.drawingData.find(function (dataItem, index, array) {
        
                return (LanguageService.stringForCurrentLanguage(dataItem.title) == title)
        
            })
    
    
            if (this.allowsCustomItem && IS_NOT(item)) {
                
                this._searchTextField.text = title
                
                this.updateFilteredData(title)
                
                item = this._customItem
                
            }
            
            
            if (IS_NOT(this.isSingleSelection)) {
                
                if (IS_NOT(this.selectedDataContains(item))) {
    
                    const selectedItemCodes = this.selectedItemCodes.copy()
    
                    selectedItemCodes.push(item.itemCode)
                    
                    this.selectedItemCodes = selectedItemCodes
                    
                    
                }
                
                return
                
            }
            
            this.selectedItemCodes = [item.itemCode]
            
            
        }
        
    }
    
    
    
    
    
    set data(data: CBDropdownDataItem<T>[]) {
        
        this._data = data
        
        this.updateFilteredData(this._searchTextField.text)
        
    }
    
    get data() {
        
        return this._data
        
    }
    
    set filteredData(data: CBDropdownDataItem<T>[]) {
        
        this._filteredData = data
        
        this._isDrawingDataValid = NO
        
    }
    
    get filteredData() {
        return this._filteredData
    }
    
    get drawingData(): CBDropdownDataItem<T>[] {
        
        if (this._isDrawingDataValid) {
            return this._drawingData
        }
    
        const result: CBDropdownDataItem<T>[] = []
    
        this._filteredData.forEach(function (this: SearchableDropdown<T>, section: CBDropdownDataItem<T>) {
            
            result.push({
                
                _id: section._id,
                title: section.title,
                
                // @ts-ignore
                rowsData: null,
                isADropdownDataSection: section.isADropdownDataSection,
                isADropdownDataRow: NO,
                
                attachedObject: section.attachedObject,
                
                itemCode: section.itemCode,
                dropdownCode: section.dropdownCode
                
                
            })
            
            if (section.rowsData) {
                
                section.rowsData.forEach(function (this: SearchableDropdown<T>, row) {
                    result.push(row)
                }.bind(this))
                
            }
            
        }.bind(this))
        
        
        this._drawingData = result
        this._isDrawingDataValid = YES
        
        return result
        
    }
    
    get selectedData() {
        
        return this._selectedData
        
    }
    
    set selectedData(selectedData: CBDropdownDataItem<T>[]) {
        
        this._selectedData = selectedData
        
    }
    
    
    
    
    
    clearSelection() {
        
        this.selectedData = []
        this.selectedIndices = []
        
        this.updateTitleWithSelection(this.selectedData)
        this.updateContentForCurrentSelection()
        this.selectionDidChange(this.selectedData)
        
    }
    
    
    
    
    
    get selectedItemCodes() {
        return this.selectedData.map(function (item) {
            return item.itemCode
        })
    }
    
    set selectedItemCodes(selectedItemCodes: string[]) {
    
        const selectedData = []
        const selectedIndices = []
    
        this._drawingData.forEach(function (item, index, array) {
            
            if (selectedItemCodes.contains(item.itemCode)) {
                
                selectedData.push(item)
                selectedIndices.push(index)
                
            }
            
        })
        
        
        this.selectedData = selectedData
        this.selectedIndices = selectedIndices
        
        this.updateTitleWithSelection(this.selectedData)
        this.updateContentForCurrentSelection()
        this.selectionDidChange(this.selectedData)
        
        
    }
    
    
    
    
    
    get selectedRowIdentifiers(): string[] {
        const result = []
        this.selectedData.forEach(function (this: SearchableDropdown<T>, selectedDatapoint: CBDropdownDataItem<T>) {
            result.push(selectedDatapoint._id)
        }.bind(this))
        return result
    }
    
    
    
    wasAddedToViewTree() {
        
        super.wasAddedToViewTree()
        
        this.setNeedsLayout()
        
        
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        super.layoutSubviews()
    
        const bounds = this.bounds
    
        const padding = RootViewController.paddingLength
        const labelHeight = padding
    
        this.updateTitleWithSelection(this.selectedData)
        
        if (this._rightImageView) {
            
            // var imageHeight = bounds.height - padding;
            // this._rightImageView.frame = new UIRectangle(bounds.width - imageHeight - padding * 0.5, padding * 0.5, imageHeight, imageHeight);
            
            this._rightImageView.frame = bounds.rectangleWithInsets(0, padding *
                0.5, 0, 0).rectangleWithWidth(24, 1).rectangleWithHeight(24, 0.5)
            
        }
        
        if (this.overflowLabel) {
            
            //this.overflowLabel.frame = bounds.rectangleWithInsets(padding, padding, 0, 0).rectangleWithInsets(0, this._rightImageView.bounds.width || 0, 0, 0).rectangleWithWidth(36, 1);
            
            //this.titleLabel.style.maxWidth = "" + (bounds.width - (bounds.width - this.titleLabel.viewHTMLElement.offsetLeft - (this.overflowLabel.frame.x || 0)))
            
            this.overflowLabel.style.width = "36px"
            this.overflowLabel.style.right = "32px"
            
            this.overflowLabel.centerYInContainer()
            
            
            this.titleLabel.style.marginRight = "60px"
            
            //this.overflowLabel.centerYInContainer();
            //this.titleLabel.frame = this.titleLabel.frame.rectangleWithWidth(bounds.width - (bounds.width - this.titleLabel.frame.x - (this.overflowLabel.frame.x || 0)))
            
        }
        
        
        
        
        
    }
    
    
    
    
    
}




























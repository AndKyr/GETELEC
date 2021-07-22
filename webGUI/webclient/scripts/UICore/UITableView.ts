/// <reference path="./UINativeScrollView.ts" />



interface UITableViewReusableViewsContainerObject {
    
    [key: string]: UIView[];
    
}


interface UITableViewReusableViewPositionObject {
    
    bottomY: number;
    topY: number;
    
    isValid: boolean;
    
}





class UITableView extends UINativeScrollView {
    
    
    
    
    allRowsHaveEqualHeight: boolean = NO
    _visibleRows: UIView[] = []
    _firstLayoutVisibleRows: UIView[] = []
    
    _rowPositions: UITableViewReusableViewPositionObject[] = []
    
    _highestValidRowPositionIndex: number = 0
    
    _reusableViews: UITableViewReusableViewsContainerObject = {}
    
    _removedReusableViews: UITableViewReusableViewsContainerObject = {}
    
    _fullHeightView: UIView
    
    _rowIDIndex: number = 0
    
    reloadsOnLanguageChange = YES
    
    sidePadding = 0
    
    _persistedData: any[] = []
    _needsDrawingOfVisibleRowsBeforeLayout = NO
    _isDrawVisibleRowsScheduled = NO
    _shouldAnimateNextLayout: boolean
    
    animationDuration = 0.25
    
    
    constructor(elementID) {
        
        super(elementID)
        
        this._class = UITableView
        this.superclass = UINativeScrollView
        
        this.scrollsX = NO
        
    }
    
    
    
    initView(elementID, viewHTMLElement) {
        
        super.initView(elementID, viewHTMLElement)
        
        this._fullHeightView = new UIView()
        this._fullHeightView.hidden = YES
        this._fullHeightView.userInteractionEnabled = NO
        this.addSubview(this._fullHeightView)
        
    }
    
    
    
    loadData() {
        
        this._persistedData = []
        
        this._calculatePositionsUntilIndex(this.numberOfRows() - 1)
        this._needsDrawingOfVisibleRowsBeforeLayout = YES
        
        this.setNeedsLayout()
        
        // this.forEachViewInSubtree(function(view) {
        
        //     view.setNeedsLayout();
        
        // })
        
    }
    
    reloadData() {
        
        this._removeVisibleRows()
        this._removeAllReusableRows()
        
        this._rowPositions = []
        this._highestValidRowPositionIndex = 0
        
        this.loadData()
        
    }
    
    
    
    highlightChanges(previousData: any[], newData: any[]) {
        
        previousData = previousData.map(function (dataPoint, index, array) {
            
            return JSON.stringify(dataPoint)
            
        })
        
        newData = newData.map(function (dataPoint, index, array) {
            
            return JSON.stringify(dataPoint)
            
        })
        
        
        const newIndexes: number[] = []
        
        newData.forEach(function (value, index, array) {
            
            if (!previousData.contains(value)) {
                
                newIndexes.push(index)
                
            }
            
        })
        
        newIndexes.forEach(function (this: UITableView, index) {
            
            if (this.isRowWithIndexVisible(index)) {
                
                this.highlightRowAsNew(this.viewForRowWithIndex(index))
                
            }
            
        }.bind(this))
        
    }
    
    
    
    highlightRowAsNew(row: UIView) {
        
        
        
    }
    
    
    
    
    
    invalidateSizeOfRowWithIndex(index: number, animateChange = NO) {
        
        if (this._rowPositions[index]) {
            this._rowPositions[index].isValid = NO
        }
        
        this._highestValidRowPositionIndex = Math.min(this._highestValidRowPositionIndex, index - 1)
        
        // if (index == 0) {
        
        //     this._highestValidRowPositionIndex = 0;
        
        //     this._rowPositions = [];
        
        // }
        
        this._needsDrawingOfVisibleRowsBeforeLayout = YES
        
        this._shouldAnimateNextLayout = animateChange
        
    }
    
    
    
    _calculateAllPositions() {
        this._calculatePositionsUntilIndex(this.numberOfRows() - 1)
    }
    
    _calculatePositionsUntilIndex(maxIndex: number) {
        
        var validPositionObject = this._rowPositions[this._highestValidRowPositionIndex]
        if (!IS(validPositionObject)) {
            validPositionObject = {
                bottomY: 0,
                topY: 0,
                isValid: YES
            }
        }
        
        var previousBottomY = validPositionObject.bottomY
        
        if (!this._rowPositions.length) {
            
            this._highestValidRowPositionIndex = -1
            
        }
        
        for (var i = this._highestValidRowPositionIndex + 1; i <= maxIndex; i++) {
            
            var height: number
            
            const rowPositionObject = this._rowPositions[i]
            
            if (IS((rowPositionObject || nil).isValid)) {
                
                height = rowPositionObject.bottomY - rowPositionObject.topY
                
            }
            else {
                
                height = this.heightForRowWithIndex(i)
                
            }
            
            
            
            
            const positionObject: UITableViewReusableViewPositionObject = {
                bottomY: previousBottomY + height,
                topY: previousBottomY,
                isValid: YES
            }
            
            if (i < this._rowPositions.length) {
                this._rowPositions[i] = positionObject
            }
            else {
                this._rowPositions.push(positionObject)
            }
            this._highestValidRowPositionIndex = i
            previousBottomY = previousBottomY + height
            
        }
        
    }
    
    
    
    indexesForVisibleRows(paddingRatio = 0.5): number[] {
        
        const firstVisibleY = this.contentOffset.y - this.bounds.height * paddingRatio
        const lastVisibleY = firstVisibleY + this.bounds.height * (1 + paddingRatio)
        
        const numberOfRows = this.numberOfRows()
        
        if (this.allRowsHaveEqualHeight) {
            
            const rowHeight = this.heightForRowWithIndex(0)
            
            var firstIndex = firstVisibleY / rowHeight
            var lastIndex = lastVisibleY / rowHeight
            
            firstIndex = Math.trunc(firstIndex)
            lastIndex = Math.trunc(lastIndex) + 1
            
            firstIndex = Math.max(firstIndex, 0)
            lastIndex = Math.min(lastIndex, numberOfRows - 1)
            
            var result = []
            for (var i = firstIndex; i < lastIndex + 1; i++) {
                result.push(i)
            }
            return result
        }
        
        var accumulatedHeight = 0
        var result = []
        
        this._calculateAllPositions()
        
        const rowPositions = this._rowPositions
        
        for (var i = 0; i < numberOfRows; i++) {
            
            const height = rowPositions[i].bottomY - rowPositions[i].topY  // this.heightForRowWithIndex(i)
            
            accumulatedHeight = accumulatedHeight + height
            if (accumulatedHeight >= firstVisibleY) {
                result.push(i)
            }
            if (accumulatedHeight >= lastVisibleY) {
                break
            }
            
        }
        
        return result
        
    }
    
    
    _removeVisibleRows() {
        
        const visibleRows = []
        this._visibleRows.forEach(function (this: UITableView, row: UIView, index: number, array: UIView[]) {
            
            this._persistedData[row._UITableViewRowIndex] = this.persistenceDataItemForRowWithIndex(
                row._UITableViewRowIndex,
                row
            )
            row.removeFromSuperview()
            this._removedReusableViews[row._UITableViewReusabilityIdentifier].push(row)
            
            
        }, this)
        this._visibleRows = visibleRows
        
    }
    
    
    
    _removeAllReusableRows() {
        this._reusableViews.forEach(function (this: UITableView, rows: UIView[]) {
            rows.forEach(function (this: UITableView, row: UIView, index: number, array: UIView[]) {
                
                this._persistedData[row._UITableViewRowIndex] = this.persistenceDataItemForRowWithIndex(
                    row._UITableViewRowIndex,
                    row
                )
                row.removeFromSuperview()
                
                this._markReusableViewAsUnused(row)
                
            }.bind(this))
        }.bind(this))
    }
    
    
    
    _markReusableViewAsUnused(row: UIView) {
        if (!this._removedReusableViews[row._UITableViewReusabilityIdentifier].contains(row)) {
            this._removedReusableViews[row._UITableViewReusabilityIdentifier].push(row)
        }
    }
    
    _drawVisibleRows() {
        
        if (!this.isMemberOfViewTree) {
            return
        }
        
        const visibleIndexes = this.indexesForVisibleRows()
        
        const minIndex = visibleIndexes[0]
        const maxIndex = visibleIndexes[visibleIndexes.length - 1]
        
        const removedViews = []
        
        const visibleRows = []
        this._visibleRows.forEach(function (this: UITableView, row: UIView, index: number, array: UIView[]) {
            if (row._UITableViewRowIndex < minIndex || row._UITableViewRowIndex > maxIndex) {
                
                //row.removeFromSuperview();
                
                this._persistedData[row._UITableViewRowIndex] = this.persistenceDataItemForRowWithIndex(
                    row._UITableViewRowIndex,
                    row
                )
                
                this._removedReusableViews[row._UITableViewReusabilityIdentifier].push(row)
                
                removedViews.push(row)
                
            }
            else {
                visibleRows.push(row)
            }
        }, this)
        this._visibleRows = visibleRows
        
        visibleIndexes.forEach(function (this: UITableView, rowIndex: number, index: number, array: number) {
            
            
            if (this.isRowWithIndexVisible(rowIndex)) {
                return
            }
            const view: UIView = this.viewForRowWithIndex(rowIndex)
            //view._UITableViewRowIndex = rowIndex;
            this._firstLayoutVisibleRows.push(view)
            this._visibleRows.push(view)
            this.addSubview(view)
            
        } as any, this)
        
        for (var i = 0; i < removedViews.length; i++) {
            
            var view: UIView = removedViews[i]
            if (this._visibleRows.indexOf(view) == -1) {
                
                //this._persistedData[view._UITableViewRowIndex] = this.persistenceDataItemForRowWithIndex(view._UITableViewRowIndex, view);
                view.removeFromSuperview()
                
                //this._removedReusableViews[view._UITableViewReusabilityIdentifier].push(view);
                
            }
            
        }
        
        //this.setNeedsLayout();
        
    }
    
    
    visibleRowWithIndex(rowIndex: number): UIView {
        for (var i = 0; i < this._visibleRows.length; i++) {
            const row = this._visibleRows[i]
            if (row._UITableViewRowIndex == rowIndex) {
                return row
            }
        }
        return nil
    }
    
    
    isRowWithIndexVisible(rowIndex: number) {
        return IS(this.visibleRowWithIndex(rowIndex))
    }
    
    
    
    reusableViewForIdentifier(identifier: string, rowIndex: number) {
        
        if (!this._removedReusableViews[identifier]) {
            this._removedReusableViews[identifier] = []
        }
        
        if (this._removedReusableViews[identifier] && this._removedReusableViews[identifier].length) {
            
            const view = this._removedReusableViews[identifier].pop()
            
            view._UITableViewRowIndex = rowIndex
            
            Object.assign(view, this._persistedData[rowIndex] || this.defaultRowPersistenceDataItem())
            
            return view
            
        }
        
        if (!this._reusableViews[identifier]) {
            this._reusableViews[identifier] = []
        }
        
        const newView = this.newReusableViewForIdentifier(identifier, this._rowIDIndex)
        this._rowIDIndex = this._rowIDIndex + 1
        
        if (this._rowIDIndex > 40) {
            
            const asd = 1
            
        }
        
        newView._UITableViewReusabilityIdentifier = identifier
        
        newView._UITableViewRowIndex = rowIndex
        
        Object.assign(newView, this._persistedData[rowIndex] || this.defaultRowPersistenceDataItem())
        this._reusableViews[identifier].push(newView)
        
        return newView
        
    }
    
    
    // Functions that should be overridden to draw the correct content START
    newReusableViewForIdentifier(identifier: string, rowIDIndex: number): UIView {
        
        const view = new UIButton(this.elementID + "Row" + rowIDIndex)
        
        view.stopsPointerEventPropagation = NO
        view.pausesPointerEvents = NO
        
        return view
        
    }
    
    heightForRowWithIndex(index: number): number {
        return 50
    }
    
    numberOfRows() {
        return 10000
    }
    
    defaultRowPersistenceDataItem() {
        
        
        
    }
    
    persistenceDataItemForRowWithIndex(rowIndex: number, row: UIView) {
        
        
        
    }
    
    viewForRowWithIndex(rowIndex: number): UIView {
        
        const row = this.reusableViewForIdentifier("Row", rowIndex);
        (row as UIButton).titleLabel.text = "Row " + rowIndex
        
        return row
        
    }
    
    // Functions that should be overridden to draw the correct content END
    
    
    // Functions that trigger redrawing of the content
    didScrollToPosition(offsetPosition: UIPoint) {
        
        super.didScrollToPosition(offsetPosition)
        
        this.forEachViewInSubtree(function (view: UIView) {
            
            view._isPointerValid = NO
            
        })
        
        if (!this._isDrawVisibleRowsScheduled) {
            
            this._isDrawVisibleRowsScheduled = YES
            
            UIView.runFunctionBeforeNextFrame(function (this: UITableView) {
                
                this._calculateAllPositions()
                
                this._drawVisibleRows()
                
                this.setNeedsLayout()
                
                this._isDrawVisibleRowsScheduled = NO
                
            }.bind(this))
            
        }
        
    }
    
    wasAddedToViewTree() {
        this.loadData()
    }
    
    setFrame(rectangle: UIRectangle, zIndex?: number, performUncheckedLayout?: boolean) {
        
        const frame = this.frame
        super.setFrame(rectangle, zIndex, performUncheckedLayout)
        if (frame.isEqualTo(rectangle) && !performUncheckedLayout) {
            return
        }
        
        this._needsDrawingOfVisibleRowsBeforeLayout = YES
        
    }
    
    
    
    
    
    didReceiveBroadcastEvent(event: UIViewBroadcastEvent) {
        
        super.didReceiveBroadcastEvent(event)
        
        if (event.name == UIView.broadcastEventName.LanguageChanged && this.reloadsOnLanguageChange) {
            
            this.reloadData()
            
        }
        
        
    }
    
    
    
    
    
    private _layoutAllRows(positions = this._rowPositions) {
        
        const bounds = this.bounds
        
        this._visibleRows.forEach(function (this: UITableView, row: UIView, index: number, array: UIView[]) {
            
            const frame = bounds.copy()
            
            const positionObject = positions[row._UITableViewRowIndex]
            frame.min.y = positionObject.topY
            frame.max.y = positionObject.bottomY
            row.frame = frame
            
            row.style.width = "" + (bounds.width - this.sidePadding * 2).integerValue + "px"
            row.style.left = "" + this.sidePadding.integerValue + "px"
            
            
        }, this)
        
        this._fullHeightView.frame = bounds.rectangleWithHeight((positions.lastElement ||
            nil).bottomY).rectangleWithWidth(bounds.width * 0.5)
        
        this._firstLayoutVisibleRows = []
        
        
    }
    
    private _animateLayoutAllRows() {
        
        UIView.animateViewOrViewsWithDurationDelayAndFunction(
            this._visibleRows,
            this.animationDuration,
            0,
            undefined,
            function (this: UITableView) {
                
                this._layoutAllRows()
                
            }.bind(this),
            function (this: UITableView) {
                
                // this._calculateAllPositions()
                // this._layoutAllRows()
                
            }.bind(this)
        )
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        const previousPositions: UITableViewReusableViewPositionObject[] = JSON.parse(JSON.stringify(this._rowPositions))
        
        const previousVisibleRowsLength = this._visibleRows.length
        
        if (this._needsDrawingOfVisibleRowsBeforeLayout) {
            
            //this._calculateAllPositions()
            
            this._drawVisibleRows()
            
            this._needsDrawingOfVisibleRowsBeforeLayout = NO
            
        }
        
        
        
        super.layoutSubviews()
        
        
        if (!this.numberOfRows() || !this.isMemberOfViewTree) {
            
            return
            
        }
        
        
        
        
        
        if (this._shouldAnimateNextLayout) {
            
            
            
            // Need to do layout with the previous positions
            
            this._layoutAllRows(previousPositions)
            
            
            if (previousVisibleRowsLength < this._visibleRows.length) {
                
                
                UIView.runFunctionBeforeNextFrame(function (this: UITableView) {
                    
                    this._animateLayoutAllRows()
                    
                }.bind(this))
                
            }
            else {
                
                this._animateLayoutAllRows()
                
            }
            
            
            this._shouldAnimateNextLayout = NO
            
        }
        else {
            
            // if (this._needsDrawingOfVisibleRowsBeforeLayout) {
            
            //     this._drawVisibleRows();
            
            //     this._needsDrawingOfVisibleRowsBeforeLayout = NO;
            
            // }
            
            this._calculateAllPositions()
            
            this._layoutAllRows()
            
            
        }
        
        
        
        
        
    }
    
    
    
    
    
    intrinsicContentHeight(constrainingWidth = 0) {
        
        
        
        var result = 10000000000000000000000000
        
        if (this._rowPositions.length) {
            
            result = this._rowPositions[this._rowPositions.length - 1].bottomY
            
        }
        
        return result
        
    }
    
    
    
    
    
}

































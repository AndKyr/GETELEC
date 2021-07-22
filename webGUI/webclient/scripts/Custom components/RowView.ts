/// <reference path="../UICore/UIView.ts" />





class RowView<CellType extends UIView = UIView> extends UIView {
    
    
    _previousLayoutBounds: UIRectangle
    _cells: CellType[]
    _cellWeights: number[]
    
    padding = 0
    _cellWidths: number[]
    _rowHeight = 50
    
    
    constructor(elementID?: string, cells: CellType[] = [], cellWidths: number[] = []) {
        
        super(elementID)
        
        this._class = RowView
        this.superclass = UIView
        
        this._cells = cells
        this._cellWeights = cellWidths
        
    }
    
    
    
    get cells() {
        return this._cells
    }
    
    
    set cells(cells: CellType[]) {
        
        const previousCells = this.cells
        
        // The cells are added to this._cells in this.addCell()
        
        previousCells.copy().forEach(function (this: RowView, cell: CellType, index: number, array: CellType[]) {
            
            if (!cells.contains(cell)) {
                
                cell.removeFromSuperview()
                this._cells.removeElement(cell)
                
            }
            
        }.bind(this))
        
        cells.copy().forEach(function (this: RowView, cell: CellType, index: number, array: CellType[]) {
            if (!IS(cell.superview)) {
                this.addCell(cell, 1, index)
            }
        }, this)
        
        this._previousLayoutBounds = nil
        this.setNeedsLayout()
        
    }
    
    
    
    removeCellAtIndex(index: number) {
        const remainingCells = this.cells.copy()
        remainingCells.removeElementAtIndex(index)
        this.cellWeights.removeElementAtIndex(index)
        this.cells = remainingCells
    }
    
    removeFirstCell() {
        this.removeCellAtIndex(0)
    }
    
    removeLastCell() {
        this.removeCellAtIndex(this.cells.length - 1)
    }
    
    
    
    addCell(cell: CellType, weight: number = 1, index = this.cells.length) {
        
        if (this.cells.contains(cell) && IS_NOT_NIL(cell)) {
            return
        }
        
        this.cells.insertElementAtIndex(index, cell)
        this.cellWeights.insertElementAtIndex(index, weight)
        this.addSubview(cell)
        this.setNeedsLayout()
        
    }
    
    
    get cellWeights() {
        return this._cellWeights
    }
    
    set cellWeights(widths: number[]) {
        this._cellWeights = widths
        this._previousLayoutBounds = nil
        this.setNeedsLayout()
    }
    
    get cellWidths() {
        return this._cellWidths
    }
    
    set cellWidths(widths: number[]) {
        this._cellWidths = widths
        this._previousLayoutBounds = nil
        this.setNeedsLayout()
    }
    
    
    
    get rowHeight() {
        
        var result = IF(this._rowHeight)(() => this._rowHeight)
        .ELSE(() => this.cells.map((value, index, array) => value.intrinsicContentHeight(this.bounds.width)).max())
        
        
        
        
        return result //this.bounds.height
        
    }
    
    
    
    
    
    layoutSubviews() {
        
        const bounds = this.bounds
        
        if (bounds.isEqualTo(this._previousLayoutBounds)) {
            
            return
            
        }
        
        super.layoutSubviews()
        
        
        this._previousLayoutBounds = bounds
        
        bounds.distributeViewsAlongWidth(this._cells, this._cellWeights, this.padding, this._cellWidths)
        
        this.cells.forEach(function (this: RowView, cell: UIView, index: number, array: UIView[]) {
            
            cell.frame = cell.frame.rectangleWithHeight(this.rowHeight)
            
        }.bind(this))
        
        
        
        
    }
    
    
    
    
    
}




























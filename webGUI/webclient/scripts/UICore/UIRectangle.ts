/// <reference path="./UIPoint.ts" />





class UIRectangle extends UIObject {
    _isBeingUpdated: boolean
    rectanglePointDidChange: (b: any) => void
    max: UIPoint
    min: UIPoint
    
    
    constructor(x: number = 0, y: number = 0, height: number = 0, width: number = 0) {
        
        super()
        
        this._class = UIRectangle
        this.superclass = UIObject
        
        
        this.min = new UIPoint(Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY)
        this.max = new UIPoint(Number.NEGATIVE_INFINITY, Number.NEGATIVE_INFINITY)
        
        this.min.didChange = this.rectanglePointDidChange
        this.max.didChange = this.rectanglePointDidChange
        
        this._isBeingUpdated = NO
        
        this.min = new UIPoint(x, y)
        this.max = new UIPoint(x + width, y + height)
        
        if (IS_NIL(height)) {
            this.max.y = height
        }
        
        if (IS_NIL(width)) {
            this.max.x = width
        }
        
        
    }
    
    
    
    
    
    copy() {
        const result = new UIRectangle(this.x, this.y, this.height, this.width)
        return result
    }
    
    isEqualTo(rectangle: UIRectangle) {
    
    
        const result = (IS(rectangle) && this.min.isEqualTo(rectangle.min) && this.max.isEqualTo(rectangle.max))
    
        return result
        
    }
    
    static zero() {
    
        const result = new UIRectangle(0, 0, 0, 0)
    
        return result
        
    }
    
    containsPoint(point: UIPoint) {
        return this.min.x <= point.x && this.min.y <= point.y &&
            point.x <= this.max.x && point.y <= this.max.y
    }
    
    updateByAddingPoint(point: UIPoint) {
        
        if (!point) {
            point = new UIPoint(0, 0)
        }
        
        this.beginUpdates()
    
        const min = this.min.copy()
        if (min.x === nil) {
            min.x = this.max.x
        }
        if (min.y === nil) {
            min.y = this.max.y
        }
    
        const max = this.max.copy()
        if (max.x === nil) {
            max.x = this.min.x
        }
        if (max.y === nil) {
            max.y = this.min.y
        }
        
        this.min.x = Math.min(min.x, point.x)
        this.min.y = Math.min(min.y, point.y)
        this.max.x = Math.max(max.x, point.x)
        this.max.y = Math.max(max.y, point.y)
        
        this.finishUpdates()
        
    }
    
    get height() {
        if (this.max.y === nil) {
            return nil
        }
        return this.max.y - this.min.y
    }
    
    set height(height: number) {
        this.max.y = this.min.y + height
    }
    
    
    
    get width() {
        if (this.max.x === nil) {
            return nil
        }
        return this.max.x - this.min.x
    }
    
    set width(width: number) {
        this.max.x = this.min.x + width
    }
    
    
    
    get x() {
        return this.min.x
    }
    
    set x(x: number) {
        
        this.beginUpdates()
    
        const width = this.width
        this.min.x = x
        this.max.x = this.min.x + width
        
        this.finishUpdates()
        
    }
    
    
    get y() {
        return this.min.y
    }
    
    
    set y(y: number) {
        
        this.beginUpdates()
    
        const height = this.height
        this.min.y = y
        this.max.y = this.min.y + height
        
        this.finishUpdates()
        
    }
    
    
    
    
    get topLeft() {
        
        return this.min.copy()
        
    }
    
    get topRight() {
        return new UIPoint(this.max.x, this.y)
    }
    
    get bottomLeft() {
        return new UIPoint(this.x, this.max.y)
    }
    
    get bottomRight() {
        
        return this.max.copy()
        
    }
    
    
    get center() {
    
        const result = this.min.copy().add(this.min.to(this.max).scale(0.5))
    
        return result
        
    }
    
    set center(center: UIPoint) {
    
        const offset = this.center.to(center)
        this.offsetByPoint(offset)
        
    }
    
    offsetByPoint(offset: UIPoint) {
        
        this.min.add(offset)
        this.max.add(offset)
        
        return this
        
    }
    
    
    
    concatenateWithRectangle(rectangle: UIRectangle) {
        
        this.updateByAddingPoint(rectangle.bottomRight)
        this.updateByAddingPoint(rectangle.topLeft)
        
        return this
        
    }
    
    
    
    
    
    intersectionRectangleWithRectangle(rectangle: UIRectangle): UIRectangle {
    
    
        const result = this.copy()
    
        result.beginUpdates()
    
        const min = result.min
        if (min.x === nil) {
            min.x = rectangle.max.x - Math.min(result.width, rectangle.width)
        }
        if (min.y === nil) {
            min.y = rectangle.max.y - Math.min(result.height, rectangle.height)
        }
    
        const max = result.max
        if (max.x === nil) {
            max.x = rectangle.min.x + Math.min(result.width, rectangle.width)
        }
        if (max.y === nil) {
            max.y = rectangle.min.y + Math.min(result.height, rectangle.height)
        }
        
        result.min.x = Math.max(result.min.x, rectangle.min.x)
        result.min.y = Math.max(result.min.y, rectangle.min.y)
        result.max.x = Math.min(result.max.x, rectangle.max.x)
        result.max.y = Math.min(result.max.y, rectangle.max.y)
        
        
        if (result.height < 0) {
    
            const averageY = (this.center.y + rectangle.center.y) * 0.5
            result.min.y = averageY
            result.max.y = averageY
            
        }
        
        if (result.width < 0) {
    
            const averageX = (this.center.x + rectangle.center.x) * 0.5
            result.min.x = averageX
            result.max.x = averageX
            
        }
        
        result.finishUpdates()
        
        return result
        
    }
    
    
    
    get area() {
        const result = this.height * this.width
        return result
    }
    
    
    intersectsWithRectangle(rectangle: UIRectangle) {
        
        return (this.intersectionRectangleWithRectangle(rectangle).area != 0)
        
    }
    
    
    
    
    
    // add some space around the rectangle
    rectangleWithInsets(left: number, right: number, bottom: number, top: number) {
        const result = this.copy()
        result.min.x = this.min.x + left
        result.max.x = this.max.x - right
        result.min.y = this.min.y + top
        result.max.y = this.max.y - bottom
        return result
    }
    
    rectangleWithInset(inset: number) {
        const result = this.rectangleWithInsets(inset, inset, inset, inset)
        return result
    }
    
    rectangleWithHeight(height: number, centeredOnPosition: number = nil) {
        
        if (isNaN(centeredOnPosition)) {
            centeredOnPosition = nil
        }
    
        const result = this.copy()
        result.height = height
        
        if (centeredOnPosition != nil) {
            const change = height - this.height
            result.offsetByPoint(new UIPoint(0, change * centeredOnPosition).scale(-1))
        }
        
        return result
        
    }
    
    rectangleWithWidth(width: number, centeredOnPosition: number = nil) {
        
        if (isNaN(centeredOnPosition)) {
            centeredOnPosition = nil
        }
    
        const result = this.copy()
        result.width = width
        
        if (centeredOnPosition != nil) {
            const change = width - this.width
            result.offsetByPoint(new UIPoint(change * centeredOnPosition, 0).scale(-1))
        }
        
        return result
        
    }
    
    rectangleWithHeightRelativeToWidth(heightRatio: number = 1, centeredOnPosition: number = nil) {
    
        const result = this.rectangleWithHeight(this.width * heightRatio, centeredOnPosition)
    
        return result
        
    }
    
    rectangleWithWidthRelativeToHeight(widthRatio: number = 1, centeredOnPosition: number = nil) {
    
        const result = this.rectangleWithWidth(this.height * widthRatio, centeredOnPosition)
    
        return result
        
    }
    
    rectangleWithX(x: number, centeredOnPosition: number = 0) {
    
        const result = this.copy()
        result.x = x - result.width * centeredOnPosition
        
        return result
        
    }
    
    rectangleWithY(y: number, centeredOnPosition: number = 0) {
    
        const result = this.copy()
        result.y = y - result.height * centeredOnPosition
        
        return result
        
    }
    
    
    
    rectangleByAddingX(x: number) {
    
        const result = this.copy()
        result.x = this.x + x
        
        return result
        
    }
    
    rectangleByAddingY(y: number) {
    
        const result = this.copy()
        result.y = this.y + y
        
        return result
        
    }
    
    
    
    
    
    rectanglesBySplittingWidth(
        weights: number[],
        paddings: number | number[] = 0,
        absoluteWidths: number | number[] = nil
    ) {
        
        if (IS_NIL(paddings)) {
            
            paddings = 1
            
        }
        
        if (!(paddings instanceof Array)) {
            
            paddings = [paddings].arrayByRepeating(weights.length - 1)
            
        }
        
        paddings = paddings.arrayByTrimmingToLengthIfLonger(weights.length - 1)
        
        if (!(absoluteWidths instanceof Array) && IS_NOT_NIL(absoluteWidths)) {
            absoluteWidths = [absoluteWidths].arrayByRepeating(weights.length)
        }
    
        const result: UIRectangle[] = []
        const sumOfWeights = weights.reduce(function (a, b, index) {
            if (IS_NOT_NIL(absoluteWidths[index])) {
                b = 0
            }
            return a + b
        }, 0)
        const sumOfPaddings = paddings.summedValue
        const sumOfAbsoluteWidths = (absoluteWidths as number[]).summedValue
        const totalRelativeWidth = this.width - sumOfPaddings - sumOfAbsoluteWidths
        var previousCellMaxX = this.x
    
        for (var i = 0; i < weights.length; i++) {
    
            var resultWidth: number
            if (IS_NOT_NIL(absoluteWidths[i])) {
                
                resultWidth = absoluteWidths[i] || 0
                
            }
            else {
                
                resultWidth = totalRelativeWidth * (weights[i] / sumOfWeights)
                
            }
    
            const rectangle = this.rectangleWithWidth(resultWidth)
    
            var padding = 0
            if (paddings.length > i && paddings[i]) {
                padding = paddings[i]
            }
            
            rectangle.x = previousCellMaxX
            previousCellMaxX = rectangle.max.x + padding
            //rectangle = rectangle.rectangleWithInsets(0, padding, 0, 0);
            result.push(rectangle)
            
        }
        
        return result
        
    }
    
    rectanglesBySplittingHeight(
        weights: number[],
        paddings: number | number[] = 0,
        absoluteHeights: number | number[] = nil
    ) {
        
        if (IS_NIL(paddings)) {
            
            paddings = 1
            
        }
        
        if (!(paddings instanceof Array)) {
            
            paddings = [paddings].arrayByRepeating(weights.length - 1)
            
        }
        
        paddings = paddings.arrayByTrimmingToLengthIfLonger(weights.length - 1)
        
        if (!(absoluteHeights instanceof Array) && IS_NOT_NIL(absoluteHeights)) {
            absoluteHeights = [absoluteHeights].arrayByRepeating(weights.length)
        }
    
        const result: UIRectangle[] = []
        const sumOfWeights = weights.reduce(function (a, b, index) {
            if (IS_NOT_NIL(absoluteHeights[index])) {
                b = 0
            }
            return a + b
        }, 0)
        const sumOfPaddings = paddings.summedValue
        const sumOfAbsoluteHeights = (absoluteHeights as number[]).summedValue
        const totalRelativeHeight = this.height - sumOfPaddings - sumOfAbsoluteHeights
        var previousCellMaxY = this.y
    
        for (var i = 0; i < weights.length; i++) {
            var resultHeight: number
            if (IS_NOT_NIL(absoluteHeights[i])) {
                
                resultHeight = absoluteHeights[i] || 0
                
            }
            else {
                
                resultHeight = totalRelativeHeight * (weights[i] / sumOfWeights)
                
            }
    
            const rectangle = this.rectangleWithHeight(resultHeight)
    
            var padding = 0
            if (paddings.length > i && paddings[i]) {
                padding = paddings[i]
            }
            
            rectangle.y = previousCellMaxY
            previousCellMaxY = rectangle.max.y + padding
            //rectangle = rectangle.rectangleWithInsets(0, 0, padding, 0);
            result.push(rectangle)
        }
        
        return result
        
    }
    
    
    
    
    
    rectanglesByEquallySplittingWidth(numberOfFrames: number, padding: number = 0) {
        const result: UIRectangle[] = []
        const totalPadding = padding * (numberOfFrames - 1)
        const resultWidth = (this.width - totalPadding) / numberOfFrames
        for (var i = 0; i < numberOfFrames; i++) {
            const rectangle = this.rectangleWithWidth(resultWidth, i / (numberOfFrames - 1))
            result.push(rectangle)
        }
        return result
    }
    
    rectanglesByEquallySplittingHeight(numberOfFrames: number, padding: number = 0) {
        const result: UIRectangle[] = []
        const totalPadding = padding * (numberOfFrames - 1)
        const resultHeight = (this.height - totalPadding) / numberOfFrames
        for (var i = 0; i < numberOfFrames; i++) {
            const rectangle = this.rectangleWithHeight(resultHeight, i / (numberOfFrames - 1))
            result.push(rectangle)
        }
        return result
    }
    
    
    
    distributeViewsAlongWidth(
        views: UIView[],
        weights: number | number[] = 1,
        paddings?: number | number[],
        absoluteWidths?: number | number[]
    ) {
        
        if (!(weights instanceof Array)) {
            weights = [weights].arrayByRepeating(views.length)
        }
    
        const frames = this.rectanglesBySplittingWidth(weights, paddings, absoluteWidths)
    
        frames.forEach(function (frame, index, array) {
            views[index].frame = frame
        })
        
        return this
        
    }
    
    distributeViewsAlongHeight(
        views: UIView[],
        weights: number | number[] = 1,
        paddings?: number | number[],
        absoluteHeights?: number | number[]
    ) {
        
        if (!(weights instanceof Array)) {
            weights = [weights].arrayByRepeating(views.length)
        }
    
        const frames = this.rectanglesBySplittingHeight(weights, paddings, absoluteHeights)
    
        frames.forEach(function (frame, index, array) {
            views[index].frame = frame
        })
        
        return this
        
    }
    
    
    distributeViewsEquallyAlongWidth(views: UIView[], padding: number) {
    
        const frames = this.rectanglesByEquallySplittingWidth(views.length, padding)
    
        frames.forEach(function (frame, index, array) {
            views[index].frame = frame
        })
        
        return this
        
    }
    
    distributeViewsEquallyAlongHeight(views: UIView[], padding: number) {
    
        const frames = this.rectanglesByEquallySplittingHeight(views.length, padding)
    
        frames.forEach(function (frame, index, array) {
            views[index].frame = frame
        })
        
        return this
        
    }
    
    
    
    rectangleForNextRow(padding: number = 0, height = this.height) {
        const result = this.rectangleWithY(this.max.y + padding)
        if (height != this.height) {
            result.height = height
        }
        return result
    }
    
    rectangleForNextColumn(padding: number = 0, width = this.width) {
        const result = this.rectangleWithX(this.max.x + padding)
        if (width != this.width) {
            result.width = width
        }
        return result
    }
    
    rectangleForPreviousRow(padding: number = 0) {
        const result = this.rectangleWithY(this.min.y - this.height - padding)
        return result
    }
    
    rectangleForPreviousColumn(padding: number = 0) {
        const result = this.rectangleWithX(this.min.x - this.width - padding)
        return result
    }
    
    
    
    // Bounding box
    static boundingBoxForPoints(points) {
        const result = new UIRectangle()
        for (var i = 0; i < points.length; i++) {
            result.updateByAddingPoint(points[i])
        }
        return result
    }
    
    
    beginUpdates() {
        this._isBeingUpdated = YES
    }
    
    finishUpdates() {
        this._isBeingUpdated = NO
        this.didChange()
    }
    
    
    didChange() {
        
        // Callback to be set by delegate
        
    }
    
    _rectanglePointDidChange() {
        
        if (!this._isBeingUpdated) {
            
            this.didChange()
            
        }
        
    }
    
    
    
    
    
}
























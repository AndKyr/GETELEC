/// <reference path="./UICore/UIViewController.ts" />
///<reference path="UICore/UITextView.ts"/>





class GETELECTableViewController extends UIViewController {
    
    
    private titleLabel: UITextView
    private inputTextArea: UITextArea
    
    
    private loadDataButton: CBButton
    chartView: UIView
    chart: any
    private resultsLabel: UIView
    
    
    
    
    
    constructor(view) {
        
        // Calling super
        super(view)
        
        // Code for further setup if necessary
        
    }
    
    
    loadIntrospectionVariables() {
        
        super.loadIntrospectionVariables()
        this.superclass = UIViewController
        
    }
    
    
    static readonly routeComponentName = "fit_funasdasd"
    
    static readonly ParameterIdentifierName = { "ID": "ID" }
    
    
    loadSubviews() {
        
        this.view.backgroundColor = UIColor.whiteColor
        
        
        this.titleLabel = new UITextView(this.view.elementID + "TitleLabel", UITextView.type.header4)
        this.titleLabel.textAlignment = UITextView.textAlignment.left
        this.titleLabel.nativeSelectionEnabled = NO
        this.titleLabel.isSingleLine = NO
        this.view.addSubview(this.titleLabel)
        
        
        this.inputTextArea = new UITextArea(this.view.elementID + "InputTextArea")
        this.inputTextArea.placeholderText = "Input your data here."
        this.inputTextArea.changesOften = YES
        this.view.addSubview(this.inputTextArea)
        
        
        this.loadDataButton = new CBButton()
        this.loadDataButton.titleLabel.text = "Load data"
        this.view.addSubview(this.loadDataButton)
        
        this.loadDataButton.enabled = NO
        
        
        
        this.resultsLabel = new UIView()
        this.resultsLabel.innerHTML = "Results"
        this.resultsLabel.hidden = YES
        this.view.addSubview(this.resultsLabel)
        
        
        this.chartView = new UIView("ResultDataChart", nil, "canvas")
        this.chartView.hidden = YES
        this.view.addSubview(this.chartView)
        
        
        this.inputTextArea.addControlEventTarget.TextChange = (sender, event) => {
            
            //console.log(sender);
            
            this.view.setNeedsLayoutUpToRootView()
            
            this.loadDataButton.enabled = IS(this.inputTextArea.text)
            
        }
        
        
        // @ts-ignore
        const ctx = this.chartView.viewHTMLElement.getContext("2d")
        
        // @ts-ignore
        this.chart = new Chart(ctx, {
            // The type of chart we want to create
            //type: "line",
            type: "scatter",
            
            // The data for our dataset
            data: {
                //labels: ['January', 'February', 'March', 'April', 'May', 'June'],
                datasets: [
                    
                    {
                        label: "Input Data",
                        fill: false,
                        lineTension: 0,
                        //backgroundColor: UIColor.transparentColor,
                        borderColor: "rgb(125, 115, 252)" //,
                        //data: [{ x: 0, y: 0 }, { x: 2, y: 10 }, { x: 5, y: 5 }, { x: 7, y: 2 }, { x: 8, y: 20 }, {
                        // x: 11, y: 45 }]
                    },
                    {
                        label: "Fitted theoretical line",
                        fill: false,
                        lineTension: 0,
                        backgroundColor: "rgb(255, 99, 132)",
                        borderColor: "rgb(255, 99, 132)" //,
                        //data: [{ x: 0, y: 0 }, { x: 2, y: 10 }, { x: 5, y: 5 }, { x: 7, y: 2 }, { x: 8, y: 20 }, {
                        // x: 11, y: 45 }]
                    }
                
                
                ]
            },
            
            // Configuration options go here
            options: {
                
                fill: false,
                lineTension: 0,
                scales: {
                    x: {
                        type: "linear",
                        position: "bottom"
                    },
                    y: {
                        type: "logarithmic",
                        position: "bottom"
                    },
                    yAxes: [
                        {
                            type: "logarithmic",
                            ticks: {
                                min: 1,
                                max: 1000000,
                                callback: function (value, index, values) {
                                    if (value === 1000000) {
                                        return "1M"
                                    }
                                    if (value === 100000) {
                                        return "100K"
                                    }
                                    if (value === 10000) {
                                        return "10K"
                                    }
                                    if (value === 1000) {
                                        return "1K"
                                    }
                                    if (value === 100) {
                                        return "100"
                                    }
                                    if (value === 10) {
                                        return "10"
                                    }
                                    if (value === 1) {
                                        return "1"
                                    }
                                    return null
                                }
                            }
                        }
                    ]
                },
                tooltips: {
                    mode: "interpolate",
                    intersect: true
                },
                plugins: {
                    crosshair: {
                        line: {
                            color: "#F66",  // crosshair line color
                            width: 1        // crosshair line width
                        },
                        sync: {
                            enabled: true,            // enable trace line syncing with other charts
                            group: 1,                 // chart group
                            suppressTooltips: false   // suppress tooltips when showing a synced tracer
                        },
                        zoom: {
                            enabled: NO,                                      // enable zooming
                            zoomboxBackgroundColor: "rgba(66,133,244,0.2)",     // background color of zoom box
                            zoomboxBorderColor: "#48F",                         // border color of zoom box
                            zoomButtonText: "Reset Zoom",                       // reset zoom button text
                            zoomButtonClass: "reset-zoom"                      // reset zoom button class
                        },
                        callbacks: {
                            beforeZoom: function (start, end) {                  // called before zoom, return false to prevent zoom
                                return true
                            },
                            afterZoom: function (start, end) {                   // called after zoom
                            }
                        }
                    }
                }
                // responsive: false,
                // maintainAspectRatio: false,
                // showScale: true
                
            }
            
        })
        
        
        
        
        
        this.loadDataButton.addControlEventTarget.EnterDown.PointerUpInside = async (sender, event) => {
            
            //console.log(sender)
            
            
            try {
                
                CBDialogViewShower.showActionIndicatorDialog("Loading.")
                
                const stringSocketClientResult = await SocketClient.PerformFitFun({
                    inputData: this.inputTextArea.text
                })
                
                if (IS(stringSocketClientResult.errorResult)) {
                    
                    console.log(stringSocketClientResult.errorResult)
                    
                    CBDialogViewShower.alert(
                        "Failed to analyze data. " + JSON.stringify(stringSocketClientResult.errorResult)
                    )
                    
                    return
                    
                }
                
                const result = JSON.parse(stringSocketClientResult.result)
                
                //console.log(result)
                
                
                this.resultsLabel.innerHTML = FIRST(result.resultHTMLString, "")
                    + IF(IS_NOT(result.resultHTMLString))(RETURNER(
                        `\r\nRadius: ${ result.Radius } \r\nbeta: ${ result.beta } \r\nsigma_Aeff: ${ result.sigma_Aeff }`
                    ))()
                
                const pointPoints = this.pointObjectsFromValues(result.xplot_mrk, result.yplot_mrk)
                const linePoints = this.pointObjectsFromValues(result.xplot_line, result.yplot_line)
                
                
                this.chart.data.datasets[0].data = pointPoints
                this.chart.data.datasets[0].showLine = NO
                
                this.chart.data.datasets[1].data = linePoints
                this.chart.data.datasets[1].showLine = YES
                this.chart.data.datasets[1].pointRadius = 0
                this.chart.data.datasets[1].interpolate = YES
                
                
                this.chart.options.scales.xAxes[0].scaleLabel.display = true
                this.chart.options.scales.xAxes[0].scaleLabel.labelString = "" + FIRST(result.xAxisUnit, "")
                
                this.chart.options.scales.yAxes[0].scaleLabel.display = true
                this.chart.options.scales.yAxes[0].scaleLabel.labelString = "" + FIRST(result.yAxisUnit, "")
                
                this.resultsLabel.hidden = NO
                this.chartView.hidden = NO
                
                this.chart.update()
                
                this.view.setNeedsLayoutUpToRootView()
                
            } catch (exception) {
                
                console.log(exception)
                
                CBDialogViewShower.alert("Failed to display results. " + JSON.stringify(exception))
                
            }
            
            CBDialogViewShower.hideActionIndicatorDialog()
            
            
        }
        
        
    }
    
    
    
    
    
    private pointObjectsFromValues(xValues, yValues) {
        
        var resultPoints: any[] = []
        
        for (var i = 0; i < xValues.length; i++) {
            
            var x = xValues[i]
            
            var y = yValues[i]
            
            resultPoints.push({ x: x, y: y })
            
        }
        
        return resultPoints
        
    }
    
    
    
    
    
    async handleRoute(route: UIRoute) {
        
        super.handleRoute(route)
        
        const inquiryComponent = route.componentWithViewController(GETELECTableViewController)
        
        this.titleLabel.text = "Enter your data and press the load button. TABLE VIEW"
        
        route.didcompleteComponent(inquiryComponent)
        
    }
    
    
    
    
    
    layoutViewsManually() {
        
        super.layoutViewsManually()
        
        const padding = RootViewController.paddingLength
        const labelHeight = padding * 1.25
        
        // View bounds
        var bounds = this.view.bounds
        
        this.view.setPaddings(0, 0, padding, 0)
        
        
        
        this.titleLabel.frame = bounds.rectangleWithInsets(padding, padding, 0, padding)
        .rectangleWithHeight(this.titleLabel.intrinsicContentHeight(bounds.width))
        
        this.inputTextArea.frame = this.titleLabel.frame.rectangleForNextRow(
            padding * 0.25,
            this.inputTextArea.intrinsicContentHeight(this.titleLabel.frame.width) + 5
        )
        
        
        this.loadDataButton.frame = this.inputTextArea.frame.rectangleForNextRow(padding, labelHeight * 2)
        
        this.loadDataButton.setMargins(0, 0, padding, 0)
        
        
        this.resultsLabel.frame = this.loadDataButton.frame.rectangleForNextRow(padding)
        
        
        this.chartView.frame = this.resultsLabel.frame.rectangleForNextRow(
            padding,
            this.loadDataButton.frame.width * 0.5
        )
        //.rectangleWithInsets(padding * 5 * 0, padding * 5, 0, 0)
        //.rectangleByAddingX(-padding)
        
        //this.chartView.setPadding(padding);
        
        this.chartView.setMaxSizes(this.chartView.frame.height, this.chartView.frame.width)
        
        
        
    }
    
    
    intrinsicViewContentHeight(constrainingWidth: number = 0): number {
        
        const padding = RootViewController.paddingLength
        const labelHeight = padding * 1.5
        
        var result = padding + this.titleLabel.intrinsicContentHeight(constrainingWidth) + padding * 0.25 +
            labelHeight * 2 + padding * 0.5 + padding +
            this.inputTextArea.intrinsicContentHeight(constrainingWidth)
        
        if (IS_NOT(this.chartView.hidden)) {
            
            result = result + padding + (constrainingWidth - padding * 2) * 0.5 + padding + labelHeight * 2
            
        }
        
        return result
        
        
    }
    
    
}


// @ts-ignore
Chart.Interaction.modes["interpolate"] = function (chart, e, options) {
    
    // This function has a separate license
    
    // MIT License
    //
    // Copyright (c) 2018 Abel Heinsbroek
    //
    // Permission is hereby granted, free of charge, to any person obtaining a copy
    // of this software and associated documentation files (the "Software"), to deal
    // in the Software without restriction, including without limitation the rights
    // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    // copies of the Software, and to permit persons to whom the Software is
    // furnished to do so, subject to the following conditions:
    //
    //     The above copyright notice and this permission notice shall be included in all
    // copies or substantial portions of the Software.
    //
    //     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    //     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    //     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    // SOFTWARE.
    
    var items = []
    
    for (var datasetIndex = 0; datasetIndex < chart.data.datasets.length; datasetIndex++) {
        
        
        // check for interpolate setting
        if (!chart.data.datasets[datasetIndex].interpolate) {
            continue
        }
        
        var meta = chart.getDatasetMeta(datasetIndex)
        // do not interpolate hidden charts
        if (meta.hidden) {
            continue
        }
        
        
        var xScale = chart.scales[meta.xAxisID]
        var yScale = chart.scales[meta.yAxisID]
        
        var xValue = xScale.getValueForPixel(e.x)
        
        
        var data = chart.data.datasets[datasetIndex].data
        
        var index = data.findIndex(function (o) {
            return o.x >= xValue
        })
        
        if (index === -1) {
            continue
        }
        
        
        // linear interpolate value
        var prev = data[index - 1]
        var next = data[index]
        
        if (prev && next) {
            var slope = (next.y - prev.y) / (next.x - prev.x)
            var interpolatedValue = prev.y + (xValue - prev.x) * slope
        }
        
        if (chart.data.datasets[datasetIndex].steppedLine && prev) {
            interpolatedValue = prev.y
        }
        
        if (isNaN(interpolatedValue)) {
            continue
        }
        
        var yPosition = yScale.getPixelForValue(interpolatedValue)
        
        // do not interpolate values outside of the axis limits
        if (isNaN(yPosition)) {
            continue
        }
        
        // create a 'fake' event point
        
        var fakePoint = {
            
            value: interpolatedValue,
            xValue: xValue,
            
            tooltipPosition: function () {
                return this._model
            },
            hasValue: function () {
                return true
            },
            _model: {
                x: e.x,
                y: yPosition
            },
            _datasetIndex: datasetIndex,
            _index: items.length,
            _xScale: {
                getLabelForIndex: function (indx) {
                    
                    let xValue: number = items[indx].xValue
                    
                    xValue = xValue.toPrecision(5) as any
                    
                    return xValue
                    
                }
            },
            _yScale: {
                getLabelForIndex: function (indx) {
                    
                    let value = items[indx].value.toPrecision(5)
                    
                    return value
                    
                }
            },
            _chart: chart
            
        }
        
        items.push(fakePoint)
        
    }
    
    
    // add other, not interpolated, items
    // @ts-ignore
    var xItems = Chart.Interaction.modes.x(chart, e, options)
    
    xItems.forEach((item, index, array) => {
        
        
        if (!chart.data.datasets[item._datasetIndex].interpolate) {
            items.push(item)
        }
        
    })
    
    
    
    return items
}








































































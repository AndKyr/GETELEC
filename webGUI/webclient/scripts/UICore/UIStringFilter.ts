///<reference path="UIObject.ts"/>










class UIStringFilter extends UIObject {
    
    
    static _sharedWebWorker = new Worker("compiledScripts/UIStringFilterWebWorker.js")
    
    static _instanceNumber = -1
    
    
    _instanceNumber: number
    
    _isThreadClosed = NO
    
    _webWorker = UIStringFilter._sharedWebWorker
    
    constructor(useSeparateWebWorker = NO) {
        
        
        super()
        
        this._class = UIStringFilter
        this.superclass = UIObject
        
        
        if (useSeparateWebWorker) {
            
            this._webWorker = new Worker("compiledScripts/UIStringFilterWebWorker.js")
            
        }
        
        UIStringFilter._instanceNumber = UIStringFilter._instanceNumber + 1
        this._instanceNumber = UIStringFilter._instanceNumber
        
        
        
        
        
    }
    
    
    
    
    
    get instanceIdentifier() {
        
        return this._instanceNumber
        
    }
    
    
    
    
    
    filterData(
        filteringString: string,
        data: string[],
        excludedData: string[],
        identifier: any,
        completion: (filteredData: string[], filteredIndexes: number[], identifier: any) => void
    ) {
        
        
        if (this._isThreadClosed) {
            
            return
            
        }
        
        //var startTime = Date.now();
    
        const instanceIdentifier = this.instanceIdentifier
    
        this._webWorker.onmessage = function (message) {
            
            if (message.data.instanceIdentifier == instanceIdentifier) {
                
                //console.log("Filtering took " + (Date.now() - startTime) + " ms");
                
                completion(message.data.filteredData, message.data.filteredIndexes, message.data.identifier)
                
            }
            
        }
        
        this._webWorker.postMessage({
            
            "filteringString": filteringString,
            "data": data,
            "excludedData": excludedData,
            "identifier": identifier,
            "instanceIdentifier": instanceIdentifier
            
        })
        
        
    }
    
    
    filteredData(
        filteringString: string,
        data: string[],
        excludedData: string[] = [],
        identifier: any = MAKE_ID()
    ) {
    
    
        const result: Promise<{
        
            filteredData: string[],
            filteredIndexes: number[],
            identifier: any
        
        }> = new Promise((resolve, reject) => {
        
            this.filterData(filteringString, data, excludedData, identifier,
                (filteredData, filteredIndexes, filteredIdentifier) => {
    
                    if (filteredIdentifier == identifier) {
    
                        resolve({
        
                            filteredData: filteredData,
                            filteredIndexes: filteredIndexes,
                            identifier: filteredIdentifier
        
                        })
        
                    }
                    
                
                }
            )
        
        
        })
    
        return result
        
        
    }
    
    
    
    
    
    closeThread() {
        
        this._isThreadClosed = YES
        
        if (this._webWorker != UIStringFilter._sharedWebWorker) {
            
            this._webWorker.terminate()
            
        }
        
        
        
    }
    
    
    
    
    
}






































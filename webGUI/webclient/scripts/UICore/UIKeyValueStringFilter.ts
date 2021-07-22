///<reference path="UIObject.ts"/>








class UIKeyValueStringFilter extends UIObject {
    
    
    static _sharedWebWorker = new Worker("compiledScripts/UIKeyValueStringFilterWebWorker.js")
    
    static _instanceNumber = -1
    
    
    _instanceNumber: number
    
    _isThreadClosed = NO
    
    _webWorker = UIKeyValueStringFilter._sharedWebWorker
    
    constructor(useSeparateWebWorker = NO) {
        
        
        super()
        
        this._class = UIKeyValueStringFilter
        this.superclass = UIObject
        
        
        if (useSeparateWebWorker) {
            
            this._webWorker = new Worker("compiledScripts/UIKeyValueStringFilterWebWorker.js")
            
        }
        
        UIKeyValueStringFilter._instanceNumber = UIKeyValueStringFilter._instanceNumber + 1
        this._instanceNumber = UIKeyValueStringFilter._instanceNumber
        
        
    }
    
    
    
    
    
    get instanceIdentifier() {
        
        return this._instanceNumber
        
    }
    
    
    
    
    
    filterData(
        filteringString: string,
        data: any[],
        excludedData: string[],
        dataKeyPath: string,
        identifier: any,
        completion: (filteredData: string[], filteredIndexes: string[], identifier: any) => void
    ) {
        
        
        if (this._isThreadClosed) {
            
            return
            
        }
    
        const startTime = Date.now()
    
        const instanceIdentifier = this.instanceIdentifier
    
        this._webWorker.onmessage = function (message) {
            
            if (message.data.instanceIdentifier == instanceIdentifier) {
                
                console.log("Filtering took " + (Date.now() - startTime) + " ms.")
                
                completion(message.data.filteredData, message.data.filteredIndexes, message.data.identifier)
                
            }
            
        }
        
        
        try {
            
            this._webWorker.postMessage({
                
                "filteringString": filteringString,
                "data": data,
                "excludedData": excludedData,
                "dataKeyPath": dataKeyPath,
                "identifier": identifier,
                "instanceIdentifier": instanceIdentifier
                
            })
            
        } catch (exception) {
            
            completion([], [], identifier)
            
        }
        
        
        
        
        
    }
    
    
    
    
    
    closeThread() {
        
        this._isThreadClosed = YES
        
        if (this._webWorker != UIKeyValueStringFilter._sharedWebWorker) {
            
            this._webWorker.terminate()
            
        }
        
        
        
    }
    
    
    
    
    
}






































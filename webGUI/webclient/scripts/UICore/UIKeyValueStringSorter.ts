///<reference path="UIObject.ts"/>











interface UIKeyValueStringSorterSortingInstruction {
    
    keyPath: string;
    
    dataType: string;
    
    direction: string;
    
    
}





class UIKeyValueStringSorter extends UIObject {
    
    
    static _sharedWebWorker = new Worker("compiledScripts/UIKeyValueStringSorterWebWorker.js")
    
    static _instanceNumber = -1
    
    
    _instanceNumber: number
    
    _isThreadClosed = NO
    
    _webWorker = UIKeyValueStringSorter._sharedWebWorker
    
    constructor(useSeparateWebWorker = NO) {
        
        
        super()
        
        this._class = UIKeyValueStringSorter
        this.superclass = UIObject
        
        
        if (useSeparateWebWorker) {
            
            this._webWorker = new Worker("compiledScripts/UIKeyValueStringSorterWebWorker.js")
            
        }
        
        UIKeyValueStringSorter._instanceNumber = UIKeyValueStringSorter._instanceNumber + 1
        this._instanceNumber = UIKeyValueStringSorter._instanceNumber
        
        
    }
    
    
    
    
    
    get instanceIdentifier() {
        
        return this._instanceNumber
        
    }
    
    
    
    
    
    static dataType = {
        
        "string": "string"
        
    }
    
    
    static direction = {
        
        "descending": "descending",
        "ascending": "ascending"
        
    }
    
    
    
    
    
    sortData<T>(
        data: T[],
        sortingInstructions: UIKeyValueStringSorterSortingInstruction[],
        identifier: any,
        completion: (sortedData: T[], sortedIndexes: number[], identifier: any) => void
    ) {
        
        
        if (this._isThreadClosed) {
            
            return
            
        }
    
        const startTime = Date.now()
    
        const instanceIdentifier = this.instanceIdentifier
    
        this._webWorker.onmessage = function (message) {
            
            if (message.data.instanceIdentifier == instanceIdentifier) {
                
                console.log("Sorting " + data.length + " items took " + (Date.now() - startTime) + " ms.")
                
                completion(message.data.sortedData, message.data.sortedIndexes, message.data.identifier)
                
            }
            
        }
        
        
        try {
            
            this._webWorker.postMessage({
                
                "data": data,
                "sortingInstructions": sortingInstructions,
                "identifier": identifier,
                "instanceIdentifier": instanceIdentifier
                
            })
            
        } catch (exception) {
            
            completion([], [], identifier)
            
        }
        
        
        
        
        
    }
    
    
    
    
    
    sortedData<T>(
        data: T[],
        sortingInstructions: UIKeyValueStringSorterSortingInstruction[],
        identifier: any = MAKE_ID()
    ) {
    
        const result: Promise<{
        
            sortedData: T[],
            sortedIndexes: number[],
            identifier: any
        
        }> = new Promise((resolve, reject) => {
        
            this.sortData(data, sortingInstructions, identifier, (sortedData, sortedIndexes, sortedIdentifier) => {
    
                if (sortedIdentifier == identifier) {
    
                    resolve({
        
                        sortedData: sortedData,
                        sortedIndexes: sortedIndexes,
                        identifier: sortedIdentifier
        
                    })
        
                }
                
                
            
            })
        
        
        
        })
    
        return result
        
    }
    
    
    
    
    
    closeThread() {
        
        this._isThreadClosed = YES
        
        if (this._webWorker != UIKeyValueStringSorter._sharedWebWorker) {
            
            this._webWorker.terminate()
            
        }
        
        
        
    }
    
    
    
    
    
}







































// @ts-check

if ("contains" in Array.prototype == false) {
    
    // @ts-ignore
    Array.prototype.contains = function (element) {
        
        var result = (this.indexOf(element) != -1);
        return result;
        
    }
    
}

if ("contains" in String.prototype == false) {
    
    // @ts-ignore
    String.prototype.contains = function (string) {
        
        var result = (this.indexOf(string) != -1);
        return result;
        
    }
    
}





onmessage = function (event) {
    
    //console.log('Message received from main script');
    var workerResult = filterData(event.data.filteringString, event.data.data, event.data.excludedData, event.data.dataKeyPath);
    
    workerResult.identifier = event.data.identifier;
    workerResult.instanceIdentifier = event.data.instanceIdentifier;
    
    
    // @ts-ignore
    postMessage(workerResult);
    
}



function valueForKeyPath(keyPath, object) {
    
    var keys = keyPath.split(".");
    var currentObject = object;
    
    keys.forEach(function (key, index, array) {
        currentObject = currentObject[key];  
    })
    
    return currentObject;
    
}





function filterData(filteringString, data, excludedData, dataKeyPath) {

    var filteredData = [];
    var filteredIndexes = [];

    if (filteringString) {

        var filteringStringWords = [];
        filteringString.split(" ").forEach(function (word, index, array) {
            if (word) {
                filteringStringWords.push(word.toLowerCase());
            }
        })
        
        data.forEach(function (dataObject, index, array) {

            var dataString = valueForKeyPath(dataKeyPath, dataObject);
            
            var lowercaseDataString = dataString.toLowerCase();
            
            // Look through all the words in the input
            var wordsFound = [];
            filteringStringWords.forEach(function (word) {
                wordsFound.push(lowercaseDataString.contains(word) && !excludedData.contains(dataString));
            })

            // Only show the dataString if it matches all of them
            // @ts-ignore
            if (wordsFound.contains(true) && !wordsFound.contains(false)) {
                
                filteredData.push(dataObject);
                filteredIndexes.push(index);
                
            }

        })

        

    }
    else if (excludedData.length) {

        
        filteredData = data.forEach(function (dataObject, index, array) {
            
            if (excludedData.indexOf(dataObject) == -1) {
                
                filteredData.push(dataObject);
                filteredIndexes.push(index);
                
            }
            
        })

    }
    else {

        filteredData = data;
        
        data.forEach(function (object, index, array) {
            
            filteredIndexes.push(index);
            
        })

    }

    
    
    return { "filteredData": filteredData, "filteredIndexes": filteredIndexes };
    
    

}

























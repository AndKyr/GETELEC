
// @ts-check









onmessage = function (event) {
    
    //console.log('Message received from main script');
    var workerResult = sortData(
        event.data.data,
        event.data.sortingInstructions
    )
    
    workerResult.identifier = event.data.identifier;
    workerResult.instanceIdentifier = event.data.instanceIdentifier;
    
    
    // @ts-ignore
    postMessage(workerResult);
    
}





function valueForKeyPath(keyPath, object) {
        
    var keys = keyPath.split(".");
    var currentObject = object;
    
    for (var i = 0; i < keys.length; i++) {
        
        var key = keys[i];
        
        if (key.substring(0, 2) == "[]") {
            
            // This next object will be an array and the rest of the keys need to be run for each of the elements
            
            currentObject = currentObject[key.substring(2)];
            
            // CurrentObject is now an array
            
            var remainingKeyPath = keys.slice(i + 1).join(".");
            
            var currentArray = currentObject;
            
            currentObject = currentArray.map(function (subObject, index, array) {
                
                var result = valueForKeyPath(remainingKeyPath, subObject);
                
                return result;
                
            })
            
            break;
            
        }
        
        currentObject = (currentObject || {})[key];
        
        
    }
    
    return currentObject;
    
}







function compare(firstObject, secondObject, sortingInstructions) {
    
    
    if (sortingInstructions.length == 0) {
        return 0;
    }
    
    
    var sortingInstruction = sortingInstructions[0];
    
    
    var directionMultiplier = 1;
    if (sortingInstruction.direction == "descending") {
        directionMultiplier = -1;
    }
    
    
    var firstDataString = firstObject[sortingInstruction.keyPath];
    
    var secondDataString = secondObject[sortingInstruction.keyPath];
    
    
    
    
    if (firstDataString < secondDataString) {
        
        return -1 * directionMultiplier;
        
    }
    
    if (firstDataString > secondDataString) {
        
        return 1 * directionMultiplier;
        
    }
    
    if (sortingInstructions.length > 1) {
        
        var remainingSortingInstructions = sortingInstructions.slice(1);
        
        
        
        return compare(firstObject, secondObject, remainingSortingInstructions);
        
        
    }
    
    return 0;
    
  }

  
  
  

function sortData(data, sortingInstructions) {
    
    
    var sortingObjects = data.map(function (dataItem, index, array) {
        
        var result = {
            
            "_UIKeyValueStringSorterWebWorkerSortingObjectIndex": index
            
        }
        
        
        sortingInstructions.forEach(function(instruction, index, instructionsArray) {
            
            result[instruction.keyPath] = JSON.stringify(valueForKeyPath(instruction.keyPath, dataItem) || {}).toLowerCase();
            
        })
        
        
        
        
        return result;
        
    
    })
    
    
    var sortedData = sortingObjects.sort(function (firstObject, secondObject) {
            
        return compare(firstObject, secondObject, sortingInstructions);

    })
    
    var sortedIndexes = sortedData.map(function (object, index, array) {
        
        var sortedIndex = object._UIKeyValueStringSorterWebWorkerSortingObjectIndex;
        
        return sortedIndex;
        
    })
    
    var result = {
        
        "sortedData": sortedIndexes.map(function (sortedIndex, index, array) {
            
            return data[sortedIndex];
            
        }),
        "sortedIndexes": sortedIndexes
        
    }
    
    
    return result;
    
    
}

























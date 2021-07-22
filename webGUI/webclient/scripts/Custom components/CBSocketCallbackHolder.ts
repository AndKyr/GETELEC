/// <reference path="../UICore/UIObject.ts" />
/// <reference path="../UICore/UICoreExtensions.ts" />





interface CBSocketCallbackHolderMessageDescriptor {
    
    key: string;
    message: {
        identifier: string;
        inResponseToIdentifier?: string;
        keepWaitingForResponses?: boolean;
    }
    
    
    sentAtTime: number;
    
    //completionTriggered: boolean;
    
    messageDataHash: string;
    
    responseDataHash?: string;
    
    mainResponseReceived: boolean;
    
    anyMainResponseReceived: boolean;
    
    completionPolicy: string;
    completionFunction: CBSocketMessagecompletionFunction;
    
}


interface CBSocketCallbackHolderStoredResponseObject {
    
    messageKey: string;
    messageData: any;
    messageDataHash: string;
    
}




class CBSocketCallbackHolder extends UIObject {
    
    
    
    
    
    messageDescriptors: {
        
        [x: string]: CBSocketCallbackHolderMessageDescriptor[]
        
    } = {}
    
    handlers: {
        [x: string]: CBSocketMessageHandlerFunction[]
    } = {}
    
    onetimeHandlers: {
        [x: string]: CBSocketMessageHandlerFunction[]
    } = {}
    
    keysForIdentifiers: {
        
        [x: string]: string
        
    } = {}
    
    
    isValid = YES
    _storeableResponseKeys: string[]
    _storedResponseHashesDictionary: {
        
        [x: string]: {
            
            hash: string,
            validityDate: number
            
            
            
            
        }
        
    }
    _verifiedResponseHashesDictionary: {
        
        [x: string]: boolean
        
    } = {}
    
    _socketClient: CBSocketClient
    
    
    
    
    
    constructor(socketClient: CBSocketClient, previousCallbackHolder?: CBSocketCallbackHolder) {
        
        super()
        
        
        this._socketClient = socketClient
        
        if (IS(previousCallbackHolder)) {
            
            this.handlers = previousCallbackHolder.handlers
            this._verifiedResponseHashesDictionary = previousCallbackHolder._verifiedResponseHashesDictionary
            
        }
        
        
        
        
        
    }
    
    
    
    
    
    triggerDisconnectHandlers() {
        
        this.messageDescriptors.forEach(function (descriptor: CBSocketCallbackHolderMessageDescriptor, key: string) {
            
            if (descriptor.mainResponseReceived) {
                
                descriptor.completionFunction(CBSocketClient.disconnectionMessage, nil)
                
            }
            
        })
        
    }
    
    
    
    
    
    registerHandler(key: string, handlerFunction: CBSocketMessageHandlerFunction) {
        
        
        if (!this.handlers[key]) {
            
            this.handlers[key] = []
            
        }
        
        this.handlers[key].push(handlerFunction)
        
        
        
    }
    
    registerOnetimeHandler(key: string, handlerFunction: CBSocketMessageHandlerFunction) {
        
        
        if (!this.onetimeHandlers[key]) {
            
            this.onetimeHandlers[key] = []
            
        }
        
        this.onetimeHandlers[key].push(handlerFunction)
        
        
        
    }
    
    
    
    
    
    get storedResponseHashesDictionary() {
        
        if (IS_NOT(this._storedResponseHashesDictionary)) {
            
            this._storedResponseHashesDictionary = JSON.parse(localStorage["CBSocketResponseHashesDictionary"] || "{}")
            
        }
        
        return this._storedResponseHashesDictionary
        
    }
    
    storedResponseHashObjectForKey(requestKey: string, requestDataHash: string) {
        
        const localStorageKey = this.keyForRequestKeyAndRequestDataHash(requestKey, requestDataHash)
        
        const hashObject = this.storedResponseHashesDictionary[localStorageKey]
        
        const result = FIRST(hashObject, {} as any)
        
        
        
        return result
        
    }
    
    storedResponseForKey(requestKey: string, requestDataHash: string) {
        
        const localStorageKey = this.keyForRequestKeyAndRequestDataHash(requestKey, requestDataHash)
        
        const storedObject = JSON.parse(localStorage[localStorageKey] || "{}")
        
        return storedObject.responseMessageData
        
    }
    
    keyForRequestKeyAndRequestDataHash(requestKey: string, requestDataHash: string) {
        
        const result = "_CBSCH_LS_key_" + requestKey + "_" + requestDataHash
        
        return result
        
    }
    
    storeResponse(
        requestKey: string,
        requestDataHash: string,
        responseMessage: CBSocketMessage<any>,
        responseDataHash: string
    ) {
        
        
        if (!responseMessage.canBeStoredAsResponse ||
            (IS_NOT(responseMessage.messageData) && IS_NOT(responseMessage.messageDataHash))) {
            
            return
            
        }
        
        
        const localStorageKey = this.keyForRequestKeyAndRequestDataHash(requestKey, requestDataHash)
        
        
        var validityDate: number
        
        if (responseMessage.responseValidityDuration) {
            
            validityDate = Date.now() + responseMessage.responseValidityDuration
            
        }
        
        const storedResponseHashesDictionary = this.storedResponseHashesDictionary
        storedResponseHashesDictionary[localStorageKey] = {
            
            hash: responseDataHash,
            validityDate: validityDate
            
        }
        
        this.saveInLocalStorage(localStorageKey, {
            
            responseMessageData: responseMessage.messageData,
            responseHash: responseDataHash
            
        })
        
        
        this.saveStoredResponseHashesDictionary(storedResponseHashesDictionary)
        
    }
    
    
    
    
    
    private saveStoredResponseHashesDictionary(storedResponseHashesDictionary: { [x: string]: { hash: string; validityDate: number; }; }) {
        
        this.saveInLocalStorage("CBSocketResponseHashesDictionary", storedResponseHashesDictionary)
        
    }
    
    saveInLocalStorage(key: string, object: any) {
        
        
        const stringToSave = JSON.stringify(object)
        
        if (stringToSave != localStorage[key]) {
            
            localStorage[key] = stringToSave
            
        }
        
        
    }
    
    
    
    
    
    socketShouldSendMessage(
        key: string,
        message: CBSocketMessage<any>,
        completionPolicy: string,
        completionFunction: CBSocketMessagecompletionFunction
    ) {
        
        
        
        var result = YES
        
        var triggerStoredResponseImmediately = NO
        
        
        // @ts-ignore
        const messageDataHash = objectHash(message.messageData || nil)
        
        const descriptorKey = "socketMessageDescriptor_" + key + messageDataHash
        
        this.messageDescriptors[descriptorKey] = (this.messageDescriptors[descriptorKey] || [])
        
        
        const hashObject = this.storedResponseHashObjectForKey(key, messageDataHash)
        message.storedResponseHash = hashObject.hash
        
        
        
        if (completionPolicy == CBSocketClient.completionPolicy.first) {
            
            const descriptorsForKey = (this.messageDescriptors[descriptorKey] || [])
            
            const matchingDescriptor = descriptorsForKey.find(function (descriptor, index, array) {
                return (descriptor.messageDataHash == messageDataHash)
            })
            
            if (matchingDescriptor) {
                
                result = NO
                
            }
            
        }
        
        if (completionPolicy == CBSocketClient.completionPolicy.storedOrFirst) {
            
            const descriptorsForKey = (this.messageDescriptors[descriptorKey] || [])
            
            const matchingDescriptor = descriptorsForKey.find(function (descriptor, index, array) {
                return (descriptor.messageDataHash == messageDataHash)
            })
            
            const storedResponse = IS(message.storedResponseHash)
            
            if (matchingDescriptor ||
                (storedResponse && this._verifiedResponseHashesDictionary[message.storedResponseHash])) {
                
                result = NO
                
                triggerStoredResponseImmediately = YES
                
            }
            
        }
        
        if (completionPolicy == CBSocketClient.completionPolicy.firstOnly) {
            
            const descriptorsForKey = (this.messageDescriptors[descriptorKey] || [])
            
            const matchingDescriptor = descriptorsForKey.find(function (descriptor, index, array) {
                return (descriptor.messageDataHash == messageDataHash)
            })
            
            if (matchingDescriptor) {
                
                return NO
                
            }
            
        }
        
        
        if (hashObject && hashObject.hash && hashObject.validityDate && message.storedResponseHash &&
            this._verifiedResponseHashesDictionary[message.storedResponseHash] && hashObject.validityDate >
            Date.now()) {
            
            result = NO
            
            triggerStoredResponseImmediately = YES
            
        }
        
        
        if (IS(completionFunction)) {
            
            this.messageDescriptors[descriptorKey].push({
                
                key: key,
                message: {
                    
                    identifier: message.identifier,
                    inResponseToIdentifier: message.inResponseToIdentifier,
                    keepWaitingForResponses: message.keepWaitingForResponses
                    
                },
                
                sentAtTime: Date.now(),
                
                //completionTriggered: NO,
                
                
                messageDataHash: messageDataHash,
                
                mainResponseReceived: NO,
                anyMainResponseReceived: NO,
                
                
                completionPolicy: completionPolicy,
                completionFunction: completionFunction
                
            })
            
            this.keysForIdentifiers[message.identifier] = descriptorKey
            
        }
        
        
        if (triggerStoredResponseImmediately) {
            
            this.socketDidReceiveMessageForKey(
                CBSocketClient.responseMessageKey,
                {
                    
                    identifier: nil,
                    messageData: nil,
                    
                    inResponseToIdentifier: message.identifier,
                    
                    useStoredResponse: YES
                    
                },
                nil
            )
            
        }
        
        
        return result
        
        
    }
    
    
    
    
    
    static defaultMultipleMessagecompletionFunction(responseMessages: any[], callcompletionFunctions: () => void) {
        callcompletionFunctions()
    }
    
    
    socketWillSendMultipleMessage(
        messageToSend: CBSocketMultipleMessage,
        completionFunction: CBSocketMultipleMessagecompletionFunction = CBSocketCallbackHolder.defaultMultipleMessagecompletionFunction
    ) {
        
        
        
        const key = CBSocketClient.multipleMessageKey
        
        
        // @ts-ignore
        const messageDataHash = objectHash(messageToSend.messageData || nil)
        
        const descriptorKey = "socketMessageDescriptor_" + key + messageDataHash
        
        this.messageDescriptors[descriptorKey] = (this.messageDescriptors[descriptorKey] || [])
        
        
        messageToSend.storedResponseHash = this.storedResponseHashObjectForKey(key, messageDataHash).hash
        
        
        
        
        
        this.messageDescriptors[descriptorKey].push({
            
            key: key,
            message: {
                
                identifier: messageToSend.identifier,
                inResponseToIdentifier: messageToSend.inResponseToIdentifier,
                keepWaitingForResponses: messageToSend.keepWaitingForResponses
                
            },
            
            sentAtTime: Date.now(),
            
            //completionTriggered: NO,
            
            
            messageDataHash: messageDataHash,
            
            mainResponseReceived: NO,
            anyMainResponseReceived: NO,
            
            
            completionPolicy: CBSocketClient.completionPolicy.directOnly,
            completionFunction: function (
                this: CBSocketCallbackHolder,
                responseMessage: CBSocketMultipleMessageObject[],
                respondWithMessage
            ) {
                
                completionFunction(
                    responseMessage.map(function (messageObject, index, array) {
                        
                        return messageObject.message.messageData
                        
                    }),
                    function (this: CBSocketCallbackHolder) {
                        
                        //console.log("Received multiple message response with length of " + responseMessage.length + ".");
                        
                        // Call all completion functions
                        responseMessage.forEach(function (
                            this: CBSocketCallbackHolder,
                            messageObject: CBSocketMultipleMessageObject,
                            index: number,
                            array: CBSocketMultipleMessageObject[]
                        ) {
                            
                            this._socketClient.didReceiveMessageForKey(messageObject.key, messageObject.message)
                            
                        }.bind(this))
                        
                    }.bind(this)
                )
                
            }.bind(this)
            
        })
        
        this.keysForIdentifiers[messageToSend.identifier] = descriptorKey
        
        
    }
    
    
    
    
    
    socketDidReceiveMessageForKey(
        key: string,
        message: CBSocketMessage<any>,
        sendResponseFunction: CBSocketMessageSendResponseFunction
    ) {
        
        
        if (!this.isValid) {
            
            return
            
        }
        
        
        // Call static handlers
        if (this.handlers[key]) {
            
            this.handlers[key].forEach(function (
                this: CBSocketCallbackHolder,
                handler: CBSocketMessageHandlerFunction,
                index,
                array
            ) {
                
                handler(message.messageData, sendResponseFunction)
                
            }.bind(this))
            
        }
        
        if (this.onetimeHandlers[key]) {
            
            this.onetimeHandlers[key].forEach(function (
                this: CBSocketCallbackHolder,
                handler: CBSocketMessageHandlerFunction,
                index,
                array
            ) {
                
                handler(message.messageData, sendResponseFunction)
                
            }.bind(this))
            
            delete this.onetimeHandlers[key]
            
        }
        
        
        
        // Temporary response handlers are evaluated here
        if (message.inResponseToIdentifier &&
            (CBSocketClient.responseMessageKey == key || CBSocketClient.multipleMessageKey == key)) {
            
            // Find descriptors for the key of the message that is being responded to
            const descriptorKey = this.keysForIdentifiers[message.inResponseToIdentifier]
            const descriptorsForKey = (this.messageDescriptors[descriptorKey] || [])
            
            // Find response data hash to check for differences
            const responseDataHash = message.messageDataHash
            
            // Remove identifier from dictionary
            if (!message.keepWaitingForResponses) {
                
                delete this.keysForIdentifiers[message.inResponseToIdentifier]
                
                delete this.messageDescriptors[descriptorKey]
                
            }
            
            
            // Function to call completion function
            const callCompletionFunction = (descriptor, storedResponseCondition = NO) => {
                
                var messageData = message.messageData
                
                if (message.useStoredResponse && storedResponseCondition) {
                    
                    messageData = this.storedResponseForKey(descriptor.key, descriptor.messageDataHash)
                    
                    const responseHash = this.storedResponseHashObjectForKey(
                        descriptor.key,
                        descriptor.messageDataHash
                    ).hash
                    
                    const localStorageKey = this.keyForRequestKeyAndRequestDataHash(
                        descriptor.key,
                        descriptor.messageDataHash
                    )
                    
                    if (message.responseValidityDuration && this.storedResponseHashesDictionary[localStorageKey]) {
                        
                        this.storedResponseHashesDictionary[localStorageKey].validityDate = Date.now() +
                            message.responseValidityDuration
                        
                        this.saveStoredResponseHashesDictionary(this.storedResponseHashesDictionary)
                        
                    }
                    
                    this._verifiedResponseHashesDictionary[responseHash] = YES
                    
                    console.log("Using stored response.")
                    
                }
                
                // Call completionFunction and set response data hash
                descriptor.completionFunction(messageData, sendResponseFunction)
                descriptor.responseDataHash = responseDataHash
                
            }
            
            
            descriptorsForKey.copy().forEach(function (
                this: CBSocketCallbackHolder,
                descriptor: CBSocketCallbackHolderMessageDescriptor,
                index: number,
                array: CBSocketCallbackHolderMessageDescriptor[]
            ) {
                
                
                if ((descriptor.completionPolicy == CBSocketClient.completionPolicy.directOnly &&
                    descriptor.message.identifier == message.inResponseToIdentifier) || descriptor.completionPolicy ==
                    CBSocketClient.completionPolicy.first || descriptor.completionPolicy ==
                    CBSocketClient.completionPolicy.firstOnly || descriptor.completionPolicy ==
                    CBSocketClient.completionPolicy.storedOrFirst) {
                    
                    // Calling completion function and removing descriptor
                    
                    if (!message.keepWaitingForResponses) {
                        
                        this.storeResponse(descriptor.key, descriptor.messageDataHash, message, responseDataHash)
                        
                        descriptorsForKey.removeElement(descriptor)
                        
                        sendResponseFunction.respondingToMainResponse = YES
                        
                    }
                    
                    callCompletionFunction(descriptor, !message.keepWaitingForResponses)
                    
                }
                else if (descriptor.completionPolicy == CBSocketClient.completionPolicy.all) {
                    
                    // Calling completion function
                    callCompletionFunction(descriptor, !message.keepWaitingForResponses)
                    
                    // Marking descriptor as having been responded to
                    if (!message.keepWaitingForResponses) {
                        
                        if (message.inResponseToIdentifier == descriptor.message.identifier) {
                            
                            sendResponseFunction.respondingToMainResponse = YES
                            descriptor.mainResponseReceived = YES
                            descriptorsForKey.removeElement(descriptor)
                            
                        }
                        
                        descriptor.anyMainResponseReceived = YES
                        
                    }
                    
                    
                    
                }
                else if (descriptor.completionPolicy == CBSocketClient.completionPolicy.allDifferent) {
                    
                    // Calling completionFunction if messageData is different from previous
                    if (descriptor.responseDataHash != responseDataHash) {
                        
                        callCompletionFunction(descriptor, !message.keepWaitingForResponses)
                        
                    }
                    
                    // Marking descriptor as having been responded to
                    if (!message.keepWaitingForResponses) {
                        
                        if (message.inResponseToIdentifier == descriptor.message.identifier) {
                            
                            sendResponseFunction.respondingToMainResponse = YES
                            descriptor.mainResponseReceived = YES
                            descriptorsForKey.removeElement(descriptor)
                            
                        }
                        
                        descriptor.anyMainResponseReceived = YES
                        
                    }
                    
                }
                else if (descriptor.completionPolicy == CBSocketClient.completionPolicy.last &&
                    descriptor.message.identifier == message.inResponseToIdentifier) {
                    
                    if (!message.keepWaitingForResponses) {
                        
                        // Marking descriptor as having been responded to
                        descriptor.mainResponseReceived = YES
                        descriptor.anyMainResponseReceived = YES
                        
                        sendResponseFunction.respondingToMainResponse = YES
                        
                    }
                    else {
                        
                        descriptor.completionFunction(message.messageData, sendResponseFunction)
                        
                    }
                    
                }
                else if (descriptor.completionPolicy == CBSocketClient.completionPolicy.firstAndLast ||
                    descriptor.completionPolicy == CBSocketClient.completionPolicy.firstAndLastIfDifferent) {
                    
                    if (!message.keepWaitingForResponses) {
                        
                        // Only calling completionFunction once as a first response call
                        if (!descriptor.anyMainResponseReceived) {
                            
                            callCompletionFunction(descriptor, !message.keepWaitingForResponses)
                            
                        }
                        
                        // Marking descriptor as having been responded to
                        if (descriptor.message.identifier == message.inResponseToIdentifier) {
                            
                            descriptor.mainResponseReceived = YES
                            sendResponseFunction.respondingToMainResponse = YES
                            
                        }
                        
                        descriptor.anyMainResponseReceived = YES
                        
                    }
                    else if (descriptor.message.identifier == message.inResponseToIdentifier &&
                        message.keepWaitingForResponses) {
                        
                        descriptor.completionFunction(message.messageData, sendResponseFunction)
                        
                    }
                    
                }
                
            }.bind(this))
            
            
            
            
            
            // Last message completion policies
            
            const allResponsesReceived = descriptorsForKey.allMatch(function (descriptorObject, index, array) {
                return descriptorObject.mainResponseReceived
            })
            
            descriptorsForKey.copy().forEach(function (
                this: CBSocketCallbackHolder,
                descriptor: CBSocketCallbackHolderMessageDescriptor,
                index: number,
                array: CBSocketCallbackHolderMessageDescriptor[]
            ) {
                
                if ((descriptor.completionPolicy == CBSocketClient.completionPolicy.last ||
                    descriptor.completionPolicy == CBSocketClient.completionPolicy.firstAndLast) &&
                    allResponsesReceived && !message.keepWaitingForResponses) {
                    
                    // Calling completionFunction
                    callCompletionFunction(descriptor, !message.keepWaitingForResponses)
                    
                    // Cleaning up
                    descriptorsForKey.removeElement(descriptor)
                    
                }
                else if (descriptor.completionPolicy == CBSocketClient.completionPolicy.firstAndLastIfDifferent &&
                    allResponsesReceived && !message.keepWaitingForResponses) {
                    
                    // Calling completionFunction if needed
                    if (descriptor.responseDataHash != responseDataHash) {
                        
                        callCompletionFunction(descriptor, !message.keepWaitingForResponses)
                        
                    }
                    
                    // Cleaning up
                    descriptorsForKey.removeElement(descriptor)
                    
                }
                
            }.bind(this))
            
            
        }
        
        
        
        
        
    }
    
    
    
    
    
}






























import * as SocketIo from "socket.io"

import { SocketController, SocketServerClient } from "./SocketController"
import objectHash = require("object-hash")
import Utils from "./Utils"





export default class SocketSession {
    
    
    
    private _socket: SocketIo.Socket
    private _userProfile: any
    
    private completionFunctions: {
        [x: string]: CBSocketMessagecompletionFunction;
    } = {}
    
    private handlers: {
        [x: string]: CBSocketMessageHandlerFunction[];
    } = {}
    
    private responseGroups: {
        
        [x: string]: {
            
            numberOfKeys: number;
            inResponseToMessage: CBSocketMessage;
            multipleMessageObjects: CBSocketMultipleMessageObject[];
            
            responseKey: string;
            
        }
        
    } = {}
    
    private _isHandlingMultipleMessages = false
    private _multipleMessageResponses: CBSocketMultipleMessageObject[] = []
    
    private responseMessageKey = "CBSocketResponseMessage"
    private multipleMessageKey = "CBSocketMultipleMessage"
    _accessToken?: string | null
    
    socketClient: SocketServerClient
    
    
    constructor(
        socket: SocketIo.Socket,
        accessToken?: string | null,
        userProfile?: any,
        
    ) {
        
        this._socket = socket
        this._accessToken = accessToken
        this._userProfile = userProfile
        
        
        const messageTargets = SocketController.sharedInstance.messageTargets as Object
        
        messageTargets.forEachValue((value, key, index) =>
            this._socket.on(key, (message: CBSocketMessage<any>) => {
                this.didReceiveMessageForKey(key, message)
            })
        )
        
        this._socket.on(this.responseMessageKey, (message: CBSocketMessage) => {
            
            this.didReceiveMessageForKey(this.responseMessageKey, message)
            
        })
        
        this._socket.on(this.multipleMessageKey, (message: CBSocketMultipleMessage) => {
            
            this._isHandlingMultipleMessages = true
            
            const responseGroup = {
                
                numberOfKeys: message.messageData.length,
                inResponseToMessage: message,
                multipleMessageObjects: [],
                
                responseKey: this.multipleMessageKey
                
            }
            
            message.messageData.forEach((messageObject, index, array) => {
                
                if (message.shouldGroupResponses) {
                    
                    this.responseGroups[messageObject.message.identifier] = responseGroup
                    
                }
                
                this.didReceiveMessageForKey(messageObject.key, messageObject.message)
                
            })
            
            if (this._multipleMessageResponses.length) {
                
                this._sendMessageForKey(this.multipleMessageKey, this._multipleMessageResponses, message)
                
                this._multipleMessageResponses = []
                
            }
            
            this._isHandlingMultipleMessages = false
            
        })
        
        
        
        const session = this
        
        this.socketClient = new Proxy({ "name": "SocketClient" }, {
            
            get(target, key) {
                
                var result = (messageData: any, resultInsteadOfError: any) => session.resultForMessageForKey(
                    key as string,
                    messageData,
                    resultInsteadOfError
                )
                
                return result
                
            }
            
        }) as any
        
        
        
        
        
    }
    
    
    
    
    
    
    get isOpen() {
        
        return this.socket.connected
        
    }
    
    get socket() {
        
        return this._socket
        
    }
    
    get accessToken() {
        
        return this._accessToken
        
    }
    
    get userProfile() {
        
        return this._userProfile
        
    }
    
    
    
    
    
    
    
    
    triggerReconnection(message?: string) {
        
        this.socket.emit("CBPerformReconnect", message)
        
        this.socket.disconnect(true)
        
        this.socketDidDisconnect()
        
    }
    
    
    
    
    
    sendMessageForKey<MessageType = any>(
        key: string,
        message: MessageType,
        completion?: CBSocketMessagecompletionFunction
    ) {
        
        this._sendMessageForKey(key, message, undefined, false, false, false, undefined, completion)
        
    }
    
    
    
    
    
    async _sendMessageForKey(
        key: string,
        message: any,
        inResponseToMessage: CBSocketMessage<any> = {} as any,
        keepMessageConnectionOpen = false,
        shouldStoreResponse = false,
        shouldUseStoredResponse = false,
        responseValidityDuration?: number,
        completion?: CBSocketMessagecompletionFunction
    ): Promise<string | undefined> {
        
        if (message instanceof Promise) {
            try {
                message = await message
            } catch (exception) {
                message = exception
            }
        }
        
        if (this.isOpen) {
            
            var identifier = Utils.makeID()
            
            if (completion) {
                this.completionFunctions[identifier] = completion
            }
            
            var messageDataHash
            
            try {
                
                messageDataHash = objectHash(JSON.parse(JSON.stringify(message)))
                
            } catch (exception) {
                
                var asd = 1
                
            }
            
            
            var messageObject: CBSocketMessage<any> = {
                
                messageData: message,
                identifier: identifier,
                inResponseToIdentifier: inResponseToMessage.identifier,
                keepWaitingForResponses: keepMessageConnectionOpen,
                
                canBeStoredAsResponse: shouldStoreResponse,
                useStoredResponse: shouldUseStoredResponse,
                
                responseValidityDuration: responseValidityDuration,
                
                messageDataHash: messageDataHash
                
            }
            
            if (inResponseToMessage.storedResponseHash == messageDataHash && messageDataHash) {
                
                messageObject = {
                    
                    messageData: "",
                    identifier: identifier,
                    inResponseToIdentifier: inResponseToMessage.identifier,
                    keepWaitingForResponses: keepMessageConnectionOpen,
                    
                    canBeStoredAsResponse: shouldStoreResponse,
                    useStoredResponse: true,
                    
                    responseValidityDuration: responseValidityDuration,
                    
                    messageDataHash: messageDataHash
                    
                }
                
                
                
            }
            
            
            const responseGroup = this.responseGroups[inResponseToMessage.identifier]
            
            if (responseGroup && !keepMessageConnectionOpen) {
                
                
                
                
                responseGroup.multipleMessageObjects.push({
                    
                    key: key,
                    message: messageObject
                    
                })
                
                
                this.removeMessageFromResponseGroups(inResponseToMessage)
                
                
                if (responseGroup.numberOfKeys == responseGroup.multipleMessageObjects.length) {
                    
                    this._sendMessageForKey(
                        responseGroup.responseKey,
                        responseGroup.multipleMessageObjects,
                        responseGroup.inResponseToMessage
                    )
                    
                }
                
                
                
            }
            else if (this._isHandlingMultipleMessages && key != this.multipleMessageKey) {
                
                this._multipleMessageResponses.push({
                    
                    key: key,
                    message: messageObject
                    
                })
                
            }
            else {
                
                this._socket.emit(key, messageObject)
                
            }
            
            
            
            return messageDataHash
            
        }
        
        
    }
    
    
    
    
    
    private removeMessageFromResponseGroups(inResponseToMessage: CBSocketMessage<any>) {
        delete this.responseGroups[inResponseToMessage.identifier]
    }
    
    private async didReceiveMessageForKey(
        key: string,
        message: CBSocketMessage<any>,
        responseGroupIdentifier?: string
    ) {
        
        // Only handle messages if connected to client
        if (!this.isOpen) {
            return
        }
        
        const startTime = Date.now()
        
        // associate responses with response group
        
        
        var clientShouldStoreResponse = false
        var clientStoredResponseTriggered = false
        var useStoredResponseWithErrorResponse = false
        
        var responseSent = false
        
        var responseValidityDuration: number | undefined
        
        const userID = "" + this.userProfile.contactInformation.email
        
        const excludedKeys = [
            "RetrieveDropdownDataForCode",
            "RetrieveNumberOfChatNotificationsForCBInquiryWithID",
            "RetrieveNumberOfOfferNotificationsForCBInquiryWithID"
        ]
        
        var excludeMessage = false
        
        var sendResponseFunction: CBSocketMessageSendResponseFunction = function (
            this: SocketSession,
            responseMessage: any,
            completion: CBSocketMessagecompletionFunction
        ): Promise<string> {
            
            if (clientStoredResponseTriggered) {
                // @ts-ignore
                return
            }
            
            if (!excludedKeys.contains(key) && !excludeMessage) {
                
                console.log("Sending response for message with key \"" + key + "\" after " + (Date.now() - startTime) +
                    "ms. User ID is " + userID + ".")
                
            }
            
            responseSent = true
            
            // @ts-ignore
            return this._sendMessageForKey(
                this.responseMessageKey,
                responseMessage,
                message,
                false,
                clientShouldStoreResponse,
                false,
                responseValidityDuration,
                completion
            )
            
        }.bind(this) as any
        
        sendResponseFunction.sendIntermediateResponse = function (
            this: SocketSession,
            updateMessage: any,
            completion: CBSocketMessagecompletionFunction
        ) {
            
            if (clientStoredResponseTriggered) {
                return
            }
            
            if (!excludedKeys.contains(key) && !excludeMessage) {
                
                console.log("Sending intermediate response for message with key \"" + key + "\" after " +
                    (Date.now() - startTime) + "ms. User ID is " + userID + ".")
                
            }
            
            responseSent = true
            
            this._sendMessageForKey(
                this.responseMessageKey,
                updateMessage,
                message,
                true,
                false,
                false,
                responseValidityDuration,
                completion
            )
            
        }.bind(this)
        
        sendResponseFunction.sendErrorResponse = function (
            this: SocketSession,
            errorMessage: any,
            completion: CBSocketMessagecompletionFunction
        ) {
            
            if (clientStoredResponseTriggered) {
                return
            }
            
            var isException = false
            
            if (errorMessage.message && errorMessage.stack) {
                
                isException = true
                
                errorMessage = {
                    
                    message: errorMessage.message,
                    stack: errorMessage.stack
                    
                }
                
            }
            
            if (!excludedKeys.contains(key) && !excludeMessage) {
                
                console.log("Sending error response for message with key \"" + key + "\" after " +
                    (Date.now() - startTime) + "ms. User ID is " + userID + ".")
                
                console.log(errorMessage)
                
                
                
            }
            
            const messageObject = {
                
                _isCBSocketErrorMessage: true,
                messageData: errorMessage
                
            }
            
            responseSent = true
            
            this._sendMessageForKey(
                this.responseMessageKey,
                messageObject,
                message,
                false,
                false,
                useStoredResponseWithErrorResponse,
                responseValidityDuration,
                completion
            )
            
        }.bind(this)
        
        sendResponseFunction.confirmStoredResponseHash = function (
            this: SocketSession,
            responseHash: any,
            completion?: CBSocketMessagecompletionFunction
        ) {
            
            clientShouldStoreResponse = true
            
            if (responseHash == message.storedResponseHash && message.storedResponseHash) {
                
                if (!excludedKeys.contains(key) && !excludeMessage) {
                    
                    console.log("Sending stored hash response for message with key \"" + key + "\" after " +
                        (Date.now() - startTime) + "ms. User ID is " + userID + ".")
                    
                    
                    
                }
                
                responseSent = true
                
                this._sendMessageForKey(
                    this.responseMessageKey,
                    "",
                    message,
                    false,
                    false,
                    true,
                    responseValidityDuration,
                    completion
                )
                
                clientStoredResponseTriggered = true
                
                return true
                
            }
            
            return false
            
        }.bind(this)
        
        
        
        sendResponseFunction.excludeMessageFromAutomaticConnectionEvents = () => {
            
            excludeMessage = true
            
        }
        
        sendResponseFunction.setResponseValidityDuration = function (duration) {
            
            responseValidityDuration = duration
            
        }
        
        sendResponseFunction.useStoredResponseWithErrorResponse = function (this: SocketSession) {
            
            useStoredResponseWithErrorResponse = true
            
        }.bind(this)
        
        if (this.handlers[key]) {
            
            await Promise.all(this.handlers[key].map((handler, index, array) => handler(
                message.messageData,
                sendResponseFunction
            )))
            
        }
        
        const messageTarget = SocketController.sharedInstance.messageTargets[key]
        
        if (messageTarget) {
            
            await messageTarget(message.messageData, this, sendResponseFunction)
            
        }
        
        if (!responseSent) {
            
            excludeMessage = true
            
            sendResponseFunction(null)
            
        }
        
        if (message.inResponseToIdentifier && this.responseMessageKey == key) {
            
            var completionFunction = this.completionFunctions[message.inResponseToIdentifier]
            
            if (completionFunction) {
                
                try {
                    
                    completionFunction(message.messageData, sendResponseFunction)
                    
                } catch (exception) {
                    
                    sendResponseFunction.sendErrorResponse(exception)
                    
                }
                
                if (!message.keepWaitingForResponses) {
                    
                    delete this.completionFunctions[message.inResponseToIdentifier]
                    
                }
                
            }
            
        }
        
    }
    
    
    
    
    
    async resultForMessageForKey<ResultType = any>(key: string, message: any, resultInsteadOfError: ResultType | null = null) {
        
        var result = resultInsteadOfError
        var errorResult
        
        var sendResponseFunction: CBSocketMessageSendResponseFunction = function (
            this: SocketSession,
            responseMessage: any,
            completion: CBSocketMessagecompletionFunction
        ) {
            
            result = responseMessage
            
        }.bind(this) as any
        
        sendResponseFunction.sendIntermediateResponse = function (
            this: SocketSession,
            updateMessage: any,
            completion: CBSocketMessagecompletionFunction
        ) {
            
            result = updateMessage
            
        }.bind(this)
        
        sendResponseFunction.sendErrorResponse = function (
            this: SocketSession,
            errorMessage: any,
            completion: CBSocketMessagecompletionFunction
        ) {
            
            errorResult = errorMessage
            
        }.bind(this)
        
        sendResponseFunction.confirmStoredResponseHash = function (
            this: SocketSession,
            responseHash: any,
            completion?: CBSocketMessagecompletionFunction
        ) {
            
            return false
            
        }.bind(this)
        
        sendResponseFunction.excludeMessageFromAutomaticConnectionEvents = Utils.nil
        
        sendResponseFunction.setResponseValidityDuration = Utils.nil
        
        sendResponseFunction.useStoredResponseWithErrorResponse = Utils.nil
        
        if (this.handlers[key]) {
            
            await Promise.all(this.handlers[key].map(async (handler, index, array) => {
                
                // Handlers can actually return promises
                await handler(message, sendResponseFunction)
                
            }))
            
        }
        
        const messageTarget = SocketController.sharedInstance.messageTargets[key]
        
        if (messageTarget) {
            
            await messageTarget(message, this, sendResponseFunction)
            
        }
        
        if (errorResult && resultInsteadOfError == null) {
            
            throw errorResult
            
        }
        
        return result
        
    }
    
    
    
    
    
    addTargetForMessagesForKey<ResponseMessageType = any>(
        key: string,
        handlerFunction: CBSocketMessageHandlerFunction<ResponseMessageType>
    ) {
        
        if (!this.handlers[key]) {
            this.handlers[key] = [] as any
        }
        
        this.handlers[key].push(handlerFunction)
        
        this._socket.on(key, (message: CBSocketMessage<any>) => {
            this.didReceiveMessageForKey(key, message)
        })
        
    }
    
    
    
    
    
    socketDidDisconnect() {
        
        // this.instanceConnection.disconnectDate = new Date()
        //
        // this.instanceConnection.save()
        //
        // this.addInstanceConnectionEvent("disconnect")
        
    }
    
    
    
    
    
}






























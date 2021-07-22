///<reference path="../webclient/scripts/SharedTypes.d.ts"/>

import SocketSession from "./SocketSession"
import Utils from "./Utils"





type ThenArg<T> = T extends PromiseLike<infer U> ? U : T;

type PropType<TObj, TProp extends keyof TObj> = TObj[TProp];

type SocketServerFunction<MessageType, ResultType> = (
    messageData: MessageType,
    socketSession: SocketSession,
    respondWithMessage: CBSocketMessageSendResponseFunction<ResultType>
) => void;

type SocketServerFunctionFromClientFunction<ClientFunctionType extends (...args: any) => any> = SocketServerFunction<Parameters<ClientFunctionType>[0], PropType<ThenArg<ReturnType<ClientFunctionType>>, "result"> | PromiseLike<PropType<ThenArg<ReturnType<ClientFunctionType>>, "result">>>;

type ServerFunctionsFromClientFunctions<T> = {
    
    [P in keyof T]: SocketServerFunctionFromClientFunction<T[P] extends (...args: any) => any ? T[P] : never>;
    
}

type SocketServer = ServerFunctionsFromClientFunctions<SocketClient>;





type SocketServerClientFunction<MessageType, ResultType> = (
    message: MessageType,
    resultInsteadOfError?: ResultType
) => ResultType;

type SocketServerClientFunctionFromClientFunction<ClientFunctionType extends (...args: any) => any> = SocketServerClientFunction<Parameters<ClientFunctionType>[0], PropType<ThenArg<ReturnType<ClientFunctionType>>, "result"> | PromiseLike<PropType<ThenArg<ReturnType<ClientFunctionType>>, "result">>>;

type ServerClientFunctionsFromClientFunctions<T> = {
    
    [P in keyof T]: SocketServerClientFunctionFromClientFunction<T[P] extends (...args: any) => any ? T[P] : never>;
    
}

export type SocketServerClient = ServerClientFunctionsFromClientFunctions<SocketClient>;

// interface SocketServer<MessageType, ResultType> {

//     [x: string]: SocketServerFunction<MessageType, ResultType>;

// }


export class SocketController {
    
    
    
    
    private static _instance: SocketController
    
    
    
    private _connectionEventTargets: ((socketSession: SocketSession) => void)[] = []
    
    messageTargets: SocketServer = {} as any
    
    private _sessions: { [x: string]: SocketSession[] } & Object = {}
    
    
    
    
    
    constructor(private _socketIO: SocketIO.Server) {
        
        
        
        
        
        _socketIO.on("connect", (socket: SocketIO.Socket) => {
            
            //console.log('Socket.io client connected. clientID = ' + socket.client.id);
            
            var currentSession: SocketSession | null = null
            
            
            socket.on(
                "CBSocketHandshakeInitMessage",
                async (message: CBSocketMessage<CBSocketHandshakeInitMessage>) => {
                    
                    
                    try {
                        
                        if (!message || !message.messageData || !message.identifier) {
                            
                            console.log("Invalid socket handshake message. " + JSON.stringify(message))
                            
                            this._socketIO.to(socket.id).emit("CBSocketHandshakeResponseMessage", {
                                
                                identifier: Utils.makeID(),
                                inResponseToIdentifier: message.identifier,
                                messageData: {
                                    
                                    accepted: false,
                                    message: "Invalid handshake message."
                                    
                                }
                                
                            })
                            
                            return
                            
                        }
                        
                        
                        var response: CBSocketMessage<CBSocketHandshakeResponseMessage> = {
                            
                            identifier: Utils.makeID(),
                            inResponseToIdentifier: message.identifier,
                            messageData: { accepted: true }
                            
                        }
                        
                        
                        var response: CBSocketMessage<CBSocketHandshakeResponseMessage> = {
                            
                            identifier: Utils.makeID(),
                            inResponseToIdentifier: message.identifier,
                            messageData: {
                                accepted: true
                            }
                            
                        }
                        
                        
                            
                            // Guest user
                            
                            console.log("Starting new socket session for a guest user.")
                            
                            currentSession = new SocketSession(
                                socket,
                                null,
                                Utils.nil
                            )
                            
                            this._addSocketSession(currentSession)
                            
                            
                            this._connectionEventTargets.forEach((target, index, array) => {
                                
                                // @ts-ignore
                                target(currentSession)
                                
                            })
                            
                            // This removes the userProfile value from the response
                            response.messageData = { accepted: true }
                            
                        
                        
                        this._socketIO.to(socket.id).emit("CBSocketHandshakeResponseMessage", response)
                        
                    } catch (exception) {
                        
                        console.log("EXCEPTION DURING SOCKET HANDSHAKE.")
                        
                        console.log(exception)
                        
                        this._socketIO.to(socket.id).emit("CBSocketHandshakeResponseMessage", {
                            
                            identifier: Utils.makeID(),
                            inResponseToIdentifier: message.identifier,
                            messageData: {
                                
                                accepted: false,
                                message: "Handshake failed with an exception."
                                
                            }
                            
                        })
                        
                    }
                    
                    
                    
                    
                    
                }
            )
            
            
            
            
            
            socket.on("disconnect", () => {
                
                console.log("Socket.io client disconnected")
                
                if (currentSession) {
                    currentSession.socketDidDisconnect()
                    
                    // Remove the current session from _sessions
                    this._removeSocketSession(currentSession)
                    
                    currentSession = null
                }
                
            })
            
        })
        
        
        
        
        
    }
    
    
    reconnectSession(session: SocketSession, message?: string) {
        
        session.triggerReconnection(message)
        
        this._removeSocketSession(session)
        
    }
    
    
    private _addSocketSession(session: SocketSession) {
        
        const userIDKey = "" + session.userProfile._id
        
        if (!this._sessions[userIDKey]) {
            this._sessions[userIDKey] = []
        }
        this._sessions[userIDKey].push(session)
        
        
        
    }
    
    
    private _removeSocketSession(session: SocketSession) {
        
        const userIDKey = "" + session.userProfile._id
        
        var userSocketSessions = this.socketSessionsForUserWithID(userIDKey)
        var index = userSocketSessions.indexOf(session)
        userSocketSessions.splice(index, 1)
        
        if (userSocketSessions.length == 0) {
            
            delete this._sessions[userIDKey]
            
            
            
        }
        
        
        
    }
    
    
    
    
    
    get connectedUserIDStrings() {
        
        return Object.keys(this._sessions)
        
    }
    
    get activeSessions() {
        
        const result: SocketSession[] = []
        
        this._sessions.forEachValue(
            (sessions: SocketSession[], key, index) => sessions.forEach(
                (session, index, array) => result.push(session)
            )
        )
        
        return result
        
    }
    
    
    
    sendMessageForKeyToUserWithID<MessageType = any>(key: string, message: MessageType, userID: string) {
        this.socketSessionsForUserWithID(userID).forEach(async (session, index, array) => {
            session.sendMessageForKey<MessageType>(key, message)
        })
    }
    
    
    
    socketSessionsForUserWithID(userID: string = "asdasd") {
        
        var result = this._sessions[userID]
        if (!result) {
            result = []
        }
       
        return result
        
    }
    
    
    
    
    addTargetForConnectionEvents(targetFunction: (socketSession: SocketSession) => void) {
        
        this._connectionEventTargets.push(targetFunction)
        
    }
    
    
    
    
    
    public static get sharedInstance() {
        return this._instance
    }
    
    public static Instance(socketIO: SocketIO.Server) {
        return this._instance || (this._instance = new this(socketIO))
    }
    
    
    
    
    
}










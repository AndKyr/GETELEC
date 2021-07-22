
interface SocketClient {
    
    [x: string]: SocketClientFunction<any, any>;
    
    PerformFitFun: SocketClientFunction<{ inputData: string }, string>;
    
    //RetrieveVehicleDetailsForCurrentUser: SocketClientNoMessageFunction<CBVehicleDetails[]>;
    
}


interface SocketClientResult<ResultType> {
    
    responseMessage: any;
    result: ResultType;
    errorResult: any;
    respondWithMessage: CBSocketMessageSendResponseFunction;
    
}


type SocketClientFunction<MessageType, ResultType> = (
    messageData: MessageType,
    completionPolicy?: string,
    isUserBound?: boolean,
    intermediateResultFunction?: any
) => Promise<SocketClientResult<ResultType>>

type SocketClientNoMessageFunction<ResultType> = (
    messageData?: null,
    completionPolicy?: string,
    isUserBound?: boolean
) => Promise<SocketClientResult<ResultType>>;





interface CBSocketMultipleMessageObject<MessageDataType = any> {
    
    key: string;
    message: CBSocketMessage<MessageDataType>;
    
}


// CBSocket communication messages
interface CBSocketMessage<MessageDataType = any> {
    
    identifier: string;
    inResponseToIdentifier?: string;
    keepWaitingForResponses?: boolean;
    
    messageData: MessageDataType;
    
    // This is sent from client to server with requests
    storedResponseHash?: string;
    
    // This is always present on messages sent from the server side
    messageDataHash?: string;
    
    // This tells the client to store this message for future use
    canBeStoredAsResponse?: boolean;
    
    // This tells the client to use the previously stored response
    useStoredResponse?: boolean;
    
    // This tells the client that the response is valid for at least this long in ms
    responseValidityDuration?: number;
    
}


interface CBSocketMultipleMessage extends CBSocketMessage<CBSocketMultipleMessageObject[]> {
    
    shouldGroupResponses: boolean;
    
}


type CBSocketMessageSendResponseFunctionBase<ResponseMessageType> = (
    responseMessage: ResponseMessageType,
    completion?: CBSocketMessagecompletionFunction
) => Promise<string>;

type CBSocketMessagecompletionFunction = (
    responseMessage: any,
    respondWithMessage: CBSocketMessageSendResponseFunction
) => void;
type CBSocketMessageHandlerFunction<ResponseMessageType = any> = (
    message: any,
    respondWithMessage: CBSocketMessageSendResponseFunction<ResponseMessageType>
) => void;

type CBSocketMultipleMessagecompletionFunction = (responseMessages: any[], callcompletionFunctions: () => void) => void;


interface CBSocketMessageSendResponseFunction<ResponseMessageType = any> extends CBSocketMessageSendResponseFunctionBase<ResponseMessageType> {
    respondingToMainResponse: boolean;
    
    excludeMessageFromAutomaticConnectionEvents: () => void;
    
    // @ts-ignore
    setResponseValidityDuration(duration: number);
    
    // @ts-ignore
    useStoredResponseWithErrorResponse();
    
    // @ts-ignore
    sendErrorResponse(message?: any, completion?: CBSocketMessagecompletionFunction);
    
    // @ts-ignore
    sendIntermediateResponse(updateMessage: any, completion?: CBSocketMessagecompletionFunction);
    
    // This tells the client to use the stored response if responseHash matches and also enables storing of responses in the client in the first place
    // Returns true if the hash matched
    confirmStoredResponseHash(responseHash: string, completion?: CBSocketMessagecompletionFunction): boolean;
    
}



// Socket handshake messages
interface CBSocketHandshakeInitMessage {
    
    accessToken?: string;
    
    instanceIdentifier: string;
    
}


interface CBSocketHandshakeResponseMessage {
    
    accepted: boolean;
    
}


type TypeWithoutKey<Type, Key> = Pick<Type, Exclude<keyof Type, Key>>;

type TypeWithoutID<Type> = TypeWithoutKey<Type, "_id">;

type Diff<T extends keyof any, U extends keyof any> =
    ({ [P in T]: P } & { [P in U]: never } & { [x: string]: never })[T];

type Overwrite<T, U> = Pick<T, Diff<keyof T, keyof U>> & U;

type RecursivePartial<T> = {
    [P in keyof T]?:
    T[P] extends (infer U)[] ? RecursivePartial<U>[] :
        T[P] extends object ? RecursivePartial<T[P]> :
            T[P];
};


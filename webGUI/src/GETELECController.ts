
import Utils from "../lib/Utils"
import { readdirSync, writeFileSync } from "fs"
import { SocketController } from "../lib/SocketController"
import * as path from "path"
import { Express } from "express"
import { Diagnostic } from "typescript"

const needle = require('needle');




export class GETELECController {
    
    private static _instance: GETELECController
    expressApplication: Express
    
    constructor(expressApplication: Express) {
        
        this.expressApplication = expressApplication
        
        const socketController = SocketController.sharedInstance
        
        // This is for testing
        expressApplication.post("/fit_fun", (request, response) => {
        
            console.log(request);
            
        })
        
        const messageTargets = socketController.messageTargets
        
        messageTargets.PerformFitFun = async (requestObject, socketSession, respondWithMessage) => {
            
            try {
                
                let result = await needle('post', 'http://172.17.207.205:8000/fit_fun', requestObject.inputData, { json: true })
                //let result = await needle('post', 'http://localhost:7005/fit_fun', requestObject.inputData, {json: true})
                
                respondWithMessage(result.body)
                
            } catch (exception) {
                
                respondWithMessage.sendErrorResponse(exception);
                
            }
            
        }
        
    }
    
    
    public static get sharedInstance() {
        return this._instance
    }
    
    // noinspection TypeScriptUMDGlobal,JSDeprecatedSymbols
    public static Instance(expressApplication: Express) {
        
        return this._instance || (this._instance = new this(expressApplication))
        
    }
    
    
}














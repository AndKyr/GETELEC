import { parseStringPromise, Builder, Parser } from "xml2js";



import objectHash = require("object-hash");
import { pbkdf2, randomBytes } from "crypto";
import { promisify } from "util";
import { request } from "http";




var YES = true;
var NO = false;





type UIIFBlockReceiver<T> = (functionToCall: () => any) => UIIFEvaluator<T>;

type UIIFEvaluatorBase<T> = () => T;

interface UIIFEvaluator<T> extends UIIFEvaluatorBase<T> {
    
    ELSE_IF: (otherValue: any) => UIIFBlockReceiver<T>;
    ELSE: (functionToCall: () => any) => T;
    
}





export default class Utils {
    
    static _builderForXMLStrings = new Builder();
    
    private static _numberOfPasswordHashIterations = 10000;
    private static _passwordKeyLength = 551;
    private static _passwordHashingMethod = "sha512";
    private static _passwordHashSeparator = "_";
    
    
    
    
    // Make an identifier with a random part and the current timestamp
    static makeID(randomPartLength = 16) {
        var result = "";
        var characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
        for (var i = 0; i < randomPartLength; i++) {
            result = result + characters.charAt(Math.floor(Math.random() * characters.length));
        }
        result = result + Date.now();
        return result;
    }
    
    
    static makeKey(randomPartLength = 16) {
        
        var result = randomBytes(randomPartLength).toString("base64");
        
        return result;
        
    }
    
    
    static async hashFromPassword(password: string) {
      
        const salt = randomBytes(128).toString("base64");
        
        const numberOfIterations = Utils._numberOfPasswordHashIterations;
        
        const actualHash = (await promisify(pbkdf2)(password, salt, numberOfIterations, Utils._passwordKeyLength, Utils._passwordHashingMethod)).toString("base64");
        
        const result = [numberOfIterations, Utils._passwordKeyLength, Utils._passwordHashingMethod, salt, actualHash].join(Utils._passwordHashSeparator);
        
        return result;
        
    }
    
    static async checkIfPasswordMatchesHash(password: string, hash: string) {
      
        const components = hash.split(Utils._passwordHashSeparator);
        
        const actualHashFromPassword = await promisify(pbkdf2)(password, components[3], parseInt(components[0], 10), parseInt(components[1], 10), components[2]);
        
        const result = (components[4] === actualHashFromPassword.toString("base64"));
        
        return result;
        
    }
    
    
    
    
    
    
    
    static async objectFromXMLString(string: string) {
        
        const result = await parseStringPromise(string, { explicitArray: false });
        
        return result;
        
    }
    
    static stringInXMLFormatFromObject(object: any) {
        
        const result = Utils._builderForXMLStrings.buildObject(object);
        
        return result;
        
    }
    
    
    static getDataFromHostAndPath(host: string, path: string) {
        
        var options = {
            host: host,
            path: path
        }
        
        var result = new Promise<string>((resolve, reject) => {
            
            var requestObject = request(options, function (res) {
            
                var data = "";
                
                res.on("data", function (chunk) {
                    
                    data += chunk;
                    
                })
                
                res.on("end", function () {
                    
                    
                    resolve(data);
                    
                    
                })
                
            })
            
            requestObject.on("error", function (error) {
                
                
                reject(error);
                
                
            })
            
            requestObject.end();
            
            
            
        })
        
        
        return result;
        
        
    }
    
    
    
    
    
    static getDataFromPostingToHostAndPath(host: string, path: string, requestData: string, port?: number | string) {
        
        var options = {
            
            method: "POST",
            host: host,
            path: path,
            port: port,
            headers: {
                
                "Content-Type": "text/plain",
                "Content-Length": Buffer.byteLength(requestData)
                
            }
            
        }
        
        var result = new Promise<string>((resolve, reject) => {
            
            
            
            var requestObject = request(options, function (res) {
            
                var data = "";
                
                res.on("data", function (chunk) {
                    
                    data += chunk;
                    
                })
                
                res.on("end", function () {
                    
                    
                    resolve(data);
                    
                    
                })
                
            })
            
            requestObject.on("error", function (error) {
                
                
                reject(error);
                
                
            })
            
            requestObject.write(requestData);
            
            requestObject.end();
            
            
            
        })
        
        
        return result;
        
        
    }
    
    
    
    
    
    
    
    
    static callFunctionEveryHour(functionToCall: Function) {
        let running = true
        let nextHour = () => {
          return 3600000 - new Date().getTime() % 3600000
        }
        let nextCall = setTimeout(() => {
          functionToCall()
          Utils.callFunctionEveryHour(functionToCall)
        }, nextHour())
        return {
          next() { return running ? nextHour() : -1 },
          exec() { functionToCall() },
          stop() {
            clearTimeout(nextCall)
            running = false
          },
          start() {
            clearTimeout(nextCall)
            nextCall = setTimeout(() => {
              functionToCall()
              Utils.callFunctionEveryHour(functionToCall)
            }, nextHour())
            running = true
          }
        }
      }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    static NilFunction() { return Utils.nil }
    
    static nil: any = new Proxy(Object.assign(Utils.NilFunction, { "className": "Nil" }), {
        
        get(target, name) {
            
            if (name == Symbol.toPrimitive) {
                
                return function (hint: string) {
                    if (hint == "number") {
                        return 0;
                    }
                    if (hint == "string") {
                        return "";
                    }
                    return false;
                }
                
            }
            
            return Utils.NilFunction();
            
        },
        
        set(target, name, value) {
            return Utils.NilFunction();
        }
        
    });
    
    
    
    
    
    
    static IS(object: any) {
        
        if (object && object !== Utils.nil) {
            
            return YES;
            
        }
        
        return NO;
        
        //return (object != nil && object);
        
    }
    
    static IS_NOT(object: any) {
        
        return !Utils.IS(object);
        
    }
    
    static IS_DEFINED(object: any) {
        
        if (object != undefined) {
            
            return YES;
            
        }
        
        return NO;
        
    }
    
    static IS_UNDEFINED(object: any) {
        
        return !Utils.IS_DEFINED(object);
        
    }
    
    static IS_NIL(object: any) {
        
        if (object === Utils.nil) {
            
            return YES;
            
        }
        
        return NO;
        
    }
    
    static IS_NOT_NIL(object: any) {
        
        return !Utils.IS_NIL(object);
        
    }
    
    
    static IS_LIKE_NULL(object: any) {
        
        return (Utils.IS_UNDEFINED(object) || Utils.IS_NIL(object) || object == null);
        
    }
    
    static IS_NOT_LIKE_NULL(object: null | undefined) {
        
        return !Utils.IS_LIKE_NULL(object);
        
    }
    
    
    static ONLY_IF_IS(value: any) {
        
        if (Utils.IS(value)) {
            
            return value;
            
        }
        
    }
    
    
    static IS_AN_EMAIL_ADDRESS(email: string) {
        var re = /\S+@\S+\.\S+/;
        return re.test(email);
    }
    
    
    static FIRST_OR_NIL<T>(...objects: (T | undefined)[]): T {
        
        var result = objects.find(function (object, index, array) {
            
            return Utils.IS(object);
            
        });
        
        return result || Utils.nil;
        
    }
    
    static FIRST<T>(...objects: T[]): T {
        
        var result = objects.find(function (object, index, array) {
            
            return Utils.IS(object);
            
        });
        
        return result || Utils.nil;
        
    }
    
    
    static RETURNER<T>(value: T) {
        
        return function (...objects: any[]) {
            
            return value;
            
        }
        
    }
    
    
    
    
    
    static IF<T = any>(value: any): UIIFBlockReceiver<T> {
        
        var thenFunction = Utils.nil;
        var elseFunction = Utils.nil;
        
        var result: any = function (functionToCall: () => T) {
            thenFunction = functionToCall;
            return result.evaluateConditions;
        }
        
        
        result.evaluateConditions = function () {
            if (Utils.IS(value)) {
                return thenFunction();
            }
            return elseFunction();
        }
        
        
        result.evaluateConditions.ELSE_IF = function (otherValue: any) {
            
            var functionResult = Utils.IF(otherValue) as (UIIFBlockReceiver<T> & { evaluateConditions: UIIFEvaluator<T> })
            elseFunction = functionResult.evaluateConditions;
            
            var functionResultEvaluateConditionsFunction: any = function () {
                return result.evaluateConditions();
            }
            functionResultEvaluateConditionsFunction.ELSE_IF = functionResult.evaluateConditions.ELSE_IF
            functionResultEvaluateConditionsFunction.ELSE = functionResult.evaluateConditions.ELSE
            
            functionResult.evaluateConditions = functionResultEvaluateConditionsFunction;
            
            return functionResult;
            
        }
        
        
        result.evaluateConditions.ELSE = function (functionToCall: () => T) {
            elseFunction = functionToCall;
            return result.evaluateConditions();
        }
        
        
        return result;
    }
    
    
    
    
    
    
    
    
    
    
    
}











/// <reference path="./Custom components/CBSocketClient.ts" />





class CBCore extends UIObject {




    private static _sharedInstance: CBCore

    _socketClient: CBSocketClient = new CBSocketClient(this)


    _functionsToCallForEachSocketClient: (() => void)[] = []





    constructor() {

        super()

        this._class = CBCore
        this.superclass = UIObject





    }





    static get sharedInstance() {
        if (!CBCore._sharedInstance) {
            CBCore._sharedInstance = new CBCore()
        }
        return CBCore._sharedInstance
    }




    broadcastMessageInRootViewTree(message: UIViewBroadcastEvent) {

        (UICore.main.rootViewController as UIViewController).view.broadcastEventInSubtree(message)

    }




    get socketClient() {
        return this._socketClient
    }






    set languageKey(languageKey: string) {

        if (IS_NOT(languageKey)) {

            localStorage.removeItem("CBLanguageKey")

        }

        localStorage.setItem("CBLanguageKey", JSON.stringify(languageKey))

        this.didSetLanguageKey()

    }

    get languageKey() {

        const result = FIRST(localStorage.getItem("CBLanguageKey"), LanguageService.defaultLanguageKey).replace(
            "\"",
            ""
        ).replace("\"", "")


        return result

    }

    didSetLanguageKey() {

        UIRoute.currentRoute.routeWithComponent(
            "settings",
            { "language": this.languageKey },
            YES
        ).applyByReplacingCurrentRouteInHistory()

    }
    
    
    readTextFileFromURL(fileURL: string) {

        var result = new Promise<string>((resolve, reject) => {

            var rawFile = new XMLHttpRequest();
            rawFile.open("GET", fileURL, false);
            rawFile.onreadystatechange = function () {

                if (rawFile.readyState === 4) {

                    if (rawFile.status === 200 || rawFile.status == 0) {

                        var allText = rawFile.responseText;

                        resolve(allText);

                        return;

                    }

                    reject(rawFile.response);

                }

            }

            rawFile.send(null);

        })
        
        return result;

    }
    
    reloadSocketConnection() {

        // @ts-ignore
        this.socketClient.socket.disconnect()



        this._socketClient = new CBSocketClient(this)
        this._socketClient._messagesToBeSent = []

        this._functionsToCallForEachSocketClient.forEach(function (functionToCall, index, array) {

            functionToCall()

        })



    }


    callFunctionForEachSocketClient(functionToCall: () => void) {
        this._functionsToCallForEachSocketClient.push(functionToCall)
        functionToCall()
    }





}





















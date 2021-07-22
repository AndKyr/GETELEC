/// <reference path="../UICore/UICore.ts" />
/// <reference path="../UICore/UIRoute.ts" />





class LanguageService implements UILanguageService {
    
    
    static _currentLanguageKey: string
    
    
    static languageValues = {
        
        "en": {
            
            "languageName": "English",
            "languageNameShort": "ENG",
            
            "topBarTitle": "Carrum",
            
            "frontPageHTML": "Asdasd",
            "bottomBarContent": "PUT BOTTOM BAR STUFF IN HERE"
            
            
        },
        "est": {
            
            "languageName": "Eesti keel",
            "languageNameShort": "EST",
            
            "topBarTitle": "Carrum"
            
            
        }
        
    }
    
    static languages = JSON.parse(JSON.stringify(LanguageService.languageValues))
    
    static useStoredLanguageValues(values = {}) {
        
        LanguageService.languages = JSON.parse(JSON.stringify(LanguageService.languageValues)).objectByCopyingValuesRecursivelyFromObject(
            values) as any
        
        LanguageService.broadcastLanguageChangeEvent()
        
    }
    
    
    
    
    
    static broadcastLanguageChangeEvent() {
        
        UICore.main.rootViewController.view.broadcastEventInSubtree({
            name: UIView.broadcastEventName.LanguageChanged,
            parameters: {}
        })
        
    }
    
    static get defaultLanguageKey() {
        
        return "en";
        
    }
    
    static get currentLanguageKey() {
        
        if (!LanguageService._currentLanguageKey) {
            
            LanguageService.updateCurrentLanguageKey()
            
        }
        
        return LanguageService._currentLanguageKey
        
    }
    
    static updateCurrentLanguageKey(route = UIRoute.currentRoute) {
    
        var result = route.componentWithName("settings").parameters.language
    
        if (IS_NOT(result)) {
            
            result = LanguageService.defaultLanguageKey
            
        }
        
        LanguageService._currentLanguageKey = result
        
    }
    
    get currentLanguageKey() {
    
        const result = LanguageService.currentLanguageKey
    
        return result
        
    }
    
    
    
    static stringForKey(
        key: string,
        languageKey: string,
        defaultString: string,
        parameters?: { [x: string]: string | UILocalizedTextObject; }
    ) {
    
        var result
    
        if (IS(key) && LanguageService.languages[languageKey] &&
            IS_DEFINED(LanguageService.languages[languageKey][key])) {
            
            result = LanguageService.languages[languageKey][key]
            
        }
        else {
            
            result = defaultString
            
        }
        
        if (IS(parameters)) {
    
            const parameterKeys = Object.keys(parameters)
    
            parameterKeys.forEach(function (key, index, array) {
    
                const keyString = "%" + key + "%"
    
                const parameter = parameters[key]
    
                var parameterString
    
                if (parameter instanceof Object) {
                    
                    parameterString = UICore.languageService.stringForCurrentLanguage(parameter as UILocalizedTextObject)
                    
                }
                else {
                    
                    parameterString = parameter
                    
                }
                
                
                result = result.replace(new RegExp(keyString, "g"), parameterString)
                
            })
            
        }
        
        return result
        
    }
    
    stringForKey(
        key: string,
        languageKey: string,
        defaultString: string,
        parameters?: { [x: string]: string | UILocalizedTextObject; }
    ) {
        
        
        return LanguageService.stringForKey(key, languageKey, defaultString, parameters)
        
        
    }
    
    
    static localizedTextObjectForKey(
        key: string,
        defaultString = key,
        parameters?: { [x: string]: string | UILocalizedTextObject; }
    ) {
    
        const result = {}
    
        LanguageService.languages.forEach(function (languageObject, languageKey) {
            
            result[languageKey] = LanguageService.stringForKey(key, languageKey, defaultString, parameters)
            
        })
        
        return result
        
    }
    
    localizedTextObjectForKey(
        key: string,
        defaultString?: string,
        parameters?: { [x: string]: string | UILocalizedTextObject; }
    ) {
    
        const result = LanguageService.localizedTextObjectForKey(key, defaultString, parameters)
    
        return result
        
    }
    
    
    static localizedTextObjectForText(text: string) {
        
        if (IS_NOT(text)) {
            
            return nil
            
        }
    
        const result = {
        
            [LanguageService.defaultLanguageKey]: text
        
        }
    
        return result
        
    }
    
    localizedTextObjectForText(text: string) {
    
        const result = LanguageService.localizedTextObjectForText(text)
    
        return result
        
    }
    
    
    static stringForCurrentLanguage(localizedTextObject: CBLocalizedTextObject) {
        
        if (!LanguageService || !localizedTextObject) {
    
            const asd = 1
    
        }
        
        localizedTextObject = FIRST_OR_NIL(localizedTextObject)
    
        var result = localizedTextObject[LanguageService.currentLanguageKey]
    
        if (IS_NOT(result)) {
            
            result = localizedTextObject[LanguageService.defaultLanguageKey]
            
        }
        
        if (IS_NOT(result)) {
            
            result = localizedTextObject["en"]
            
        }
        
        if (IS_NOT(result)) {
            
            result = ""
            
        }
        
        return result
        
    }
    
    stringForCurrentLanguage(localizedTextObject: CBLocalizedTextObject) {
        
        return LanguageService.stringForCurrentLanguage(localizedTextObject)
        
    }
    
    
    
    
}





UICore.languageService = LanguageService




























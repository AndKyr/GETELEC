type ValueOf<T> = T[keyof T];
type PropType<TObj, TProp extends keyof TObj> = TObj[TProp];

type UIRouteParameters<T = any> = {
    
    [key: string]: string;
    
} | T;


interface UIRouteComponent<T = any> {
    
    name: string;
    parameters: UIRouteParameters<T>;
    
}


// @ts-ignore
class UIRoute extends Array<UIRouteComponent> {
    
    _isHandled: boolean = NO
    completedComponents: UIRouteComponent[] = []
    
    parentRoute: UIRoute
    
    
    constructor(hash?: string) {
        
        super()
        
        if (!hash || !hash.startsWith) {
            
            return
            
        }
        
        if (hash.startsWith("#")) {
            hash = hash.slice(1)
        }
        
        hash = decodeURIComponent(hash)
        
        const components = hash.split("]")
        components.forEach(function (this: UIRoute, component: string, index: number, array: string[]) {
            
            const componentName = component.split("[")[0]
            const parameters = {}
            
            if (!componentName) {
                
                return
                
            }
            
            const parametersString = component.split("[")[1] || ""
            const parameterPairStrings = parametersString.split(",") || []
            
            parameterPairStrings.forEach(function (pairString, index, array) {
                
                const keyAndValueArray = pairString.split(":")
                const key = decodeURIComponent(keyAndValueArray[0])
                const value = decodeURIComponent(keyAndValueArray[1])
                
                if (key) {
                    
                    parameters[key] = value
                    
                }
                
                
                
            })
            
            
            
            this.push({
                name: componentName,
                parameters: parameters
            })
            
            
            
            
            
        }, this)
        
        
    }
    
    
    
    static get currentRoute() {
        
        return new UIRoute(window.location.hash)
        
    }
    
    
    
    
    
    apply() {
        
        window.location.hash = this.stringRepresentation
        
    }
    
    
    applyByReplacingCurrentRouteInHistory() {
        
        window.location.replace(this.linkRepresentation)
        
    }
    
    
    
    copy() {
        var result = new UIRoute()
        result = Object.assign(result, this)
        return result
    }
    
    
    childRoute() {
        
        var result = this.copy()
        
        result.completedComponents.forEach(function (component, index, completedComponents) {
            
            var indexInResult = result.indexOf(component)
            
            if (indexInResult > -1) {
                
                result.splice(indexInResult, 1)
                
            }
            
        })
        
        result.completedComponents = []
        
        result.parentRoute = this
        return result
        
    }
    
    
    
    routeByRemovingComponentsOtherThanOnesNamed(componentNames: string[]) {
        const result = this.copy()
        const indexesToRemove: number[] = []
        result.forEach(function (component, index, array) {
            if (!componentNames.contains(component.name)) {
                indexesToRemove.push(index)
            }
        })
        indexesToRemove.forEach(function (indexToRemove, index, array) {
            result.removeElementAtIndex(indexToRemove)
        })
        return result
    }
    
    
    
    routeByRemovingComponentNamed(componentName: string) {
        const result = this.copy()
        const componentIndex = result.findIndex(function (component, index) {
            return (component.name == componentName)
        })
        if (componentIndex != -1) {
            result.splice(componentIndex, 1)
        }
        return result
    }
    
    
    routeByRemovingParameterInComponent(componentName: string, parameterName: string, removeComponentIfEmpty = NO) {
        var result = this.copy()
        var parameters = result.componentWithName(componentName).parameters
        if (IS_NOT(parameters)) {
            parameters = {}
        }
        delete parameters[parameterName]
        result = result.routeWithComponent(componentName, parameters)
        if (removeComponentIfEmpty && Object.keys(parameters).length == 0) {
            result = result.routeByRemovingComponentNamed(componentName)
        }
        return result
    }
    
    routeBySettingParameterInComponent(componentName: string, parameterName: string, valueToSet: string) {
        var result = this.copy()
        if (IS_NIL(valueToSet) || IS_NIL(parameterName)) {
            return result
        }
        var parameters = result.componentWithName(componentName).parameters
        if (IS_NOT(parameters)) {
            parameters = {}
        }
        parameters[parameterName] = valueToSet
        result = result.routeWithComponent(componentName, parameters)
        return result
    }
    
    
    routeWithViewControllerComponent<T extends typeof UIViewController>(
        viewController: T,
        parameters: UIRouteParameters<{ [P in keyof T["ParameterIdentifierName"]]: string }>,
        extendParameters: boolean = NO
    ) {
        
        return this.routeWithComponent(viewController.routeComponentName, parameters, extendParameters)
        
    }
    
    routeWithComponent(name: string, parameters: UIRouteParameters, extendParameters: boolean = NO) {
        
        const result = this.copy()
        var component = result.componentWithName(name)
        if (IS_NOT(component)) {
            component = {
                name: name,
                parameters: {}
            }
            result.push(component)
        }
        
        if (IS_NOT(parameters)) {
            
            parameters = {}
            
        }
        
        if (extendParameters) {
            component.parameters = Object.assign(component.parameters, parameters)
        }
        else {
            component.parameters = parameters
        }
        
        return result
        
    }
    
    navigateBySettingComponent(name: string, parameters: UIRouteParameters, extendParameters: boolean = NO) {
        
        this.routeWithComponent(name, parameters, extendParameters).apply()
        
    }
    
    
    
    componentWithViewController<T extends typeof UIViewController>(viewController: T): UIRouteComponent<{ [P in keyof T["ParameterIdentifierName"]]: string }> {
        
        return this.componentWithName(viewController.routeComponentName)
        
    }
    
    componentWithName(name: string): UIRouteComponent {
        var result = nil
        this.forEach(function (component, index, self) {
            if (component.name == name) {
                result = component
            }
        })
        return result
    }
    
    
    
    didcompleteComponent(component: UIRouteComponent) {
        
        const self: UIRoute = this
        const index = self.indexOf(component, 0)
        if (index > -1) {
            
            
            self.completedComponents.push(self.splice(index, 1)[0])
            
            //self.completedComponents.push(component);
            
        }
        
    }
    
    set isHandled(isHandled: boolean) {
        
        this._isHandled = isHandled
        
    }
    
    get isHandled() {
        
        return (this._isHandled || (this.length == 0 && this.completedComponents.length != 0))
        
    }
    
    
    get linkRepresentation() {
        return "#" + this.stringRepresentation
    }
    
    
    get stringRepresentation() {
        var result = ""
        this.completedComponents.forEach(function (component, index, self) {
            result = result + component.name
            const parameters = component.parameters
            result = result + "["
            Object.keys(parameters).forEach(function (key, index, keys) {
                if (index) {
                    result = result + ","
                }
                result = result + encodeURIComponent(key) + ":" + encodeURIComponent(parameters[key])
            })
            result = result + "]"
        })
        this.forEach(function (component, index, self) {
            result = result + component.name
            const parameters = component.parameters
            result = result + "["
            Object.keys(parameters).forEach(function (key, index, keys) {
                if (index) {
                    result = result + ","
                }
                result = result + encodeURIComponent(key) + ":" + encodeURIComponent(parameters[key])
            })
            result = result + "]"
        })
        return result
    }
    
    
    
    
    
}


















































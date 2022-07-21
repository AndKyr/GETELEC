const addressLines = {

	mainPage: "https://getelec.org/index.html",
	ivCalculationPage: "https://getelec.org/ivCalculations.html",
	threeDGraphing: "https://getelec.org/3dGraphing.html",
	developersPage: "https://getelec.org/developers.html",
	documentationPage: "https://getelec.org/documentation.html",
	briefDocumentationPage: "https://getelec.org/BriefDoc.html"

}

const typeDict = {
    0: "String", 1: "Array",
    99: "Unknown"
}

const separatorDict = {
    0: ".", 1: ",",
    2: " ", 99: "Unknown"
}

const multDict = {
    0: "*", 1: "e", 99: "Unknown"
}

const commaDict = {
    0: ",", 1: ".", 99: "Unknown"
}

export { addressLines, typeDict, separatorDict, multDict, commaDict };
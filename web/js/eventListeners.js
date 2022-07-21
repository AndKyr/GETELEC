export function getButtonsForMain() {

	let ivDataBtn = document.getElementById("ivDataBtn");
	let threeDGraphingBtn = document.getElementById("3dGraphingBtn");
	let documentationBtn = document.getElementById("documentationBtn");

	ivDataBtn.addEventListener("click", goToIvCalculationPage);
	threeDGraphingBtn.addEventListener("click", goToThreeDGraphingPage);
	documentationBtn.addEventListener("click", goToDocumentationPage);
}

export function getButtonsForIvCalculations() {
	let logo = document.getElementById("logo");

	logo.addEventListener("click", goToMainPage());
}

export function getButtonsForDocumentration() {
	let logo = document.getElementById("logo");

	logo.addEventListener("click", goToMainPage());
}

function goToMainPage() {
	window.location.href = "./index.html";
}

function goToIvCalculationPage() {
	window.location.href = "./ivCalculations.html";
}

function goToThreeDGraphingPage() {
	window.location.href = "./3dGraphing.html";
}

function goToDocumentationPage() {
	window.location.href = "./BriefDoc.html";
}
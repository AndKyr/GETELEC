function main() {

    loadEventListeners();

}

main();


function loadEventListeners() {

	let ivDataBtn = document.getElementById("ivDataBtn");
	let threeDGraphingBtn = document.getElementById("emissionDataBtn");
	let documentationBtn = document.getElementById("documentationBtn");

	ivDataBtn.addEventListener("click", goToIvCalculationPage);
	threeDGraphingBtn.addEventListener("click", goToThreeDGraphingPage);
	documentationBtn.addEventListener("click", goToDocumentationPage);
}

function getButtonsForIvCalculations() {
	let logo = document.getElementById("logo");

	logo.addEventListener("click", window.location.href="./index.html");
}

function getButtonsForDocumentration() {
	let logo = document.getElementById("logo");

	logo.addEventListener("click", window.location.href = "./index.html");
}

function goToIvCalculationPage() {
	window.location.href = "./ivCalculations.html";
}

function goToThreeDGraphingPage() {
	window.location.href = "./emissionCalc.html";
}

function goToDocumentationPage() {
	window.location.href = "./briefDoc.html";
}
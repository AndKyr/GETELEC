import * as path from "path"
import "./lib/Extensions"
import express from "express"
import { SocketController } from "./lib/SocketController"
import { GETELECController } from "./src/GETELECController"
import SocketIO = require("socket.io")


var expressApplication = express()

expressApplication.use(express.json())

//expressApp.set('view engine', 'ejs');// tell Express we're using EJS
expressApplication.set("views", path.join(__dirname + "/webclient"))// set path to *.ejs files

expressApplication.get("/", (request, response) => {
    
    response.sendFile(__dirname + "/webclient/Test.html")
    
})

expressApplication.get("/frontpage.html", (request, response) => {
    
    response.sendFile(__dirname + "/webclient/Frontpage.html")
    
})

expressApplication.get("/styles.css", (request, response) => {
    
    response.sendFile(__dirname + "/webclient/styles.css")
    
})
expressApplication.get("/compiledScripts/:filename", (request, response) => {
    
    response.sendFile(__dirname + "/webclient/compiledScripts/" + request.params.filename)
    
})

expressApplication.get("/images/:filename", (request, response) => {
    
    response.sendFile(__dirname + "/webclient/images/" + request.params.filename)
    
})

let server = expressApplication.listen(3000)

var socketIO = SocketIO(server, { cookie: false, pingTimeout: 30000 })

SocketController.Instance(socketIO)

let asd = GETELECController.Instance(expressApplication)












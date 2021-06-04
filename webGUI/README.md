# edit-qd-script

This program makes it easier to use GETELEC.

# Description

Start the application.
Enter your data.
Press the load button.

# Prerequisites

1. Install Typescript and TS-node
```npm install -g typescript``` and ```npm install -g ts-node```
2. To be able to distribute, install PKG
```npm install -g pkg```

# Installation

1. Clone this repo
2. Initialize with ```npm install```
3. Run ```tsc``` in the ```webclient``` folder
4. Run ```tsc``` in the root of the project

# Running

Either run using VSCode by pressing the run button or run ```ts-node -r ./app.ts``` .
Then open a browser and go to localhost:7005

# Running the getelec python server
Run: hug -f API.py

# Packaging

Package with ```pkg .``` in the root folder.





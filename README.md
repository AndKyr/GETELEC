# GETELEC
General Tool for Electron Emission Calculations - A computational tool for calculating thermal-field electron emission current density and Nottingham effect heat for metallic and semiconducting emitters.

For details, see the associated publication http://dx.doi.org/10.1016/j.commatsci.2016.11.010, https://arxiv.org/pdf/1609.02364. If you use the code, please cite these papers.

# HOW TO RUN

## Prepare library for use. Cook tables. Setting up GETELEC

- [ ] Make sure to be in the main directory

```
make clean-all
make
```

## Start server for web client

Default port used by node server: 3000

- [ ] Make sure to be in /web directory
- [ ] Launch node server using following command: 

```
node server
```

## Server commands

- [ ] /q : Stop server (Alternative: /stop, /exit, q, stop, exit)
- [ ] /kickAll : Kick all clients (WIP)

## Required packages for node

- [ ] npm install winston
- [ ] npm install socket.io

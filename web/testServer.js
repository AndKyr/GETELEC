const { PythonShell } = require('python-shell');

const params = JSON.stringify({  
     
    "materialType":[1],"sweepParam":[2],"field":[2,2.48,2.95,3.43,3.9,4.38,4.86,5.33,5.81,6.29,6.76,7.24,7.71,8.19,8.67,9.14,9.62,10.1,10.57,11.05,11.52,12],"radius":[50],"work_function":[4.5],"temperature":[300],"ec":[4.05],"ef":[4.61],"eg":[1.12],"gammaMetal":[10],"gammaSemi":[10],"me":[0.98],"mp":[0.5],"calculateEC":[1],"calculateNH":[0],"calculateES":[0]
  });

  const options = {

    scriptPath: "",
    args: [params]

  }

  console.log(`${updateTime()} Got a job: ${params}`);

  PythonShell.run("python/returnHW.py", options, (err, res) => {

    if(err) console.log(err);

    console.log(res)

  })

  PythonShell.run("python/DELETEME.py", options, (err, res) => {

    if(err) console.log(err);

    console.log(res)

  })




  function updateTime() {

    let date_ob = new Date();
    let hours = date_ob.getHours();
    let minutes = date_ob.getMinutes();
    let seconds = date_ob.getSeconds();
    time = "\x1b[36m[" + hours + ":" + minutes + ":" + seconds + " INFO]:\x1b[0m ";
  
    return time;
  
  }
  
  function updateTimeError() {
  
    let date_ob = new Date();
    let hours = date_ob.getHours();
    let minutes = date_ob.getMinutes();
    let seconds = date_ob.getSeconds();
    time = "\x1b[31m[" + hours + ":" + minutes + ":" + seconds + " INFO]:\x1b[0m ";
  
    return time;
  
  }
  

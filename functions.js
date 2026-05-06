function arrayToCSV(data){    
    var csv = [];
    for (var i=0; i < data[0].length; i++){
        row = []
        for (var j=0; j < data.length; j++){
            row.push(data[j][i]);
        }
        csv.push(row);
    }
    return csv.map(row =>
        row
        .map(String)
        .map(v => v.replaceAll('"','""'))
        .join(',')
        ).join('\r\n');
}

function nestedArrayToCSV(data){
    var csv = [];
    for(var a=0; a<data.length;a++){
        data[a]=data[a].flat();
    }
    // flatten all the subarrays in each array of data
    for (var i=0; i < data[0].flat().length; i++){
        row=[];
        for (var j=0; j < data.length; j++){
            row.push(data[j][i]);
            }
        csv.push(row);           
    }
        
    return csv.map(row =>
        row
        .map(String)
        .map(v => v.replaceAll('"','""'))
        .join(',')
        ).join('\r\n');
}

function downloadBlob(content,buttonname,filename,contentType){
    var blob = new Blob([content], {type: contentType});
    var url=URL.createObjectURL(blob);
    var pom = document.getElementById(buttonname);
    pom.setAttribute("href",url);
    pom.setAttribute('download',filename);
    return(pom)
}

function getBessel(x){
    let bessel=[besselK(0,x),besselK(1,x)]
    return bessel
}

function calcSmin(){
    var d0=document.getElementById("d0").value;
    var l=document.getElementById("l").value;
    var r0=document.getElementById("r0").value;
    var ka=document.getElementById("Ka").value;
    ka = 4*ka;
    var kc=document.getElementById("Kc").value;
    var s=document.getElementById("s").value;
    var c0=document.getElementById("c0").value;

    var u0=div(sub(d0,l),2);

    landa = pow((d0*d0*kc)/ka,0.25);
    alpha = 2*Math.PI*r0*kc;
    r0landa = r0/landa;
    multFact = 2*Math.PI*r0landa;
    ipsqrt=sqrt(complex(0,1));
    insqrt=sqrt(complex(0,-1));
    argp=div(mul(ipsqrt,r0),landa);
    argn=div(mul(insqrt,r0),landa);
    Kargp=getBessel(argp);
    Kargn=getBessel(argn);

    // Eq.25
    temp1=mul(Kargp[0],Kargn[1]);
    temp2=mul(mul(Kargn[0],Kargp[1]),complex(0,1));
    // divi is the denominator for Eqs.25, 26, and 27
    divi=sub(temp1,temp2);
    f1=div(mul(mul(insqrt,Kargp[1]),Kargn[1]),divi);
    f1=mul(f1,(multFact*r0landa*r0landa));
    // Eq.26
    // temp1 and temp2 here used for numerator
    temp1=mul(Kargn[0],Kargp[1]);
    temp2=mul(mul(Kargp[0],Kargn[1]),complex(0,1));
    f2=div(sub(temp1,temp2),divi);
    f2=mul(f2,multFact*r0landa);
    // Eq.27
    f3=div(mul(insqrt,mul(Kargp[0],Kargn[0])),divi);
    f3=mul(f3,multFact);

    a1=mul(f1,div(kc,r0*r0));
    a2=mul(f2,div(kc,r0));
    a3=mul(f3,kc);

    temp1=add(mul(a2,u0),alpha*c0);
    temp2=mul(a3,-2.0);
    comSmin=div(temp1,temp2)
    return(comSmin.re)
}

function rStopMin(){
    var r0=document.getElementById("r0");
    var rStop=document.getElementById("rstop");
    rStop.min=r0.value;
}

function calcU(r){
    var changeParams=4.114;
    var imagErrorThreshold=Math.pow(10,-10);

    var d0=document.getElementById("d0").value;
    var l=document.getElementById("l").value;
    var r0=document.getElementById("r0").value;
    var ka=document.getElementById("Ka").value;
    var kc=document.getElementById("Kc").value;
    var sigma=document.getElementById("sigma").value;
    var s=document.getElementById("s").value;
    if(document.getElementById("sRelaxed_button").checked){
        s=calcSmin();
    }
    else if(document.getElementById("sConstrained_button").checked){
        s=0;
    }
    var c0=document.getElementById("c0").value;

    var u0=div(sub(d0,l),2);

    var gamma=sigma/kc;
    var beta=(4*ka)/(d0*d0*kc);

    var temp=sqrt(complex((gamma*gamma)-(4*beta),0));
    var kp2=div(add(complex(gamma,0),temp),2);
    var kn2=div(sub(complex(gamma,0),temp),2);
    var kp=sqrt(kp2);
    var kn=sqrt(kn2);

    // Kkpr and Kknr are arrays of length 2, based on the input r;
    // Kkpr0 and Kknr0 are arrays of length 2, based on r0;
    // array[0] is an object with the real (.re) and imaginary (.im) output of besselK order 0 
    // array [1] is an object with the real (.re) and imaginary (.im) of besselK order 1 
    var Kkpr0=getBessel(mul(kp,r0));
    var Kknr0=getBessel(mul(kn,r0));
    var Kkpr=getBessel(mul(kp,r));
    var Kknr=getBessel(mul(kn,r));
    var divi= sub(mul(kn,(mul(Kkpr0[0],Kknr0[1]))),mul(kp,(mul(Kknr0[0],Kkpr0[1]))));
    var Ap=div(add(mul(kn,mul(Kknr0[1],u0)),mul(Kknr0[0],s)),divi);
    var An=div(sub(mul(kp,mul(Kkpr0[1],-u0)),mul(Kkpr0[0],s)),divi);
    var res=add(mul(Ap,Kkpr[0]),mul(An,Kknr[0]));

    return(res.re)
}

function calcdU(r){
    var changeParams=4.114;
    var imagErrorThreshold=Math.pow(10,-10)

    var d0=document.getElementById("d0").value;
    var l=document.getElementById("l").value;
    var r0=document.getElementById("r0").value;
    if(document.getElementById("r0_checkbox").checked){
            //update for the case of l1 -> l2
            r0=r0;
        }
    var ka=document.getElementById("Ka").value;
    var kc=document.getElementById("Kc").value;
    var alpha=document.getElementById("alpha").value;
    var s=document.getElementById("s").value;
    if(document.getElementById("sRelaxed_button").checked){
        s=calcSmin();
    }
    else if(document.getElementById("sConstrained_button").checked){
        s=0;
    }
    var c0=document.getElementById("c0").value;

    var u0=div(sub(d0,l),2);

    var u0=div(sub(d0,l),2);

    var gamma=alpha/kc;
    var beta=(4*ka)/(d0*d0*kc);

    var temp=sqrt(complex((gamma*gamma)-(4*beta),0));
    var kp2=div(add(complex(gamma,0),temp),2);
    var kn2=div(sub(complex(gamma,0),temp),2);
    var kp=sqrt(kp2);
    var kn=sqrt(kn2);

    // Kkpr and Kknr are arrays of length 2, based on the input r;
    // Kkpr0 and Kknr0 are arrays of length 2, based on r0;
    // array[0] is an object with the real (.re) and imaginary (.im) output of besselK order 0 
    // array [1] is an object with the real (.re) and imaginary (.im) of besselK order 1 
    var Kkpr=getBessel(mul(kp,r));
    var Kknr=getBessel(mul(kn,r));
    var Kkpr0=getBessel(mul(kp,r0));
    var Kknr0=getBessel(mul(kn,r0));
    var divi= sub(mul(kn,(mul(Kkpr0[0],Kknr0[1]))),mul(kp,(mul(Kknr0[0],Kkpr0[1]))));
    var Ap=div(add(mul(kn,mul(Kknr0[1],u0)),mul(Kknr0[0],s)),divi);
    var An=div(sub(mul(kp,mul(Kkpr0[1],-u0)),mul(Kkpr0[0],s)),divi);
    var res=add(mul(mul(Ap,kp),Kkpr[1]),mul(mul(An,kn),Kknr[1]));
    return(-1*res.re)
}

function getUforRange(startR,stopR,interval){
    var i=0;
    // due to the floating numbers, this isn't always an integer
    var length=(stopR-startR)/interval;
    //NOTE: division of (1.4-1.0)/0.1 gives 3.999, not 4.0
    var rArray=[];
    var uArray=[];
    while(i<length){
        var currR;
        if(i==0){
            currR=startR;
            rArray.push(currR);
        }
        else{
            currR=Number(rArray[i-1])+Number(interval);
            rArray.push(currR);
        }
        var currU=calcU(currR);
        uArray.push(currU);
        i++;
    }
    return([rArray,uArray])
}

function formatPertPlotData(getURangeOutput){
    var d0=document.getElementById("d0").value;
    var monolayerWidth=Number(d0)/2;
    let plotData = getURangeOutput;
    let coords = plotData[0].map( (r,i) => ({x: r, y: Number(Number(monolayerWidth) - plotData[1][i])}));
    return(coords)
}

function fillPertCoordData(getURangeOutput){
    for(var i=0; i < getURangeOutput[0].length; i++){
        document.getElementById("pertCoordTable").innerHTML+="<tr><td>"+Number(getURangeOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getURangeOutput[1][i]).toPrecision(5)+"</td></tr>"
    }
}

function fillPertCoordDataFast(getURangeOutput){
    var pertCoords;
    for(var i=0; i < getURangeOutput[0].length; i++){
        if(pertCoords == undefined){
            pertCoords="<tr><td>"+Number(getURangeOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getURangeOutput[1][i]).toPrecision(5)+"</td></tr>"
        }
        else{
            pertCoords+="<tr><td>"+Number(getURangeOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getURangeOutput[1][i]).toPrecision(5)+"</td></tr>"
        }
        document.getElementById("pertCoordTable").innerHTML+=pertCoords;
    }
}

function getDeltaGce(r){
    var changeParams=4.114;
    var d0=document.getElementById("d0").value;
    var ka=document.getElementById("Ka").value;
    var u=calcU(r);
    return((u*u*Math.PI*r*4*ka)/(d0*d0)/changeParams)
}

function getDeltaGsd(r){
    var changeParams=4.114;
    var imagErrorThreshold=Math.pow(10,-10)

    var d0=document.getElementById("d0").value;
    var l=document.getElementById("l").value;
    var r0=document.getElementById("r0").value;
    var ka=document.getElementById("Ka").value;
    var kc=document.getElementById("Kc").value;
    var alpha=document.getElementById("alpha").value;
    var s=document.getElementById("s").value;
    if(document.getElementById("sRelaxed_button").checked){
        s=calcSmin();
    }
    else if(document.getElementById("sConstrained_button").checked){
        s=0;
    }
    var c0=document.getElementById("c0").value;

    var u0=div(sub(d0,l),2);

    var u0=div(sub(d0,l),2);

    var gamma=alpha/kc;
    var beta=(4*ka)/(d0*d0*kc);

    var temp=sqrt(complex((gamma*gamma)-(4*beta),0));
    var kp2=div(add(complex(gamma,0),temp),2);
    var kn2=div(sub(complex(gamma,0),temp),2);
    var kp=sqrt(kp2);
    var kn=sqrt(kn2);

    // Kkpr0 and Kknr0 are arrays of length 2;
    // array[0] is an object with the real (.re) and imaginary (.im) output of besselK order 0 
    // array [1] is an object with the real (.re) and imaginary (.im) of besselK order 1 
    var Kkpr=getBessel(mul(kp,r));
    var Kknr=getBessel(mul(kn,r));
    var Kkpr0=getBessel(mul(kp,r0));
    var Kknr0=getBessel(mul(kn,r0));
    var divi= sub(mul(kn,(mul(Kkpr0[0],Kknr0[1]))),mul(kp,(mul(Kknr0[0],Kkpr0[1]))));
    var Ap=div(add(mul(kn,mul(Kknr0[1],u0)),mul(Kknr0[0],s)),divi);
    var An=div(sub(mul(kp,mul(Kkpr0[1],-u0)),mul(Kkpr0[0],s)),divi);
    var temp1sd=mul(kp2,mul(Ap,Kkpr[0]));
    var temp2sd=mul(kn2,mul(An,Kknr[0]));
    var res=mul(pow(add(temp1sd,temp2sd),2.0),(Math.PI*r*kc));
    return(res.re/changeParams)
}

function getDeltaGmec(r){
    var changeParams=4.114;
    var imagErrorThreshold=Math.pow(10,-10)

    var d0=document.getElementById("d0").value;
    var l=document.getElementById("l").value;
    var r0=document.getElementById("r0").value;
    var ka=document.getElementById("Ka").value;
    var kc=document.getElementById("Kc").value;
    var alpha=document.getElementById("alpha").value;
    var s=document.getElementById("s").value;
    if(document.getElementById("sRelaxed_button").checked){
        s=calcSmin();
    }
    else if(document.getElementById("sConstrained_button").checked){
        s=0;
    }
    var c0=document.getElementById("c0").value;

    var u0=div(sub(d0,l),2);

    var u0=div(sub(d0,l),2);

    var gamma=alpha/kc;
    var beta=(4*ka)/(d0*d0*kc);

    var temp=sqrt(complex((gamma*gamma)-(4*beta),0));
    var kp2=div(add(complex(gamma,0),temp),2);
    var kn2=div(sub(complex(gamma,0),temp),2);
    var kp=sqrt(kp2);
    var kn=sqrt(kn2);

    // Kkpr0 and Kknr0 are arrays of length 2;
    // array[0] is an object with the real (.re) and imaginary (.im) output of besselK order 0 
    // array [1] is an object with the real (.re) and imaginary (.im) of besselK order 1 
    var Kkpr=getBessel(mul(kp,r));
    var Kknr=getBessel(mul(kn,r));
    var Kkpr0=getBessel(mul(kp,r0));
    var Kknr0=getBessel(mul(kn,r0));
    var divi= sub(mul(kn,(mul(Kkpr0[0],Kknr0[1]))),mul(kp,(mul(Kknr0[0],Kkpr0[1]))));
    var Ap=div(add(mul(kn,mul(Kknr0[1],u0)),mul(Kknr0[0],s)),divi);
    var An=div(sub(mul(kp,mul(Kkpr0[1],-u0)),mul(Kkpr0[0],s)),divi);
    var temp1mec=mul(kp2,mul(Ap,Kkpr[0]));
    var temp2mec=mul(kn2,mul(An,Kknr[0]));
    var res=mul(add(temp1mec,temp2mec),(-2*Math.PI*kc*c0*r));
    return(res.re/changeParams)
}

function getDeltaGst(r){
    var changeParams=4.114;
    var imagErrorThreshold=Math.pow(10,-10)

    var d0=document.getElementById("d0").value;
    var l=document.getElementById("l").value;
    var r0=document.getElementById("r0").value;
    var ka=document.getElementById("Ka").value;
    var kc=document.getElementById("Kc").value;
    var alpha=document.getElementById("alpha").value;
    var s=document.getElementById("s").value;
    if(document.getElementById("sRelaxed_button").checked){
        s=calcSmin();
    }
    else if(document.getElementById("sConstrained_button").checked){
        s=0;
    }
    var c0=document.getElementById("c0").value;

    var u0=div(sub(d0,l),2);

    var u0=div(sub(d0,l),2);

    var gamma=alpha/kc;
    var beta=(4*ka)/(d0*d0*kc);

    var temp=sqrt(complex((gamma*gamma)-(4*beta),0));
    var kp2=div(add(complex(gamma,0),temp),2);
    var kn2=div(sub(complex(gamma,0),temp),2);
    var kp=sqrt(kp2);
    var kn=sqrt(kn2);

    // Kkpr0 and Kknr0 are arrays of length 2;
    // array[0] is an object with the real (.re) and imaginary (.im) output of besselK order 0 
    // array [1] is an object with the real (.re) and imaginary (.im) of besselK order 1 
    var Kkpr=getBessel(mul(kp,r));
    var Kknr=getBessel(mul(kn,r));
    var Kkpr0=getBessel(mul(kp,r0));
    var Kknr0=getBessel(mul(kn,r0));
    var divi= sub(mul(kn,(mul(Kkpr0[0],Kknr0[1]))),mul(kp,(mul(Kknr0[0],Kkpr0[1]))));
    var Ap=div(add(mul(kn,mul(Kknr0[1],u0)),mul(Kknr0[0],s)),divi);
    var An=div(sub(mul(kp,mul(Kkpr0[1],-u0)),mul(Kkpr0[0],s)),divi);
    var temp1st=mul(kp,mul(Ap,Kkpr[1]));
    var temp2st=mul(kn,mul(An,Kknr[1]));
    var res=mul(pow(add(temp1st,temp2st),2),(Math.PI*r*alpha));
    return(res.re/changeParams)
}

function getEnergyDecompforRange(startR,stopR,interval){
    var i=0;
    // due to the floating numbers, this isn't always an integer
    var length=(stopR-startR)/interval;
    //NOTE: division of (1.4-1.0)/0.1 gives 3.999, not 4.0
    var rArray=[];
    var deltaGceArray=[];
    var deltaGsdArray=[];
    var deltaGmecArray=[];
    var deltaGstArray=[];
    while(i<length){
        var currR;
        if(i==0){
            currR=startR;
            rArray.push(currR);
        }
        else{
            currR=Number(rArray[i-1])+Number(interval);
            rArray.push(currR);
        }
        var currDeltaGce=getDeltaGce(currR);
        var currDeltaGsd=getDeltaGsd(currR);
        var currDeltaGmec=getDeltaGmec(currR);
        var currDeltaGst=getDeltaGst(currR);
        deltaGceArray.push(currDeltaGce);
        deltaGsdArray.push(currDeltaGsd);
        deltaGmecArray.push(currDeltaGmec);
        deltaGstArray.push(currDeltaGst);
        i++;
    }
    return([rArray,deltaGceArray,deltaGsdArray,deltaGmecArray,deltaGstArray])
}

function formatEnergyPlotData(getEnergyDecompOutput_r,getEnergyDecompOutputDeltaG){
    let rData = getEnergyDecompOutput_r;
    let deltaData = getEnergyDecompOutputDeltaG;
    let coords = rData.map( (r,i) => ({x: r, y: Number(deltaData[i])}));
    return(coords)
}

function fillEnergyCoordData(getEnergyDecompOutput){
    for(var i=0; i < getEnergyDecompOutput[0].length; i++){
        document.getElementById("energyCoordTable").innerHTML+="<tr><td>"+Number(getEnergyDecompOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getEnergyDecompOutput[1][i]).toPrecision(5)+"</td><td>"+Number(getEnergyDecompOutput[2][i]).toPrecision(5)+"</td><td>"+Number(getEnergyDecompOutput[3][i]).toPrecision(5)+"</td></tr>"
    }
}

function fillEnergyCoordDataFast(getEnergyDecompOutput){
    var energyCoords;
    for(var i=0; i < getEnergyDecompOutput[0].length; i++){
        if(energyCoords == undefined){
            energyCoords="<tr><td>"+Number(getEnergyDecompOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getEnergyDecompOutput[1][i]).toPrecision(5)+"</td><td>"+Number(getEnergyDecompOutput[2][i]).toPrecision(5)+"</td><td>"+Number(getEnergyDecompOutput[3][i]).toPrecision(5)+"</td><td>"+Number(getEnergyDecompOutput[4][i]).toPrecision(5)+"</td></tr>"
        }
        else{
            energyCoords+="<tr><td>"+Number(getEnergyDecompOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getEnergyDecompOutput[1][i]).toPrecision(5)+"</td><td>"+Number(getEnergyDecompOutput[2][i]).toPrecision(5)+"</td><td>"+Number(getEnergyDecompOutput[3][i]).toPrecision(5)+"</td><td>"+Number(getEnergyDecompOutput[4][i]).toPrecision(5)+"</td></tr>"
        }
    }
    document.getElementById("energyCoordTable").innerHTML+=energyCoords;
}

function getCurvature(r){
    var changeParams=4.114;
    var imagErrorThreshold=Math.pow(10,-10)

    var d0=document.getElementById("d0").value;
    var l=document.getElementById("l").value;
    var r0=document.getElementById("r0").value;
    var ka=document.getElementById("Ka").value;
    var kc=document.getElementById("Kc").value;
    var alpha=document.getElementById("alpha").value;
    var s=document.getElementById("s").value;
    if(document.getElementById("sRelaxed_button").checked){
        s=calcSmin();
    }
    else if(document.getElementById("sConstrained_button").checked){
        s=0;
    }
    var c0=document.getElementById("c0").value;

    var u0=div(sub(d0,l),2);

    var u0=div(sub(d0,l),2);

    var gamma=alpha/kc;
    var beta=(4*ka)/(d0*d0*kc);

    var temp=sqrt(complex((gamma*gamma)-(4*beta),0));
    var kp2=div(add(complex(gamma,0),temp),2);
    var kn2=div(sub(complex(gamma,0),temp),2);
    var kp=sqrt(kp2);
    var kn=sqrt(kn2);

    // Kkpr0 and Kknr0 are arrays of length 2;
    // array[0] is an object with the real (.re) and imaginary (.im) output of besselK order 0 
    // array [1] is an object with the real (.re) and imaginary (.im) of besselK order 1 
    var Kkpr=getBessel(mul(kp,r));
    var Kknr=getBessel(mul(kn,r));
    var Kkpr0=getBessel(mul(kp,r0));
    var Kknr0=getBessel(mul(kn,r0));
    var divi= sub(mul(kn,(mul(Kkpr0[0],Kknr0[1]))),mul(kp,(mul(Kknr0[0],Kkpr0[1]))));
    var Ap=div(add(mul(kn,mul(Kknr0[1],u0)),mul(Kknr0[0],s)),divi);
    var An=div(sub(mul(kp,mul(Kkpr0[1],-u0)),mul(Kkpr0[0],s)),divi);
    var temp1curv=mul(kp2,mul(Ap,Kkpr[0]));
    var temp2curv=mul(kn2,mul(An,Kknr[0]));
    var res=add(temp1curv,temp2curv);
    return(res.re/changeParams)
}

function getCurvatureforRange(startR,stopR,interval){
    var i=0;
    // due to the floating numbers, this isn't always an integer
    var length=(stopR-startR)/interval;
    //NOTE: division of (1.4-1.0)/0.1 gives 3.999, not 4.0
    var rArray=[];
    var curvatureArray=[];
    while(i<length){
        var currR;
        if(i==0){
            currR=startR;
            rArray.push(currR);
        }
        else{
            currR=Number(rArray[i-1])+Number(interval);
            rArray.push(currR);
        }
        var currCurvature=getCurvature(currR);
        curvatureArray.push(currCurvature);
        i++;
    }
    return([rArray,curvatureArray])
}

function fillCurvatureCoordData(getCurvatureOutput){
    for(var i=0; i < getCurvatureOutput[0].length; i++){
        document.getElementById("curvatureTable").innerHTML+="<tr><td>"+Number(getCurvatureOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getCurvatureOutput[1][i]).toPrecision(5)+"</td></tr>"
    }
}

function fillCurvatureCoordDataFast(getCurvatureOutput){
    var curvatureCoords;
    for(var i=0; i < getCurvatureOutput[0].length; i++){
        if(curvatureCoords == undefined){
            curvatureCoords="<tr><td>"+Number(getCurvatureOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getCurvatureOutput[1][i]).toPrecision(5)+"</td></tr>"           
        }
        else{
            curvatureCoords+="<tr><td>"+Number(getCurvatureOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getCurvatureOutput[1][i]).toPrecision(5)+"</td></tr>"
        }
    }
    document.getElementById("curvatureTable").innerHTML+=curvatureCoords;
}

function curvature_Analytical(){
var d0=document.getElementById("d0").value;
var l=document.getElementById("l").value;
var r0=document.getElementById("r0").value;
var ka=document.getElementById("Ka").value;
ka = 4*ka;
var kc=document.getElementById("Kc").value;
var s=document.getElementById("s").value;
if(document.getElementById("sRelaxed_button").checked){
    s=calcSmin();
}
else if(document.getElementById("sConstrained_button").checked){
    s=0;
}
var c0=document.getElementById("c0").value;

var u0=div(sub(d0,l),2);

landa = pow((d0*d0*kc)/ka,0.25);
alpha = 2*Math.PI*r0*kc;
r0landa = r0/landa;
multFact = 2*Math.PI*r0landa;
ipsqrt=sqrt(complex(0,1));
insqrt=sqrt(complex(0,-1));
argp=div(mul(ipsqrt,r0),landa);
argn=div(mul(insqrt,r0),landa);
Kargp=getBessel(argp);
Kargn=getBessel(argn);

// Eq.25
temp1=mul(Kargp[0],Kargn[1]);
temp2=mul(mul(Kargn[0],Kargp[1]),complex(0,1));
// divi is the denominator for Eqs.25, 26, and 27
divi=sub(temp1,temp2);
f1=div(mul(mul(insqrt,Kargp[1]),Kargn[1]),divi);
f1=mul(f1,(multFact*r0landa*r0landa));
// Eq.26
// temp1 and temp2 here used for numerator
temp1=mul(Kargn[0],Kargp[1]);
temp2=mul(mul(Kargp[0],Kargn[1]),complex(0,1));
f2=div(sub(temp1,temp2),divi);
f2=mul(f2,multFact*r0landa);
// Eq.27
f3=div(mul(insqrt,mul(Kargp[0],Kargn[0])),divi);
f3=mul(f3,multFact);

a1=mul(f1,div(kc,r0*r0));
a2=mul(f2,div(kc,r0));
a3=mul(f3,kc);

var changeParams=4.114;

Hb=sub(a1,div(mul(a2,a2),mul(a3,4)));
Hb=div(Hb,changeParams);
Hx=div(mul(a2,(Math.PI*kc*r0*-1)),a3);
Hx=div(Hx,changeParams);
Hc=div(pow((Math.PI*kc*r0),2)*-1,a3);
Hc=div(Hc,changeParams);

// calculate smin analytically
temp1=add(mul(a2,u0),alpha*c0);
temp2=mul(a3,-2);
smin=div(temp1,temp2);

//calculate deltaGdef analytical 
temp1=mul(a1,u0*u0);
temp2=mul(a2,s*u0);
temp3=mul(a3,s*s);
res=add(add(add(temp1,temp2),temp3),alpha*s*c0);
res=div(res,changeParams);
return([a1,a2,a3,Hb,Hx,Hc,smin,res]);
}

function fillBiquadraticCoordData(getBiquadraticOutput){
    for(var i=0; i < getBiquadraticOutput[0].length; i++){
        var currArray = getBiquadraticOutput[0][i];
        for(var j=0; j<currArray.length;j++){
            document.getElementById("bqCoordTable").innerHTML+="<tr><td>"+Number(getBiquadraticOutput[0][i][j]).toPrecision(3)+"</td><td>"+Number(getBiquadraticOutput[1][i][j]).toPrecision(3)+"</td><td>"+Number(getBiquadraticOutput[2][i][j]).toPrecision(5)+"</td></tr>"
        }        
    }
} 

function fillBiquadraticCoordDataFast(getBiquadraticOutput){
    var biquadraticCoords;
    for(var i=0; i < getBiquadraticOutput[0].length; i++){
        var currArray = getBiquadraticOutput[0][i];
        for(var j=0; j<currArray.length;j++){
            if(biquadraticCoords == undefined){
                biquadraticCoords="<tr><td>"+Number(getBiquadraticOutput[0][i][j]).toPrecision(3)+"</td><td>"+Number(getBiquadraticOutput[1][i][j]).toPrecision(3)+"</td><td>"+Number(getBiquadraticOutput[2][i][j]).toPrecision(5)+"</td></tr>";
            }
            else{
                biquadraticCoords+="<tr><td>"+Number(getBiquadraticOutput[0][i][j]).toPrecision(3)+"</td><td>"+Number(getBiquadraticOutput[1][i][j]).toPrecision(3)+"</td><td>"+Number(getBiquadraticOutput[2][i][j]).toPrecision(5)+"</td></tr>";
            }
        }
    }
    document.getElementById("bqCoordTable").innerHTML+=biquadraticCoords;
}

function fillQuadraticCoordData(getURangeOutput,tableName){
    for(var i=0; i < getURangeOutput[0].length; i++){
        document.getElementById(tableName).innerHTML+="<tr><td>"+Number(getURangeOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getURangeOutput[1][i]).toPrecision(5)+"</td></tr>"
    }
}

function fillQuadraticCoordDataFast(getURangeOutput,tableName){
    var quadraticCoords;
    for(var i=0; i < getURangeOutput[0].length; i++){
        if(quadraticCoords == undefined){
            quadraticCoords="<tr><td>"+Number(getURangeOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getURangeOutput[1][i]).toPrecision(5)+"</td></tr>"           
        }
        else{
            quadraticCoords+="<tr><td>"+Number(getURangeOutput[0][i]).toPrecision(3)+"</td><td>"+Number(getURangeOutput[1][i]).toPrecision(5)+"</td></tr>"
        }
    }
    document.getElementById(tableName).innerHTML+=quadraticCoords;
}

function fillFdisCoordDataFast(getFdisOutput){
    var FdisCoords;
    for(var i=0; i < getFdisOutput[0].length; i++){
        var currArray = getFdisOutput[0][i];
        for(var j=0; j<currArray.length;j++){
            if(FdisCoords == undefined){
                FdisCoords="<tr><td>"+Number(getFdisOutput[0][i][j]).toPrecision(3)+"</td><td>"+Number(getFdisOutput[1][i][j]).toPrecision(3)+"</td><td>"+Number(getFdisOutput[2][i][j]).toPrecision(5)+"</td></tr>";
            }
            else{
                FdisCoords+="<tr><td>"+Number(getFdisOutput[0][i][j]).toPrecision(3)+"</td><td>"+Number(getFdisOutput[1][i][j]).toPrecision(3)+"</td><td>"+Number(getFdisOutput[2][i][j]).toPrecision(5)+"</td></tr>";
            }
        }
    }
    document.getElementById("fdisCoordTable").innerHTML+=FdisCoords;
}

//console.log("a1:",curvature_Analytical()[0]);
//console.log("a2:",curvature_Analytical()[1]);
//console.log("a3:",curvature_Analytical()[2]);
//console.log("Hb:",curvature_Analytical()[3]);
//console.log("Hx:",curvature_Analytical()[4]);
//console.log("Hc:",curvature_Analytical()[5]);
//console.log("smin:",curvature_Analytical()[6]);
//console.log("deltaGdef_analytical:",curvature_Analytical()[7]);
//console.log("biquadratic equation Gdef_RB:",biquadratic());
//console.log("u0_p for Gdef_RB=0:",biquadratic()[1]);
//console.log("u0_n for Gdef_RB=0:",biquadratic()[2]);
//console.log("besselK 0th order, x=4",besselK(0,4));
//console.log("besselK 1st order, x=4",besselK(1,4));
//console.log("besselK 0th order, x=4+0i",besselK(0,complex(4,0)));
//console.log("besselK 1st order, x=4+0i",besselK(1,complex(4,0)));
//console.log("besselK 0th order, x=4+1i",besselK(0,complex(4,1)));
//console.log("besselK 1st order, x=4+1i",besselK(1,complex(4,1)));

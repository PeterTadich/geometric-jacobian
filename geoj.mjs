//geoj = geometric jacobian

//ECMAScript module

//To do:
//    - remove svdClean()

import * as hlao from 'matrix-computations';
import * as ludcmp from 'lu-decomposition';
import * as svdcmp from 'singular-value-decomposition';
import * as mcht from 'homogeneous-transformations';

/*
%MATLAB
mdl_puma560
q = [0.6,0.19,0.75,0.91,0.11,0.71];
p560.jacob0(q)
%ans =
%    0.0745   -0.2910   -0.2237         0         0         0
%    0.1568   -0.1991   -0.1531         0         0         0
%   -0.0000    0.0873   -0.3367         0         0         0
%    0.0000    0.5646    0.5646   -0.6665    0.7309   -0.6463
%   -0.0000   -0.8253   -0.8253   -0.4560   -0.2436   -0.5472
%    1.0000    0.0000    0.0000    0.5898    0.6376    0.5318
*/

/*
%MATLAB
mdl_stanford
q = [0.6,0.19,0.75,0.91,0.11,0.71];
stanf.jacob0(q)
%ans =
%   -0.2330    0.8162    0.1559    0.0015    0.2534         0
%    0.0992    0.5584    0.1066    0.0286   -0.0211         0
%         0   -0.2134    0.9820   -0.0033   -0.0673         0
%   -0.0000   -0.5646         0    0.1559    0.0516    0.2632
%         0    0.8253         0    0.1066    0.9919    0.0984
%    1.0000    0.0000         0    0.9820   -0.1159    0.9597
*/

//see also:
//DUMMYgeometricJacobian(HTM) in I:\code\spatial_v2\js\graph\transformGraph.js
//geometricJacobian.js
function geometricJacobian(mua){
    /*
    //example of usage: three-link planar arm
    //example from - REF: Robotics Modelling, Planning and Control, Page 143 (section 3.7.5)
    var ai = 500.0; //units mm
    var vi1 = Math.PI; //units radians
    var vi2 = -1.0*Math.PI/2.0;
    var vi3 = -1.0*Math.PI/2.0;
    var Links = [
        [ai, 0.0, 0.0, vi1],
        [ai, 0.0, 0.0, vi2],
        [ai, 0.0, 0.0, vi3]
    ];
    var Jaco = geometricJacobian(Links);
    console.log(Jaco);

    Jaco = 
    [
       [[-500.0], [-500.0],   [0.0]],
       [[   0.0], [ 500.0], [500.0]],
       [[   0.0], [   0.0], [  0.0]],
       [[   0.0], [   0.0], [  0.0]],
       [[   0.0], [   0.0], [  0.0]],
       [[   1.0], [   1.0], [  1.0]]
    ]
    */
    
    // geometric Jacobian depends on:
    // - the manipulator configuration
    // - the frame in which the end-effector velocity is expressed
    
    /*
            -    -
    -   -   |zi-1| for a prismatic joint
    |JPi|   |   0|
    |   |   -    -
    |   | =                                               Equ. 3.30
    |   |   -                  -
    |JOi|   |zi-1 x (pe - pi-1)| for a revolute joint
    -   -   |              zi-1|
            -                  -
    
    where:
    zi-1 = R0,1(q1)...Ri-2,i-1(qi-1)z0, z0 = [0 0 1]T     Equ. 3.31
    pe~ = A0,1(q1)...An-1,n(qn)p0~, p0~ = [0 0 0 1]T      Equ. 3.32
    pi-1~ = A0,1(q1)...Ai-2,i-1(qi-1)P0~                  Equ. 3.33
    
    These equations are for the computation of the geometric Jacobian with respect to the base frame.
    
    REF: Robotics Modelling, Planning and Control, Page 112
    */
    
    var debug = 0;
    
    var Li = mua.DH;
    var n = Li.length;
    var jt = []; //joint type
    for(var i=1;i<=n;i=i+1) jt[i] = mua.jt[i-1];
    
    if(debug){
        var Ri_neg1i = "R0," + n + " = "; //Ri-1,i for i = 1 o n, transforms a vector from frame {i} to frame {i-1} (base frame)
        for(var i=1;i<=n;i=i+1){
            Ri_neg1i = Ri_neg1i + "R" + (i-1) + "," + i;
            if(i != n) Ri_neg1i = Ri_neg1i  + " * ";
        }
        console.log(Ri_neg1i);
    }
    
    var J = []; //the Jacobian matrix
    var z0 = [[0.0],[0.0],[1.0]]; //Allows selection of the third column - (equ. 3.31).
    var p0 = [[0.0],[0.0],[0.0],[1.0]]; //Allows the selection of the fourth column - (equ. 3.32).
                                        //REF: Robotics Modelling, Planning and Control, page 112 .
    
    //calc. pe (Equ. 3.32 REF: Robotics Modelling, Planning and Control, Page 112)
    for(var i=0;i<n;i=i+1){ //for all links
        if(i == 0) var T0n = mcht.Aij(Li[i]);
        else T0n = hlao.matrix_multiplication(T0n,mcht.Aij(Li[i])); //Aij - transforms a vector from frame {j} to frame {i}
    }                                                     //T0n - transforms a vector from frame {n} (link n) to frame {0} (link 0)
    var pe = hlao.matrix_multiplication(T0n,p0);
    //only select the first three
    pe.pop();
    if(debug) console.log(pe);
    //console.log(pe);
    
    for(var j=1;j<=n;j=j+1){ //for all links (start at link i+1)
        //calc. Ai-1,i as the following is based on this calculation:
        //  - zi-1
        //  - pi-1
        
        var R0ineg1 = hlao.identity_matrix(3);
        var A0ineg1 = hlao.identity_matrix(4);
        
        for(var i=0;i<j;i=i+1){ //now look at link i
            if(i == 0) var T0n = mcht.Aij(Li[i]);
            else T0n = hlao.matrix_multiplication(T0n,mcht.Aij(Li[i])); //Aij - transforms a vector from frame {j} to frame {i}
            if(i == (j - 2)){                                           //T0n - transforms a vector from frame {n} (link n) to frame {0} (link 0)
                var R0ineg1 = [ //R0,i-1
                    [T0n[0][0],T0n[0][1],T0n[0][2]],
                    [T0n[1][0],T0n[1][1],T0n[1][2]],
                    [T0n[2][0],T0n[2][1],T0n[2][2]]
                ];
                var A0ineg1 = T0n; //A0,i-1
            }
        }        
        
        var zineg1 = hlao.matrix_multiplication(R0ineg1,z0); //zi-1 (Equ. 3.31 REF: Robotics Modelling, Planning and Control, Page 112)
        var pineg1 = hlao.matrix_multiplication(A0ineg1,p0); //pi-1 (Equ. 3.33 REF: Robotics Modelling, Planning and Control, Page 113)
        //console.log(j);
        //console.log(T0n);
        
        //only select the first three
        pineg1.pop();
        
        if(debug){
            console.log(zineg1);
            console.log(pineg1);
        }
        
        //   - revolute
        if(jt[j].includes('R')){
            J[j-1] = hlao.vector_cross(zineg1,hlao.matrix_arithmetic(pe,pineg1,'-'));
            for(var i=0;i<3;i=i+1) J[j-1].push([zineg1[i][0]]);
        }
        //   - prismatic
        if(jt[j].includes('P')){
            J[j-1] = zineg1;
            for(var i=0;i<3;i=i+1) J[j-1].push([0.0]);
        }
        if(debug) console.log(J[j-1]);
    }
    
    J = hlao.matrix_transpose(J);
    //if(debug) console.log(size(J));
    return J;
}

//see also:
//I:\code\spatial_v2\js\RMC\RMC_torso.js "gJ()"
//I:\code\spatial_v2\js\inverseKinematics\geometricJacobian.js
//I:\code\spatial_v2\js\graph\transformGraph.js "DUMMYgeometricJacobian(HTM)"
//pass in the transforms defined in the inertial frame
/*
//TEST:
//   - see Robotics page 113
//   - I:\code\spatial_v2\js\inverseKinematics\geometricJacobian.js
//      - 'threeLinkPlanarArm()' returns Jacobian for three Link planar arm
//initialization
var qa = [];
//   - select a robot initial pose
qn.forEach(q => {qa.push([q])});
Links.forEach((q, i) => {Links[i][3] = qa[i][0]});

var JtoFix = geometricJacobian(Links); //ref: I:\code\spatial_v2\js\inverseKinematics\geometricJacobian.js
var J = [];
var Jstr = "";
for(var i=0;i<6;i=i+1){
    var row = [];
    for(var j=0;j<6;j=j+1){
        Jstr = Jstr + JtoFix[i][j][0].toFixed(4) + ',';
        row.push(JtoFix[i][j][0]);
    }
    Jstr = Jstr + '\n';
    J.push(row);
}
console.log(Jstr);

var TOn = CTdIF(Links);
var pe = directKinematicsDH(Links); //"directKinematicsDH()" from I:\code\spatial_v2\js\DenavitHartenberg\DenavitHartenberg.js use "Tb0" (base) and "Tne" (tool)
var pe = [[pe[0][3]],[pe[1][3]],[pe[2][3]]];
var Jcheck = gJ(TOn,pe);
printJacobian(Jcheck);
*/
function gJ(T0,pe,jt){ //T0 is "T" "zero"
    var n = T0.length;
    
    //check for homogeneous representation
    if(pe.length === 4) var pe = [[pe[0][0]],[pe[1][0]],[pe[2][0]]];
    
    //revolute equ. (3.30) (REF: Robotics Modelling, Planning and Control, Page 112 (section 3.1.3))
    var JP = []; //contribution to the linear velocity
    var JO = []; //contribution to the angular velocity
    var J = []; //Jacobian
    var z = []; //unit vector of the joint axis
    var p = []; //position vector of the joint
    
    //var pe = [[openLink.endEffector.T0[0][3]], [openLink.endEffector.T0[1][3]], [openLink.endEffector.T0[2][3]]];
    
    //base frame
    T0.unshift(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0]
        ]
    );
    
    //IMPORTANT: push tool tip
    //T0.push(tooltip);
    
    /*
    for(var i=0;i<n;i=i+1){
        //print_homogenous_transform(T0[i]);
        p.push([[T0[i][0][3]], [T0[i][1][3]], [T0[i][2][3]]]); //matrix_multiplication(T0[i],p0) where p0 = [[0.0],[0.0],[0.0],[1.0]], see (Equ. 3.33)
        z.push(                                                //matrix_multiplication(T0[i],z0) where z0 = [[0.0],[0.0],[1.0]], see (equ. 3.31)
            matrix_multiplication( //Not needed. Just select [T0[i][0][2], T0[i][1][2], T0[i][2][2]]
                    [
                        [[T0[i][0][0]], [T0[i][0][1]], [T0[i][0][2]]],
                        [[T0[i][1][0]], [T0[i][1][1]], [T0[i][1][2]]],
                        [[T0[i][2][0]], [T0[i][2][1]], [T0[i][2][2]]]
                    ],
                    [[0.0],[0.0],[1.0]]
                )
            );
    }
    */
    
    for(var i=0;i<n;i=i+1){
        p.push([[T0[i][0][3]], [T0[i][1][3]], [T0[i][2][3]]]); //matrix_multiplication(T0[i],p0) where p0 = [[0.0],[0.0],[0.0],[1.0]], see (Equ. 3.33)
        z.push([[T0[i][0][2]], [T0[i][1][2]], [T0[i][2][2]]]); //matrix_multiplication(T0[i],z0) where z0 = [[0.0],[0.0],[1.0]], see (equ. 3.31)
    }
    
    T0.shift();
    
    for(var i=1;i<=n;i=i+1){
        //   - revolute
        if(jt[i-1].includes('R')){
            //   - linear velocity
            JP[i] = hlao.vector_cross(z[i-1],hlao.matrix_arithmetic(pe,p[i-1],'-'));
            //   - angular velocity
            JO[i] = z[i-1];
        }
        //   - prismatic
        if(jt[i-1].includes('P')){
            //   - linear velocity
            JP[i] = z[i-1];
            //   - angular velocity
            JO[i] = [[0.0],[0.0],[0.0]];
        }
        
        J.push([JP[i][0][0],JP[i][1][0],JP[i][2][0],JO[i][0][0],JO[i][1][0],JO[i][2][0]]); //row vector
    }
    
    //console.log(J);
    J = hlao.matrix_transpose(J); //convert each row vector into a column vector
    //console.log(J);
    
    return J;
}

//Jacobian inversion (returns J^-1)
//IMPORTANT: better to use svdcmp() to inspect singular values
//see also I:\code\spatial_v2\js\RMC\RMC_torso.js
function JacobianInverse(J){
    var m; var n;
    [m,n] = [J.length,J[0].length];
    //var dim = size(J); //get the dimension of the Jacobian matrix
    //var m = dim[0]; //number of rows
    
    //   - adjust matrix for inversion (add dummy zeroes to the start of each row)
    for(var j=0;j<m;j=j+1){ //For each row of J[] add an extra element '0.0' to the beginning of the array (beginning of the row).
        J[j].unshift(0.0);
    }
    //   - create an array of zeroes (row vector) to the J[]
    var offsetRow = [];
    for(var j=0;j<(m+1);j=j+1){ //only 6 + 1 elements as J[] is 6 x 6 matrix or 6 x 7 matrix with padded zeroes
        offsetRow.push(0.0);
    }
    J.unshift(offsetRow); //add the row vector of zeroes to the beginning of J[] array - now a 7 x 7 matrix
    //console.log('Adjusted Jacobian (padded with zeroes):');
    //printJacobian(J);
    
    //   - calc. J^-1
    var Jinverse = ludcmp.matrixInverseLU(J,m);
    
    //   - re-adjust the matrix
    Jinverse.shift();
    
    //   - quick print
    //for(var j=0;j<m;j=j+1){
    //    console.log(Jinverse[j]);
    //}
    
    for(var j=0;j<m;j=j+1){
        Jinverse[j].shift(); //For each row of Jinverse[] remove the first element.
        //console.log(Jinverse[j]);
    }
    //console.log(Jinverse);
    
    return Jinverse;
}

//see also I:\code\spatial_v2\js\RMC\RMC_torso.js
function JacobianInverse_svdcmp(J){
    var debug = 0;
    
    //var dim = size(J);
    //var m = dim[0]; //matrix rank - number of independent rows (number of rows)
    //var n = dim[1]; //number of columns
    var m; var n;
    [m,n] = [J.length,J[0].length];
    if(debug) console.log('m: ' + m + ', n: ' + n);

    var Jnull = [];
    for(var j=0;j<m;j=j+1){
        var rowData = [];
            for(var k=0;k<n;k=k+1){
                rowData.push(J[j][k]);
            }
        Jnull.push(rowData);
    }
    if(debug) console.log('Jnull:');
    if(debug) printJacobian(Jnull);
    
    //   - adjust matrix for inversion
    for(var j=0;j<m;j=j+1){ //For each row of Jnull[] add an extra element '0.0' to the beginning of the array (beginning of the row).
        Jnull[j].unshift(0.0);
    }
    //   - create an array of zeroes (row vector)
    var offsetRow = [];
    for(var j=0;j<(n+1);j=j+1){
        offsetRow.push(0.0);
    }
    Jnull.unshift(offsetRow); //add the row vector of zeroes to the beginning of Jnull[] array - now a 7 x 7 matrix

    //if(debug) console.log(size(Jnull));
    if(debug) console.log('Jnull adjusted:');
    if(debug) printJacobian(Jnull);

    var w = hlao.zeros_vector((n+1),'row'); //row vector where index = 0 is undefined
    var v = hlao.zeros_matrix((n+1),(n+1)); //matrix
    var uwv = svdcmp.svdcmp(Jnull, m, n, w, v);
    //if(debug) console.log(uwv);

    var U = svdcmp.svdClean(uwv[0]); //Drop the first element in the array as it is zero.
    var S = svdcmp.svdClean(uwv[1]); //W
    var V = svdcmp.svdClean(uwv[2]);
    if(debug){
        console.log('U:');
        print_multi_array(U);
        console.log('S (' + S.length + '):');
        print_multi_array(S);
        console.log('V:');
        print_multi_array(V);
    }
    
    //run diagnostics
    //   - Get the singular values (positive or zero).
    var w = uwv[1]; //IMPORTANT: this also alters 'S'
    //w.shift(); // Drop the first element in the array as it is zero. IMPORTANT: not required due to 'S' using 'svdClean()'
    if(debug) console.log('W = ' + w);

    //   - Get the rank.
    var rank = hlao.matrix_rank(w);
    if(debug) console.log('Matrix rank = ' + rank);

    //   - Get the condition number of the matrix.
    var w_min = +1e6;
    var w_max = -1e6;
    for(var j=0;j<w.length;j=j+1){
        if(w_max < w[j]) w_max = w[j];
        if((w_min > w[j]) && (w[j] > 1e-6)) w_min = w[j]; // w_min must be non zero.
    }
    var cond_numb = w_max/w_min;
    if(debug) console.log('sigma 1 = ' + w_max);
    if(debug) console.log('sigma r = ' + w_min);
    if(debug) console.log('Matrix condition number = ' + cond_numb);
    if(cond_numb > 900.0) alert('Condition number: ' + cond_numb);
    
    //Calc. inverse
    var diag = [];
    var TOL = 1e-12; //Tolerance for the evaluation of 'S'.
    if(debug) console.log('S (' + S.length + '):');
    for(var i=0;i<S.length;i=i+1){ //IMPORTANT: Magic number.
        if(Math.abs(S[i]) > (0.0 + TOL)) diag[i] = 1.0/S[i];
        else diag[i] = 0.0;
    }
    if(debug) console.log(diag);
    var Sinv = [ //IMPORTANT: fix magic number.
        [diag[0],     0.0,     0.0,     0.0,     0.0,     0.0],
        [    0.0, diag[1],     0.0,     0.0,     0.0,     0.0],
        [    0.0,     0.0, diag[2],     0.0,     0.0,     0.0],
        [    0.0,     0.0,     0.0, diag[3],     0.0,     0.0],
        [    0.0,     0.0,     0.0,     0.0, diag[4],     0.0],
        [    0.0,     0.0,     0.0,     0.0,     0.0, diag[5]]
    ];
    if(debug) console.log(Sinv);
    
    var JJTinv = hlao.matrix_multiplication(
            hlao.matrix_multiplication(
                V,Sinv
            ),
            hlao.matrix_transpose(U)
        );
    if(debug){
        console.log('(J x JT)-1:');
        print_multi_array(JJTinv);
    }
    
    //return matrix_multiplication(U,matrix_transpose(V));
    return JJTinv;
}

/*
//see I:\code\spatial_v2\js\FIT5147\rot3dfit.js
function svdClean(A){
    var dim = size(A);
    var m = dim[0]; //Rows.
    var n = dim[0]; //Columns.
    //console.log('m: ' + m + ', n:' + n);
    
    if(m > 1){
        for(var i=0;i<m;i=i+1){ //For each row of 'A' drop the extra element '0.0' at the beginning of the array (beginning of the row).
            A[i].shift();
        }
        A.shift(); //Drop the first row.
    } else {
        A.shift(); //It is a row vector hence just drop the first element.
    }
    
    return A;
}
*/

//print the Jacobian
//see also I:\code\spatial_v2\js\RMC\RMC_torso.js
function printJacobian(J){
    //var dim = size(J); //get the dimension of the Jacobian matrix
    //var m = dim[0]; //number of rows
    //var n = dim[1]; //number of columns
    var m; var n;
    [m,n] = [J.length,J[0].length];
    
    var Jstr = "";
    for(var j=0;j<m;j=j+1){ //row
        for(var k=0;k<n;k=k+1){ //col
            Jstr = Jstr + J[j][k].toFixed(4) + ' ';
        }
        Jstr = Jstr + '\n';
    }
    console.log(Jstr);
}

//see I:\code\spatial_v2\js\FischerTechnik\industry_robots_rob3.js
function print_multi_array(A){
    var m = A.length;
    var n = A[0].length;
    var str = "";
    
    if(typeof n == 'undefined'){
        //row vector
        m = 1;
        n = A.length;
        str = str + '[';
        for(var j=0;j<n;j=j+1){ //column
            str = str + A[j];
            if(j < n-1) str = str + ',';
        }
        str = str + '];';
    } else {
        str = str + '[\n';
        for(var i=0;i<m;i=i+1){ //row
            str = str + '    [';
            for(var j=0;j<n;j=j+1){ //column
                str = str + A[i][j][0].toFixed(4);
                if(j < n-1) str = str + ',';
            }
            str = str + ']';
            if(i < m-1) str = str + ',';
            str = str + '\n';
        }
        str = str + '];';
    }
    
    console.log(str);
}

export {
    geometricJacobian,
    gJ,
    JacobianInverse,
    JacobianInverse_svdcmp,
    printJacobian,
    print_multi_array
};
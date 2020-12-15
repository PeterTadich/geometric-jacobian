# geometric-jacobian
Geometric Jacobian

## Dependencies

There are 4 dependencies 'matrix-computations', 'singular-value-decomposition', 'lu-decomposition' and 'homogeneous-transformations'.

```bash
https://github.com/PeterTadich/matrix-computations
https://github.com/PeterTadich/singular-value-decomposition
https://github.com/PeterTadich/lu-decomposition
https://github.com/PeterTadich/homogeneous-transformations
```

## Installation

### Node.js

```bash
npm install https://github.com/PeterTadich/geometric-jacobian
```

### Google Chrome Web browser

No installation required for the Google Chrome Web browser.

## How to use

### Node.js

```js
import * as geoj from 'geometric-jacobian';
```

### Google Chrome Web browser

```js
import * as geoj from './geoj.mjs';
```

## Examples

### Node.js (server side)

#### geoj.geometricJacobian()

Get the Jacobian matrix - extended version. Copy the following code to index.mjs

```js
//npm install https://github.com/PeterTadich/industrial-robots https://github.com/PeterTadich/geometric-jacobian

import * as mcir from 'industrial-robots';
import * as geoj from 'geometric-jacobian';

var mua = mcir.puma560(); //get the robot
var J = geoj.geometricJacobian(mua); //calc geometric Jacobian with default q's
geoj.print_multi_array(J); //multi-dimensional array
```

Then run:

```bash
npm init -y
npm install https://github.com/PeterTadich/industrial-robots https://github.com/PeterTadich/geometric-jacobian
node index.mjs
```

If the above does not work modify the package.json file as follows:
Helpful ref: [https://stackoverflow.com/questions/45854169/how-can-i-use-an-es6-import-in-node-js](https://stackoverflow.com/questions/45854169/how-can-i-use-an-es6-import-in-node-js)

```js
"scripts": {
    "test": "echo \"Error: no test specified\" && exit 1",
    "start": "node --experimental-modules index.mjs"
  },
"type": "module",
```

```bash
npm start
```

Result:

```js
[
    [0.0745,-0.2910,-0.2237, 0.0000, 0.0000, 0.0000],
    [0.1567,-0.1991,-0.1531, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.0873,-0.3367, 0.0000, 0.0000, 0.0000],
    [0.0000, 0.5646, 0.5646,-0.6665, 0.7309,-0.6463],
    [0.0000,-0.8253,-0.8253,-0.4560,-0.2436,-0.5472],
    [1.0000, 0.0000, 0.0000, 0.5898, 0.6376, 0.5318]
];
```

#### geoj.gJ()

Get the Jacobian matrix (using T0, pe). Copy the following code to index.mjs

```js
//npm install https://github.com/PeterTadich/industrial-robots https://github.com/PeterTadich/geometric-jacobian https://github.com/PeterTadich/homogeneous-transformations

import * as mcir from 'industrial-robots';
import * as geoj from 'geometric-jacobian';
import * as mcht from 'homogeneous-transformations';

var T0n = mcht.CTdIF(mua.DH); //link transforms w.r.t inertial frame
var pe = mcht.directKinematicsDH(mua.DH); //calc position of end-effector in inertial frame
var pe = [[pe[0][3]],[pe[1][3]],[pe[2][3]]];
var J = geoj.gJ(T0n,pe,mua.jt); //jt - joint types
geoj.printJacobian(J); //two-dimensional array
```

Result:

```js
    0.0745 -0.2910 -0.2237  0.0000  0.0000  0.0000
    0.1567 -0.1991 -0.1531  0.0000  0.0000  0.0000
    0.0000  0.0873 -0.3367  0.0000  0.0000  0.0000
    0.0000  0.5646  0.5646 -0.6665  0.7309 -0.6463
    0.0000 -0.8253 -0.8253 -0.4560 -0.2436 -0.5472
    1.0000  0.0000  0.0000  0.5898  0.6376  0.5318
```

#### geoj.JacobianInverse_svdcmp()

Run diagnostic's. Get the rank and condition number of the Jacobian matrix. Copy the following code to index.mjs

```js
//npm install https://github.com/PeterTadich/industrial-robots https://github.com/PeterTadich/geometric-jacobian

import * as mcir from 'industrial-robots';
import * as geoj from 'geometric-jacobian';

var mua = mcir.puma560(); //get the robot
var J = geoj.geometricJacobian(mua);

//convert multi to 2 dimensional array
var m; var n;
[m,n] = [J.length,J[0].length];
var Jc = [];
for(var i=0;i<m;i=i+1){
    Jc[i] = [];
    for(var j=0;j<n;j=j+1){
        Jc[i][j] = J[i][j][0];
    }
}
var Jinv = geoj.JacobianInverse_svdcmp(Jc,1); //set 'Dx' = 1 (diagnostic)
geoj.printJacobian(Jinv);
```

Result:

```js
Matrix rank = 6
sigma 1 = 1.6725
sigma r = 0.0231
Matrix condition number = 72.3347
Jacobian inverse:
    -6.4679   9.4542  0.0000  0.0000 -0.0000 -0.0000
    -4.2459   2.0182  1.9035  0.0000 -0.0000 -0.0000
    -1.1008   0.5232 -2.4762 -0.0000 -0.0000 -0.0000
    -5.3793 -29.8352 -4.0937 -1.9978  7.2944  5.0774
     7.4052  -7.5875  0.3515  0.7309 -0.2436  0.6376
     9.2500  24.4067  4.1186  1.3394 -7.7975 -4.5149
```
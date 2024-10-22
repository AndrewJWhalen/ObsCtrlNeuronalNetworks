(* ::Package:: *)

(* NonLinear Controllability for a 3 Node FN Network with Sigmoidal \
input activation, Measuring each node and statistical parameter \
averaging *)

(*from the paper: \
 Observability and Controllability of Nonlinear Networks: The Role of Symmetry \
 Andrew J. Whalen, Sean N. Brennan, Timothy D. Sauer, and Steven J. Schiff \
 Physical Review X, 2014 \
Please cite the paper if you make use of or modify this code. \
11/18/2014 \
 Andrew Whalen \
 Center for Neural Engineering \
 Penn State University \
 email: andrew.whalen@yale.edu (Updated email as of 1/1/2023) *)

IDcoupling=0;
ICflag=1;

For[k = 1, k < 4, k++,

For[j = 1, j < 21, j++, 

  (* Import the Matlab sim data *)
  If [IDcoupling == 1,
  X1 = Import["/MatlabData/ChaosHetIC1IdK3NodeFNCoupling" <> ToString[j] <> ".mat"];
  X2 = Import["/MatlabData/LCHetIC1IdK3NodeFNCoupling" <> ToString[j] <> ".mat"];
  X3 = Import["/MatlabData/CHetIC1IdK3NodeFNCoupling" <> ToString[j] <> ".mat"];
  ];

  If [IDcoupling == 0,
  X1 = Import["/MatlabData/ChaosHetIC1HetK3NodeFNCoupling" <> ToString[j] <> ".mat"];
  X2 = Import["/MatlabData/LCHetIC1HetK3NodeFNCoupling" <> ToString[j] <> ".mat"];
  X3 = Import["/MatlabData/CHetIC1HetK3NodeFNCoupling" <> ToString[j] <> ".mat"];
  ];

  (* Identical Nodal dynamics *)
  p1 = 0.7; p2 = 0.8; p3 = 10.0;
  q1 = 0.7; q2 = 0.8; q3 = 10.0;
  r1 = 0.7; r2 = 0.8; r3 = 10.0; 

If [IDcoupling == 1,
  If[j == 1, d12 = d21 = d23 = d32 = d31 = d13 = 0.004000000000000;];
  If[j == 2, d12 = d21 = d23 = d32 = d31 = d13 = 0.0564157894736842;];
  If[j == 3, d12 = d21 = d23 = d32 = d31 = d13 = 0.108831578947368;];
  If[j == 4, d12 = d21 = d23 = d32 = d31 = d13 = 0.161247368421053;];
  If[j == 5, d12 = d21 = d23 = d32 = d31 = d13 = 0.213663157894737;];
  If[j == 6, d12 = d21 = d23 = d32 = d31 = d13 = 0.266078947368421;];
  If[j == 7, d12 = d21 = d23 = d32 = d31 = d13 = 0.318494736842105;];
  If[j == 8, d12 = d21 = d23 = d32 = d31 = d13 = 0.370910526315789;];
  If[j == 9, d12 = d21 = d23 = d32 = d31 = d13 = 0.423326315789474;];
  If[j == 10, d12 = d21 = d23 = d32 = d31 = d13 = 0.475742105263158;];
  If[j == 11, d12 = d21 = d23 = d32 = d31 = d13 = 0.528157894736842;];
  If[j == 12, d12 = d21 = d23 = d32 = d31 = d13 = 0.580573684210526;];
  If[j == 13, d12 = d21 = d23 = d32 = d31 = d13 = 0.632989473684211;];
  If[j == 14, d12 = d21 = d23 = d32 = d31 = d13 = 0.685405263157895;];
  If[j == 15, d12 = d21 = d23 = d32 = d31 = d13 = 0.737821052631579;];
  If[j == 16, d12 = d21 = d23 = d32 = d31 = d13 = 0.790236842105263;];
  If[j == 17, d12 = d21 = d23 = d32 = d31 = d13 = 0.842652631578947;];
  If[j == 18, d12 = d21 = d23 = d32 = d31 = d13 = 0.895068421052632;];
  If[j == 19, d12 = d21 = d23 = d32 = d31 = d13 = 0.947484210526316;];
  If[j == 20, d12 = d21 = d23 = d32 = d31 = d13 = 0.999900000000000;];
];

If [IDcoupling == 0,
  rrr={0.10,-0.10,0.20,-0.20,0.30};
  If[j == 1, d32 = 0.004000000000000;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 2, d32 = 0.0564157894736842;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 3, d32 = 0.108831578947368;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 4, d32 = 0.161247368421053;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 5, d32 = 0.213663157894737;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 6, d32 = 0.266078947368421;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 7, d32 = 0.318494736842105;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 8, d32 = 0.370910526315789;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 9, d32 = 0.423326315789474;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 10, d32 = 0.475742105263158;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 11, d32 = 0.528157894736842;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 12, d32 = 0.580573684210526;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 13, d32 = 0.632989473684211;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 14, d32 = 0.685405263157895;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 15, d32 = 0.737821052631579;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 16, d32 = 0.790236842105263;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 17, d32 = 0.842652631578947;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 18, d32 = 0.895068421052632;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 19, d32 = 0.947484210526316;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
  If[j == 20, d32 = 0.999900000000000;d12=d32+rrr[[1]]*d32;d21=d32+rrr[[2]]*d32;d13=d32+rrr[[3]]*d32;d31=d32+rrr[[4]]*d32;d23=d32+rrr[[5]]*d32;];
];

(* Threshold Function - Hyperbolic Tangent out[0,1], in[-2,2] *)
  \[Kappa] = 1;
  h = 0;
  w = 1/4;
  
  (* Control Variables *)
  contr = k; (* measured node *)
  
  (************************************************************)
  (************************ MOTIF 1 ***************************)
  (************************************************************)
  
  (* NonLinear Vector Field *)
  (* Standard Form Equations *)
  f1 = p3*(v1 - v1^3/3 + w1 + d21*\[Kappa]/2*(Tanh[(v2 - h)/(2*w)] + 1) + d31*\[Kappa]/2*(Tanh[(v3 - h)/(2*w)] + 1));
  f2 = -(v1 - p1 + p2*w1)/p3;
  f3 = q3*(v2 - v2^3/3 + w2 + d12*\[Kappa]/2*(Tanh[(v1 - h)/(2*w)] + 1) + d32*\[Kappa]/2*(Tanh[(v3 - h)/(2*w)] + 1));
  f4 = -(v2 - q1 + q2*w2)/q3;
  f5 = r3*(v3 - v3^3/3 + w3 + d13*\[Kappa]/2*(Tanh[(v1 - h)/(2*w)] + 1) + d23*\[Kappa]/2*(Tanh[(v2 - h)/(2*w)] + 1));
  f6 = -(v3 - r1 + r2*w3)/r3;

 (* Chaos Form Equations *)
  f1c = p3*(v1 - v1^3/3 - w1 + d21*\[Kappa]/2*(Tanh[(v2 - h)/(2*w)] + 1) + d31*\[Kappa]/2*(Tanh[(v3 - h)/(2*w)] + 1));
  f2c = v1 - p2*w1 + p1;
  f3c = q3*(v2 - v2^3/3 - w2 + d12*\[Kappa]/2*(Tanh[(v1 - h)/(2*w)] + 1) + d32*\[Kappa]/2*(Tanh[(v3 - h)/(2*w)] + 1));
  f4c = v2 - q2*w2 + q1;
  f5c = r3*(v3 - v3^3/3 - w3 + d13*\[Kappa]/2*(Tanh[(v1 - h)/(2*w)] + 1) + d23*\[Kappa]/2*(Tanh[(v2 - h)/(2*w)] + 1));
  f6c = v3 - r2*w3 + r1;
 
F = {f1, f2, f3, f4, f5, f6};
Fc = {f1c, f2c, f3c, f4c, f5c, f6c};

(* Input Vector *)
If [contr == 1, g1 = 1; g2 = 0; g3 = 0; g4 = 0; g5 = 0; g6 = 0;];
If [contr == 2, g1 = 0; g2 = 0; g3 = 1; g4 = 0; g5 = 0; g6 = 0;];
If [contr == 3, g1 = 0; g2 = 0; g3 = 0; g4 = 0; g5 = 1; g6 = 0;];
G = {g1, g2, g3, g4, g5, g6};

(* Calculate the Jacobian of the Nonlinear Vector Field and input vector *)
DFdx = {{D[f1, v1], D[f1, w1], D[f1, v2], D[f1, w2], D[f1, v3], 
    D[f1, w3]}, {D[f2, v1], D[f2, w1], D[f2, v2], D[f2, w2], 
    D[f2, v3], D[f2, w3]}, {D[f3, v1], D[f3, w1], D[f3, v2], 
    D[f3, w2], D[f3, v3], D[f3, w3]}, {D[f4, v1], D[f4, w1], 
    D[f4, v2], D[f4, w2], D[f4, v3], D[f4, w3]}, {D[f5, v1], 
    D[f5, w1], D[f5, v2], D[f5, w2], D[f5, v3], 
    D[f5, w3]}, {D[f6, v1], D[f6, w1], D[f6, v2], D[f6, w2], 
    D[f6, v3], D[f6, w3]}};

DFcdx = {{D[f1c, v1], D[f1c, w1], D[f1c, v2], D[f1c, w2], D[f1c, v3], 
    D[f1c, w3]}, {D[f2c, v1], D[f2c, w1], D[f2c, v2], D[f2c, w2], 
    D[f2c, v3], D[f2c, w3]}, {D[f3c, v1], D[f3c, w1], D[f3c, v2], 
    D[f3c, w2], D[f3c, v3], D[f3c, w3]}, {D[f4c, v1], D[f4c, w1], 
    D[f4c, v2], D[f4c, w2], D[f4c, v3], D[f4c, w3]}, {D[f5c, v1], 
    D[f5c, w1], D[f5c, v2], D[f5c, w2], D[f5c, v3], 
    D[f5c, w3]}, {D[f6c, v1], D[f6c, w1], D[f6c, v2], D[f6c, w2], 
    D[f6c, v3], D[f6c, w3]}};

Dgdx = {{D[g1, v1], D[g1, w1], D[g1, v2], D[g1, w2], D[g1, v3], 
    D[g1, w3]}, {D[g2, v1], D[g2, w1], D[g2, v2], D[g2, w2], 
    D[g2, v3], D[g2, w3]}, {D[g3, v1], D[g3, w1], D[g3, v2], 
    D[g3, w2], D[g3, v3], D[g3, w3]}, {D[g4, v1], D[g4, w1], 
    D[g4, v2], D[g4, w2], D[g4, v3], D[g4, w3]}, {D[g5, v1], 
    D[g5, w1], D[g5, v2], D[g5, w2], D[g5, v3], 
    D[g5, w3]}, {D[g6, v1], D[g6, w1], D[g6, v2], D[g6, w2], 
    D[g6, v3], D[g6, w3]}};

(* Lie Brackets *)
LB1 = Dgdx.F - DFdx.G;
LB1c = Dgdx.Fc - DFcdx.G;

DLB1 = {{D[LB1[[1]], v1], D[LB1[[1]], w1], D[LB1[[1]], v2], 
    D[LB1[[1]], w2], D[LB1[[1]], v3], D[LB1[[1]], w3]},
   {D[LB1[[2]], v1], D[LB1[[2]], w1], D[LB1[[2]], v2], 
    D[LB1[[2]], w2], D[LB1[[2]], v3], 
    D[LB1[[2]], w3]}, {D[LB1[[3]], v1], D[LB1[[3]], w1], 
    D[LB1[[3]], v2], D[LB1[[3]], w2], D[LB1[[3]], v3], 
    D[LB1[[3]], w3]}, {D[LB1[[4]], v1], D[LB1[[4]], w1], 
    D[LB1[[4]], v2], D[LB1[[4]], w2], D[LB1[[4]], v3], 
    D[LB1[[4]], w3]}, {D[LB1[[5]], v1], D[LB1[[5]], w1], 
    D[LB1[[5]], v2], D[LB1[[5]], w2], D[LB1[[5]], v3], 
    D[LB1[[5]], w3]}, {D[LB1[[6]], v1], D[LB1[[6]], w1], 
    D[LB1[[6]], v2], D[LB1[[6]], w2], D[LB1[[6]], v3], 
    D[LB1[[6]], w3]}};

DLB1c = {{D[LB1c[[1]], v1], D[LB1c[[1]], w1], D[LB1c[[1]], v2], 
    D[LB1c[[1]], w2], D[LB1c[[1]], v3], D[LB1c[[1]], w3]},
   {D[LB1c[[2]], v1], D[LB1c[[2]], w1], D[LB1c[[2]], v2], 
    D[LB1c[[2]], w2], D[LB1c[[2]], v3], 
    D[LB1c[[2]], w3]}, {D[LB1c[[3]], v1], D[LB1c[[3]], w1], 
    D[LB1c[[3]], v2], D[LB1c[[3]], w2], D[LB1c[[3]], v3], 
    D[LB1c[[3]], w3]}, {D[LB1c[[4]], v1], D[LB1c[[4]], w1], 
    D[LB1c[[4]], v2], D[LB1c[[4]], w2], D[LB1c[[4]], v3], 
    D[LB1c[[4]], w3]}, {D[LB1c[[5]], v1], D[LB1c[[5]], w1], 
    D[LB1c[[5]], v2], D[LB1c[[5]], w2], D[LB1c[[5]], v3], 
    D[LB1c[[5]], w3]}, {D[LB1c[[6]], v1], D[LB1c[[6]], w1], 
    D[LB1c[[6]], v2], D[LB1c[[6]], w2], D[LB1c[[6]], v3], 
    D[LB1c[[6]], w3]}};

LB2 = DLB1.F - DFdx.LB1;
LB2c = DLB1c.Fc - DFcdx.LB1c;

DLB2 = {{D[LB2[[1]], v1], D[LB2[[1]], w1], D[LB2[[1]], v2], 
    D[LB2[[1]], w2], D[LB2[[1]], v3], D[LB2[[1]], w3]},
   {D[LB2[[2]], v1], D[LB2[[2]], w1], D[LB2[[2]], v2], 
    D[LB2[[2]], w2], D[LB2[[2]], v3], 
    D[LB2[[2]], w3]}, {D[LB2[[3]], v1], D[LB2[[3]], w1], 
    D[LB2[[3]], v2], D[LB2[[3]], w2], D[LB2[[3]], v3], 
    D[LB2[[3]], w3]}, {D[LB2[[4]], v1], D[LB2[[4]], w1], 
    D[LB2[[4]], v2], D[LB2[[4]], w2], D[LB2[[4]], v3], 
    D[LB2[[4]], w3]}, {D[LB2[[5]], v1], D[LB2[[5]], w1], 
    D[LB2[[5]], v2], D[LB2[[5]], w2], D[LB2[[5]], v3], 
    D[LB2[[5]], w3]}, {D[LB2[[6]], v1], D[LB2[[6]], w1], 
    D[LB2[[6]], v2], D[LB2[[6]], w2], D[LB2[[6]], v3], 
    D[LB2[[6]], w3]}};

DLB2c = {{D[LB2c[[1]], v1], D[LB2c[[1]], w1], D[LB2c[[1]], v2], 
    D[LB2c[[1]], w2], D[LB2c[[1]], v3], D[LB2c[[1]], w3]},
   {D[LB2c[[2]], v1], D[LB2c[[2]], w1], D[LB2c[[2]], v2], 
    D[LB2c[[2]], w2], D[LB2c[[2]], v3], 
    D[LB2c[[2]], w3]}, {D[LB2c[[3]], v1], D[LB2c[[3]], w1], 
    D[LB2c[[3]], v2], D[LB2c[[3]], w2], D[LB2c[[3]], v3], 
    D[LB2c[[3]], w3]}, {D[LB2c[[4]], v1], D[LB2c[[4]], w1], 
    D[LB2c[[4]], v2], D[LB2c[[4]], w2], D[LB2c[[4]], v3], 
    D[LB2c[[4]], w3]}, {D[LB2c[[5]], v1], D[LB2c[[5]], w1], 
    D[LB2c[[5]], v2], D[LB2c[[5]], w2], D[LB2c[[5]], v3], 
    D[LB2c[[5]], w3]}, {D[LB2c[[6]], v1], D[LB2c[[6]], w1], 
    D[LB2c[[6]], v2], D[LB2c[[6]], w2], D[LB2c[[6]], v3], 
    D[LB2c[[6]], w3]}};

LB3 = DLB2.F - DFdx.LB2;
LB3c = DLB2c.Fc - DFcdx.LB2c;

DLB3 = {{D[LB3[[1]], v1], D[LB3[[1]], w1], D[LB3[[1]], v2], 
    D[LB3[[1]], w2], D[LB3[[1]], v3], D[LB3[[1]], w3]},
   {D[LB3[[2]], v1], D[LB3[[2]], w1], D[LB3[[2]], v2], 
    D[LB3[[2]], w2], D[LB3[[2]], v3], 
    D[LB3[[2]], w3]}, {D[LB3[[3]], v1], D[LB3[[3]], w1], 
    D[LB3[[3]], v2], D[LB3[[3]], w2], D[LB3[[3]], v3], 
    D[LB3[[3]], w3]}, {D[LB3[[4]], v1], D[LB3[[4]], w1], 
    D[LB3[[4]], v2], D[LB3[[4]], w2], D[LB3[[4]], v3], 
    D[LB3[[4]], w3]}, {D[LB3[[5]], v1], D[LB3[[5]], w1], 
    D[LB3[[5]], v2], D[LB3[[5]], w2], D[LB3[[5]], v3], 
    D[LB3[[5]], w3]}, {D[LB3[[6]], v1], D[LB3[[6]], w1], 
    D[LB3[[6]], v2], D[LB3[[6]], w2], D[LB3[[6]], v3], 
    D[LB3[[6]], w3]}};

DLB3c = {{D[LB3c[[1]], v1], D[LB3c[[1]], w1], D[LB3c[[1]], v2], 
    D[LB3c[[1]], w2], D[LB3c[[1]], v3], D[LB3c[[1]], w3]},
   {D[LB3c[[2]], v1], D[LB3c[[2]], w1], D[LB3c[[2]], v2], 
    D[LB3c[[2]], w2], D[LB3c[[2]], v3], 
    D[LB3c[[2]], w3]}, {D[LB3c[[3]], v1], D[LB3c[[3]], w1], 
    D[LB3c[[3]], v2], D[LB3c[[3]], w2], D[LB3c[[3]], v3], 
    D[LB3c[[3]], w3]}, {D[LB3c[[4]], v1], D[LB3c[[4]], w1], 
    D[LB3c[[4]], v2], D[LB3c[[4]], w2], D[LB3c[[4]], v3], 
    D[LB3c[[4]], w3]}, {D[LB3c[[5]], v1], D[LB3c[[5]], w1], 
    D[LB3c[[5]], v2], D[LB3c[[5]], w2], D[LB3c[[5]], v3], 
    D[LB3c[[5]], w3]}, {D[LB3c[[6]], v1], D[LB3c[[6]], w1], 
    D[LB3c[[6]], v2], D[LB3c[[6]], w2], D[LB3c[[6]], v3], 
    D[LB3c[[6]], w3]}};

LB4 = DLB3.F - DFdx.LB3;
LB4c = DLB3c.Fc - DFcdx.LB3c;

DLB4 = {{D[LB4[[1]], v1], D[LB4[[1]], w1], D[LB4[[1]], v2], 
    D[LB4[[1]], w2], D[LB4[[1]], v3], D[LB4[[1]], w3]},
   {D[LB4[[2]], v1], D[LB4[[2]], w1], D[LB4[[2]], v2], 
    D[LB4[[2]], w2], D[LB4[[2]], v3], 
    D[LB4[[2]], w3]}, {D[LB4[[3]], v1], D[LB4[[3]], w1], 
    D[LB4[[3]], v2], D[LB4[[3]], w2], D[LB4[[3]], v3], 
    D[LB4[[3]], w3]}, {D[LB4[[4]], v1], D[LB4[[4]], w1], 
    D[LB4[[4]], v2], D[LB4[[4]], w2], D[LB4[[4]], v3], 
    D[LB4[[4]], w3]}, {D[LB4[[5]], v1], D[LB4[[5]], w1], 
    D[LB4[[5]], v2], D[LB4[[5]], w2], D[LB4[[5]], v3], 
    D[LB4[[5]], w3]}, {D[LB4[[6]], v1], D[LB4[[6]], w1], 
    D[LB4[[6]], v2], D[LB4[[6]], w2], D[LB4[[6]], v3], 
    D[LB4[[6]], w3]}};

DLB4c = {{D[LB4c[[1]], v1], D[LB4c[[1]], w1], D[LB4c[[1]], v2], 
    D[LB4c[[1]], w2], D[LB4c[[1]], v3], D[LB4c[[1]], w3]},
   {D[LB4c[[2]], v1], D[LB4c[[2]], w1], D[LB4c[[2]], v2], 
    D[LB4c[[2]], w2], D[LB4c[[2]], v3], 
    D[LB4c[[2]], w3]}, {D[LB4c[[3]], v1], D[LB4c[[3]], w1], 
    D[LB4c[[3]], v2], D[LB4c[[3]], w2], D[LB4c[[3]], v3], 
    D[LB4c[[3]], w3]}, {D[LB4c[[4]], v1], D[LB4c[[4]], w1], 
    D[LB4c[[4]], v2], D[LB4c[[4]], w2], D[LB4c[[4]], v3], 
    D[LB4c[[4]], w3]}, {D[LB4c[[5]], v1], D[LB4c[[5]], w1], 
    D[LB4c[[5]], v2], D[LB4c[[5]], w2], D[LB4c[[5]], v3], 
    D[LB4c[[5]], w3]}, {D[LB4c[[6]], v1], D[LB4c[[6]], w1], 
    D[LB4c[[6]], v2], D[LB4c[[6]], w2], D[LB4c[[6]], v3], 
    D[LB4c[[6]], w3]}};

LB5 = DLB4.F - DFdx.LB4;
LB5c = DLB4c.Fc - DFcdx.LB4c;

(* Controllability Matrix *)
Ctrlb = {G, LB1, LB2, LB3, LB4, LB5};
Ctrlbc = {G, LB1c, LB2c, LB3c, LB4c, LB5c};
  
  (* Calculate the Observability at each point along the 12,000 point state trajectory *)
  For[i = 1, i < 12001, i++, 
    
    Ctrlb1ins = Ctrlbc /. {v1 -> X1[[1, 1, i]], w1 -> X1[[1, 25, i]], v2 -> X1[[1, 2, i]], w2 -> X1[[1, 26, i]], v3 -> X1[[1, 3, i]], w3 -> X1[[1, 27, i]]}; 
    sv1 = SingularValueList[Ctrlb1ins,6];
    If[i == 1, \[Delta]1 = {Abs[sv1[[-1]]]/Abs[sv1[[1]]]}, AppendTo[\[Delta]1, Abs[sv1[[-1]]]/Abs[sv1[[1]]]], Print[fail]];
    
    Ctrlb2ins = Ctrlbc /. {v1 -> X2[[1, 1, i]], w1 -> X2[[1, 25, i]], v2 -> X2[[1, 2, i]], w2 -> X2[[1, 26, i]], v3 -> X2[[1, 3, i]], w3 -> X2[[1, 27, i]]};
    sv2 = SingularValueList[Ctrlb2ins,6];
    If[i == 1, \[Delta]2 = {Abs[sv2[[-1]]]/Abs[sv2[[1]]]}, AppendTo[\[Delta]2, Abs[sv2[[-1]]]/Abs[sv2[[1]]]], Print[fail]];
    
    Ctrlb3ins = Ctrlb /. {v1 -> X3[[1, 1, i]], w1 -> X3[[1, 25, i]], v2 -> X3[[1, 2, i]], w2 -> X3[[1, 26, i]], v3 -> X3[[1, 3, i]], w3 -> X3[[1, 27, i]]};
    sv3 = SingularValueList[Ctrlb1ins,6];
    If[i == 1, \[Delta]3 = {Abs[sv3[[-1]]]/Abs[sv3[[1]]]}, AppendTo[\[Delta]3, Abs[sv3[[-1]]]/Abs[sv3[[1]]]], Print[fail]];

  ];
  
  (* File save *)

  If [IDcoupling == 1,
    Export[ToString["FNCtrlbDELTAIdN"]<>ToString[k]<>ToString["M1ChaosK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]1,"CSV"];
    Export[ToString["FNCtrlbDELTAIdN"]<>ToString[k]<>ToString["M1LCK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]2, "CSV"];
    Export[ToString["FNCtrlbDELTAIdN"]<>ToString[k]<>ToString["M1ConsK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]3,"CSV"];  
  ];

  If [IDcoupling == 0,
    Export[ToString["FNCtrlbDELTAHetN"]<>ToString[k]<>ToString["M1ChaosK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]1,"CSV"];
    Export[ToString["FNCtrlbDELTAHetN"]<>ToString[k]<>ToString["M1LCK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]2, "CSV"];
    Export[ToString["FNCtrlbDELTAHetN"]<>ToString[k]<>ToString["M1ConsK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]3,"CSV"];
  ];


  ];

];

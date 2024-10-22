(* ::Package:: *)

(* NonLinear Observability for a 3 Node FN Network with Sigmoidal \
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
  meas = k; (* measured node *)
  
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
  f2c = v1 - p2*w1+p1;
  f3c = q3*(v2 - v2^3/3 - w2 + d12*\[Kappa]/2*(Tanh[(v1 - h)/(2*w)] + 1) + d32*\[Kappa]/2*(Tanh[(v3 - h)/(2*w)] + 1));
  f4c = v2 - q2*w2 + q1;
  f5c = r3*(v3 - v3^3/3 - w3 + d13*\[Kappa]/2*(Tanh[(v1 - h)/(2*w)] + 1) + d23*\[Kappa]/2*(Tanh[(v2 - h)/(2*w)] + 1));
  f6c = v3 - r2*w3 +r1;

  (* Differential Coordinate Transform Map *)
  If[meas == 1, X = v1; Y = f1; Yc = f1c;];
  If[meas == 2, X = v2; Y = f3; Yc = f3c;];
  If[meas == 3, X = v3; Y = f5; Yc = f5c;];
  Z = D[Y, v1]*f1 + D[Y, w1]*f2 + D[Y, v2]*f3 + D[Y, w2]*f4 + D[Y, v3]*f5 + D[Y, w3]*f6;
  W = D[Z, v1]*f1 + D[Z, w1]*f2 + D[Z, v2]*f3 + D[Z, w2]*f4 + D[Z, v3]*f5 + D[Z, w3]*f6;
  U = D[W, v1]*f1 + D[W, w1]*f2 + D[W, v2]*f3 + D[W, w2]*f4 + D[W, v3]*f5 + D[W, w3]*f6;
  V = D[U, v1]*f1 + D[U, w1]*f2 + D[U, v2]*f3 + D[U, w2]*f4 + D[U, v3]*f5 + D[U, w3]*f6;

  Zc = D[Yc, v1]*f1c + D[Yc, w1]*f2c + D[Yc, v2]*f3c + D[Yc, w2]*f4c + D[Yc, v3]*f5c + D[Yc, w3]*f6c;
  Wc = D[Zc, v1]*f1c + D[Zc, w1]*f2c + D[Zc, v2]*f3c + D[Zc, w2]*f4c + D[Zc, v3]*f5c + D[Zc, w3]*f6c;
  Uc = D[Wc, v1]*f1c + D[Wc, w1]*f2c + D[Wc, v2]*f3c + D[Wc, w2]*f4c + D[Wc, v3]*f5c + D[Wc, w3]*f6c;
  Vc = D[Uc, v1]*f1c + D[Uc, w1]*f2c + D[Uc, v2]*f3c + D[Uc, w2]*f4c + D[Uc, v3]*f5c + D[Uc, w3]*f6c;
  
  (* Calculate the Jacobian of the map \[Phi] *)
  J\[Phi] = {{D[X, v1], D[X, w1], D[X, v2], D[X, w2], D[X, v3], D[X, w3]},
		{D[Y, v1], D[Y, w1], D[Y, v2], D[Y, w2], D[Y, v3], D[Y, w3]},
		{D[Z, v1], D[Z, w1], D[Z, v2], D[Z, w2], D[Z, v3], D[Z, w3]},
		{D[W, v1], D[W, w1], D[W, v2], D[W, w2], D[W, v3], D[W, w3]},
		{D[U, v1], D[U, w1], D[U, v2], D[U, w2], D[U, v3], D[U, w3]},
		{D[V, v1], D[V, w1], D[V, v2], D[V, w2], D[V, v3], D[V, w3]}};

  J\[Phi]c = {{D[X, v1], D[X, w1], D[X, v2], D[X, w2], D[X, v3], D[X, w3]},
		{D[Yc, v1], D[Yc, w1], D[Yc, v2], D[Yc, w2], D[Yc, v3], D[Yc, w3]},
		{D[Zc, v1], D[Zc, w1], D[Zc, v2], D[Zc, w2], D[Zc, v3], D[Zc, w3]},
		{D[Wc, v1], D[Wc, w1], D[Wc, v2], D[Wc, w2], D[Wc, v3], D[Wc, w3]},
		{D[Uc, v1], D[Uc, w1], D[Uc, v2], D[Uc, w2], D[Uc, v3], D[Uc, w3]},
		{D[Vc, v1], D[Vc, w1], D[Vc, v2], D[Vc, w2], D[Vc, v3], D[Vc, w3]}};
 
  (* Calculate the Observability at each point along the 12,000 point state trajectory *)
  For[i = 1, i < 12001, i++, 
   
   J\[Phi]1ins = J\[Phi]c /. {v1 -> X1[[1, 1, i]], w1 -> X1[[1, 25, i]], v2 -> X1[[1, 2, i]], w2 -> X1[[1, 26, i]], v3 -> X1[[1, 3, i]], w3 -> X1[[1, 27, i]]};
   sv1 = SingularValueList[J\[Phi]1ins,6];
   If[i == 1, \[Delta]1 = {Abs[sv1[[-1]]]/Abs[sv1[[1]]]}, AppendTo[\[Delta]1, Abs[sv1[[-1]]]/Abs[sv1[[1]]]], Print[fail]];

   J\[Phi]2ins = J\[Phi]c /. {v1 -> X2[[1, 1, i]], w1 -> X2[[1, 25, i]], v2 -> X2[[1, 2, i]], w2 -> X2[[1, 26, i]], v3 -> X2[[1, 3, i]], w3 -> X2[[1, 27, i]]};
   sv2 = SingularValueList[J\[Phi]2ins,6];
   If[i == 1, \[Delta]2 = {Abs[sv2[[-1]]]/Abs[sv2[[1]]]}, AppendTo[\[Delta]2, Abs[sv2[[-1]]]/Abs[sv2[[1]]]], Print[fail]];

   J\[Phi]3ins = J\[Phi] /. {v1 -> X3[[1, 1, i]], w1 -> X3[[1, 25, i]], v2 -> X3[[1, 2, i]], w2 -> X3[[1, 26, i]], v3 -> X3[[1, 3, i]], w3 -> X3[[1, 27, i]]};
   sv3 = SingularValueList[J\[Phi]3ins,6];
   If[i == 1, \[Delta]3 = {Abs[sv3[[-1]]]/Abs[sv3[[1]]]}, AppendTo[\[Delta]3, Abs[sv3[[-1]]]/Abs[sv3[[1]]]], Print[fail]];

  ];
  
  (* File save *)

  If [IDcoupling == 1,
    Export[ToString["FNObsDELTAIdN"]<>ToString[k]<>ToString["M1ChaosK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]1,"CSV"];
    Export[ToString["FNObsDELTAIdN"]<>ToString[k]<>ToString["M1LCK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]2, "CSV"];
    Export[ToString["FNObsDELTAIdN"]<>ToString[k]<>ToString["M1ConsK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]3,"CSV"];
  ];

  If [IDcoupling == 0,
    Export[ToString["FNObsDELTAHetN"]<>ToString[k]<>ToString["M1ChaosK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]1,"CSV"];
    Export[ToString["FNObsDELTAHetN"]<>ToString[k]<>ToString["M1LCK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]2, "CSV"];
    Export[ToString["FNObsDELTAHetN"]<>ToString[k]<>ToString["M1ConsK"]<>ToString[j]<>"IC"<>ToString[ICflag], \[Delta]3,"CSV"];
  ];

  ];

];

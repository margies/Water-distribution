set NODE;
set EDGE within {NODE, NODE};
set DiameterSet;

#---------------------------
# Junctions aren't reservoirs
param elev {NODE} >= 0; # Physical elevation of junction i, i.e., the height of junction i ([m])
param dem {NODE} >= 0; # Water demand at junction i ([m^3/s])

# Pipes
param length {EDGE} >= 0; # Length of pipe
# param v_max 1; # Upper bound on the velocity of water in pipe e ([m/s])
# k
# param roughness 100; # Physical constant depending on the roughness of pipe e
# C
param cost {DiameterSet};   # Cost of the r-th diameter that can be chosen for pipe e ([EUR/m])
# D
param diameter {DiameterSet} >= 0; # Diameter of pipe e

# Decision Variables of he Problem
# Q
var flow {EDGE}; # Flow in pipe e ([m3/s]), i.e., the volume of water which passes through pipe per unit time

# H
var hydraulicH {NODE} >= 0; # Fixed hydraulic head at junction i
# X
var select {EDGE, DiameterSet} binary;
# A
var A {EDGE};

#---------------------------
# Objection function
minimize Total_Cost:
	sum{(i, j) in EDGE, k in DiameterSet} select[i, j, k] * length[i, j] * cost[k];

#subject to Reservior:
#	hydraulicH[37] = 121;

subject to Flow {(i, j) in EDGE}:
	- (3.14 * ((sum{c in DiameterSet}select[i,j,c]*cost[c])/2)^2) <= flow[i, j] and flow[i, j] <= 3.14 * ((sum{c in DiameterSet}select[i,j,c]*cost[c])/2)^2;

subject to PH {i in NODE}:
	hydraulicH[i] - elev[i] >=40;

subject to PH2 {i in NODE}:
	hydraulicH[i] - elev[i] <= 121;

subject to HLoss {(i, j) in EDGE}:
	hydraulicH[i] - hydraulicH[j] = 
		(if flow[i, j] < 0 then -(-flow[i, j])^1.852 else flow[i, j]^1.852) * 10.7 * length[i, j] * (100^(-1.852))/(3.14 * ((sum{c in DiameterSet}select[i,j,c]*cost[c])/2)^2)^2.435;

subject to oneDia {(i, j) in EDGE}:
	sum{c in DiameterSet} select[i, j, c] = 1;

subject to c3 {i in NODE}:
	sum{(j, k) in EDGE: j = i} flow[j, k] - sum{(j, k) in EDGE: k = i} flow[j, k] = dem[i];

#---------------------------
# Data
data;

set NODE :=
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37;

set EDGE:=
	(1, 17) (17, 2) (2, 3) (3, 4) (4, 5) (5, 6) (6, 7) (7, 24) (24, 8) (8, 28) 
	(28, 9) (9, 36) (36, 1) (1, 31) (31, 10) (10, 11) (11, 19) (19, 12) (12, 4) 
	(2, 18) (18, 10) (10, 32) (32, 27) (27, 16) (16, 25) (25, 8) (3, 11) (11, 26) 
	(26, 15) (15, 22) (22, 7) (5, 13) (13, 14) (14, 20) (20, 15) (15, 16) (16, 29) 
	(29, 30) (30, 9) (17, 18) (12, 13) (19, 20) (14, 21) (21, 6) (21, 22) (22, 23)
	(24, 23) (23, 25) (26, 27) (28, 29) (29, 33) (32, 33) (33, 34) (31, 34) (34, 35)
	(30, 35) (35, 36) (37, 1);

set DiameterSet :=
r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13;

param: 
dem elev := 
	1	65.16	0.00049
	2	64.4	0.00104
	3	63.35	0.00102
	4	62.5	0.00081
	5	61.24	0.00063
	6	65.4	0.00079
	7	67.9	0.00026
	8	66.5	0.00058
	9	66	0.00054
	10	64.17	0.00111
	11	63.7	0.00175
	12	62.64	0.00091
	13	61.9	0.00116
	14	62.6	0.00054
	15	63.5	0.0011
	16	64.3	0.00121
	17	65.5	0.00127
	18	64.1	0.00202
	19	62.9	0.00188
	20	62.83	0.00093
	21	62.8	0.00096
	22	63.9	0.00097
	23	64.2	0.00086
	24	67.5	0.00067
	25	64.4	0.00077
	26	63.4	0.00169
	27	63.9	0.00142
	28	65.65	0.0003
	29	64.5	0.00062
	30	64.1	0.00054
	31	64.4	0.0009
	32	64.2	0.00103
	33	64.6	0.00077
	34	64.7	0.00074
	35	65.43	0.00116
	36	65.9	0.00047
	37	0	0;

param: diameter cost:=
    r1   0.06   19.8
    r2   0.08   24.5
    r3   0.1   27.2
    r4   0.125   37
    r5   0.15   39.4
    r6   0.2   54.4
    r7   0.25   72.9
    r8   0.3   90.7
    r9   0.35   119.5
    r10   0.4   139.1
    r11   0.45   164.4
    r12   0.5   186
    r13   0.6   241.3;

param: length:=
	1 17 132.76
	17 2 374.68
	2 3 119.74
	3 4 312.72
	4 5 289.09
	5 6 336.33
	6 7 135.81
	7 24 201.26
	24 8 132.53
	8 28 144.66
	28 9 175.72
	9 36 112.17
	36 1 210.74
	1 31 75.41
	31 10 181.42
	10 11 146.96
	11 19 162.69
	19 12 99.64
	12 4 52.98
	2 18 162.97
	18 10 83.96
	10 32 49.82
	32 27 78.5
	27 16 99.27
	16 25 82.29
	25 8 147.49
	3 11 197.32
	11 26 83.3
	26 15 113.8
	15 22 80.82
	22 7 340.97
	5 13 77.39
	13 14 112.37
	14 20 37.34
	20 15 108.85
	15 16 182.82
	16 29 136.01
	29 30 56.7
	30 9 124.08
	17 18 234.6
	12 13 203.83
	19 20 248.05
	14 21 65.19
	21 6 210.09
	21 22 147.57
	22 23 103.8
	24 23 210.95
	23 25 75.08
	26 27 180.29
	28 29 149.05
	29 33 215.05
	32 33 144.44
	33 34 34.74
	31 34 59.93
	34 35 165.67
	30 35 119.97
	35 36 83.17
	37 1 1;

option solver bonmin;
solve;
printf "result:/n";
display select;
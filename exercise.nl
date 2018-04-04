set NODE;
set EDGE within {NODE, NODE};
set DiameterSet;

#---------------------------
# Junctions aren't reservoirs
param elev {NODE}; # Physical elevation of junction i, i.e., the height of junction i ([m])
param dem {NODE} >= 0; # Water demand at junction i ([m^3/s])
# param ph_min {NODE} >= 0; # Lower bound on pressure head at junction i ([m])
param ph_max {NODE} >= 0; # Upper bound on pressure head at junction i ([m])

# Pipes
param length {EDGE} >= 0; # Length of pipe
param v_max {EDGE}; # Upper bound on the velocity of water in pipe e ([m/s])
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

subject to Reservior:
	hydraulicH[1] = 100;

subject to Flow {(i, j) in EDGE}:
	- (3.14 * v_max[i, j] * ((sum{c in DiameterSet}select[i,j,c]*cost[c])/2)^2) <= flow[i, j] and flow[i, j] <= 3.14 * v_max[i, j] * ((sum{c in DiameterSet}select[i,j,c]*cost[c])/2)^2;

subject to PH {i in NODE}:
	30 <= hydraulicH[i] - elev[i];

subject to PH2 {i in NODE}:
	hydraulicH[i] - elev[i] <= 100;

subject to HLoss {(i, j) in EDGE}:
	hydraulicH[i] - hydraulicH[j] = 
		(if flow[i, j] < 0 then -(-flow[i, j])^1.852 else flow[i, j]^1.852) * 10.7 * length[i, j] * (130^(-1.852))/(3.14 * ((sum{c in DiameterSet}select[i,j,c]*cost[c])/2)^2)^2.435;

subject to oneDia {(i, j) in EDGE}:
	sum{c in DiameterSet} select[i, j, c] = 1;

subject to c3 {i in NODE}:
	sum{(j, k) in EDGE: j = i} flow[j, k] - sum{(j, k) in EDGE: k = i} flow[j, k] = dem[i];

#---------------------------
# Data
data;

set NODE :=
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32;

set EDGE:=
	(1,2) (2,3) (3,4) (4,5) (5,6) (6,7) (7,8) (8,9) (9,10) (10,11) (11,12) (12,13) (10,14) (14,15) (15,16) (16,17) (17,18) (18,19) (19,3) (3,20) (20,21) (21,22) (20,23) (23,24) (24,25) (25,26) (26,27) (27,16) (23,28) (28,29) (29,30) (30,31) (31,32) (32,25);

set DiameterSet :=
r1 r2 r3 r4 r5 r6;

param: 
dem elev := 
	1	0	0
	2	0.24722	0
	3	0.23611	0
	4	0.03611	0
	5	0.20139	0
	6	0.27917	0
	7	0.375	0
	8	0.15278	0
	9	0.14583	0
	10	0.14583	0
	11	0.13889	0
	12	0.15556	0
	13	0.26111	0
	14	0.17083	0
	15	0.07778	0
	16	0.08611	0
	17	0.24028	0
	18	0.37361	0
	19	0.01667	0
	20	0.35417	0
	21	0.25833	0
	22	0.13472	0
	23	0.29028	0
	24	0.22778	0
	25	0.04722	0
	26	0.25	0
	27	0.10278	0
	28	0.08056	0
	29	0.1	0
	30	0.1	0
	31	0.02917	0
	32	0.22361	0;

param: diameter cost:=
    r1   0.3048   45.73
    r2   0.4064   70.40
    r3   0.5080   98.39
    r4   0.6096   129.33
    r5   0.7620   180.75
    r6   1.0160   278.28;

param: length v_max:=
	1	2	100	7
	2	3	1350	7
	3	4	900	3
	4	5	1150	3
	5	6	1450	2.5
	6	7	450	2.5
	7	8	850	2
	8	9	850	2
	9	10	800	2
	10	11	950	2
	11	12	1200	2
	12	13	3500	2
	10	14	800	2
	14	15	500	2
	15	16	550	2
	16	17	2730	2
	17	18	1750	2
	18	19	800	3.5
	19	3	400	3.5
	3	20	2200	3
	20	21	1500	2
	21	22	500	2
	20	23	2650	2
	23	24	1230	3
	24	25	1300	2
	25	26	850	2
	26	27	300	2
	27	16	750	2
	23	28	1500	2
	28	29	2000	2
	29	30	1600	2
	30	31	150	2
	31	32	860	2
	32	25	950	2;


option solver bonmin;
solve;
printf "result:/n";
display select;
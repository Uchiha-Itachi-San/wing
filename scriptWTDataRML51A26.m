% NACA RM-L51A26 Wind Tunnel Data
% Inboard Flap
% delta = 0deg
dataRML51A26(1).CL = [-5.273724	-3.766917	-1.849465	6.84E-02	2.054399	4.040422	6.026639	8.012467	9.930114	11.98432	14.03911	15.88818	18.01076	19.92802	21.98203	23.96727	26.02109	27.93776	30.12735	31.97448	33.20544; ...
                      -0.2164683	-0.1539399	-8.17E-02	2.83E-05	7.70E-02	0.1540533	0.2358277	0.3080782	0.3850624	0.4573413	0.5439059	0.6160998	0.688407	0.7558674	0.8233843	0.8813492	0.9441043	0.9972789	1.036281	1.060856	1.066128]';
                 
                 
dataRML51A26(1).CD = [2.25E-02	1.74E-02	1.42E-02	1.16E-02	1.02E-02	1.08E-02	0.0113903	1.39E-02	1.70E-02	0.0214168	2.71E-02	3.28E-02	4.04E-02	4.92E-02	6.00E-02	7.15E-02	8.42E-02	9.75E-02	0.1108706	0.1242266	0.1420447	0.1579581	0.1725966	0.18915	0.2063408	0.226079	0.251549	0.2706522	0.2929376	0.3152255	0.3432428	0.3642609	0.3846441	0.4075719	0.4349645	0.4540779	0.4821028	0.5082256	0.5286087	0.5489995	0.5706627	0.5929582; ...
                     -0.2189142	-0.1694827	-0.1200274	-7.53E-02	-2.58E-02	2.37E-02	7.79E-02	0.1321611	0.1840526	0.238317	0.2996678	0.353948	0.4011813	0.4507873	0.4933465	0.5359136	0.5737827	0.6045891	0.6353955	0.6638449	0.69235	0.7137605	0.7351552	0.7542167	0.7732862	0.7947443	0.8209876	0.840081	0.8639277	0.8854176	0.9140496	0.9308098	0.9452052	0.9643461	0.9811856	0.9908513	1.012413	1.024522	1.038918	1.046242	1.05594	1.070359]';

dataRML51A26(1).Cm = [-0.2164948	-0.1202749	-5.50E-02	0	0.1030928	0.1993127	0.3058419	0.4020618	0.5051546	0.604811	0.6872852	0.7594502	0.8419244	0.90378	0.9690722	1.037801; ...
                      3.05E-02	1.81E-02	1.05E-02	5.71E-03	-7.62E-03	-1.71E-02	-2.86E-02	-3.81E-02	-5.14E-02	-6.86E-02	-8.19E-02	-9.71E-02	-0.1133333	-0.1257143	-0.14	-0.1552381]';

dataRML51A26(1).IB(1).deltas_deg = -30:5:30;
dataRML51A26(1).IB(1).alpha_deg = 0;
dataRML51A26(1).IB(1).deltaCL = [-0.1142857	-0.1	-0.0809524	-0.0619048	-0.0428571	-0.0190476	0	0.0333333	0.052381	0.0761905	0.1047619	0.1333333	0.1571428];
dataRML51A26(1).IB(1).deltaCD = [0.03264	0.02736	0.02144128	0.01744	0.01472128	0.01264	0.0092787	0.00912	0.0115187	0.0139213	0.01696	0.02256	0.0281613];


dataRML51A26(1).IB(2).deltas_deg = -30:5:30;
dataRML51A26(1).IB(2).alpha_deg = 8;
dataRML51A26(1).IB(2).deltaCL = [-0.1190476	-0.1	-0.0904762	-0.07619046	-0.05238094	-0.0285714	0	0.0285714	0.052381	0.0666667	0.1047619	0.1333333	0.1476193];

dataRML51A26(1).IB(3).deltas_deg = [0.00	3.87	7.82	13.48	15.91	19.87	23.73	29.93];
dataRML51A26(1).IB(3).alpha_deg = 0;
dataRML51A26(1).IB(3).deltaCL = [0	1.59E-02	3.17E-02	5.08E-02	6.35E-02	8.89E-02	0.1079365	0.1333333];
dataRML51A26(1).IB(3).deltaCm = [4.44E-03	-1.71E-03	-9.15E-03	-0.0202812	-2.46E-02	-3.20E-02	-3.76E-02	-5.06E-02];

dataRML51A26(1).IB(4).deltas_deg = [0	3.89836	7.715153	13.39261	15.92085	19.69014	23.60285	30.04221];
dataRML51A26(1).IB(4).Cl = [8.9294E-04	3.0849E-03	4.5162E-03	6.6644E-03	7.9374E-03	9.1789E-03	1.0226E-02	1.3883E-02] - 8.9294E-04;

dataRML51A26(1).IB(5).deltas_deg = [0	3.84	7.71	12.39	15.85	19.77	23.65	29.98];
dataRML51A26(1).IB(5).Cn = [-0.000588	-0.000852	-0.001264	-0.001554	-0.001659	-0.001779	-0.002339	-0.003267] + 0.000588;

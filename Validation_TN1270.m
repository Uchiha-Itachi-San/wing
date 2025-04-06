%% Experimental setup
nu = 1.46 * 10^-5;
rho = 1.225;
ReNumber    = 4.32e6;
MaNumber    = 0.13;
V_A = MaNumber * 340;
TaperRatio  = 2.5; % NACA convention
AspectRatio = 8;
alphas_deg = -5:2:21; % [°]

%% Calculate wing geometry

c_mean = ReNumber*nu / V_A;
c_tip  = 3/2*c_mean*(TaperRatio+1)/(TaperRatio^2+TaperRatio+1);
c_root = TaperRatio*c_tip;
span = AspectRatio*c_tip*(TaperRatio+1)/2;
nPanelsY = 20;

%% Calculate user-specified reference values
GEO.Ref.q_Ref = V_A^2*rho/2;
GEO.Ref.S_Ref = span*c_tip*(1+TaperRatio)/2;
GEO.Ref.r_RefRef_Ref = [0.015; 0; 0];


%% Calculate VLM
nPanelsX = 10;
airfoil1 = Airfoil('NACA', '4416');
airfoil2 = Airfoil('NACA', '4412');
wing_VLM = TaperedWing(GEO.Ref.S_Ref, [airfoil1, airfoil2], AspectRatio, 1/TaperRatio);
% wing_VLM.setRefPointMom([0; 0; 0]);
wing_VLM.addAirfoilZone(0, 1, 1, 2, 'linear');
wing_VLM.addTwistTransition(0, 1, -4.5*pi/180, 'linear');

results = wing_VLM.calculateVLM(alphas_deg /180*pi, 0, nPanelsX, nPanelsY);
VLM.CL = [results().c_L];
VLM.CD = [results().c_D] + wing_VLM.estimateZeroLiftDrag(V_A / nu);
VLM.Cm = [results().c_m];

%% Plotting results
load ValidationData_TN1270;

figureFlap = 'TN-1270 - Aerodynamics over alpha';
handleCoeff = findobj('type', 'figure', 'Name', figureFlap);
if isempty(handleCoeff)
    handleCoeff = figure('Name', figureFlap);
end
figure(handleCoeff);
clf

ColorIt = {'red', 'green', 'blue', 'magenta'};

subplot(221);

plot(alphas_deg, VLM.CL, 'LineWidth', 2, 'Color', ColorIt{2});
hold on;
plot(ValidationData1270.CL(:,1), ValidationData1270.CL(:,2), 'LineWidth', 2, 'LineStyle', ':', 'Color', ColorIt{3});
xlabel('\alpha [°]');
ylabel('C_L [-]');
legend('VLM', 'Exp. NACA-TN-1270', 'Location', 'SouthEast');
grid on;

subplot(222);
plot(alphas_deg, VLM.CD, 'LineWidth', 2, 'Color', ColorIt{2});
hold on;
xlabel('\alpha [°]');
ylabel('C_D [-]');
grid on;

subplot(223);
plot(VLM.CD, VLM.CL, 'LineWidth', 2, 'Color', ColorIt{2});
hold on;
plot(ValidationData1270.CD(:,1), ValidationData1270.CD(:,2), 'LineWidth', 2, 'LineStyle', ':', 'Color', ColorIt{3});
xlabel('C_D [-]');
ylabel('C_L [-]');
grid on;

subplot(224);
plot(VLM.Cm, VLM.CL, 'LineWidth', 2, 'Color', ColorIt{2});
hold on;
plot(ValidationData1270.Cm(:,1), ValidationData1270.Cm(:,2), 'LineWidth', 2, 'LineStyle', ':', 'Color', ColorIt{3});
% plot(C_L_PAWAT, C_m_PAWAT, 'LineWidth', 2, 'LineStyle', '--', 'Color', ColorIt{1});
% plot(VOLANOS.CL(:,3), VoLaNoSData1270.Cm(:,3), 'LineWidth', 2, 'LineStyle', '--', 'Color', ColorIt{2});
xlabel('C_m [-]');
ylabel('C_L [-]');
%legend('Fixed-Wake Lifting Line', 'VLM Tornado', 'NACA-TN-1270', 'Fixed-Wake Lifting Line mod. Ref', 'VLM Tornado mod. Ref');
grid on;


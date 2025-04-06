% Validation NACA TN-2080
CaseID = 1;
alphas_deg = 0;%[4:4:8,10:2:24];
deltas_deg = 0:5:60;

deltas = deltas_deg / 180*pi;
% alphas_deg = unique(sort([alphas_deg, 4]));
alphas = alphas_deg /180*pi;
i_Alpha_spanLoading = find(alphas_deg == 4);

nPanelsY = 20;

%% Geometry
S_ref = 1.8;
AR = 3.1;
TR = 1;

%% Create Flight State object

V_A           = 92;
nu = 1.46*10^-5;

%% Load airfoil database
% load('N:\TUM-PC\Dokumente\MATLAB\Finite Wing\Input\NACA64A012\Airfoil_NACA64A012_v4.mat');
% airfoil.importCoordinates('NACA64A010.dat');
% airfoilDB = airfoil;
% airfoilDB.FlapMethod = 'correctionDATCOM';
% airfoilDB.extendAlphaRange();

%% Calculation VOLANOS

nPanelsX = 10;
clear('wing_VLM_TN2080');
clear('VLM_TN2080');

for i = 1:length(deltas)
    wing_VLM_TN2080(i) = TaperedWing(S_ref, [], AR, TR);
    wing_VLM_TN2080(i).setRefPointMom([0; 0; 0]);
    % Inboard flap
    wing_VLM_TN2080(i).addFlaps([0, 0.484], [0.625, 0.625]*0.3048, 2, 1);
%     wing_VoLaNoS_TN2080(i).addFlaps([0, 1], [0.625, 0.625]*0.3048, 2, 1);
    % Outboard flap
%     wing_VoLaNoS(i).addFlaps([0.484, 1], [0.625, 0.625]*0.3048, 2, 1);
    wing_VLM_TN2080(i).setFlapDeflection(2, deltas(i));
    
    fprintf('Calculating flap deflection angle %.1f°\n', deltas_deg(i));
    results = wing_VLM_TN2080(i).calculateVLM(alphas, 0, nPanelsX, nPanelsY);
    VLM_TN2080.CL(i) = [results.c_L];
    VLM_TN2080.CD(i) = [results.c_D] + wing_VLM_TN2080(i).estimateZeroLiftDrag(V_A / nu);
    VLM_TN2080.Cm(i) = [results.c_m];
%     VLM_TN2080.Cl(i) = -[results.c_l];
%     VLM_TN2080.Cn(i) = -[results.c_n];

%     panelData_VOLANOS_TN2080(i) = wing_VLM_TN2080(i).aeroProp.panelData;
end

wing_VLM_TN2080(i).plotVLMLattice(nPanelsX, nPanelsY);

%% Plot results

run scriptWTDataTN2080;

figureFlap = 'TN2080: Aerodynamics over delta';
handleCoeff = findobj('type', 'figure', 'Name', figureFlap);
if isempty(handleCoeff)
    handleCoeff = figure('Name', figureFlap);
end
figure(handleCoeff);
clf

Colors = {'red', 'green', 'blue', 'magenta'};

subplot(221);
plot(deltas_deg, VLM_TN2080.CL, 'LineWidth', 2, 'Color', Colors{2});
hold on;
plot(dataTN2080(CaseID).IB(1).deltas_deg, dataTN2080(CaseID).IB(1).deltaCL, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('\delta [°]');
ylabel('C_L [-]');
legend('VLM', 'Exp. NACA-TN-2080', 'Location', 'SouthEast');
grid on;

subplot(222);
plot(deltas_deg, VLM_TN2080.CD, 'LineWidth', 2, 'Color', Colors{2});
hold on
plot(dataTN2080(CaseID).IB(1).deltas_deg, dataTN2080(CaseID).IB(1).CD, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('\delta [°]');
ylabel('C_D [-]');
% legend('VLM', 'Exp. NACA-TN-2080', 'Location', 'SouthEast');
grid on;

subplot(223);
plot(VLM_TN2080.CD, VLM_TN2080.CL, 'LineWidth', 2, 'Color', Colors{2});
hold on;
plot(dataTN2080(CaseID).IB(1).CD, dataTN2080(CaseID).IB(1).deltaCL, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('C_D [-]');
ylabel('C_L [-]');
grid on;

subplot(224);
plot(VLM_TN2080.CL, VLM_TN2080.Cm, 'LineWidth', 2, 'Color', Colors{2});
hold on;
plot(dataTN2080(CaseID).IB(1).deltaCL, dataTN2080(CaseID).IB(1).Cm, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('C_L [-]');
ylabel('C_m [-]');
grid on;

% %% Plot span distribution
% i_wing = length(deltas);
% i_alpha = 1; % length(alphas_deg);
% 
% figureFlap = ['TN2080: Span distributions at alpha = ', num2str(alphas_deg(i_alpha)), ', delta = ', num2str(deltas_deg(i_wing))];
% handleCoeff = findobj('type', 'figure', 'Name', figureFlap);
% if isempty(handleCoeff)
%     handleCoeff = figure('Name', figureFlap);
% end
% figure(handleCoeff);
% clf;
% 
% yLabels = {'CL [-]', 'CD [-]', 'Cm [-]', 'Re [-]', '\alpha [°]', 'CL * c [m]'};
% % legends = {{'CL', 'CL_{max}'}, 'CD', 'Cm', 'Re', {'alpha_{local}', 'alpha_{CLmax}'}, {'CL * c', 'CL_{max} * c'}};
% legends = {{'CL', 'CL_{visc}', '\DeltaCL_{corr}'}, {'CD', 'CD_{visc}', 'CD_i'}, {'Cm', 'Cm_{visc}', '\DeltaCm_{corr}'}, 'Re', {'alpha_{local}', '\delta_{1}', '\delta_{2}'}, {'CL * c', 'CL_{visc} * c'}};
% data1 = {'%CL', '%CD', '%Cm', '%Re', '(%alpha_local) /pi*180', '%CL .* %c'};
% data2 = {'%CL_visc', '%CD_visc', '%Cm_visc', [], '%delta1 /pi*180', '%CL_visc .* %c'};
% data3 = {'%deltaCl_corr', '%CDi', '%deltaCm_corr', [], '%delta2 /pi*180', []};
% 
% for i = 1:length(yLabels)
%     subplot(2,3,i);
%     hold off
%     try
%         evalc(['data_temp3 = ', strrep(data1{i}, '%', ['wing_VoLaNoS_TN2080(i_wing).aeroProp.panelData(', num2str(i_alpha), ').'])]);
%         plot(wing_VLM_TN2080(i_wing).aeroProp.panelData(i_alpha).y, data_temp3, 'Color', 'g', 'Marker', '+', 'LineWidth', 2);
%         hold on
%         if ~isempty(data2{i})
%             evalc(['data_temp3 = ', strrep(data2{i}, '%', ['wing_VoLaNoS_TN2080(i_wing).aeroProp.panelData(', num2str(i_alpha), ').'])]);
%             plot(wing_VLM_TN2080(i_wing).aeroProp.panelData(i_alpha).y, data_temp3, 'Color', 'b', 'Marker', 'none', 'LineWidth', 2);
%         end
%         if ~isempty(data3{i})
%             evalc(['data_temp3 = ', strrep(data3{i}, '%', ['wing_VoLaNoS_TN2080(i_wing).aeroProp.panelData(', num2str(i_alpha), ').'])]);
%             plot(wing_VLM_TN2080(i_wing).aeroProp.panelData(i_alpha).y, data_temp3, 'Color', 'r', 'Marker', 'none', 'LineWidth', 2);
%         end
%     end
%     xlabel('y [m]');
%     ylabel(yLabels{i});
%     legend(legends{i}, 'Location', 'southeast');
%     grid on
% end
% 
% %% Plot local parameters
% y_target = -0.5 * (wing_VLM_TN2080(1).flapBreaks.relYPos1 + wing_VLM_TN2080(1).flapBreaks.relYPos2) * 0.5 * wing_VLM_TN2080(1).getSpan();
% i_y = find(wing_VLM_TN2080(1).aeroProp.panelData(1).y > y_target, 1, 'first');
% 
% figureFlap = 'TN2080: Plot of local parameters';
% handleCoeff = findobj('type', 'figure', 'Name', figureFlap);
% if isempty(handleCoeff)
%     handleCoeff = figure('Name', figureFlap);
% end
% figure(handleCoeff);
% clf;
% 
% CL_local = zeros(length(deltas_deg), 1);
% CL_local_visc = zeros(length(deltas_deg), 1);
% Cm_local = zeros(length(deltas_deg), 1);
% Cm_local_visc = zeros(length(deltas_deg), 1);
% delta1 = zeros(length(deltas_deg), 1);
% delta2 = zeros(length(deltas_deg), 1);
% alpha_loc = zeros(length(deltas_deg), 1);
% 
% for j = 1:length(deltas_deg)
%     
%     CL_local(j) = wing_VLM_TN2080(j).aeroProp.panelData(1).CL(i_y);
%     CL_local_visc(j) = wing_VLM_TN2080(j).aeroProp.panelData(1).CL_visc(i_y);
%     Cm_local(j) = wing_VLM_TN2080(j).aeroProp.panelData(1).Cm(i_y);
%     Cm_local_visc(j) = wing_VLM_TN2080(j).aeroProp.panelData(1).Cm_visc(i_y);
%     delta1(j) = wing_VLM_TN2080(j).aeroProp.panelData(1).delta1(i_y) / pi*180;
%     delta2(j) = wing_VLM_TN2080(j).aeroProp.panelData(1).delta2(i_y) / pi*180;
%     alpha_loc(j) = wing_VLM_TN2080(j).aeroProp.panelData(1).alpha_local(i_y) / pi*180;
% end
% 
% subplot(1,3,1)
% plot(deltas_deg, CL_local, 'LineWidth', 2, 'Color', 'g');
% hold on;
% plot(deltas_deg, CL_local_visc, 'LineWidth', 2, 'Color', 'b');
% xlabel('\delta [°]');
% ylabel('C_L [-]');
% legend('CL_{local}', 'CL_{local,visc}', 'Location', 'SouthEast');
% grid on;
% 
% subplot(1,3,2)
% plot(deltas_deg, Cm_local, 'LineWidth', 2, 'Color', 'g');
% hold on;
% plot(deltas_deg, Cm_local_visc, 'LineWidth', 2, 'Color', 'b');
% xlabel('\delta [°]');
% ylabel('C_m [-]');
% legend('Cm_{local}', 'Cm_{local,visc}', 'Location', 'SouthEast');
% grid on;
% 
% subplot(1,3,3)
% plot(deltas_deg, delta1, 'LineWidth', 2, 'Color', 'b');
% hold on;
% plot(deltas_deg, delta2, 'LineWidth', 2, 'Color', 'r');
% plot(deltas_deg, alpha_loc, 'LineWidth', 2, 'Color', 'g');
% xlabel('\delta [°]');
% ylabel('angle [°]');
% legend('\delta_1', '\delta_2', '\alpha_{loc}', 'Location', 'SouthEast');
% grid on;
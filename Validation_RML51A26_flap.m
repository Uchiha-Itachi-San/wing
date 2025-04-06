
CaseID = 1;
alphas_deg = 0;%[4:4:8,10:2:24];
deltas_deg = 0:0.5:30;

deltas = deltas_deg / 180*pi;
% alphas_deg = unique(sort([alphas_deg, 4]));
alphas = alphas_deg /180*pi;
betas = 0;
i_Alpha_spanLoading = find(alphas_deg == 4);

nPanelsY = 15;

%% Test case setup
phi_0  = 60*pi/180; % Sweep angle of leading edge
TR     = 0.001;       % Taper ratio
S_ref = 21.46; % [m²]
ReNumber = 6.0e6;


%% Calculate wing geometry
span = sqrt(4 * S_ref / tan(phi_0));
c_root = tan(phi_0) * 0.5 * span;
c_tip = 0;
phi_25 = atan(0.75 * c_root / (0.5 * span));
% c_mean = ReNumber*State.KinViscosity/State.V_A;
% lambda = 1/TR(CaseID);
% c_tip  = 3/2*c_mean*(lambda+1)/(lambda^2+lambda+1);
% c_root = lambda*c_tip;
% span   = AR(CaseID)*c_tip*(lambda+1)/2;
twist_tip = 0;
% S_ref = 0.5 * (c_root + c_tip) * span;
c_avg = S_ref / span;
AR    = span^2 / S_ref;        % Aspect ratio
mac = 2 / 3 * c_root * (1 + TR + TR^2) / (1 + TR);
y_mac = 1/3 * ((1 + 2 * TR) / (1 + TR)) * 0.5 * span;
deltaX_mac = tan(phi_0) * y_mac;

%% Calculation VLM

nPanelsX = 8;
clear('wing_RML51A26');
clear('VLM_RML51A26');

for i = 1:length(deltas)
    wing_RML51A26(i) = TaperedWing(S_ref, [], AR, TR);
    wing_RML51A26(i).setRefPointMom([deltaX_mac + 0.25 * mac - 0.25 * c_root; 0; 0]);
    wing_RML51A26(i).addSweepKink(0, phi_25);
    wing_RML51A26(i).addFlaps([0.0288, 0.5], [0.762, 0.762], 2, 0);
%     wing_VoLaNoS_RML51A26(i).addFlaps([0.5, 1], [0.762, 0.762], 2, 0);
    wing_RML51A26(i).setFlapDeflection(2, deltas(i));
    % wing_VoLaNoS.addAirfoilZone(0, 1, 1, 2, 'linear');
    % wing_VoLaNoS.addTwistTransition(0, 1, twist_tip, 'linear', 0);
    
    fprintf('Calculating flap deflection angle %.1f°\n', deltas_deg(i));
    results = wing_RML51A26(i).calculateVLM(alphas, betas, nPanelsX, nPanelsY);
    VLM_RML51A26.CL(i) = [results.c_L];
    VLM_RML51A26.CD(i) = [results.c_D] + wing_RML51A26(i).estimateZeroLiftDrag(ReNumber / mac);
    VLM_RML51A26.Cm(i) = [results.c_m];
    VLM_RML51A26.Cl(i) = [results.c_l];
    VLM_RML51A26.Cn(i) = [results.c_n];

%     panelData_VOLANOS_RML51A26(i) = wing_RML51A26(i).aeroProp.panelData;
end
wing_RML51A26(i).plotVLMLattice(nPanelsX, nPanelsY);

%% Plot results

run scriptWTDataRML51A26;

figureFlap = 'RML51A26: Aerodynamics over delta';
handleCoeff = findobj('type', 'figure', 'Name', figureFlap);
if isempty(handleCoeff)
    handleCoeff = figure('Name', figureFlap);
end
figure(handleCoeff);
clf

Colors = {'red', 'green', 'blue', 'magenta'};

subplot(221);
plot(deltas_deg, VLM_RML51A26.CL, 'LineWidth', 2, 'Color', Colors{2});
hold on;
plot(dataRML51A26(CaseID).IB(1).deltas_deg, dataRML51A26(CaseID).IB(1).deltaCL, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('\delta [°]');
ylabel('C_L [-]');
% legend('PAWAT', 'VoLaNoS', 'Exp. NACA-TN-1270', 'Location', 'SouthEast');
grid on;

subplot(222);
plot(deltas_deg, VLM_RML51A26.Cm, 'LineWidth', 2, 'Color', Colors{2});
hold on
plot(dataRML51A26(CaseID).IB(3).deltas_deg, dataRML51A26(CaseID).IB(3).deltaCm - dataRML51A26(CaseID).IB(3).deltaCm(1), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('\delta [°]');
ylabel('C_m [-]');
grid on;

subplot(223);
plot(VLM_RML51A26.CD, VLM_RML51A26.CL, 'LineWidth', 2, 'Color', Colors{2});
hold on;
plot(dataRML51A26(CaseID).IB(1).deltaCD(dataRML51A26(CaseID).IB(1).deltas_deg >= 0), dataRML51A26(CaseID).IB(1).deltaCL(dataRML51A26(CaseID).IB(1).deltas_deg >= 0), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('C_D [-]');
ylabel('C_L [-]');
grid on;

subplot(224);
plot(VLM_RML51A26.Cm, VLM_RML51A26.CL, 'LineWidth', 2, 'Color', Colors{2});
hold on;
plot(dataRML51A26(CaseID).IB(3).deltaCm, dataRML51A26(CaseID).IB(3).deltaCL, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('C_m [-]');
ylabel('C_L [-]');
grid on;

% %% Plot span distribution
% i_wing = 2; length(deltas);
% i_alpha = 1; %length(alphas_deg);
% 
% figureFlap = ['RML51A26: Span distributions at alpha = ', num2str(alphas_deg(i_alpha)), ', delta = ', num2str(deltas_deg(i_wing))];
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
%         evalc(['data_temp3 = ', strrep(data1{i}, '%', ['wing_VoLaNoS_RML51A26(i_wing).aeroProp.panelData(', num2str(i_alpha), ').'])]);
%         plot(wing_RML51A26(i_wing).aeroProp.panelData(i_alpha).y, data_temp3, 'Color', 'g', 'Marker', '+', 'LineWidth', 2);
%         hold on
%         if ~isempty(data2{i})
%             evalc(['data_temp3 = ', strrep(data2{i}, '%', ['wing_VoLaNoS_RML51A26(i_wing).aeroProp.panelData(', num2str(i_alpha), ').'])]);
%             plot(wing_RML51A26(i_wing).aeroProp.panelData(i_alpha).y, data_temp3, 'Color', 'b', 'Marker', 'none', 'LineWidth', 2);
%         end
%         if ~isempty(data3{i})
%             evalc(['data_temp3 = ', strrep(data3{i}, '%', ['wing_VoLaNoS_RML51A26(i_wing).aeroProp.panelData(', num2str(i_alpha), ').'])]);
%             plot(wing_RML51A26(i_wing).aeroProp.panelData(i_alpha).y, data_temp3, 'Color', 'r', 'Marker', 'none', 'LineWidth', 2);
%         end
%     end
%     xlabel('y [m]');
%     ylabel(yLabels{i});
%     legend(legends{i}, 'Location', 'southeast');
%     grid on
% end
%     
% %% Plot local parameters
% y_target = -0.5 * (0.0288 + 0.5) * 0.5 * wing_RML51A26(1).getSpan();
% i_y = find(wing_RML51A26(1).aeroProp.panelData(1).y > y_target, 1, 'first');
% 
% figureFlap = 'RML51A26: Plot of local parameters';
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
%     CL_local(j) = wing_RML51A26(j).aeroProp.panelData(1).CL(i_y);
%     CL_local_visc(j) = wing_RML51A26(j).aeroProp.panelData(1).CL_visc(i_y);
%     Cm_local(j) = wing_RML51A26(j).aeroProp.panelData(1).Cm(i_y);
%     Cm_local_visc(j) = wing_RML51A26(j).aeroProp.panelData(1).Cm_visc(i_y);
%     delta1(j) = wing_RML51A26(j).aeroProp.panelData(1).delta1(i_y) / pi*180;
%     delta2(j) = wing_RML51A26(j).aeroProp.panelData(1).delta2(i_y) / pi*180;
%     alpha_loc(j) = wing_RML51A26(j).aeroProp.panelData(1).alpha_local(i_y) / pi*180;
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

%% Control moments
figureFlap = 'RML51A26: Moments over delta';
handleCoeff = findobj('type', 'figure', 'Name', figureFlap);
if isempty(handleCoeff)
    handleCoeff = figure('Name', figureFlap);
end
figure(handleCoeff);
clf

Colors = {'red', 'green', 'blue', 'magenta'};

subplot(121);
plot(deltas_deg, VLM_RML51A26.Cl, 'LineWidth', 2, 'Color', Colors{2});
hold on;
plot(dataRML51A26(CaseID).IB(4).deltas_deg, dataRML51A26(CaseID).IB(4).Cl, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('\delta [°]');
ylabel('C_l [-]');
% legend('PAWAT', 'VoLaNoS', 'Exp. NACA-TN-1270', 'Location', 'SouthEast');
grid on;

subplot(122);
plot(deltas_deg, VLM_RML51A26.Cn, 'LineWidth', 2, 'Color', Colors{2});
hold on;
plot(dataRML51A26(CaseID).IB(5).deltas_deg, dataRML51A26(CaseID).IB(5).Cn, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('\delta [°]');
ylabel('C_n [-]');
% legend('PAWAT', 'VoLaNoS', 'Exp. NACA-TN-1270', 'Location', 'SouthEast');
grid on;

CaseID = 1;
alphas_deg = 0:2:24;%[4:4:8,10:2:24];

% alphas_deg = unique(sort([alphas_deg, 4]));
alphas = alphas_deg /180*pi;
i_Alpha_spanLoading = find(alphas_deg == 4);
betas = 0;

nPanelsY = 20;

%% Test case setup
phi_0  = 60*pi/180; % Sweep angle of leading edge
TR     = 0;       % Taper ratio
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

nPanelsX = 10;

wing_RML51A26 = TaperedWing(S_ref, [], AR(CaseID), TR(CaseID));
wing_RML51A26.setRefPointMom([deltaX_mac + 0.25 * mac - 0.25 * c_root; 0; 0]);
% wing_VoLaNoS.addAirfoilZone(0, 1, 1, 2, 'linear');
% wing_VoLaNoS.addTwistTransition(0, 1, twist_tip, 'linear', 0);
wing_RML51A26.addSweepKink(0, phi_25(CaseID));

results = wing_RML51A26.calculateVLM(alphas, betas, nPanelsX, nPanelsY);
VLM.CL = [results().c_L];
VLM.CD = [results().c_D] + wing_RML51A26.estimateZeroLiftDrag(ReNumber / mac);
VLM.Cm = [results().c_m];

% panelData_VOLANOS = wing_RML51A26.aeroProp.panelData(i_Alpha_spanLoading);

%% Plot results

% figureFlap = 'Span Loading';
% handleCoeff = findobj('type', 'figure', 'Name', figureFlap);
% if isempty(handleCoeff)
%     handleCoeff = figure('Name', figureFlap);
% end
% figure(handleCoeff);
% clf
% 
% ind_VOLANOS = ceil(0.5 * length(panelData_VOLANOS.y)):length(panelData_VOLANOS.y);
% 
% subplot(131)
% 
% subplot(132)
% plot(panelData_PAWAT(:,1)/(span/2),panelData_PAWAT(:,3)/C_L_PAWAT(i_Alpha_spanLoading),'b');
% hold on
% plot(panelData_VOLANOS.y(ind_VOLANOS)/(span/2), panelData_VOLANOS.CL(ind_VOLANOS) / VOLANOS.CL(i_Alpha_spanLoading),'r');
% % plot(ExpData{CaseID}.LiftCoeff.eta,ExpData{CaseID}.LiftCoeff.Val,'ko');
% xlabel('\eta [-]')
% ylabel('C_l/C_L')
% set(gca,'XLim',[0,1.1]);
% grid on
% legend('PAWAT', 'VoLaNoS', 'Location', 'South');
% subplot(133)
% plot(panelData_PAWAT(:,1) / (span/2), panelData_PAWAT(:,3) .* panelData_PAWAT(:,2) / C_L_PAWAT(i_Alpha_spanLoading) / c_avg,'b');
% hold on
% plot(panelData_VOLANOS.y(ind_VOLANOS) / (span/2),panelData_VOLANOS.CL(ind_VOLANOS) .* panelData_VOLANOS.c(ind_VOLANOS) / VOLANOS.CL(i_Alpha_spanLoading) / c_avg,'r');
% % plot(ExpData{CaseID}.LoadCoeff.eta,ExpData{CaseID}.LoadCoeff.Val,'ko');
% xlabel('\eta [-]')
% ylabel('C_lc/C_Lc_a_v')
% set(gca,'XLim',[0,1.1]);
% grid on
% legend('PAWAT',  'VoLaNoS', 'Location', 'South');

%%

run scriptWTDataRML51A26;

figureFlap = 'Aerodynamics over alpha';
handleCoeff = findobj('type', 'figure', 'Name', figureFlap);
if isempty(handleCoeff)
    handleCoeff = figure('Name', figureFlap);
end
figure(handleCoeff);
clf

Colors = {'red', 'green', 'blue', 'magenta'};

subplot(221);
plot(alphas_deg, VLM.CL, 'LineWidth', 2, 'Color', Colors{2});
hold on;
plot(dataRML51A26(CaseID).CL(:,1), dataRML51A26(CaseID).CL(:,2), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('\alpha [°]');
ylabel('C_L [-]');
% legend('PAWAT', 'VoLaNoS', 'Exp. NACA-TN-1270', 'Location', 'SouthEast');
grid on;

subplot(222);
plot(alphas_deg, VLM.CD, 'LineWidth', 2, 'Color', Colors{2});
hold on;

xlabel('\alpha [°]');
ylabel('C_D [-]');
grid on;

subplot(223);
plot(VLM.CD, VLM.CL, 'LineWidth', 2, 'Color', Colors{2});
hold on;
plot(dataRML51A26(CaseID).CD(:,1), dataRML51A26(CaseID).CD(:,2), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('C_D [-]');
ylabel('C_L [-]');
grid on;

subplot(224);
plot(VLM.Cm, VLM.CL, 'LineWidth', 2, 'Color', Colors{2});
hold on;
plot(dataRML51A26(CaseID).Cm(:,2), dataRML51A26(CaseID).Cm(:,1), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
xlabel('C_m [-]');
ylabel('C_L [-]');
grid on;

%% Plot span distribution

% i_alpha = length(alphas_deg)-1;
% 
% figureFlap = ['Span distributions at alpha = ', num2str(alphas_deg(i_alpha))];
% handleCoeff = findobj('type', 'figure', 'Name', figureFlap);
% if isempty(handleCoeff)
%     handleCoeff = figure('Name', figureFlap);
% end
% figure(handleCoeff);
% clf;
% 
% yLabels = {'CL [-]', 'CD [-]', 'Cm [-]', 'Re [-]', '\alpha [°]', 'CL * c [m]'};
% % legends = {{'CL', 'CL_{max}'}, 'CD', 'Cm', 'Re', {'alpha_{local}', 'alpha_{CLmax}'}, {'CL * c', 'CL_{max} * c'}};
% legends = {{'CL', 'CL_{visc}'}, {'CD', 'CD_{visc}', 'CD_i'}, 'Cm', 'Re', {'alpha_{local}', '\delta_{1}', '\delta_{2}'}, {'CL * c', 'CL_{visc} * c'}};
% data1 = {'%CL', '%CD', '%Cm', '%Re', '(%alpha_local) /pi*180', '%CL .* %c'};
% data2 = {'%CL_visc', '%CD_visc', '%Cm_visc', [], '%delta1 /pi*180', '%CL_visc .* %c'};
% data3 = {[], '%CDi', [], [], '%delta2 /pi*180', []};
% 
% for i = 1:length(yLabels)
%     subplot(2,3,i);
%     hold off
%     try
%         evalc(['data_temp3 = ', strrep(data1{i}, '%', ['wing_VoLaNoS.aeroProp.panelData(', num2str(i_alpha), ').'])]);
%         plot(wing_RML51A26.aeroProp.panelData(i_alpha).y, data_temp3, 'Color', 'g', 'Marker', '+', 'LineWidth', 2);
%         hold on
%         if ~isempty(data2{i})
%             evalc(['data_temp3 = ', strrep(data2{i}, '%', ['wing_VoLaNoS.aeroProp.panelData(', num2str(i_alpha), ').'])]);
%             plot(wing_RML51A26.aeroProp.panelData(i_alpha).y, data_temp3, 'Color', 'b', 'Marker', 'none', 'LineWidth', 2);
%         end
%         if ~isempty(data3{i})
%             evalc(['data_temp3 = ', strrep(data3{i}, '%', ['wing_VoLaNoS.aeroProp.panelData(', num2str(i_alpha), ').'])]);
%             plot(wing_RML51A26.aeroProp.panelData(i_alpha).y, data_temp3, 'Color', 'r', 'Marker', 'none', 'LineWidth', 2);
%         end
%     end
%     xlabel('y [m]');
%     ylabel(yLabels{i});
%     legend(legends{i}, 'Location', 'southeast');
%     grid on
%     end
    
    
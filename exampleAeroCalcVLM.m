% Script for aerodynamic analysis of one wing by use of VLM
addpath("classes\");
%% Create Wing Geometry
% Baseline wing planform data
AR = 11; % aspect ratio (for TaperedWing object only)
TR = 1; % taper ratio TR = c_t/c_r (for TaperedWing object only)
sweepAngle = 4.5 /180*pi; % sweep of quarter-chord line in rad
% sweepAngle = atan(tan(60/180*pi) - 4/AR*(0.25-0) * (1-TR)/(1+TR));
% dihedralAngle = 0 /180*pi; % diheadral angle in rad
S_ref = 40; % Wing reference area in m² (for TaperedWing object only)

% State
delta_flap = 15 /180*pi; % flap deflection in rad

% Create airfoil object
airfoil1 = Airfoil('Custom'); % Airfoil object type 'Custom' or 'NACA'
% DOA5 Dornier Airfoil 5
airfoil1.importAirfoilCoordinates('DO-A-5.dat'); % import airfoil coordinates for 'Custom' Airfoil object -> relevant if the airfoil is cambered
airfoil1.plotAirfoil(); % plot of airfoil geometry
% NACA 4/5 digit airfoil
airfoil2 = Airfoil('NACA', '0012'); % Defintion of airfoil object by use of a 4/5-digit NACA code

% Create multi-partition wing object from class "MultiPartitionWing
chordTable = [0,  2.7, 10.5; ... % y-station on half wing in [m]
              2.7,  2.5, 1];     % chord length at y-station in [m]

wingObject = MultiPartitionWing(chordTable, [airfoil1, airfoil2]); % MultiPartitionWing object (airfoil IDs are counted in the given order)

% For a straight tapered wing alternatively create a tapered wing object from class "TaperedWing"
% wingObject = TaperedWing(S_ref, [airfoil1, airfoil2], AR, TR); % TaperedWing object (airfoil IDs are counted in the given order)

% Usage of different airfoils
wingObject.addAirfoilZone(0, 1, 1, 1); % addAirfoilZone(relSpanPos1, relSpanPos2, airfoil ID1, airfoil ID2); 

% Sweep
wingObject.addSweepKink(0, sweepAngle); % The quarter-chord (25%) line is changed at y/b=0 (root) to have the given sweep angle up to the tip or next sweep kink

% Dihedral
% wingObject.addDihedralKink(0, dihedralAngle); % The quarter-chord (25%) line is changed at y/b=0 (root) to have the given dihedral/anhedral angle up to the tip or next sweep kink

% Twist
% wingObject.addTwistTransition(0.25, 0.85, -5/180*pi, 'linear'); addTwistTransition(relSpanPos1, relSpanPos2, twist angle in rad, twist type: 'linear' or 'uniform')

% Flap Geometry
wingObject.addFlaps([0.7, 0.95], 0.1 * chordTable(2,2) * [1,1], 1, -1); % addFlaps(relSpanPos, absFlapDepths in m, flapID, flapControlMode: (1)symmetric (-1)mirrored (0)antisymmetric)

% Flap Deflections
wingObject.setFlapDeflection(1, delta_flap); % setFlapDeflection(flap ID, flap deflection in rad)

% Plot wing planform geometry
wingObject.plotGeometry(); % Simplified visualization of the wing

%% Calculate aerodynamics (VLM)
alphas_deg = (-10:5:30); % angle of attack in deg
alphas = alphas_deg /180*pi; % angle of attack in rad
betas = (0:5:10) /180*pi; % sideslip angle in rad
nPanelsX = 10; % chordwise number of panels
nPanelsY = 20; % spanwise number of panels for one half span

% Moment coefficient reference point
[relPosition, c_mac, x_mac] = wingObject.getACPosition(); % Get x-offset of aerodynamic center w.r.t. wing apex (by use of wing geometry)
% Note that the origin of the wing's coordinate system is located in the 25% root chord point.
wingObject.setRefPointMom([relPosition(1) - 0.25 * wingObject.getRootChord(); 0; 0]); % Set the moment coefficient reference point to the estimated aerodynamic center
% Note that by default the moment coefficients are established about the 25% root chord point, i.e. (0|0|0) in the wing objects
% coordinate system. For instance, to move the reference point to the wing apex use wingObject.setRefPointMom([-0.25 * wingObject.getRootChord(); 0; 0]);

wingObject.plotVLMLattice(nPanelsX, nPanelsY); % plot of generated VLM lattice of given panel numbers
results = wingObject.calculateVLM(alphas, betas, nPanelsX, nPanelsY); % calculation of VLM

% Zero-lift drag estimation
V = 150; % flight speed in m/s
nu = 1.46 * 10^-5; % kinematic viscosity
CD0 = wingObject.estimateZeroLiftDrag(V / nu); % Estimation of zero lift drag coefficient by use of empirical relations

% Results of aerodynamic calculation (static coefficients in aerodynamic coordinate system)
VLM.alphas = alphas';
VLM.betas = betas;
VLM.CL = reshape([results().c_L], length(alphas), length(betas));           % Lift coefficient
VLM.CS = reshape([results().c_S], length(alphas), length(betas));           % Side force coefficient
VLM.CD = reshape([results().c_D], length(alphas), length(betas)) + CD0;     % Drag coefficient (VLM method only provides induced drag)
VLM.Cm = reshape([results().c_m], length(alphas), length(betas));           % Pitching moment coefficient based on mean aerodynamic chord
VLM.Cn = reshape([results().c_n], length(alphas), length(betas));           % Yawing moment coefficient based on span (not halfspan)
VLM.Cl = reshape([results().c_l], length(alphas), length(betas));           % Rolling moment coefficient based on span (not halfspan)

% Calculation of dynamic coefficients
rot_X0 = [-0.25 * wingObject.getRootChord(); 0; 0]; % Rotation axis reference point in body fixed coordinate system (here wing apex as an example)

% Roll damping
rot_dir_axis = [1;0;0]; % Direction of rotation axis in body fixed coordinate system
rot_omega = 5 / 180*pi; % roll rate p [rad/s] 
results_0 = wingObject.calculateVLM(0, 0, nPanelsX, nPanelsY); % reference values
results = wingObject.calculateVLM(0, 0, nPanelsX, nPanelsY, rot_omega, rot_X0, rot_dir_axis); % rotation
VLM.Cl_p = (results.c_l - results_0.c_l) / (rot_omega * 0.5 * wingObject.getSpan()); % Note that the VLM method is independet from freestream velocity, i.e. V_inf = 1 m/s for the nondimensional roll rate.

% Pitch damping
rot_dir_axis = [0;1;0]; % Direction of rotation axis in body fixed coordinate system
rot_omega = 5 / 180*pi; % pitch rate q [rad/s] 
results_0 = wingObject.calculateVLM(0, 0, nPanelsX, nPanelsY); % reference values
results = wingObject.calculateVLM(0, 0, nPanelsX, nPanelsY, rot_omega, rot_X0, rot_dir_axis); % rotation
VLM.Cm_q = (results.c_m - results_0.c_m) / (rot_omega * 0.5 * c_mac); % Note that the VLM method is independet from freestream velocity, i.e. V_inf = 1 m/s for the nondimensional roll rate.

% Yaw damping
rot_dir_axis = [0;0;1]; % Direction of rotation axis in body fixed coordinate system
rot_omega = 5 / 180*pi; % yaw rate r [rad/s] 
results_0 = wingObject.calculateVLM(0, 0, nPanelsX, nPanelsY); % reference values
results = wingObject.calculateVLM(0, 0, nPanelsX, nPanelsY, rot_omega, rot_X0, rot_dir_axis); % rotation
VLM.Cn_r = (results.c_n - results_0.c_n) / (rot_omega * 0.5 * wingObject.getSpan()); % Note that the VLM method is independet from freestream velocity, i.e. V_inf = 1 m/s for the nondimensional roll rate.

%% Plot of calculated data
figureFlap = 'VLM aerodynamics';
handleCoeff = findobj('type', 'figure', 'Name', figureFlap);
if isempty(handleCoeff)
    handleCoeff = figure('Name', figureFlap);
end
figure(handleCoeff);
clf

subplot(161);
plot(alphas_deg, VLM.CL, 'LineWidth', 2);
xlabel('\alpha [°]');
ylabel('C_L [-]');
legend(compose('beta = %g°', betas /pi*180), 'Location', 'SouthEast');
grid on;

subplot(162);
plot(VLM.CD, VLM.CL, 'LineWidth', 2);
xlabel('C_D [-]');
ylabel('C_L [-]');
grid on;

subplot(163);
plot(VLM.CS, VLM.CL, 'LineWidth', 2);
xlabel('C_S [-]');
ylabel('C_L [-]');
grid on;

subplot(164);
plot(VLM.Cm, VLM.CL, 'LineWidth', 2);
xlabel('C_m [-]');
ylabel('C_L [-]');
grid on;

subplot(165);
plot(VLM.Cn, VLM.CL, 'LineWidth', 2);
xlabel('C_n [-]');
ylabel('C_L [-]');
grid on;

subplot(166);
plot(VLM.Cl, VLM.CL, 'LineWidth', 2);
xlabel('C_l [-]');
ylabel('C_L [-]');
grid on;
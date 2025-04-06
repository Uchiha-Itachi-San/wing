%
%****************************************************************************************************************
%*                                                                                                              *
%*                               F A C H B E R E I C H  L U F T F A H R T T E C H N I K                         *
%*                           T E C H N I S C H E  H O C H S C H U L E  I N G O L S T A D T                      *
%*                                                                                                              *
%*                                 A E R O N A U T I C A L  E N G I N E E R I N G                               *
%*                            U N I V E R S I T Y  O F  A P P L I E D  S C I E N C E S                          *
%*                                                                                                              *
%*                                     Technische Hochschule Ingolstadt (THI)                                   *
%*                                   Esplanade 10, D-85049 Ingolstadt - Germany                                 *
%*                    Phone: +49 841 9348-4416, eMail: korbinian.stadlberger@thi.de, Web:www.thi.de             *
%*                                                                                                              *
%*   (c) 2024 by Technische Hochschule Ingolstadt (University of Applied Sciences), All Rights Reserved         *
%*                                                                                                              *
%****************************************************************************************************************
%*                                                                                                              *
%*                                                   Wing		                                                *
%*                                                                                                              *
%****************************************************************************************************************
%*                                                                                                              *
%* Type        :   Abstract Class Definition                                                                    *
%*                                                                                                              *
%* Circulation :                                                                    							*
%*                                                                                                              *
%* Purpose     :   Definition of abstract class Wing                                                            *
%*                                                                                                              *
%* Version     :   1.0                                                                                          *
%*                                                                                                              *
%****************************************************************************************************************
%*                                                                                                              *
%* Remarks  :                                                                      								*
%*                                                                                                              *
%****************************************************************************************************************
%*  Author                  *     Date    *     Description                                                     *
%****************************************************************************************************************
%*  Stadlberger, Korbinian  * 2010-SEP-13 *     Initial Release (LLS/TUM)
%*  Stadlberger, Korbinian  * 2024-NOV-25 *     Adaption for THI lecture
%*  Stadlberger, Korbinian  * 2025-FEB-07 *     VLM calculation method added for THI lecture
%****************************************************************************************************************

classdef Wing < AerodynamicGeometry
   
   
%%  %*                                         Properties                                                      *
                                                                     %                                  *
    %*  Property               Description                                   Type             Dim         Unit	* 
    properties (SetAccess = protected)                                                                    
        
        symmetry               % Indicates whether wing is symmetric (1)     boolean          (1x1)       [-]
                               % or not (0)
        aspectRatio            % Aspect ratio                                double           (1x1)       [-]           
        taperRatio             % Taper ratio / meaning depends	 double           (1x1)       [-] 
        					   % on the class inherited from this class                
        nPartitions            % Number of partitions used for the halfspan  int              (1x1)       [-]  
        chordKinks             % Collection of chord kinks                   struct           (1xm)       [-]
        relSpanPositions       % Span positions relative to the halfspan     double           (1xn+1)     [-]             
        chordLengths           % Chord lengths at relative span positions    double           (1xn+1)     [m]              
        twistAngles            % Twist angles for each partition	         double           (2xn)       [rad]
        twistRotPoints         % Twist rotation points for each partition    double           (2xn)       [-]
        dihedralAngles         % Dihedral angles for each span position      double           (1x1)       [rad]       
        sweepAngles            % Sweep angles of quarter chord line for each span position         double           (1x1)       [rad]       
        incidenceAngle = 0     % Incidence angle of entire wing              double           (1x1)       [rad]       
        airfoilDB              % Airfoil database                            Airfoil          (1xn)       [-]       
        airfoilsUsed           % Airfoil IDs for each span position          double           (2xn)       [-]       
        airfoilBlending        % Fraction of inboard airfoil                 double           (2xn)       [-]       
        flapChordLengths       % Relative flap chord length for each         double           (1xn)       [-]       
                               % partition                                                    		       
        flapDeflections        % Flap deflection for each partition          double           (1x1)       [rad]       
        deflectionModes        % Deflection Modes for each partition         int              (1x1)       [-]   
        					   % (1)symmetric (-1)mirrored (0)antisymmetric                               		       
        flapControlIDs		   % Flap Control IDs for each partition         int              (1x1)       [-]   
        hingeLines             % Flap hinge line points for each partition   double           (2xn)       [-]   
        dihedralKinks          % Collection of dihedral kinks                struct           (1x1)       [-]      
        sweepKinks             % Collection of sweep kinks                   struct           (1x1)       [-]      
        twistTransitions       % Collection of twist transitions             struct           (1x1)       [-]      
        airfoilZones           % Collection of airfoil zone definitions      struct           (1x1)       [-]      
        flapBreaks             % Collection of flap breaks                   struct           (1x1)       [-]        
        t_c_max_global         % Maximum thickness ratio                     double           (1x1)       [-]
        devices                % Device objects to account for actuators etc Device           (1xn)       [-]
        k_retractable          % Indicates mass supplement if retractable    double           (1x1)       [-]
        k_FRP                  % Indicates mass reduction due to use of composite materials
        roughness_surface
    end
    
    methods (Abstract)
        
        calculateChordLengths(this)
        % This function is meant to calculate the chord lengths at the relative span positions 
        % depending on the wing's specific planform and its properties (changing during runtime.)
        
% %         % refreshWing(this)
        % Refreshes all important properties which are dependent of other properties having changed 
        % during runtime, which define the wing shape.
                
    end %methods abstract
    
    methods
        
        
        function refreshMass(this, refreshOption)
            % Refreshes the wing mass (property mass) as well as the relative c.g. position of the 
            % wing.
            % Called functions:
            %	getCm0() of class AerodynamicProperties
            %	getRootChord() of this class
            %	getSpan() of this class
            
            if nargin == 1
                refreshOption = 0;
            end
            
            if refreshOption ~= 0
            
            % MASS ESTIMATION TO BE ADDED !
            
            % Devices
            nbDevices = size(this.devices,2);
            sum_mass_devices = 0;
            moment_devices = [0; 0; 0];

            for i = 1:1:nbDevices
                moment_devices = moment_devices + this.devices(i).mass * (this.devices(i).refPoint + this.devices(i).relCG);
                sum_mass_devices = sum_mass_devices + this.devices(i).mass;
            end
            relCG_devices = moment_devices / sum_mass_devices;
            
            this.mass = m_wing_readytofly + sum_mass_devices;
            this.relCG = (sum_mass_devices * relCG_devices...
                + m_wing_readytofly * relCG_wing)...
                / (this.mass);    
            end
        end
        
        function refreshRelCG(this, refreshOption)
            % Refreshes the property relCG. This might be done within the function refreshMass().
            % 
            % Called functions:
            %	refreshMass() of this class
            
            if nargin == 1
                refreshOption = 0;
            end
            
            if refreshOption ~= 0
                this.refreshMass();
                
                % Kundu
                ac = this.getACPosition();
                this.relCG = [ac(1) + ((- 0.25 + 0.30) * this.getRootChord());0;0];
            end
        end
        
        function setIncidenceAngle(this, newAngle)
            % Setter for wing incidence angle (property incidenceAngle).
            %
            % Input variables:
            %	newAngle: incidence angle in radians
            
            this.incidenceAngle = newAngle;
        end
        
        function refreshMCGList(this)
           
            this.mcg_display.relevantObjects = num2cell(this.devices)';
            this.mcg_display.indentString = '   ';
        end
        
        function addDevices(this, newDevices)
            % Adds new devices to the array already existing (property devices).
            %
            % Input variables:
            %	newDevices: array of Device objects to be added, Device (1xn)
            
            if isempty(this.devices)
                this.devices = newDevices;
            else
                this.devices = [this.devices, newDevices];
            end
            
            this.refreshMCGList();
        end
        
        function removeAllDevices(this)
            
            this.devices = [];
            this.refreshMCGList();
        end
        
%         function setSref(this, newSref)
%             % Setter for wing reference area (property Sref).
%             %
%             % Input variables:
%             %	newSref: Reference area in m²
%             %
%             % Called functions:
%             %	refreshWing (abstract) of any class inherited from this class
%             
%             this.Sref = newSref;
%             this.refreshWing();
%             
%         end
%         
%         %{
%         function setRefPoint(this, newRefPoint)
%             % This function only exists because this class is inherited from the abstract class 
%             % 'Geometry'. This function is unnecessary because the wing reference point is always 
%             % [0;0;0] and is the same as the root chord nose. This reference point is the origin of 
%             % the global coordinate system.
%             
%             disp('Wing.setRefPoint(newPoint): Reference point for wing [0 0 0] not modifiable.')
%         end
%         %}
%         
%         function setAR(this, newAR)
%             % Setter for wing aspect ratio (property aspectRatio).
%             %
%             % Input variables:
%             %	newAR: wing aspect ratio
%             %
%             % Called functions:
%             %	refreshWing (abstract) of any class inherited from this class
%             
%             this.aspectRatio = newAR;
%             this.refreshWing();
%         end
%         
%         function setETR(this, newETR)
%             % Setter for wing equivalent taper ratio (property taperRatio).
%             %
%             % Input variables:
%             %	newETR: wing equivalent taper ratio
%             %
%             % Called functions:
%             %	refreshWing (abstract) of any class inherited from this class
%             
%             this.taperRatio = newETR;
%             this.refreshWing();
%         end
%         
%         function setARSTR(this, newAR, newS, newTR)
%             
%             this.aspectRatio = newAR;
%             this.Sref = newS;
%             this.taperRatio = newTR;
%             this.refreshWing();
%         end
        
        function refreshSref(this)
            % This function only exists because this class is inherited from the abstract class 
            % 'AerodynamicGeometry'. This function is unnecessary because the wing reference area is 
            % a design variable and therefore does not need to be refreshed.
        end
        
        function refreshLref(this)
            % Refreshes the reference length of the wing by assigning the mean  aerodynamic chord (MAC).

            [~, this.Lref, ~] = this.getACPosition();
        end

        function [c_MAC, x_MAC] = getMAC(this) % !!!! Check redundancy to function getACPosition !!!
            %Calculation of the mean aerodynamic chord (MAC)
            % Returns the mean aerodynamic chord (MAC).
            % Called functions:
            %	getSpan() of this class

            span = this.getSpan();
            sumS = 0;
            sumSc = 0;
            sumSx = 0;
            sumSy = 0;
            
            for i = 1:1:this.nPartitions
                
                rootChord_partition = this.chordLengths(i);
                tipChord_partition = this.chordLengths(i+1);
                taperratio_partition = tipChord_partition / rootChord_partition;
                span_partition = span * (this.relSpanPositions(i+1) - this.relSpanPositions(i));
                S_partition = span_partition / 2 * (rootChord_partition + tipChord_partition) / 2;
                sumS = sumS + S_partition;
                c_MAC_partition = rootChord_partition * 2/3 * (( 1 + taperratio_partition + taperratio_partition^2 ) / ( 1 + taperratio_partition ));
                sumSc = sumSc + c_MAC_partition * S_partition;
                sweepAngle_0 = this.getSweepAngle(0.5 * (this.relSpanPositions(i+1) - this.relSpanPositions(i)), 0);
                x_MAC_partition = span_partition / 2 * tan(sweepAngle_0) * (2 * taperratio_partition + 1) / (3 * taperratio_partition + 3);
                sumSx = sumSx + x_MAC_partition * S_partition;
            end

            c_MAC = sumSc / sumS;
            x_MAC = sumSx / sumS;
        end
        
        function J = getInertia(this)
           
            % This function calculates and returns the moments of inertia around the c.g. of
            % the geometry and overloads the function inherited by the class 'Geometry'. 
            % Called functions:
				%	refreshRelCG: refreshing the c.g. position of this object
            %    
            % Output:
                %   J: vector (3x1) containing the moments of inertia around the three axis
            this.refreshSize(); 
            
            J = zeros(3,1);
            
            % Devices
            for i = 1:1:size(this.devices, 2)
                delta = this.relCG + this.devices(i).refPoint + this.devices(i).relCG;
                dist = [sqrt(delta(2)^2 + delta(3)^2);
                       sqrt(delta(1)^2 + delta(3)^2);
                       sqrt(delta(1)^2 + delta(2)^2)];
                J = J + this.devices(i).getInertia() + this.devices(i).mass * dist.^2;
            end
            
            % Wing panels
            % Approximation by homogenious cube!!!
            
            J(1) = J(1) + (this.size(2)^2 + this.size(3)^2) * this.mass / 12;
            J(2) = J(2) + (this.size(1)^2 + this.size(3)^2) * this.mass / 12;
            J(3) = J(3) + (this.size(1)^2 + this.size(2)^2) * this.mass / 12;
        end
        
        function span = getSpan(this)
            % Returns the span of the wing using its aspect ratio and reference area.
            %
            % Output variables:
            %	span: span of the wing in m
            
            if this.symmetry == 1
                span = sqrt(this.Sref * this.aspectRatio);
            else
                span = sqrt(2 * this.Sref * this.aspectRatio);
            end
        end
        
        function rootChord = getRootChord(this)
            % Returns the root chord length of the wing.
            % 
            % Output variables:
            %	rootChord: wing root chord in m
            
            rootChord = this.chordLengths(1);
        end

        function tipChord = getTipChord(this)
            % Returns the tip chord length of the wing.
            % 
            % Output variables:
            %	tipChord: wing tip chord in m
            
            tipChord = this.chordLengths(end);
        end

        function twistAngle = getTwistAngle(this, yPosRel)
            % Returns the twist angle at station with relative position yPosRel
            % 
            % Output variables:
            %	twistAngle: section twist angle in rad
            
            twistAngle = interp1(this.relSpanPositions, this.twistAngles, yPosRel);
        end

        function dihedralAngle = getDihedralAngle(this, yPosRel)
            % Returns the dihedral angle at station with relative position yPosRel
            % 
            % Output variables:
            %	dihedralAngle: section dihedral angle in rad
            
            dihedralAngle = interp1(this.relSpanPositions, [this.dihedralAngles, this.dihedralAngles(end)], yPosRel, 'previous');
        end

        function sweepAngle = getSweepAngle(this, yPosRel, relXPos)
            % Returns the sweep angle of the relXPos chord line at relative spanwise position yPosRel
            % 
            % Output variables:
            %	sweepAngle: section sweep angle in rad
            
            sweepAngle = interp1(this.relSpanPositions, [this.sweepAngles, this.sweepAngles(end)], yPosRel, 'previous');
            i_yPos = interp1(this.relSpanPositions, 1:length(this.relSpanPositions), yPosRel, 'previous');

            if nargin > 2
                rootChord_partition = this.chordLengths(i_yPos);
                tipChord_partition = this.chordLengths(i_yPos+1);
                taperratio_partition = tipChord_partition / rootChord_partition;
                span_partition = this.getSpan() * (this.relSpanPositions(i_yPos+1) - this.relSpanPositions(i_yPos));
                AR_partition = 2 * span_partition / (rootChord_partition * (1 + taperratio_partition));
                sweepAngle = atan(tan(sweepAngle) - 4 / AR_partition * (relXPos - 0.25) * (1-taperratio_partition)/(1+taperratio_partition));
            end
        end
        
        function setNPanels(this, globalNPanels)
            % Setter for the number of panels each partition will be devided into (property 
            % nPanels).
            %
            % Input variables:
            %	globalNPanels: number of panels for each partition (1x1 double)
            
            this.nPanels = ones(1,this.nPartitions) * globalNPanels;
        end
        
        function setFlapProperties(this, flapChords, flapIDs, flapControlModes) % OBSOLETE ??????
            % Setter for the definition of flap configuration. It defines the relative lengths of 
            % flap chords, the flap identification numbers, and the control modes (if deflected 
            % symetrically or antisymetrically) for each partition.
            %
            % Input variables:
            %	flapChords: relative chord lengths of the flaps, double (1xn)
            %	flapIDs: flap IDs, int (1xn)
            %	flapControlMOdes: control mode for each control mode (1 symetric, -1 antisymetric), 
            %	int (1xn)
            
            if (size(flapChords,2) ~= this.nPartitions)
                disp(strcat(['Dimension of input vector for flap chords must have dim 1x', num2str(this.nPartitions)]));
            elseif (size(flapIDs,2) ~= this.nPartitions)
                disp(strcat(['Dimension of input vector for flap IDs must have dim 1x', num2str(this.nPartitions)]));
            elseif (size(flapControlModes,2) ~= this.nPartitions)
                disp(strcat(['Dimension of input vector for flap control modes (symmetric=1, antisymmetric=-1) must have dim 1x', num2str(this.nPartitions)]));
            else
                this.flapChordLengths = flapChords;
                this.flapControlIDs = flapIDs;
                this.deflectionModes = flapControlModes;
            end
        end
        
        function addFlaps(this, relSpanPos, absFlapDepths, flapID, flapControlMode, absHingeLineShift)
           
            if nargin == 5
                absHingeLineShift = [0,0];
            end
            
            [relSpanPos, order] = sort(relSpanPos);
            absFlapDepths = absFlapDepths(order);
            
            if relSpanPos(1) >= 1
                error('Flap inboard position exceeds wing span.');
            end
            if relSpanPos(2) > 1
                relSpanPos(2) = 1;
                disp('WARNING: Flap outboard position exceeds wing span. Outboard position will be set to wing tip.');
            end
            
            checkOK = 1;
            for i = 1:size(this.flapBreaks,2)
                if (relSpanPos(1) > this.flapBreaks(i).relYPos1 && relSpanPos(1) < this.flapBreaks(i).relYPos2) || ...
                   (relSpanPos(2) > this.flapBreaks(i).relYPos1 && relSpanPos(2) < this.flapBreaks(i).relYPos2)
                    checkOK = 0;
                end
            end
            
            
            if checkOK == 1

                nK = size(this.flapBreaks,2);
                this.flapBreaks(nK+1).relYPos1 = relSpanPos(1);
                this.flapBreaks(nK+1).relYPos2 = relSpanPos(2);
                this.flapBreaks(nK+1).flapID = flapID;
                this.flapBreaks(nK+1).flapControlMode = flapControlMode;
%                 this.flapBreaks(nK+1).hingeLinePoints = absFlapDepths + absHingeLineShift;

                % Sort flaps along span
                [~, order] = sort([this.flapBreaks(:).relYPos1]);
                this.flapBreaks = this.flapBreaks(order);
            
                % Recalculate the span positions of sections, their chord lengths
                this.calculateRelSpanPositions();
                this.calculateChordLengths();
                
                % Calculate relative flap depths
                fractions = absFlapDepths ./ this.getChordInterpolated(relSpanPos);
                if fractions(2) > 1 && fractions(1) < 1
                    nStripes = 100;
                    yPos = relSpanPos(1) + (0:nStripes)/nStripes * (relSpanPos(2) - relSpanPos(1));
                    wingChords = this.getChordInterpolated(yPos);
                    flapChords = interp1(relSpanPos, absFlapDepths, yPos);
                    [~, index] = min(abs(wingChords - flapChords));
                    this.flapBreaks(nK+1).relYPos2 = yPos(index);
                    this.flapBreaks(nK+1).fractions = [fractions(1), 1];
                    this.flapBreaks(nK+2).relYPos1 = this.flapBreaks(nK+1).relYPos2;
                    this.flapBreaks(nK+2).relYPos2 = relSpanPos(2);
                    this.flapBreaks(nK+2).fractions = [1,1];
                    this.flapBreaks(nK+2).flapID = flapID;
                    this.flapBreaks(nK+2).flapControlMode = flapControlMode;
%                     this.flapBreaks(nK+2).hingeLinePoints = absFlapDepths + absHingeLineShift;
                    
                    % Recalculate the span positions of sections, their chord lengths
                    this.calculateRelSpanPositions();
                    this.calculateChordLengths();
                elseif fractions(1) > 1 && fractions(2) > 1
                    this.flapBreaks(nK+1).fractions = [1,1];
                else
                    this.flapBreaks(nK+1).fractions = fractions;
                end
                this.update4newNPartitions();
            else
                
                disp('Flap could not be added because of collision with an already existing flap at desired span position');
            end
        end
        
        function removeAllFlaps(this)
            
            this.flapBreaks = struct('relYPos1', {}, 'relYPos2', {}, 'fraction', {}, 'flapID', {}, 'flapControlMode', {});
            this.initialiseFlapChordLengths();
            this.initialiseFlapDeflections();
            this.initialiseDeflectionModes();
            this.initialiseFlapControlID();
        end
        
        function chordLength = getChordInterpolated(this, relYPos)
           
            chordLength = interp1(this.relSpanPositions, this.chordLengths, relYPos);
        end
        
        function fraction = getFlapChordFraction(this, relYPos)
           
            fraction = zeros(1,length(relYPos(:)));
            for j = 1:length(relYPos)
                index = this.findIndicesPanels(relYPos);
                fraction(j) = interp1(this.relSpanPositions(index:index+1), ...
                                   [this.chordLengths(index) * this.flapChordLengths(1,index), this.chordLengths(index+1) * this.flapChordLengths(2,index)], ...
                                   relYPos) / this.getChordInterpolated(relYPos);
                if isnan(fraction)
                    disp('stop');
                end
            end
        end
        
        function setFlapDeflection(this, ID, deflectionAngle)
            % Sets the current flap deflection for all flaps identified by the given ID to the given 
            % deflection angle.
            % 
            % Input variables:
            % 	ID: identification number of flaps to be deflected, int (1x1)
            % 	deflectionAngle: deflection angle in radians for identified flaps, double (1x1)
            
            this.flapDeflections(this.flapControlIDs == ID) = deflectionAngle;
            
%             for i = 1:this.nPartitions ----> obsolete
%                 
%                 if this.flapControlIDs(i) == ID
%                     
%                     this.flapDeflections(i) = deflectionAngle;
%                 end
%             end
        end
        
        function setGlobalAirfoilID(this, newAirfoilID)
            % Assigns the same airfoil using its airfoil ID to every partition, so that the same 
            % airfoil is used globally for the entire wing.
            % 
            % Input variables:
            % 	newAirfoilID: airfoil ID number assigned to every partition, int (1x1)
            
            if (newAirfoilID <= size(this.airfoilDB,2)) || isempty(this.airfoilDB)
                this.airfoilsUsed = ones(2, this.nPartitions) * newAirfoilID;
                this.airfoilBlending = ones(2, this.nPartitions);
            else
                error('Airfoil data bank too small to assign desired airfoil ID!');
            end
        end
        
        function setNPartitions(this, newNPartitions)
            % Setter for the number of partitions used to model the wing (property nPartitions).                                                                     
				%                                                                               
				% Input variables:                                                              
				%	globalNPanels: number of panels for each partition (1x1 double)  
				%
				% Called functions:
				%	refreshWing() (abstract) of this class            

            this.nPartitions = newNPartitions;
            
            refreshWing(this);
        end
        
        
        function addDihedralKink(this, relYPos, dihedralAngle)
            % Adds a dihedral kink at a specific relative span positon. The information is stored in 
            % a struct for each kink and added to the struct array of property dihedralKinks.
            % 
            % Input variables:
            % 	relYPos: relative span position of dihedral kink, double (1x1)
            % 	dihedralAngle: dihedral angle in radians, double (1x1)
            % 
            % Called functions:
            % 	calculateRelSpanPositions() of this class
            % 	calculateChordLengths() (abstract) of this class
            % 	setDihedralAngles() of this class
            
            if relYPos < 1
                nK = size(this.dihedralKinks,2);
                this.dihedralKinks(nK+1).relYPos = relYPos;
                this.dihedralKinks(nK+1).angle = dihedralAngle;

                % Sort kinks along span
                [~, order] = sort([this.dihedralKinks(:).relYPos]);
                this.dihedralKinks = this.dihedralKinks(order);

                % Recalculate the span positions of sections, their chord lengths and the dihedral angles of the partitions
                this.calculateRelSpanPositions();
                this.calculateChordLengths();
                this.update4newNPartitions();
            else
                error('Relative span position must be smaller than 1!');
            end
        end
        
        function deleteAllDihedralKinks(this)
            % Deletes all existing dihedral kinks.
            %
            % Called functions:
            % 	calculateRelSpanPositions() of this class        
				% 	calculateChordLengths() (abstract) of this class 
				% 	setDihedralAngles() of this class                
              
            % Create an empty struct                                                 
            this.dihedralKinks = struct('relYPos', {}, 'angle', {}, 'affectedPanels', {});
            
            % Recalculate the span positions of sections, their chord lengths and the dihedral angles of the partitions
            this.calculateRelSpanPositions();
            this.calculateChordLengths();
            this.update4newNPartitions();
        end
        
        function deleteAllSweepKinks(this)
            % Deletes all existing sweep kinks.
            %
            % Called functions:
            % 	calculateRelSpanPositions() of this class        
				% 	calculateChordLengths() (abstract) of this class 
				% 	setSweepAngles() of this class
				
				% Create an empty struct
            this.sweepKinks = struct('relYPos', {}, 'angle', {}, 'affectedPanels', {});
            
            % Recalculate the span positions of sections, their chord lengths and the sweep angles of the partitions
            this.calculateRelSpanPositions();
            this.calculateChordLengths();
            this.update4newNPartitions();
        end
        
        function addSweepKink(this, relYPos, sweepAngle)
            % Adds a sweep kink at a specific relative span positon. The information is stored in 
            % a struct for each kink and added to the struct array of property sweepKinks.
            % 
            % Input variables:
            % 	relYPos: relative span position of sweep kink, double (1x1)
            % 	sweepAngle: sweep angle in radians, double (1x1)
            % 
            % Called functions:
            % 	calculateRelSpanPositions() of this class
            % 	calculateChordLengths() (abstract) of this class
            % 	setSweepAngles() of this class
            
            if relYPos < 1
                nK = size(this.sweepKinks,2);
                this.sweepKinks(nK+1).relYPos = relYPos;
                this.sweepKinks(nK+1).angle = sweepAngle;

                % Sort kinks along span
                [~, order] = sort([this.sweepKinks(:).relYPos]);
                this.sweepKinks = this.sweepKinks(order);

                % Recalculate the span positions of sections, their chord lengths and the sweep angles of the partitions
                this.calculateRelSpanPositions();
                this.calculateChordLengths();
                this.update4newNPartitions();
            else
                error('Relative span position must be smaller than 1!');
            end
        end

        function update4newNPartitions(this)

            this.setSweepAngles();
            this.setDihedralAngles();
            this.setTwistAngles();
            this.setAirfoilIDsAndBlending();
            this.setFlapChordLengths();
            this.initialiseFlapDeflections;
        end
        
        function addTwistTransition(this, relYPos1, relYPos2, twistAngle, transitionType, rotPoint, shapeMatrix)
            % Adds a transition zone between two specific relative span positons. The information is stored in 
            % a struct for each zone and added to the struct array of property twistTransitions.
            % 
            % Input variables:
            % 	relYPos1: relative span position of the beginning of the transition, double (1x1) (unused for "singular")
            %   relYPos2: relative span position of the end of the transition, double (1x1)
            % 	twistAngle: twist angle at the end of the transition in radians, double (1x1)
            % 	transitionType: type of twist transition ('singular', 'linear', 'custom'), string (1x1) 
            %   shapeMatrix: matrix containing points representing a custom
            %   shape for the transition length, double (2xn)
            % 
            % Called functions:
            % 	calculateRelSpanPositions() of this class
            % 	calculateChordLengths() (abstract) of this class
            % 	
            
            if nargin == 5
                shapeMatrix = [];
                rotPoint = 0.25;
            elseif nargin == 6
                shapeMatrix = [];
            end
            
            if relYPos2 == 0
                disp('Input of twist transition ignored. Transition cannot be set in wing root. Use incidence angle instead.');
            elseif sum(find([this.twistTransitions.relYPos2] == relYPos2)) ~= 0
                fprintf('Input of twist transition ignored. Endpoint of transition already exists: y = %f (relative)\n', relYPos2);
            else
                
                nK = size(this.twistTransitions,2);
                if strcmp(transitionType, 'singular')
                    this.twistTransitions(nK+1).relYPos1 = [];
                else
                    this.twistTransitions(nK+1).relYPos1 = relYPos1;
                end
                this.twistTransitions(nK+1).relYPos2 = relYPos2;
                this.twistTransitions(nK+1).angle = twistAngle;
                this.twistTransitions(nK+1).rotPoint = rotPoint;
                this.twistTransitions(nK+1).transitionType = transitionType;
                this.twistTransitions(nK+1).shapeMatrix = shapeMatrix;

                % Sort transitions along span
                [~, order] = sort([this.twistTransitions(:).relYPos2]);
                this.twistTransitions = this.twistTransitions(order);

                % Recalculate the span positions of sections, their chord
                % lengths and the twist angles of the partition borders
                this.calculateRelSpanPositions();
                this.calculateChordLengths();
            end
                this.update4newNPartitions();
        end
        
        function addAirfoilZone(this, relYPos1, relYPos2, airfoilID1, airfoilID2, zoneType)
            % Adds a transition zone between two specific relative span positons. The information is stored in 
            % a struct for each zone and added to the struct array of property airfoilZones.
            % 
            % Input variables:
            % 	relYPos1: relative span position of the beginning of the transition, double (1x1)
            %   relYPos2: relative span position of the end of the transition, double (1x1)
            %   airfoilID1: airfoil ID at relYPos1, int (1x1)
            %   airfoilID2: airfoil ID at relYPos2, int (1x1)
            % 	zoneType: type of airfoil zone definition ('uniform', 'linear'), string (1x1)
            % 
            % Called functions:
            % 	calculateRelSpanPositions() of this class
            % 	calculateChordLengths() (abstract) of this class
            % 	
            
            if nargin < 6
                zoneType = 'linear';
            end

            if airfoilID1 > length(this.airfoilDB) || airfoilID2 > length(this.airfoilDB)
                error('Airfoil ID exceeds number of airfoils N = %i in airfoil data bank "Wing.airfoilDB".', length(this.airfoilDB));
            end

            nK = size(this.airfoilZones,2);
            for i = 1:nK
                if (this.airfoilZones(i).relYPos1 < relYPos1 && this.airfoilZones(i).relYPos2 > relYPos1) ...
                || (this.airfoilZones(i).relYPos1 < relYPos2 && this.airfoilZones(i).relYPos2 > relYPos2)
                    error('Airfoil transition zones may not overlap!');
                end
            end

            this.airfoilZones(nK+1).relYPos1 = relYPos1;
            this.airfoilZones(nK+1).relYPos2 = relYPos2;
            this.airfoilZones(nK+1).airfoilID1 = airfoilID1;
            
            this.airfoilZones(nK+1).zoneType = zoneType;

            if strcmp(zoneType, 'uniform')
                fprintf('Uniform airfoil zone definition: second airfoil ID ignored\n');
                this.airfoilZones(nK+1).airfoilID2 = [];
            else
                this.airfoilZones(nK+1).airfoilID2 = airfoilID2;
            end
            % Sort transitions along span
            [~, order] = sort([this.airfoilZones(:).relYPos2]);
            this.airfoilZones = this.airfoilZones(order);

            % Recalculate the span positions of sections, their chord
            % lengths and the twist angles of the partition borders
            this.calculateRelSpanPositions();
            this.calculateChordLengths();
            this.update4newNPartitions();
        end
        
        function deleteAllTwistTransitions(this)
            % Deletes all existing twist transition zones.
            %
            % Called functions:
            % 	calculateRelSpanPositions() of this class        
				% 	calculateChordLengths() (abstract) of this class 
				% 	setTwistAngles() of this class
				
				% Create an empty struct
            this.twistTransitions = struct('relYPos1', {}, 'relYPos2', {}, 'angle', {}, 'transitionType', {});
            
            % Recalculate the span positions of sections, their chord lengths and the sweep angles of the partitions
            this.calculateRelSpanPositions();
            this.calculateChordLengths();
            this.update4newNPartitions();
        end
        
        function [posVectors, sweepAngles, chordLengths] = getCorrectedLiftingLine(this)
            
            % Correct lifting line position (Ref. Barnes, J.P.), i.e.
            % vectors and chord lengths
            posVectors = this.getPositionVectors();
            % Correction for sweep and taper
            sweepAngles_extrap = [this.sweepAngles(1), 0.5 * (this.sweepAngles(1:end-1) + this.sweepAngles(2:end)), this.sweepAngles(end)];
            deltaXi_sweep_taper = -0.026 * exp(-this.aspectRatio / 7) + 0.12 * sign(sweepAngles_extrap) * exp(-this.aspectRatio / 3);
            % Correction for root and tip effects
            mu = 0.12 * (this.aspectRatio / 5)^0.7;
            deltaXi_root_tip = mu * tan(sweepAngles_extrap) .* (exp(-posVectors(2,:) / (mu * this.chordLengths(1))) - exp((posVectors(2,:) - posVectors(2,end)) / (mu * this.chordLengths(end))));
            deltaXi_total = deltaXi_sweep_taper + deltaXi_root_tip;
            deltaXi_total(isnan(deltaXi_total)) = 0;
            posVectors(1,:) =  posVectors(1,:) + deltaXi_total .* this.chordLengths;
            sweepAngles = atan((posVectors(1,2:end) - posVectors(1,1:end-1)) ./ (posVectors(2,2:end) - posVectors(2,1:end-1)));
            chordLengths = this.chordLengths .* (1 - 2 * deltaXi_total);
        end
        
        function handleFigure = plotGeometry(this, offsetXYZ)
            % Plots the wing geometry. By now, this function is incomplete. 
            % Optionally a translation vector can be given to shift the geometry to another point 
            % (only for visualisation). Otherwise the reference point is situated at [0;0;0].
            % This function is inherited from the abstract class 'Geometry'.
            % 
            % Input variables:
            % 	translation: Column vector to shift the geometries reference point (only for 
            % 	visualisation), double (3x1)
            % 	
            % Called functions:
				%	
            
            if nargin == 2
                shiftVector = offsetXYZ;
            else
                shiftVector = [0;0;0];
            end
            
            % PLOT FUNCTION TO BE TESTED/EXTENDED !
            % Get coordinates of quarter chord points
            positionVectors = this.getPositionVectors();

            % Half wing
            for j = 1:this.nPartitions

                % Leading edge and trailing edge
%                 x_LE(i) = shiftVector(1) + positionVectors(i,1) + (this.twistRotPoints(i) - 0.25) * this.chordLengths(i) - cos(this.twistAngles(i)) * this.twistRotPoints(i) * this.chordLengths(i);
%                 x_TE(i) = x_LE(i) + cos(this.twistAngles(i)) * this.chordLengths(i);
%                 z_LE(i) = shiftVector(3) + positionVectors(i,3) + sin(this.twistAngles(i)) * this.twistRotPoints(i) * this.chordLengths(i);
%                 z_TE(i) = z_LE(i) - sin(this.twistAngles(i)) * this.chordLengths(i);
                x_LE_inboard(j) = shiftVector(1) + positionVectors(1,j) + (this.twistRotPoints(1,j) - 0.25) * this.chordLengths(j) - cos(this.twistAngles(1,j)) * this.twistRotPoints(1,j) * this.chordLengths(j);
                x_TE_inboard(j) = x_LE_inboard(j) + cos(this.twistAngles(1,j)) * this.chordLengths(j);
                z_LE_inboard(j) = positionVectors(3,j) + sin(this.twistAngles(1,j)) * this.twistRotPoints(1,j) * this.chordLengths(j);
                z_TE_inboard(j) = z_LE_inboard(j) - sin(this.twistAngles(1,j)) * this.chordLengths(j);
                x_LE_outboard(j) = positionVectors(1,j+1) + (this.twistRotPoints(2,j) - 0.25) * this.chordLengths(j+1) - cos(this.twistAngles(2,j)) * this.twistRotPoints(2,j) * this.chordLengths(j+1);
                x_TE_outboard(j) = x_LE_outboard(j) + cos(this.twistAngles(2,j)) * this.chordLengths(j+1);
                z_LE_outboard(j) = shiftVector(3) + positionVectors(3,j+1) + sin(this.twistAngles(2,j)) * this.twistRotPoints(2,j) * this.chordLengths(j+1);
                z_TE_outboard(j) = z_LE_outboard(j) - sin(this.twistAngles(2,j)) * this.chordLengths(j+1);
                y_LE_inboard(j) = shiftVector(2) + positionVectors(2,j);
                y_LE_outboard(j) = shiftVector(2) + positionVectors(2,j+1);
            end
            y_sections = shiftVector(2) + positionVectors(2,:);

            % Flap hinge lines
            x_hinge = [];
            y_hinge = [];
            z_hinge = [];

            for j = 1:this.nPartitions
                
                x_hinge(end+1) = x_LE_inboard(j) + (1 - this.flapChordLengths(1,j)) * (x_TE_inboard(j) - x_LE_inboard(j));
                x_hinge(end+1) = x_LE_outboard(j) + (1 - this.flapChordLengths(2,j)) * (x_TE_outboard(j) - x_LE_outboard(j));
                z_hinge(end+1) = z_LE_inboard(j) + (1 - this.flapChordLengths(1,j)) * (z_TE_inboard(j) - z_LE_inboard(j));
                z_hinge(end+1) = z_LE_outboard(j) + (1 - this.flapChordLengths(2,j)) * (z_TE_outboard(j) - z_LE_outboard(j));
                y_hinge(end+1) = y_LE_inboard(j);
                y_hinge(end+1) = y_LE_outboard(j);
            end

            

            % Plot of Wing Planform Geometry
            figureA = 'Plot of Wing Planform Geometry';
            handleFigure = findobj('type', 'figure', 'Name', figureA);
            if isempty(handleFigure)
                handleFigure = figure('Name', figureA);
            end
            figure(handleFigure);
            clf

            % General plotting settings
            color_LE_TE = 'b';
            lineWidth_LE_TE = 4;
            color_hinge = 'c';
            lineWidth_hinge = 2;

            % Mirror Half Wing Geometry for Plot
            x2plot_LE = [x_LE_inboard;x_LE_outboard];
            x2plot_LE = [fliplr(x2plot_LE(:)'),x2plot_LE(:)'];
            y2plot = [y_LE_inboard;y_LE_outboard];
            y2plot = [fliplr(y2plot(:)'),-y2plot(:)'];
            z2plot_LE = [z_LE_inboard;z_LE_outboard];
            z2plot_LE = [fliplr(z2plot_LE(:)'),z2plot_LE(:)'];
            x2plot_TE = [x_TE_inboard;x_TE_outboard];
            x2plot_TE = [fliplr(x2plot_TE(:)'),x2plot_TE(:)'];
            z2plot_TE = [z_TE_inboard;z_TE_outboard];
            z2plot_TE = [fliplr(z2plot_TE(:)'),z2plot_TE(:)'];
            x2plot_hinge = [fliplr(x_hinge),x_hinge];
            y2plot_hinge = [fliplr(y_hinge),-y_hinge];
            z2plot_hinge = [fliplr(z_hinge),z_hinge];

            % Leading Edge Line
            plot3(x2plot_LE, y2plot, z2plot_LE, 'Color', color_LE_TE, 'LineWidth', lineWidth_LE_TE);
            hold on
            % Flap Hinge Lines
            plot3(x2plot_hinge, y2plot_hinge, z2plot_hinge, 'Color', color_hinge, 'LineWidth', lineWidth_hinge);
            % Trailing Edge Line
            plot3(x2plot_TE, y2plot, z2plot_TE, 'Color', color_LE_TE, 'LineWidth', lineWidth_LE_TE);
            % Quarter-chord Line
            plot3([fliplr(positionVectors(1,:)),positionVectors(1,:)], [-fliplr(positionVectors(2,:)),positionVectors(2,:)], [fliplr(positionVectors(3,:)),positionVectors(3,:)], 'Color', 'red', 'LineWidth', 2);
            % Section Lines
            for j = 1:length(x2plot_LE)
                if j == 1 || j == length(x2plot_LE) % Tip sections
                    color_sectionLine = color_LE_TE;
                    lineWidth_SectionLine = lineWidth_LE_TE;
                else % All other sections
                    color_sectionLine = 'black';
                    lineWidth_SectionLine = 2;
                end
                plot3([x2plot_LE(j), x2plot_TE(j)], [y2plot(j), y2plot(j)], [z2plot_LE(j), z2plot_TE(j)],  'Color', color_sectionLine, 'LineWidth', lineWidth_SectionLine);
            end

            % Flaps




            axis equal
        end
        
        function [relPosition, c_mac, x_mac] = getACPosition(this)
            % Calculates the relative position of the aerodynamic center of the wing in x-direction. 
            % The coordinates of the y- and z-direction are set to zero. For the calculation a 
            % method by XXXXX???? was applied.
            % 
            % Output variables:
            % 	relPosition: column vector containing the relative a.c. position, double (3x1)
            % 	
            % Called functions:
            % 	getPositionVectors() of this class
            
            quarterlinePositions = this.getPositionVectors();
            sumS = 0;
            sum = 0;
            sumle = 0;
            refX = quarterlinePositions(1,1) - this.getRootChord()/4;
            
            for i = 1:1:this.nPartitions
            
                rootChord = this.chordLengths(i);
                tipChord = this.chordLengths(i+1);
                spanPanel = this.getSpan() / 2 * (this.relSpanPositions(i+1) - this.relSpanPositions(i));
                SPartition = spanPanel * (rootChord + tipChord) / 2;
                sumS = sumS + SPartition;
                
                x_here = quarterlinePositions(1,i) - this.chordLengths(i) / 4 - refX;
                x_next = quarterlinePositions(1,i+1) - this.chordLengths(i+1) / 4 - refX;
                
                sum = sum + spanPanel * (x_here * (2 * this.chordLengths(i) + this.chordLengths(i+1)) + x_next * (this.chordLengths(i) + 2 * this.chordLengths(i+1)));
                 
                sumle = sumle + spanPanel * (this.chordLengths(i)^2 + this.chordLengths(i) * this.chordLengths(i+1) + this.chordLengths(i+1)^2);
            end
            
            c_mac = sumle / (3 * sumS);
            x_mac = sum / (3 * 2 * sumS);
            
            xPos = x_mac + c_mac / 4;
            
            relPosition =  [xPos; 0; 0];
        end
        
        function area = getVerticalProjectedArea(this)
            % Calculates the vertical projected area of the wing by summing up the projected area of 
            % each partition.
            % 
            % Output variables:
            % 	area: projected area in m², double (1x1)
            % 	
            % Called functions:
            % 	getSpan() of this class
            
            % Sum up projected area for one halfspan
            sumS = 0;
            for i = 1:1:this.nPartitions
            
                rootChord = this.chordLengths(i);
                tipChord = this.chordLengths(i+1);
                spanPanel = this.getSpan() / 2 * (this.relSpanPositions(i+1) - this.relSpanPositions(i));
                SPartition = sin(this.dihedralAngles(i)) * spanPanel * (rootChord + tipChord) / 2;
                sumS = sumS + SPartition;
             end
             
            area = sumS;
        end
        
        function area = getHorizontalProjectedArea(this)
            % Calculates the horizontal projected area of the wing by summing up the projected area of 
            % each partition.
            % 
            % Output variables:
            % 	area: projected area in m², double (1x1)
            % 	
            % Called functions:
            % 	getSpan() of this class
            
            % Sum up projected area for one halfspan
            sumS = 0;
            for i = 1:1:this.nPartitions
            
                rootChord = this.chordLengths(i);
                tipChord = this.chordLengths(i+1);
                spanPanel = this.getSpan() / 2 * (this.relSpanPositions(i+1) - this.relSpanPositions(i));
                SPartition = cos(this.dihedralAngles(i)) * spanPanel * (rootChord + tipChord) / 2;
                sumS = sumS + SPartition;
             end
            
            % Double area to get projected area of entire wing          
            area = 2 * sumS;
        end
        
        function refreshSize(this)
            % Refreshes the property size by searching the minimum and maximum position of the 
            % quarterline points
            %
            % Called functions:
            %	getSpan() of this class
           
            positionVectors = this.getPositionVectors();
            
            % X-direction: quarter chord points
            [xMin, iXMin] = min(positionVectors(1,:));
            [xMax, iXMax] = max(positionVectors(1,:));
            this.size(1) = xMax + this.chordLengths(iXMax) * 3/4 - xMin - this.chordLengths(iXMin) / 4;
            
            % Y-direction
            this.size(2) = this.getSpan();
            
            % Z-direction: quarter chord points
            zMin = min(positionVectors(3,:));
            zMax = max(positionVectors(3,:));
            this.size(3) = zMax  - zMin;
        end
        
        function positionVectors = getPositionVectors(this)
            % Calculates and returns the position vectors of the quarterline point of each section 
            % point of the halfspan where all dihedral and sweep kinks are taken into account.
            % This function still has to be validated whether the combination of dihedral and sweep kinks 
            % works properly.
            %
            % Output variables:
            %	positionVectors: 	matrix containing the absolute positions of quaterline section 
            % 							point, double (3xn+1)
            %
            % Called functions:
            %	getSpan() of this class
            
            nP = this.nPartitions;
            span = this.getSpan();
            positionVectors = zeros(3,nP+1);
            positionVectors(:,1) = [0;0;0];
            
            for i = 2:1:nP+1
                % Calculate the distance difference in y-direction between the sections
                deltaY = (this.relSpanPositions(i) - this.relSpanPositions(i-1)) * span/2;
                
                % x-Position
                positionVectors(1,i) = positionVectors(1,i-1) + deltaY * tan(this.sweepAngles(i-1));
                
                % y-Position
                positionVectors(2,i) = positionVectors(2,i-1) + cos(this.dihedralAngles(i-1)) * deltaY;
                
                % z-Position
                positionVectors(3,i) = positionVectors(3,i-1) + deltaY * sin(this.dihedralAngles(i-1));
            end
        end
        
        function setGlobalThicknessRatio(this, ratio)
            
            this.t_c_max_global = ratio;
        end
        
        function setKretractable(this, k)
            
            this.k_retractable = k;
        end
        
        function angle = getEquivalentSweepAngle(this)
            
            sumS = 0;
            sumSweep = 0;
            for i = 1:1:this.nPartitions
            
                rootChord = this.chordLengths(i);
                tipChord = this.chordLengths(i+1);
                spanPanel = this.getSpan() / 2 * (this.relSpanPositions(i+1) - this.relSpanPositions(i));
                SPartition = spanPanel * (rootChord + tipChord) / 2;
                sumS = sumS + SPartition;
                sumSweep = sumSweep + abs(this.sweepAngles(i)) * SPartition;
            end
             
            angle = sumSweep / sumS;
        end
        
        function angle = getEquivalentDihedralAngle(this)
            
            sumS = 0;
            sumDihedral = 0;
            for i = 1:1:this.nPartitions
            
                rootChord = this.chordLengths(i);
                tipChord = this.chordLengths(i+1);
                spanPanel = this.getSpan() / 2 * (this.relSpanPositions(i+1) - this.relSpanPositions(i));
                SPartition = spanPanel * (rootChord + tipChord) / 2;
                sumS = sumS + SPartition;
                sumDihedral = sumDihedral + abs(this.dihedralAngles(i)) * SPartition;
            end
             
            angle = sumDihedral / sumS;
        end
        
        function setSymmetry(this, newValue)
            
            if newValue == 0 || newValue == 1
                this.symmetry = newValue;
            else
                fprintf('Unknown value ''%i'' to define wing symmetry! Must be boolean (0 or 1).', newValue);
            end
        end
        
        function dataLine = getDataLine(this)
           
            airfoils = unique(this.airfoilsUsed(:));
            airfoilString = '';
            for i = 1:size(airfoils, 2)
                airfoilString = [airfoilString, this.airfoilDB(i).Name, ' '];
            end
            
            dataLine = {'Wing type', class(this);...
                        'Symmetry', this.symmetry;...
                        'Reference area [m²]', this.Sref;...
                        'Aspect ratio [-]', this.aspectRatio;...
                        '(Equivalent) Taper ratio [-]', this.taperRatio;...
                        'Span [m]', this.getSpan();...
                        'Root chord length [m]', this.getRootChord();...
                        'Mean sweep angle [°]', this.getEquivalentSweepAngle() /pi*180;...
                        'Mean dihedral angle [°]', this.getEquivalentDihedralAngle() /pi*180;...
                        'Used airfoils', airfoilString;...
                        'Wing mass [kg]', this.mass;...
                        'Wing retractable?', this.k_retractable};
                    
            if size(this.devices,2) > 0
                devicesString = '';
                for i = 1:size(this.devices, 2)
                    devicesString = [devicesString, this.devices(i).name, ' '];
                end
                dataLine = vertcat(dataLine, {'Devices', devicesString});
            end
        end
            
%         function setAirfoil(this, relSpanPositionBeginn, relSpanPositionEnd, airfoil)
%             OBSOLETE --> Use this.addAirfoilTransition instead
%             indices = this.findIndicesPanels(relSpanPositionBeginn, relSpanPositionEnd);
%             nAirfoils = size(this.airfoilDB, 2);
%             this.airfoilDB = [this.airfoilDB, airfoil];
%             for i = indices(1:end-1)
%                 this.airfoilsUsed(:,i) = (nAirfoils + 1) * [1;1];
%                 this.airfoilBlending(:,i) = [1;1];
%             end
%         end
        
%         function reprocessAeroData(this)
%            
%             aeroData = this.aeroProp.customData;
%             nPanelsX = this.aeroProp.customData.nPanelsX;
%             nPanelsY = this.aeroProp.customData.nPanelsY;
%             
%             for k = 1:length(aeroData.dataPoints)
%                 CDi = 0;
%                 Cn = 0;
%                 for j = 1:nPanelsY
%                     for i = 1:nPanelsX
%                         CXcorrected = min(1, max(0, aeroData.dataPoints(k).panelData(i,j).Cx));
%                         CDi = CDi + CXcorrected * aeroData.dataPoints(k).panelData(i,j).Sref;
%                         Cn = Cn - CXcorrected * aeroData.dataPoints(k).panelData(i,j).Sref * aeroData.dataPoints(k).panelData(i,j).PA(2);
%                     end
%                 end
%                 aeroData.dataPoints(k).CDi = CDi / this.aeroProp.customData.Sref;
%                 aeroData.dataPoints(k).Cn = Cn / (this.aeroProp.customData.Sref * this.aeroProp.customData.Lref_Cn);
%             end
%             this.aeroProp.setCustomData(aeroData);
%         end
%         
%         function aeroData = calculateFlapViscEffects(this, nStripes)
%            
%             if nargin == 1
%                 nStripes = 50;
%             end
%             
%             span = this.getSpan();
%             sumCD = 0;
%             sumCn = 0;
%             for j = 1:this.nPartitions        
%                 if (this.flapChordLengths(1,j) > 0 || this.flapChordLengths(2,j) > 0) %&& this.flapDeflections(j) > 0
%                     for i = 1:nStripes
%                         deltaRelY = (this.relSpanPositions(j+1) - this.relSpanPositions(j)) / nStripes;
%                         relYPos = this.relSpanPositions(j) + (i-1 + 0.5) * deltaRelY;
%                         relFlapChordLength = this.getFlapChordFraction(relYPos);
%                         CdFlap = Airfoil.getCdflapSectionDATCOM(relFlapChordLength, abs(this.flapDeflections(j)));
% 
%                         deltaCD = CdFlap * this.getChordInterpolated(relYPos) * deltaRelY * span / this.Sref;
%                         sumCD = sumCD + deltaCD;
%                         sumCn = sumCn + deltaCD * relYPos; % * span / span
%                     end
%                 end
%             end
%             aeroData.CDFlapVisc = sumCD;
%             aeroData.CnFlapVisc = sumCn;  
%         end

        function plotVLMLattice(this, nPanelsX, nPanelsY)
        
            geo = this.getVLMLattice(nPanelsX, nPanelsY);

            % Plot of Wing Planform Geometry
            figureA = 'Plot of Vortex Lattice Geometry';
            handleFigure = findobj('type', 'figure', 'Name', figureA);
            if isempty(handleFigure)
                handleFigure = figure('Name', figureA);
            end
            figure(handleFigure);
            clf

            % General plotting settings
            color_panels = 'black';
            color_vortex_filaments = 'b';
            color_collocPoints = 'r';
            lineWidth_LE_TE = 4;
            color_hinge = 'c';
            lineWidth_hinge = 2;

            for j = 1:size(geo.X_corners_inboard,2)
                
                plot3(geo.X_corners_inboard(:,j), geo.Y_corners_inboard(:,j), geo.Z_corners_inboard(:,j), 'Color', color_panels);
                hold on
                plot3(geo.X_corners_outboard(:,j), geo.Y_corners_outboard(:,j), geo.Z_corners_outboard(:,j), 'Color', color_panels);
                for i = 1:size(geo.X_corners_inboard,1)
                    plot3([geo.X_corners_inboard(i,j), geo.X_corners_outboard(i,j)], ...
                          [geo.Y_corners_inboard(i,j), geo.Y_corners_outboard(i,j)], ...
                          [geo.Z_corners_inboard(i,j), geo.Z_corners_outboard(i,j)], 'Color', color_panels);
                end
            end

            for i = 1:length(geo.X_A_vortex(:))

                plot3([geo.X_A_vortex(i), geo.X_B_vortex(i)], [geo.Y_A_vortex(i), geo.Y_B_vortex(i)], [geo.Z_A_vortex(i), geo.Z_B_vortex(i)], 'Color', color_vortex_filaments); % quarter line vortex
                plot3([geo.X_A_vortex(i), geo.X_D_vortex(i)], [geo.Y_A_vortex(i), geo.Y_D_vortex(i)], [geo.Z_A_vortex(i), geo.Z_D_vortex(i)], 'Color', color_vortex_filaments); % inboard quarter-to-hinge-line vortex
                plot3([geo.X_F_vortex(i), geo.X_D_vortex(i)], [geo.Y_F_vortex(i), geo.Y_D_vortex(i)], [geo.Z_F_vortex(i), geo.Z_D_vortex(i)], 'Color', color_vortex_filaments); % inboard hinge-line-to-trailing-edge vortex
                plot3([geo.X_E_vortex(i), geo.X_B_vortex(i)], [geo.Y_E_vortex(i), geo.Y_B_vortex(i)], [geo.Z_E_vortex(i), geo.Z_B_vortex(i)], 'Color', color_vortex_filaments); % outboard quarter-to-hinge-line vortex
                plot3([geo.X_E_vortex(i), geo.X_G_vortex(i)], [geo.Y_E_vortex(i), geo.Y_G_vortex(i)], [geo.Z_E_vortex(i), geo.Z_G_vortex(i)], 'Color', color_vortex_filaments); % outboard hinge-line-to-trailing-edge vortex
                plot3(geo.X_C_colloc(i), geo.Y_C_colloc(i), geo.Z_C_colloc(i), 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10, 'Color', color_collocPoints); % collocation points
            end

            axis equal
        end
        

        function [results, panelData] = calculateVLM(this, alphas, betas, nPanelsX, nPanelsY, rot_omega, rot_X0, rot_dir_axis)

            
            % Generate Vortex Lattice
            geo = this.getVLMLattice(nPanelsX, nPanelsY);

            % Reference values
            ref.b_ref = this.getSpan();
            ref.S_ref = this.Sref;
            ref.c_ref = this.Lref;
            ref.X_mom_ref = this.refPoint_mom(1);
            ref.Y_mom_ref = this.refPoint_mom(2);
            ref.Z_mom_ref = this.refPoint_mom(3);

            % Aerodynamic flow conditions (state)
            state.alpha = alphas;
            state.beta = betas; 

            if nargin > 5
                state.rot_omega = rot_omega;
                state.rot_X0 = rot_X0;
                state.rot_axis = rot_dir_axis;
            end
                
            % Solve VLM
            results = this.solveVLM7Segment(geo, state, ref);

            % Postprocessing for chordwise/spanwise data
            for j = 1:length(betas)
                for i = 1:length(alphas)
%                     
                    % TODO: extract panel data

%                     sinAlpha = sin(alphas(i));
%                     cosAlpha = cos(alphas(i));
%                     sinBeta = sin(betas(j));
%                     cosBeta = cos(betas(j));
% 
%                     t_1_1 = cosBeta * cosAlpha;
%                     t_2_1 = cosAlpha * sinBeta;
%                     t_3_1 = -sinAlpha;
%                     t_1_2 = -sinBeta;
%                     t_2_2 = cosBeta;
%                     t_3_2 = 0;
%                     t_1_3 = cosBeta * sinAlpha;
%                     t_2_3 = sinBeta * sinAlpha;
%                     t_3_3 = cosAlpha;
% 
%                     panelData(i,j).alpha_global_deg = alphas(i) /pi*180;
%                     panelData(i,j).beta_global_deg = betas(j)  /pi*180;
%                     panelData(i,j).y = geo.Y_C_colloc(1,:);
%                     panelData(i,j).c = [fliplr(c_sec_mid), c_sec_mid];
%                     CDi_sec = sum(t_1_1 * F_X + t_1_2 * F_Y + t_1_3 * F_Z) ./ qS_sec;
%                     CC_sec = sum(t_2_1 * F_X + t_2_2 * F_Y + t_2_3 * F_Z) ./ qS_sec;
%                     panelData(i,j).CL = sum(t_3_1 * F_X + t_3_2 * F_Y + t_3_3 * F_Z) ./ qS_sec;
%                     panelData(i,j).CC = CC_sec;
%                     panelData(i,j).CDi = CDi_sec;
%                     panelData(i,j).S_sections = S_sec;
%                     CDi(i,j) = results.c_D;
%                     CC(i,j) = results.c_Y;
%                     CL(i,j) = results.c_L;
%                     CD(i,j) = CDi(i,j);
%                     Cl(i,j) = results.c_l;
%                     Cm(i,j) = results.c_m;
%                     Cn(i,j) = results.c_n;
                end
            end
        end

        function CD0 = estimateZeroLiftDrag(this, Re_x)

%             S_partition = 0.5 * (this.chordLengths(1:end-1) + this.chordLengths(2:end)) .* (this.relSpanPositions(2:end) - this.relSpanPositions(1:end-1)) * 0.5 * this.getSpan();
            Re_section = Re_x * [this.chordLengths(1:end-1); this.chordLengths(2:end)];
            C_f_section = 0.074 ./ (Re_section).^0.2;% - 1742 ./ Re_section;
            if ~isempty(this.airfoilDB)
                L_section = (reshape([this.airfoilDB(this.airfoilsUsed).posMaxThickness], [2, this.nPartitions]) < 0.3) * 2 + (reshape([this.airfoilDB(this.airfoilsUsed).posMaxThickness], [2, this.nPartitions]) >= 0.3) * 1.2;
                FF = 1 + L_section .* reshape([this.airfoilDB(this.airfoilsUsed).maxThickness], [2, this.nPartitions]) + 100 * reshape([this.airfoilDB(this.airfoilsUsed).maxThickness], [2, this.nPartitions]).^4;
            else
                L_section = 2;
                FF = 1 + L_section * 0.1 + 100 * 0.1^4;
            end
%             CD0 = 2 * C_f_section .* FF .* S_partition / sum(S_partition);
            CD0 = sum(sum(C_f_section .* FF) .* (this.relSpanPositions(2:end) - this.relSpanPositions(1:end-1))) * (this.getSpan() / this.Sref);

        end
            
    end % methods public
    
    methods (Access = protected)
        
        function this = Wing()
            % "Constructor" for abstract class 'Wing' which should be called by every inherited 
            % class to initialise certain properties.
            % 
            % Called functions:
            % 	AerodynamicGeometry() of class AerodynamicGeometry (Constructor of abstract class)
            % 	AerodynamicProperties() of class AerodynamicProperties (Constructor)
            	
            this = this@AerodynamicGeometry();
            this.symmetry = 1;
            this.aeroProp = AerodynamicProperties();
            this.dihedralKinks = struct('relYPos', {}, 'angle', {}, 'affectedPanels', {});
            this.sweepKinks = struct('relYPos', {}, 'angle', {}, 'affectedPanels', {});
            this.twistTransitions = struct('relYPos1', {}, 'relYPos2', {}, 'angle', {}, 'transitionType', {});
            this.airfoilZones = struct('relYPos1', {}, 'relYPos2', {}, 'airfoilID1', {}, 'airfoilID2', {}, 'zoneType', {});
            this.flapBreaks = struct('relYPos1', {}, 'relYPos2', {}, 'fraction', {}, 'flapID', {}, 'flapControlMode', {});
            this.refPoint = [0; 0; 0];
            this.type = 'cubic';
            this.t_c_max_global = 0.10;
            this.k_retractable = 1.0;
            this.k_FRP = 0.20;
            this.roughness_surface = 0.25 * 10^-3;
            this.refreshMCGList();
            
        end

        function calculateRelSpanPositions(this)
            % Accumulates the relative span positions of the section points where the wing is 
            % modelised taking into account all kinks stored in the properties.
            
            % At first create a set of all kink positions
            yAllKinks = [this.dihedralKinks.relYPos, ...
                         this.sweepKinks.relYPos, ...
                         this.twistTransitions.relYPos1, ...
                         this.twistTransitions.relYPos2, ...
                         this.airfoilZones.relYPos1, ...
                         this.airfoilZones.relYPos2, ...
                         this.flapBreaks.relYPos1, ...
                         this.flapBreaks.relYPos2];
                     
            if ~isempty(this.chordKinks)
                yAllKinks = [this.chordKinks.relYPos, yAllKinks];
            end

            % Delete all double positions (if there are discontinuities at the same position)
            yPos = unique([0 sort(yAllKinks) 1]);
            this.relSpanPositions = yPos;
            this.nPartitions = size(yPos, 2) - 1;
        end

        function relSpanPositions = getPanelPositions(nPanel)
            % Calculates the relative span positions of the section points where the wing is 
            % modelised taking into account all kinks and the number of input panels.
            % The method tries to refine the partition sizes towards the wing tip by using 
            % the sinus function. If the global number of partitions is too small according to the 
            % number of kinks an error may occur.

            % Initialisation
            nPglobal = nPanel;
            
            % At first create a set of all kink positions
            yAllKinks = [this.dihedralKinks.relYPos, ...
                         this.sweepKinks.relYPos, ...
                         this.twistTransitions.relYPos1, ...
                         this.twistTransitions.relYPos2, ...
                         this.airfoilZones.relYPos1, ...
                         this.airfoilZones.relYPos2, ...
                         this.flapBreaks.relYPos1, ...
                         this.flapBreaks.relYPos2];
                     
            if ~isempty(this.chordKinks)
                yAllKinks = [this.chordKinks.relYPos, yAllKinks];
            end
                     
            % Delete all double positions (if there are discontinuities at the same position)
            yAllKinks = unique(yAllKinks);
            nK = size(yAllKinks,2) * ~isempty(yAllKinks);
            
            % if there are kinks to be considered
            if (nK ~= 0)
                
                yAllKinks = sort(yAllKinks);

					 % The kinks divide the halfspan into parts which a certain number of partitions is assigned
                nPlocal = zeros(1,nK);			% Vector containing the number of partitions for each part of the wing
                
                % Initialise the position of section points
                yPosKinks = [0 yAllKinks 1];
                % Delete all double positions (if there is a kink at the root chord)
                yPosKinks = unique(yPosKinks);
                yPos = zeros(1,nPglobal+1);
                indices = ones(size(yPosKinks));
                for i = 2:length(yPosKinks)
                    nPcurr = sum(nPlocal);
%                     if length(yPosKinks) == 2 && all(yPosKinks == [0,1])
% %                         yPos(nPcurr+2:end) = yPosKinks(i-1) + (1 - yPosKinks(i-1)) * 0.5 * (sin(-0.5*pi + (1:nPglobal-nPcurr) / (nPglobal-nPcurr) * pi) + 1);
%                         yPos(nPcurr+2:end) = yPosKinks(i-1) + (1 - yPosKinks(i-1)) * sin((1:nPglobal-nPcurr) / (nPglobal-nPcurr) * 0.5*pi);
%                     else
%                         yPos(nPcurr+2:end) = yPosKinks(i-1) + (1 - yPosKinks(i-1)) * sin((1:nPglobal-nPcurr) / (nPglobal-nPcurr) * 0.5*pi);
%                     end
                    yPos(nPcurr+2:end) = yPosKinks(i-1) + (1 - yPosKinks(i-1)) * sin((1:nPglobal-nPcurr) / (nPglobal-nPcurr) * 0.5*pi);

                    indices(i) = find(yPos >= yPosKinks(i), 1, 'first');
                    %if abs(yPosKinks(i) - yPos(indices(i))) < abs(yPosKinks(i) - yPos(indices(i)-1))
                        yPos(indices(i)) = yPosKinks(i);
                        nPlocal(i-1) = indices(i) - indices(i-1);
%                     else
%                         yPos(indices(i)-1) = yPosKinks(i);
%                         nPlocal(i-1) = indices(i)-1 - indices(i-1);
%                     end
                end
                %{
					 % For each initial part assign the number of partitions according to the part's ratio of occupied
					 % length to total length
                for i = 2:1:size(yPos,2)

                    nPlocal(i-1) = round(nPglobal * (yPos(i) - yPos(i-1)));
                    % one kink partition must consist of at least one partition
                    if (nPlocal(i-1) < 1) 
                        nPlocal(i-1) = 1;
                    end
                end
					 
					 % If there is still partition available to be assigned due to the round operation 
					 % assign it to the partitions beginning at the wing tip
                while (sum(nPlocal) < nPglobal)
                    nPlocal(size(nPlocal,2)) = nPlocal(size(nPlocal,2)) + 1;
                end

					 % If there are more partitions assigned than allowed reduce the number of partitions
					 % beginning at the wing root
                i = 1;
                while (sum(nPlocal) > nPglobal)
                    if (i <= nPglobal)
                        if (nPlocal(i) > 1)
                            nPlocal(i) = nPlocal(i) - 1;
                        else
                            i = i+1;
                        end
                    else
                        error('Number of wing partitions too small!');
                    end
                end
                %}
					 % Calculate the final relative y positions of each section point
                
                for i = 1:1:(size(yPosKinks,2)-1)
                    for j = 1:1:nPlocal(i)-1

                        if (i == size(yPosKinks,2)-1)
%                             if length(yPosKinks) == 2 && all(yPosKinks == [0,1])
% %                                 yPos(sum(nPlocal(1:i-1))+j+1) = yPosKinks(i) + (yPosKinks(i+1) - yPosKinks(i)) * (0.5 * (sin(-0.5*pi + j * pi() / nPlocal(i)) + 1));
%                                 yPos(sum(nPlocal(1:i-1))+j+1) = yPosKinks(i) + (yPosKinks(i+1) - yPosKinks(i)) * (sin(j * (pi()/2) / nPlocal(i)));
%                             else
%                                 % Use of the sinus function to refine the partitions towards the wing tip
%                                 yPos(sum(nPlocal(1:i-1))+j+1) = yPosKinks(i) + (yPosKinks(i+1) - yPosKinks(i)) * (sin(j * (pi()/2) / nPlocal(i)));
%                             end
                            yPos(sum(nPlocal(1:i-1))+j+1) = yPosKinks(i) + (yPosKinks(i+1) - yPosKinks(i)) * (sin(j * (pi()/2) / nPlocal(i)));
                        else
                            yPos(sum(nPlocal(1:i-1))+j+1) = yPosKinks(i) + (yPosKinks(i+1) - yPosKinks(i)) * (j / nPlocal(i));
                        end
                    end
                end
            
            % if there are no kinks assigned
            else
            	 % Calculate the final relative y positions of each section point
                yPos = zeros(1,nPglobal+1);
                
                for i = 1:1:nPglobal
		                % Use of the sinus function to refine the partitions towards the wing tip
		                yPos(i+1) = sin(i * (pi()/2) / nPglobal);
                end
            end

            relSpanPositions = yPos;
        end
        
        function setDihedralAngles(this)
            % Refreshes the dihedral angles stored in the property dihedralAngles according to the 
            % kink information stored in the property dihedralKinks.

            this.dihedralAngles = zeros(1, this.nPartitions);
            
            for i = 1:1:size(this.dihedralKinks,2)
                for j = 1:1:(this.nPartitions+1)
                
                    if (this.relSpanPositions(j) > this.dihedralKinks(i).relYPos)
                        this.dihedralAngles(j-1) = this.dihedralKinks(i).angle;
                    end
                end
            end
        end
        
        function setSweepAngles(this)
           % Refreshes the sweep angles stored in the property sweepAngles according to the 
           % kink information stored in the property sweepKinks.
           
           this.sweepAngles = zeros(1, this.nPartitions);

           for i = 1:1:size(this.sweepKinks,2)
                for j = 1:1:(this.nPartitions+1)
                
                    if (this.relSpanPositions(j) > this.sweepKinks(i).relYPos)
                        this.sweepAngles(j-1) = this.sweepKinks(i).angle;
                    end
                end
            end
        end
        
        function setAirfoilIDsAndBlending(this)
            
            if ~isempty(this.airfoilZones)
                [airfoilsUsed, airfoilBlending] = this.getAirfoilData(this.relSpanPositions);
    
                this.airfoilsUsed = airfoilsUsed;
                this.airfoilBlending = airfoilBlending;
            else
                this.setGlobalAirfoilID(1);
            end
        end

        function [airfoilsUsed, airfoilBlending] = getAirfoilData(this, relSpanPositions)

            airfoilsUsed = ones(2,length(relSpanPositions)-1);
            airfoilBlending = ones(2,length(relSpanPositions)-1);

            if ~isempty(this.airfoilZones)
                for j = 1:length(this.airfoilZones)

                    index_start = find(relSpanPositions >= this.airfoilZones(j).relYPos1, 1, 'first');

                    if this.airfoilZones(j).relYPos2 < 1
                        index_end = find(relSpanPositions >= this.airfoilZones(j).relYPos2, 1, 'first') - 1;
                    else
                        index_end = size(relSpanPositions, 2) - 1;
                    end
                    
                    switch this.airfoilZones.zoneType
                        case 'uniform'

                            airfoilsUsed(:,index_start:index_end) = this.airfoilZones(j).airfoilID1;
                            airfoilBlending(:,index_start:index_end) = 1;
                            
                        case 'linear'
                            
                            airfoilsUsed(1,index_start:index_end) = this.airfoilZones(j).airfoilID1;
                            airfoilsUsed(2,index_start:index_end) = this.airfoilZones(j).airfoilID2;
                            fractions = 1 - (relSpanPositions(index_start:index_end+1) - relSpanPositions(index_start)) / (relSpanPositions(index_end+1) - relSpanPositions(index_start));
                            airfoilBlending(1,index_start:index_end) = fractions(1:end-1); % inboard weight value of airfoil 1
                            airfoilBlending(2,index_start:index_end) = fractions(2:end); % outboard weight value of airfoil 1
                        otherwise
                            error('Unknown transition type %s. Known types are "uniform" and "linear".', this.airfoilZone.zoneType);
                    end
                end
            end
        end
        
        function setTwistAngles(this)
           % Refreshes the twist angles according to the 
           % information stored in the property twistTransitions.

           [twistAngles, twistRotPoints] = this.getTwistData(this.relSpanPositions);
            this.twistAngles = twistAngles;
            this.twistRotPoints = twistRotPoints;

%            this.twistAngles = zeros(2, this.nPartitions);
%            this.twistRotPoints = ones(2,this.nPartitions) * 0.25;
           
%            if ~isempty(this.twistTransitions)
%            startPoints = {[], this.twistTransitions(:).relYPos1};
%            transitionTypes = {'root', this.twistTransitions(:).transitionType};
%            endPoints = [0, this.twistTransitions(:).relYPos2];
%            angles = [0, this.twistTransitions(:).angle];
%            rotPoints = [0.25, this.twistTransitions(:).rotPoint];
%            
%            if endPoints(end) ~= 1
%                endPoints = [endPoints, 1];
%                startPoints = [startPoints, {[]}];
%                transitionTypes = [transitionTypes, {'singular'}];
%                angles = [angles, this.twistTransitions(end).angle];
%                rotPoints = [rotPoints, this.twistTransitions(end).rotPoint];
%            else
%                transitionTypes = {transitionTypes{:}, 'singular'};
%                endPoints = [endPoints, 1];
%                angles = [angles, angles(end)];
%                rotPoints = [rotPoints, rotPoints(end)];
%            end
%            
%            indexSingularPoints = find(ismember(transitionTypes, 'singular') == 1);
%            
%            for j = 1:length(indexSingularPoints)
%            
%                if j == 1
%                     coordY = endPoints(1:indexSingularPoints(j));
%                     coordTw = [angles(1:indexSingularPoints(j)-1), angles(indexSingularPoints(j)-1)];
%                     coordRotP = [rotPoints(1:indexSingularPoints(j)-1), rotPoints(indexSingularPoints(j)-1)];
%                     indexBegin = 1;
%                else
%                     coordY = endPoints(indexSingularPoints(max(1,j-1)):indexSingularPoints(j));
%                     coordTw = [angles(indexSingularPoints(max(1,j-1)):indexSingularPoints(j)-1), angles(indexSingularPoints(j)-1)];
%                     coordRotP = [rotPoints(indexSingularPoints(max(1,j-1)):indexSingularPoints(j)-1), rotPoints(indexSingularPoints(j)-1)];
%                     indexBegin = indexSingularPoints(j-1);
%                end
%            
%            for i = 1:(indexSingularPoints(j) - indexBegin)
% 
%                if ~isempty(startPoints{indexBegin-1+i}) && startPoints{indexBegin-1+i} > endPoints(indexBegin-1+i-1)
%                    if strcmp(transitionTypes{indexBegin-1+i}, 'custom')
%                        
%                        shapeMatrix = this.twistTransitions(indexBegin-1+i-1).shapeMatrix;
%                        shapeY = shapeMatrix(1,:) - shapeMatrix(1,1);
%                        shapeY = startPoints{indexBegin-1+i} + (endPoints(indexBegin-1+i) - startPoints{indexBegin-1+i}) * shapeY / max(shapeY);
%                        shapeTw = shapeMatrix(2,:) - shapeMatrix(2,1);
%                        shapeTw = angles(indexBegin-1+i-1) + angles(indexBegin-1+i) * shapeTw / shapeTw(end);
%                        
%                        coordY = [coordY, shapeY(1:end-1)];
%                        coordTw = [coordTw, shapeTw(1:end-1)];
%                    else
%                        coordY = [coordY, startPoints{indexBegin-1+i}];
%                        coordTw = [coordTw, angles(indexBegin-1+i-1)];
%                        coordRotP = [coordRotP, rotPoints(indexBegin-1+i-1)];
%                    end
%                end
%            end
%            
%            % Sort unique points along span
%             %[~, order] = sort(coordY);
%             [coordY,order,~] = unique(coordY, 'first');
%             %coordY = coordY(order);
%             coordTw = coordTw(order);
%             coordRotP = coordRotP(order);
%             
%            % Set twist angles and rotation Points for every section
%            for i = 1:this.nPartitions
%                if this.relSpanPositions(i) >= coordY(1) && ...
%                   this.relSpanPositions(i) < coordY(end)
% 
%                     this.twistAngles(:,i) = [interp1(coordY, coordTw, this.relSpanPositions(i)); ...
%                                         interp1(coordY, coordTw, this.relSpanPositions(i+1))];
%                     this.twistRotPoints(:,i) = [interp1(coordY, coordRotP, this.relSpanPositions(i)); ...
%                                         interp1(coordY, coordRotP, this.relSpanPositions(i+1))];                
%                end
%            end
%            end
%            end
           
        end

        function [twistAngles, twistRotPoints] = getTwistData(this, relSpanPositions)

           twistAngles = zeros(2, length(relSpanPositions) - 1);
           twistRotPoints = ones(2,length(relSpanPositions) - 1) * 0.25;
           
           if ~isempty(this.twistTransitions)
           startPoints = {[], this.twistTransitions(:).relYPos1};
           transitionTypes = {'root', this.twistTransitions(:).transitionType};
           endPoints = [0, this.twistTransitions(:).relYPos2];
           angles = [0, this.twistTransitions(:).angle];
           rotPoints = [0.25, this.twistTransitions(:).rotPoint];
           
           if endPoints(end) ~= 1
               endPoints = [endPoints, 1];
               startPoints = [startPoints, {[]}];
               transitionTypes = [transitionTypes, {'singular'}];
               angles = [angles, this.twistTransitions(end).angle];
               rotPoints = [rotPoints, this.twistTransitions(end).rotPoint];
           else
               transitionTypes = {transitionTypes{:}, 'singular'};
               endPoints = [endPoints, 1];
               angles = [angles, angles(end)];
               rotPoints = [rotPoints, rotPoints(end)];
           end
           
           indexSingularPoints = find(ismember(transitionTypes, 'singular') == 1);
           
           for j = 1:length(indexSingularPoints)
           
               if j == 1
                    coordY = endPoints(1:indexSingularPoints(j));
                    coordTw = [angles(1:indexSingularPoints(j)-1), angles(indexSingularPoints(j)-1)];
                    coordRotP = [rotPoints(1:indexSingularPoints(j)-1), rotPoints(indexSingularPoints(j)-1)];
                    indexBegin = 1;
               else
                    coordY = endPoints(indexSingularPoints(max(1,j-1)):indexSingularPoints(j));
                    coordTw = [angles(indexSingularPoints(max(1,j-1)):indexSingularPoints(j)-1), angles(indexSingularPoints(j)-1)];
                    coordRotP = [rotPoints(indexSingularPoints(max(1,j-1)):indexSingularPoints(j)-1), rotPoints(indexSingularPoints(j)-1)];
                    indexBegin = indexSingularPoints(j-1);
               end
           
           for i = 1:(indexSingularPoints(j) - indexBegin)

               if ~isempty(startPoints{indexBegin-1+i}) && startPoints{indexBegin-1+i} > endPoints(indexBegin-1+i-1)
                   if strcmp(transitionTypes{indexBegin-1+i}, 'custom')
                       
                       shapeMatrix = this.twistTransitions(indexBegin-1+i-1).shapeMatrix;
                       shapeY = shapeMatrix(1,:) - shapeMatrix(1,1);
                       shapeY = startPoints{indexBegin-1+i} + (endPoints(indexBegin-1+i) - startPoints{indexBegin-1+i}) * shapeY / max(shapeY);
                       shapeTw = shapeMatrix(2,:) - shapeMatrix(2,1);
                       shapeTw = angles(indexBegin-1+i-1) + angles(indexBegin-1+i) * shapeTw / shapeTw(end);
                       
                       coordY = [coordY, shapeY(1:end-1)];
                       coordTw = [coordTw, shapeTw(1:end-1)];
                   else
                       coordY = [coordY, startPoints{indexBegin-1+i}];
                       coordTw = [coordTw, angles(indexBegin-1+i-1)];
                       coordRotP = [coordRotP, rotPoints(indexBegin-1+i-1)];
                   end
               end
           end
           
           % Sort unique points along span
            %[~, order] = sort(coordY);
            [coordY,order,~] = unique(coordY, 'first');
            %coordY = coordY(order);
            coordTw = coordTw(order);
            coordRotP = coordRotP(order);
            
           % Set twist angles and rotation Points for every section
           for i = 1:(length(relSpanPositions) - 1)
               if relSpanPositions(i) >= coordY(1) && ...
                  relSpanPositions(i) < coordY(end)

                    twistAngles(:,i) = [interp1(coordY, coordTw, relSpanPositions(i)); ...
                                        interp1(coordY, coordTw, relSpanPositions(i+1))];
                    twistRotPoints(:,i) = [interp1(coordY, coordRotP, relSpanPositions(i)); ...
                                        interp1(coordY, coordRotP, relSpanPositions(i+1))];                
               end
           end
           end
           end
        end
        
        function setFlapChordLengths(this)
            
            [flapChordLengths, flapControlIDs, deflectionModes] = this.getFlapData(this.relSpanPositions);

            this.flapControlIDs = flapControlIDs;
            this.deflectionModes = deflectionModes;
            this.flapChordLengths = flapChordLengths;
        end

        function [flapChordLengths, flapControlIDs, deflectionModes] = getFlapData(this, relSpanPositions)

            chordLengths = this.getChordInterpolated(relSpanPositions);
            flapControlIDs = nan(1,length(relSpanPositions)-1);
            deflectionModes = nan(1,length(relSpanPositions)-1);
            flapChordLengths = zeros(2,length(relSpanPositions)-1);

            for j = 1:size(this.flapBreaks,2)
            
                relSpanPos = [this.flapBreaks(j).relYPos1, this.flapBreaks(j).relYPos2];
                
                index_start = find(relSpanPositions >= relSpanPos(1), 1, 'first');
                
                if relSpanPos(2) < 1
                    index_end = find(relSpanPositions >= relSpanPos(2), 1, 'first') - 1;
                else
                    index_end = size(relSpanPositions, 2) - 1;
                end

                for i = index_start:index_end
%                     this.flapChordLengths(1,i) = interp1(relSpanPos, this.flapBreaks(j).fractions .* this.getChordInterpolated(relSpanPos), this.relSpanPositions(i)) / this.chordLengths(i);
%                     this.flapChordLengths(2,i) = interp1(relSpanPos, this.flapBreaks(j).fractions .* this.chordLengths(i+1), this.relSpanPositions(i+1)) / this.chordLengths(i+1);
                    flapChordLengths(:,i) = (roundn(interp1(relSpanPos, this.flapBreaks(j).fractions .* this.getChordInterpolated(relSpanPos), relSpanPositions(i:i+1)) ./ chordLengths(i:i+1), -6))';
                    flapControlIDs(i) = this.flapBreaks(j).flapID;
                    deflectionModes(i) = this.flapBreaks(j).flapControlMode;
%                     this.hingeLines(:,i) = this.flapBreaks(j).hingeLinePoints';
                end
            end
        end
        
        function bool = isFlapBreak(this, YPos, option)
            
            if nargin == 2
                option = 'rel';
            end
            
            positions = [this.flapBreaks.relYPos1, this.flapBreaks.relYPos2];
            
            if strcmp(option, 'rel')
                bool = sum(positions == YPos) > 0;
            elseif strcmp(option, 'abs')
                bool = sum((positions * 0.5*this.getSpan()) == YPos) > 0;
            else
                error('Unknown option flag');
            end
        end
        
        function bool = isTwistBreak(this, YPos, option)
            
            if nargin == 2
                option = 'rel';
            end
            
            positions = [this.twistTransitions.relYPos2];
            types = {this.twistTransitions.transitionType};
            
            if strcmp(option, 'rel')
                found = find(positions == YPos);
            elseif strcmp(option, 'abs')
                found = find((positions * 0.5*this.getSpan()) == YPos);
            else
                error('Unknown option flag');
            end
            
            if (sum(found) > 0)
                bool = sum(strcmp('singular', types(found))) > 0;
            else
                bool = 0;
            end
            
        end
        
        function initialiseDihedralAngles(this)
            % Initialises the property dihedralAngles and sets the dihedral angle to zero for every 
            % partition
            
            this.dihedralAngles = zeros(1,this.nPartitions);
        end
        
        function initialiseSweepAngles(this)
            % Initialises the property sweepAngles and sets the sweep angle to zero for every 
            % partition
            
            this.sweepAngles = zeros(1,this.nPartitions);
        end
        
        function initialiseTwistAngles(this)
            % Initialises the property twistAngles and sets the twist angle to zero for every 
            % partition
            
            this.twistAngles = zeros(2,this.nPartitions);
            this.twistRotPoints = ones(2,this.nPartitions) * 0.25;
        end
            
            
        function initialiseAirfoils(this)
            % Initialises the airfoil used for each partition where the first Airfoil in database 
            % is used as default airfoil for all partitions
            
            this.airfoilsUsed = ones(2,this.nPartitions);
            this.airfoilBlending = ones(2,this.nPartitions);
        end
            
        function initialiseFlapChordLengths(this)
            % Initialises the property flapChordLengths and sets the chord lengths to zero for 
            % every partition
            
            this.flapChordLengths = zeros(2,this.nPartitions);
        end
        
        function initialiseFlapDeflections(this)
            % Initialises the property flapDeflections and sets the flap deflection angle to zero 
            % for every partition
            
            this.flapDeflections = zeros(1,this.nPartitions);
        end
        
        function initialiseDeflectionModes(this)
            % Initialises the property deflectionModes and sets the mode to one i.e. symetric 
            % deflection for every partition
            
            this.deflectionModes = ones(1,this.nPartitions);
        end
        
        function initialiseFlapControlID(this)
            % Initialises the property flapControlIDs and sets the ID to one
            % for every partition
            
            this.flapControlIDs = ones(1,this.nPartitions);
        end
        
        function indices = findIndicesPanels(this, relYPosBeginn, relYPosEnd)
           
            if nargin == 2
                indices = find(this.relSpanPositions <= relYPosBeginn, 1, 'last');
            elseif nargin == 3
                iStart = find(this.relSpanPositions <= relYPosBeginn, 1, 'last');
                if nargin < 3
                    indices = iStart;
                else
                    iEnd = find(this.relSpanPositions <= relYPosEnd, 1, 'last');
                    indices = iStart:1:iEnd;
                end
            end
        end
       
        function geo = getVLMLattice(this, nPanelsX, nPanelsY)

            % Distribute number of spanwise panels
            if nPanelsY < this.nPartitions
                
                warning('The input value for the number of spanwise panels is lower than the minimum required number.');
                nPanelsY = this.nPartitions;
                relSpanPositions = this.relSpanPositions;
            else
                relSpanPositions = this.getPanelDistribution(this.relSpanPositions, nPanelsY); % TODO: think about non-linear distribution towards wing tip
            end

            % Initialise arrays
            fractionsX_inboard = zeros(nPanelsX+1, nPanelsY);
            fractionsX_outboard = zeros(nPanelsX+1, nPanelsY);
            x_LE_inboard = zeros(1, nPanelsY);
            z_LE_inboard = zeros(1, nPanelsY);
            x_LE_outboard = zeros(1, nPanelsY);
            z_LE_outboard = zeros(1, nPanelsY);
            x_TE_inboard = zeros(1, nPanelsY);
            z_TE_inboard = zeros(1, nPanelsY);
            x_TE_outboard = zeros(1, nPanelsY);
            z_TE_outboard = zeros(1, nPanelsY);
            geo.X_corners_inboard = zeros(nPanelsX+1, 2 * nPanelsY);
            geo.X_corners_outboard = zeros(nPanelsX+1, 2 * nPanelsY);
            geo.Y_corners_inboard = zeros(nPanelsX+1, 2 * nPanelsY);
            geo.Y_corners_outboard = zeros(nPanelsX+1, 2 * nPanelsY);
            geo.Z_corners_inboard = zeros(nPanelsX+1, 2 * nPanelsY);
            geo.Z_corners_outboard = zeros(nPanelsX+1, 2 * nPanelsY);

            for s = {'A_vortex', 'B_vortex', 'D_vortex', 'E_vortex', 'F_vortex', 'G_vortex', 'C_colloc', 'n'}
                geo.(['X_',s{:}]) = zeros(nPanelsX, 2 * nPanelsY);
                geo.(['Y_',s{:}]) = zeros(nPanelsX, 2 * nPanelsY);
                geo.(['Z_',s{:}]) = zeros(nPanelsX, 2 * nPanelsY);
            end

            % Get the interpolated spanwise panel data
            positionVectors = interp1(this.relSpanPositions, this.getPositionVectors().', relSpanPositions).';
            chordLengths = this.getChordInterpolated(relSpanPositions);
            [flapChordLengths, flapControlIDs, deflectionModes] = this.getFlapData(relSpanPositions);
            flapDeflectionAngles = interp1(this.relSpanPositions, [0,this.flapDeflections], relSpanPositions(2:end), "next");
            [twistAngles, twistRotPoints] = this.getTwistData(relSpanPositions);
            [airfoilsUsed, airfoilBlending] = this.getAirfoilData(relSpanPositions);

            % Wing installation incidence angle
            twistAngles = twistAngles + this.incidenceAngle;

            % Prepare airfoil camber data
            if ~isempty(this.airfoilDB)
                usedAirfoilIDs = unique(this.airfoilsUsed(:));
                for k = 1:length(usedAirfoilIDs)
                    dx_camber = diff(this.airfoilDB(usedAirfoilIDs(k)).coordinates.camberLinePointArray(:,1));
                    dy_camber = diff(this.airfoilDB(usedAirfoilIDs(k)).coordinates.camberLinePointArray(:,2));
                    camberData(k).x_theta = [(this.airfoilDB(usedAirfoilIDs(k)).coordinates.camberLinePointArray(1:end-1,1) + 0.5 * dx_camber) / max(this.airfoilDB(usedAirfoilIDs(k)).coordinates.camberLinePointArray(:,1)); 1];
                    camberData(k).theta = atan2(dy_camber, dx_camber);
                    camberData(k).theta = [camberData(k).theta; camberData(k).theta(end)];
                end
            end

            for j = 1:nPanelsY

                % Determine chordwise panel fractions respecting flaps
                if any(flapChordLengths(:,j) > 0)

                    fractionsX_inboard(:,j) = this.getPanelDistribution([0, 1 - flapChordLengths(1,j), 1], nPanelsX); % TODO: think about non-linear distribution
                    i_hinge(j) = find(fractionsX_inboard(:,j) >= (1 - flapChordLengths(1,j)), 1, "first");
                    fractionsX_outboard(1:i_hinge(j),j) = linspace(0, 1 - flapChordLengths(2,j), i_hinge(j))';
                    fractionsX_outboard(i_hinge(j):end,j) = linspace(1 - flapChordLengths(2,j), 1, nPanelsX - i_hinge(j) + 2)'; % linear distribution
                else
                    i_hinge(j) = nPanelsX+1;
                    fractionsX_inboard(:,j) = linspace(0, 1, nPanelsX+1);
                    fractionsX_outboard(:,j) = fractionsX_inboard(:,j);
                end

                x_LE_inboard(j) = positionVectors(1,j) + (twistRotPoints(1,j) - 0.25) * chordLengths(j) - cos(twistAngles(1,j)) * twistRotPoints(1,j) * chordLengths(j);
                x_TE_inboard(j) = x_LE_inboard(j) + cos(twistAngles(1,j)) * chordLengths(j);
                z_LE_inboard(j) = positionVectors(3,j) + sin(twistAngles(1,j)) * twistRotPoints(1,j) * chordLengths(j);
                z_TE_inboard(j) = z_LE_inboard(j) - sin(twistAngles(1,j)) * chordLengths(j);
                x_LE_outboard(j) = positionVectors(1,j+1) + (twistRotPoints(2,j) - 0.25) * chordLengths(j+1) - cos(twistAngles(2,j)) * twistRotPoints(2,j) * chordLengths(j+1);
                x_TE_outboard(j) = x_LE_outboard(j) + cos(twistAngles(2,j)) * chordLengths(j+1);
                z_LE_outboard(j) = positionVectors(3,j+1) + sin(twistAngles(2,j)) * twistRotPoints(2,j) * chordLengths(j+1);
                z_TE_outboard(j) = z_LE_outboard(j) - sin(twistAngles(2,j)) * chordLengths(j+1);
    
                % Panel edge points, vortex filament edge points
                geo.X_corners_inboard(:,j) = x_LE_inboard(j) + fractionsX_inboard(:,j) * (x_TE_inboard(j) - x_LE_inboard(j));
                geo.X_corners_outboard(:,j) = x_LE_outboard(j) + fractionsX_outboard(:,j) * (x_TE_outboard(j) - x_LE_outboard(j));
                geo.X_A_vortex(:,j) = geo.X_corners_inboard(1:end-1,j) + 0.25 * diff(geo.X_corners_inboard(:,j));
                geo.X_B_vortex(:,j) = geo.X_corners_outboard(1:end-1,j) + 0.25 * diff(geo.X_corners_outboard(:,j));
                geo.X_F_vortex(:,j) = repmat(x_TE_inboard(j), nPanelsX, 1);
                geo.X_D_vortex(:,j) = [repmat(geo.X_corners_inboard(i_hinge(j),j), i_hinge(j)-1, 1); geo.X_F_vortex(i_hinge(j):end,j)];
                geo.X_G_vortex(:,j) = repmat(x_TE_outboard(j), nPanelsX, 1);
                geo.X_E_vortex(:,j) = [repmat(geo.X_corners_outboard(i_hinge(j),j), i_hinge(j)-1, 1); geo.X_G_vortex(i_hinge(j):end,j)];
                geo.Y_corners_inboard(:,j) = repmat(positionVectors(2,j), nPanelsX+1, 1);
                geo.Y_corners_outboard(:,j) = repmat(positionVectors(2,j+1), nPanelsX+1, 1);
                geo.Y_A_vortex(:,j) = geo.Y_corners_inboard(1:end-1,j);
                geo.Y_B_vortex(:,j) = geo.Y_corners_outboard(1:end-1,j);
                geo.Y_F_vortex(:,j) = geo.Y_A_vortex(:,j);
                geo.Y_D_vortex(:,j) = geo.Y_A_vortex(:,j);
                geo.Y_G_vortex(:,j) = geo.Y_B_vortex(:,j);
                geo.Y_E_vortex(:,j) = geo.Y_B_vortex(:,j);
                geo.Z_corners_inboard(:,j) = z_LE_inboard(j) + fractionsX_inboard(:,j) * (z_TE_inboard(j) - z_LE_inboard(j));
                geo.Z_corners_outboard(:,j) = z_LE_outboard(j) + fractionsX_outboard(:,j) * (z_TE_outboard(j) - z_LE_outboard(j));
                geo.Z_A_vortex(:,j) = geo.Z_corners_inboard(1:end-1,j) + 0.25 * diff(geo.Z_corners_inboard(:,j));
                geo.Z_B_vortex(:,j) = geo.Z_corners_outboard(1:end-1,j) + 0.25 * diff(geo.Z_corners_outboard(:,j));
                geo.Z_F_vortex(:,j) = repmat(z_TE_inboard(j), nPanelsX, 1);
                geo.Z_D_vortex(:,j) = [repmat(geo.Z_corners_inboard(i_hinge(j),j), i_hinge(j)-1, 1); geo.Z_F_vortex(i_hinge(j):end,j)];
                geo.Z_G_vortex(:,j) = repmat(z_TE_outboard(j), nPanelsX, 1);
                geo.Z_E_vortex(:,j) = [repmat(geo.Z_corners_outboard(i_hinge(j),j), i_hinge(j)-1, 1); geo.Z_G_vortex(i_hinge(j):end,j)];

                % Collocation points
                geo.X_C_colloc(:,j) = 0.5 * (geo.X_corners_inboard(1:end-1,j) + 0.75 * diff(geo.X_corners_inboard(:,j)) + geo.X_corners_outboard(1:end-1,j) + 0.75 * diff(geo.X_corners_outboard(:,j)));
                geo.Y_C_colloc(:,j) = repmat(0.5 * (positionVectors(2,j) + positionVectors(2,j+1)), nPanelsX, 1);
                geo.Z_C_colloc(:,j) = 0.5 * (geo.Z_corners_inboard(1:end-1,j) + 0.75 * diff(geo.Z_corners_inboard(:,j)) + geo.Z_corners_outboard(1:end-1,j) + 0.75 * diff(geo.Z_corners_outboard(:,j)));
                
                % Normal vectors on collocation points
                Xdir1 = geo.X_A_vortex(:,j) - geo.X_C_colloc(:,j);
                Ydir1 = geo.Y_A_vortex(:,j) - geo.Y_C_colloc(:,j);
                Zdir1 = geo.Z_A_vortex(:,j) - geo.Z_C_colloc(:,j);
                Xdir2 = geo.X_B_vortex(:,j) - geo.X_C_colloc(:,j);
                Ydir2 = geo.Y_B_vortex(:,j) - geo.Y_C_colloc(:,j);
                Zdir2 = geo.Z_B_vortex(:,j) - geo.Z_C_colloc(:,j);
                X_N = Ydir1 .* Zdir2 - Zdir1 .* Ydir2;
                Y_N = Zdir1 .* Xdir2 - Xdir1 .* Zdir2;
                Z_N = Xdir1 .* Ydir2 - Ydir1 .* Xdir2;
                r_N = sqrt(X_N.^2 + Y_N.^2 + Z_N.^2);
                geo.X_n(:,j) = X_N ./ r_N;
                geo.Y_n(:,j) = Y_N ./ r_N;
                geo.Z_n(:,j) = Z_N ./ r_N;

                % Adaption for airfoil camber line (rotate normal vectors)
                if ~isempty(this.airfoilDB)
                    k_inboard = find(usedAirfoilIDs == airfoilsUsed(1,j));
                    k_outboard = find(usedAirfoilIDs == airfoilsUsed(2,j));
                    thetas_airfoil_inboard = interp1(camberData(k_inboard).x_theta, camberData(k_inboard).theta, fractionsX_inboard(1:end-1,j) + 0.75 * diff(fractionsX_inboard(:,j)), "pchip");
                    if k_inboard == k_outboard
                        rotAngle = thetas_airfoil_inboard;
                    else
                        thetas_airfoil2_outboard = interp1(camberData(k_outboard).x_theta, camberData(k_outboard).theta, fractionsX_inboard(1:end-1,j) + 0.75 * diff(fractionsX_inboard(:,j)), "pchip");
                        rotAngle_inboard = airfoilBlending(1,j) * thetas_airfoil_inboard + (1 - airfoilBlending(1,j)) * thetas_airfoil2_outboard;
                        rotAngle_outboard = airfoilBlending(2,j) * thetas_airfoil_inboard + (1 - airfoilBlending(2,j)) * thetas_airfoil2_outboard;
                        rotAngle = 0.5 * (rotAngle_inboard + rotAngle_outboard);
                    end
                    geo.X_n(:,j) = cos(rotAngle) .* geo.X_n(:,j) - sin(rotAngle) .* geo.Z_n(:,j);
                    geo.Z_n(:,j) = sin(rotAngle) .* geo.X_n(:,j) + cos(rotAngle) .* geo.Z_n(:,j);
                end
            end

            % Adaption for flap deflection (rotate trailing edge panels) on half wing (for symmetric or mirrored flap
            % deflection)
            for j = 1:nPanelsY
                if i_hinge(j) <= nPanelsX && flapDeflectionAngles(j) ~= 0 && (deflectionModes(j) == 1 || deflectionModes(j) == -1)
                    
                    geo = rotateFlapPanelPoints(geo, flapDeflectionAngles(j)/pi*180, ...
                            [geo.X_corners_outboard(i_hinge(j),j) - geo.X_corners_inboard(i_hinge(j),j); ...
                             geo.Y_corners_outboard(i_hinge(j),j) - geo.Y_corners_inboard(i_hinge(j),j); ...
                             geo.Z_corners_outboard(i_hinge(j),j) - geo.Z_corners_inboard(i_hinge(j),j)], ...
                            [geo.X_corners_inboard(i_hinge(j),j); geo.Y_corners_inboard(i_hinge(j),j); geo.Z_corners_inboard(i_hinge(j),j)]);
                end
            end

            % Mirror for full span wing (Note that panel data matrices are simply concatenated)
            for s = {'A_vortex', 'B_vortex', 'D_vortex', 'E_vortex', 'F_vortex', 'G_vortex', 'C_colloc', 'n', 'corners_inboard', 'corners_outboard'}
                geo.(['X_',s{:}])(:,nPanelsY+1:end) =  geo.(['X_',s{:}])(:,1:nPanelsY);
                geo.(['Y_',s{:}])(:,nPanelsY+1:end) = -geo.(['Y_',s{:}])(:,1:nPanelsY);
                geo.(['Z_',s{:}])(:,nPanelsY+1:end) =  geo.(['Z_',s{:}])(:,1:nPanelsY);
            end

            % Rotate flap panels back on one half wing if there is an asymetric flap deflection or rotate one side if
            % there is an antisymetric deflection
            for j = 1:nPanelsY
                if i_hinge(j) <= nPanelsX && flapDeflectionAngles(j) ~= 0 && (deflectionModes(j) == -1 || deflectionModes(j) == 0)
                    
                    if deflectionModes(j) == -1
                        rotAngle = -2 * flapDeflectionAngles(j)/pi*180;
                    elseif deflectionModes(j) == 0
                        rotAngle = flapDeflectionAngles(j)/pi*180;
                    end

                    geo = rotateFlapPanelPoints(geo, rotAngle, ...
                            [geo.X_corners_outboard(i_hinge(j),j) - geo.X_corners_inboard(i_hinge(j),j); ...
                             geo.Y_corners_outboard(i_hinge(j),j) - geo.Y_corners_inboard(i_hinge(j),j); ...
                             geo.Z_corners_outboard(i_hinge(j),j) - geo.Z_corners_inboard(i_hinge(j),j)], ...
                            [geo.X_corners_inboard(i_hinge(j),j); geo.Y_corners_inboard(i_hinge(j),j); geo.Z_corners_inboard(i_hinge(j),j)]);
                end
            end


            function geoStruct = rotateFlapPanelPoints(geoStruct, angle_rotAxis_deg, dir_rotAxis, x0_rotAxis)
                % Nested function to rotate VLM points due to flap deflection
                [R,t] =  this.axelRot(angle_rotAxis_deg, dir_rotAxis, x0_rotAxis); 

                % Panels aft of flap hinge line
                for s = {'A_vortex', 'B_vortex', 'D_vortex', 'E_vortex', 'F_vortex', 'G_vortex', 'C_colloc', 'corners_inboard', 'corners_outboard'}
                    XYZ_rotated = bsxfun(@plus, R * [geoStruct.(['X_',s{:}])(i_hinge(j):end,j)'; geoStruct.(['Y_',s{:}])(i_hinge(j):end,j)'; geoStruct.(['Z_',s{:}])(i_hinge(j):end,j)'], t)';
                    geoStruct.(['X_',s{:}])(i_hinge(j):end,j) = XYZ_rotated(:,1);
                    geoStruct.(['Y_',s{:}])(i_hinge(j):end,j) = XYZ_rotated(:,2);
                    geoStruct.(['Z_',s{:}])(i_hinge(j):end,j) = XYZ_rotated(:,3);
                end
                N_rotated = (R * [geoStruct.X_n(i_hinge(j):end,j)'; geoStruct.Y_n(i_hinge(j):end,j)'; geoStruct.Z_n(i_hinge(j):end,j)'])';
                geoStruct.X_n(i_hinge(j):end,j) = N_rotated(:,1);
                geoStruct.Y_n(i_hinge(j):end,j) = N_rotated(:,2);
                geoStruct.Z_n(i_hinge(j):end,j) = N_rotated(:,3);

                % Panels in front of hinge line
                for s = {'F_vortex', 'G_vortex'}
                    XYZ_rotated = bsxfun(@plus, R * [geoStruct.(['X_',s{:}])(1:i_hinge(j)-1,j)'; geoStruct.(['Y_',s{:}])(1:i_hinge(j)-1,j)'; geoStruct.(['Z_',s{:}])(1:i_hinge(j)-1,j)'], t)';
                    geoStruct.(['X_',s{:}])(1:i_hinge(j)-1,j) = XYZ_rotated(:,1);
                    geoStruct.(['Y_',s{:}])(1:i_hinge(j)-1,j) = XYZ_rotated(:,2);
                    geoStruct.(['Z_',s{:}])(1:i_hinge(j)-1,j) = XYZ_rotated(:,3);
                end
            end
        end
    end %methods protected
    
    methods (Static)

        function fractionsOutput = getPanelDistribution(fractionsInput, nPanels)

            nZonesInput = length(fractionsInput) - 1;
            nPanelsPerPartition = ones(1, nZonesInput); % Ensure each partition gets at least one spanwise panel
            proportionalPanels = diff(fractionsInput) * (nPanels - nZonesInput); % Calculate initial proportional panel allocation % TODO: think about refinement towards wing tips
            nPanelsPerPartition = nPanelsPerPartition + floor(proportionalPanels); % Assign integer part first
            remainingPanels = nPanels - nZonesInput - sum(floor(proportionalPanels));
            remainderPanels = proportionalPanels - floor(proportionalPanels);
            % Assign remaining seats based on highest remainders
            [~, idx] = sort(remainderPanels, 'descend');
            for i = 1:remainingPanels
                nPanelsPerPartition(idx(i)) = nPanelsPerPartition(idx(i)) + 1; % Distribute remaining seats
            end
            fractionsOutput = 0;
            for j = 1:nZonesInput
                frac_temp = linspace(fractionsInput(j), fractionsInput(j+1), nPanelsPerPartition(j)+1);
                fractionsOutput = [fractionsOutput, frac_temp(2:end)];
            end
        end

        function results = solveVLM7Segment(geo, state, ref)
            % VLM solver for sling horseshoe vortices

            % Coordinates in global COS
            X_A_vortex = geo.X_A_vortex(:);
            Y_A_vortex = geo.Y_A_vortex(:);
            Z_A_vortex = geo.Z_A_vortex(:);
            X_B_vortex = geo.X_B_vortex(:);
            Y_B_vortex = geo.Y_B_vortex(:);
            Z_B_vortex = geo.Z_B_vortex(:);
            X_D_vortex = geo.X_D_vortex(:);
            Y_D_vortex = geo.Y_D_vortex(:);
            Z_D_vortex = geo.Z_D_vortex(:);
            X_E_vortex = geo.X_E_vortex(:);
            Y_E_vortex = geo.Y_E_vortex(:);
            Z_E_vortex = geo.Z_E_vortex(:);
            X_F_vortex = geo.X_F_vortex(:);
            Y_F_vortex = geo.Y_F_vortex(:);
            Z_F_vortex = geo.Z_F_vortex(:);
            X_G_vortex = geo.X_G_vortex(:);
            Y_G_vortex = geo.Y_G_vortex(:);
            Z_G_vortex = geo.Z_G_vortex(:);
            
            if isfield(geo, 'dx_AB') && ~isfield(geo, 'dy_AB') && ~isfield(geo, 'dz_AB')
                dx_AB = geo.dx_AB(:);
                dy_AB = geo.dy_AB(:);
                dz_AB = geo.dz_AB(:);
            else
                dx_AB = X_B_vortex - X_A_vortex;
                dy_AB = Y_B_vortex - Y_A_vortex;
                dz_AB = Z_B_vortex - Z_A_vortex;
            end
            
            deltaX_wake = 100 * (max(X_A_vortex) - min(X_A_vortex));
            X_F_inf_vortex = geo.X_F_vortex(:) + deltaX_wake; % Note that the trailing (semi-infinite) trailing vortex filament are not aligned with the freestream direction !
            Y_F_inf_vortex = geo.Y_F_vortex(:);
            Z_F_inf_vortex = geo.Z_F_vortex(:);
            X_G_inf_vortex = geo.X_G_vortex(:) + deltaX_wake;
            Y_G_inf_vortex = geo.Y_G_vortex(:);
            Z_G_inf_vortex = geo.Z_G_vortex(:);
            X_C = geo.X_C_colloc(:);
            Y_C = geo.Y_C_colloc(:);
            Z_C = geo.Z_C_colloc(:);
            X_n = geo.X_n(:);
            Y_n = geo.Y_n(:);
            Z_n = geo.Z_n(:);
            
            if isfield(geo, 'X_F') && ~isfield(geo, 'Y_F') && ~isfield(geo, 'Z_F')
                X_F = geo.X_F(:);
                Y_F = geo.Y_F(:);
                Z_F = geo.Z_F(:);
            else
                X_F = 0.5 * (X_A_vortex + X_B_vortex);
                Y_F = 0.5 * (Y_A_vortex + Y_B_vortex);
                Z_F = 0.5 * (Z_A_vortex + Z_B_vortex);
            end
            
            nPanels = length(X_A_vortex);
            
            % Left hand side (LHS) of system of linear equations (SLE)
            f_0 = 1 / (4 * pi);
            lhs = f_0 * (getCoeffMatrixFromVortexFilament(X_A_vortex, Y_A_vortex, Z_A_vortex, X_B_vortex, Y_B_vortex, Z_B_vortex, X_C, Y_C, Z_C) ...             % quarter line vortex
                       + getCoeffMatrixFromVortexFilament(X_D_vortex, Y_D_vortex, Z_D_vortex, X_A_vortex, Y_A_vortex, Z_A_vortex, X_C, Y_C, Z_C) ...             % left hinge-line-to-quarter-line
                       + getCoeffMatrixFromVortexFilament(X_B_vortex, Y_B_vortex, Z_B_vortex, X_E_vortex, Y_E_vortex, Z_E_vortex, X_C, Y_C, Z_C) ...             % right quarter-line-to-hinge-line
                       + getCoeffMatrixFromVortexFilament(X_F_vortex, Y_F_vortex, Z_F_vortex, X_D_vortex, Y_D_vortex, Z_D_vortex, X_C, Y_C, Z_C) ...             % left trailing-edge-to-hinge-line
                       + getCoeffMatrixFromVortexFilament(X_E_vortex, Y_E_vortex, Z_E_vortex, X_G_vortex, Y_G_vortex, Z_G_vortex, X_C, Y_C, Z_C) ...             % right hinge-line-to-trailing-edge
                       + getCoeffMatrixFromVortexFilament(X_F_inf_vortex, Y_F_inf_vortex, Z_F_inf_vortex, X_F_vortex, Y_F_vortex, Z_F_vortex, X_C, Y_C, Z_C) ... % left wake vortex
                       + getCoeffMatrixFromVortexFilament(X_G_vortex, Y_G_vortex, Z_G_vortex, X_G_inf_vortex, Y_G_inf_vortex, Z_G_inf_vortex, X_C, Y_C, Z_C));   % right wake vortex
            
            alphas = state.alpha(:);
            betas = state.beta(:);

            % Rotation of entire wing around an arbitrary axis
            if ~isfield(state, 'rot_omega')
                state.rot_omega = 0;
                omega_vec = [0;0;0];
                state.rot_X0 = [0;0;0];
            else
                if ~isfield(state, 'rot_X0') || ~isfield(state, 'rot_axis')

                    warning('Rotation axis for dynamic wing must be given by the fields rot_x0 (X|Y|Z) and rot_axis (X|Y|Z) in the state struct. Rotation frequency ignored.');
                    state.rot_omega = 0;
                    omega_vec = [0;0;0];
                    state.rot_X0 = [0;0;0];
                end
            end
            
            for j = 1:length(betas)
                for i = 1:length(alphas)
            
                    fprintf('VLM: Calculating beta = %.1f°, alpha = %.1f° ... \n', betas(j)/pi*180, alphas(i)/pi*180);
            
                    sinAlpha = sin(alphas(i));
                    cosAlpha = cos(alphas(i));
                    sinBeta = sin(betas(j));
                    cosBeta = cos(betas(j));
                    
                    % T = [ cos(beta)*cos(alpha), -sin(beta), cos(beta)*sin(alpha);...
                    %       cos(alpha)*sin(beta),  cos(beta), sin(beta)*sin(alpha);...
                    %      -sin(alpha),            0,         cos(alpha)];
            
                    t_1_1 = cosBeta * cosAlpha;
                    t_2_1 = cosAlpha * sinBeta;
                    t_3_1 = -sinAlpha;
                    t_1_2 = -sinBeta;
                    t_2_2 = cosBeta;
                    t_3_2 = 0;
                    t_1_3 = cosBeta * sinAlpha;
                    t_2_3 = sinBeta * sinAlpha;
                    t_3_3 = cosAlpha;
                           
                    V_inf_X =  t_1_1;
                    V_inf_Y = -t_2_1;
                    V_inf_Z = -t_3_1;
                    
                    % Right hand side (RHS) of system of linear equations (SLE)
                    if state.rot_omega == 0
                        rhs = -(V_inf_X * X_n + V_inf_Y * Y_n + V_inf_Z * Z_n); % no rotation
                    else
                        omega_vec = state.rot_omega * state.rot_axis / norm(state.rot_axis); % rotation vector (omega)
                        % Distance of rotation axis to collocation points (r)
                        deltaX_rot = X_C - state.rot_X0(1); 
                        deltaY_rot = Y_C - state.rot_X0(2);
                        deltaZ_rot = Z_C - state.rot_X0(3);
                        % Velocities in collocation points due to rotation (vector cross product: omega x r)
                        V_rot_X = omega_vec(2) * deltaZ_rot - omega_vec(3) * deltaY_rot;
                        V_rot_Y = omega_vec(3) * deltaX_rot - omega_vec(1) * deltaZ_rot;
                        V_rot_Z = omega_vec(1) * deltaY_rot - omega_vec(2) * deltaX_rot;
                        
                        % Right hand side (RHS) of system of linear equations (SLE) including rotation
                        rhs = -((V_inf_X - V_rot_X) .* X_n + (V_inf_Y - V_rot_Y) .* Y_n + (V_inf_Z - V_rot_Z) .* Z_n);
                    end
                    
                    % solve SLE
                    gamma = lhs \ rhs;
                    
                    % Aerodynamic forces
                    V_ind = f_0 * (getInducedVelocitiesFromVortexFilament(X_A_vortex, Y_A_vortex, Z_A_vortex, X_B_vortex, Y_B_vortex, Z_B_vortex, X_F, Y_F, Z_F, gamma, true) ...              % quarter line vortex
                                 + getInducedVelocitiesFromVortexFilament(X_D_vortex, Y_D_vortex, Z_D_vortex, X_A_vortex, Y_A_vortex, Z_A_vortex, X_F, Y_F, Z_F, gamma, false) ...             % left hinge-line-to-quarter-line
                                 + getInducedVelocitiesFromVortexFilament(X_B_vortex, Y_B_vortex, Z_B_vortex, X_E_vortex, Y_E_vortex, Z_E_vortex, X_F, Y_F, Z_F, gamma, false) ...             % right quarter-line-to-hinge-line
                                 + getInducedVelocitiesFromVortexFilament(X_F_vortex, Y_F_vortex, Z_F_vortex, X_D_vortex, Y_D_vortex, Z_D_vortex, X_F, Y_F, Z_F, gamma, false) ...             % left trailing-edge-to-hinge-line
                                 + getInducedVelocitiesFromVortexFilament(X_E_vortex, Y_E_vortex, Z_E_vortex, X_G_vortex, Y_G_vortex, Z_G_vortex, X_F, Y_F, Z_F, gamma, false) ...             % right hinge-line-to-trailing-edge
                                 + getInducedVelocitiesFromVortexFilament(X_F_inf_vortex, Y_F_inf_vortex, Z_F_inf_vortex, X_F_vortex, Y_F_vortex, Z_F_vortex, X_F, Y_F, Z_F, gamma, false) ... % left wake vortex
                                 + getInducedVelocitiesFromVortexFilament(X_G_vortex, Y_G_vortex, Z_G_vortex, X_G_inf_vortex, Y_G_inf_vortex, Z_G_inf_vortex, X_F, Y_F, Z_F, gamma, false));   % right wake vortex
                    
                    Veff_X = V_ind(:,1) + V_inf_X;
                    Veff_Y = V_ind(:,2) + V_inf_Y;
                    Veff_Z = V_ind(:,3) + V_inf_Z;
                             
                    F_X = gamma .* (Veff_Y .* dz_AB - Veff_Z .* dy_AB);
                    F_Y = gamma .* (Veff_Z .* dx_AB - Veff_X .* dz_AB);
                    F_Z = gamma .* (Veff_X .* dy_AB - Veff_Y .* dx_AB);
                    
                    results(i,j).c_D = 2 * sum(t_1_1 * F_X + t_1_2 * F_Y + t_1_3 * F_Z) / ref.S_ref;
                    results(i,j).c_S = 2 * sum(t_2_1 * F_X + t_2_2 * F_Y + t_2_3 * F_Z) / ref.S_ref;
                    results(i,j).c_L = 2 * sum(t_3_1 * F_X + t_3_2 * F_Y + t_3_3 * F_Z) / ref.S_ref;
                    results(i,j).c_X = 2 * sum(F_X) / ref.S_ref;
                    results(i,j).c_Y = 2 * sum(F_X) / ref.S_ref;
                    results(i,j).c_Z = 2 * sum(F_Z) / ref.S_ref;
                    
                    deltaX_ref = X_F - ref.X_mom_ref;
                    deltaY_ref = Y_F - ref.Y_mom_ref;
                    deltaZ_ref = Z_F - ref.Z_mom_ref;
                    L_body_panels = -F_Y .* deltaZ_ref + F_Z .* deltaY_ref;
                    M_body_panels = -F_Z .* deltaX_ref + F_X .* deltaZ_ref;
                    N_body_panels = -F_X .* deltaY_ref + F_Y .* deltaX_ref;
                    results(i,j).L_body = sum(L_body_panels);
                    results(i,j).M_body = sum(M_body_panels);
                    results(i,j).N_body = sum(N_body_panels);
                    results(i,j).L_aero = t_1_1 * results(i,j).L_body + t_1_2 * results(i,j).M_body + t_1_3 * results(i,j).N_body;
                    results(i,j).M_aero = t_2_1 * results(i,j).L_body + t_2_2 * results(i,j).M_body + t_2_3 * results(i,j).N_body;
                    results(i,j).N_aero = t_3_1 * results(i,j).L_body + t_3_2 * results(i,j).M_body + t_3_3 * results(i,j).N_body;
                    results(i,j).c_l = 2 * sum(results(i,j).L_aero) / ref.S_ref / ref.b_ref;
                    results(i,j).c_m = 2 * sum(results(i,j).M_aero) / ref.S_ref / ref.c_ref;
                    results(i,j).c_n = 2 * sum(results(i,j).N_aero) / ref.S_ref / ref.b_ref;
                    results(i,j).S_ref = ref.S_ref;
                    results(i,j).c_ref = ref.c_ref;
                    results(i,j).b_ref = ref.b_ref;
                    results(i,j).p_ref_mom = [ref.X_mom_ref; ref.Y_mom_ref; ref.Z_mom_ref];
                    results(i,j).alpha = alphas(i);
                    results(i,j).beta = betas(j);
                    results(i,j).omega_rot = omega_vec;
                    results(i,j).P0_rot = state.rot_X0;
                    results(i,j).X_F = X_F;
                    results(i,j).Y_F = Y_F;
                    results(i,j).Z_F = Z_F;
                    results(i,j).F_X = F_X;
                    results(i,j).F_Y = F_Y;
                    results(i,j).F_Z = F_Z;
                    results(i,j).Vind = V_ind;
                    results(i,j).Veff = [Veff_X, Veff_Y, Veff_Z];
                    results(i,j).T = [t_1_1, t_1_2, t_1_3;...
                                      t_2_1, t_2_2, t_2_3;...
                                      t_3_1, t_3_2, t_3_3];
                    
                    fprintf('CD: %.3f | CS: %.4f | CL: %.3f \nCl: %.3f | Cm: %.4f | Cn: %.3f \n', results(i,j).c_D, results(i,j).c_S, results(i,j).c_L, results(i,j).c_l, results(i,j).c_m, results(i,j).c_n);
                end
            end
            
            fprintf('Reference point (X|Y|Z) of moment coefficients:  (%.1f|%.1f|%.1f) \n', ref.X_mom_ref, ref.Y_mom_ref, ref.Z_mom_ref);

            if state.rot_omega ~= 0
                fprintf('Imposed rotation of %.1f rad/s about axis (%.1f|%.1f|%.1f) in point (%.1f|%.1f|%.1f) \n', state.rot_omega, omega_vec(1) / norm(omega_vec), omega_vec(2) / norm(omega_vec), omega_vec(3) / norm(omega_vec), state.rot_X0(1), state.rot_X0(2), state.rot_X0(3));
            end
            
            function A = getCoeffMatrixFromVortexFilament(X_A, Y_A, Z_A, X_B, Y_B, Z_B, X_P, Y_P, Z_P)
            
                [f_4, dx_AP, dy_AP, dz_AP, dx_BP, dy_BP, dz_BP] = getInfluenceFactorFromVortexFilament(X_A, Y_A, Z_A, X_B, Y_B, Z_B, X_P, Y_P, Z_P);
            
                A = f_4 .* (repmat(X_n, 1, nPanels) .* (dy_AP .* dz_BP - dz_AP .* dy_BP) ...
                          + repmat(Y_n, 1, nPanels) .* (dz_AP .* dx_BP - dx_AP .* dz_BP) ...
                          + repmat(Z_n, 1, nPanels) .* (dx_AP .* dy_BP - dy_AP .* dx_BP));
            end
            
            function [f_4, dx_AP, dy_AP, dz_AP, dx_BP, dy_BP, dz_BP] = getInfluenceFactorFromVortexFilament(X_A, Y_A, Z_A, X_B, Y_B, Z_B, X_P, Y_P, Z_P)
            
                dx_AP = bsxfun(@plus, X_P, -X_A');
                dy_AP = bsxfun(@plus, Y_P, -Y_A');
                dz_AP = bsxfun(@plus, Z_P, -Z_A');
                dr_AP = sqrt(dx_AP.^2 + dy_AP.^2 + dz_AP.^2);
                dx_BP = bsxfun(@plus, X_P, -X_B');
                dy_BP = bsxfun(@plus, Y_P, -Y_B');
                dz_BP = bsxfun(@plus, Z_P, -Z_B');
                dr_BP = sqrt(dx_BP.^2 + dy_BP.^2 + dz_BP.^2);
            
                f_1 = dr_AP + dr_BP;
                f_2 = dr_AP .* dr_BP;
                f_3 = dx_AP .* dx_BP + dy_AP .* dy_BP + dz_AP .* dz_BP;
                f_4 = f_1 ./ (f_2 .* (f_2 + f_3));
            end
            
            function Vind = getInducedVelocitiesFromVortexFilament(X_A, Y_A, Z_A, X_B, Y_B, Z_B, X_P, Y_P, Z_P, gamma, isPonAB)
            
                [f_4, dx_AP, dy_AP, dz_AP, dx_BP, dy_BP, dz_BP] = getInfluenceFactorFromVortexFilament(X_A, Y_A, Z_A, X_B, Y_B, Z_B, X_P, Y_P, Z_P);
                
                if isPonAB % if point P is located on vortex filament AB
                    i_diag = logical(eye(nPanels));
                    f_4(i_diag) = 0;
                end
            
                Vind = [(f_4 .* (dy_AP .* dz_BP - dz_AP .* dy_BP)) * gamma, ...
                        (f_4 .* (dz_AP .* dx_BP - dx_AP .* dz_BP)) * gamma, ...
                        (f_4 .* (dx_AP .* dy_BP - dy_AP .* dx_BP)) * gamma];
            end
            end
    end
end
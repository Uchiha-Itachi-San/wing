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
%*  Stadlberger, Korbinian  * 2012-AUG-06 *     Initial Release (LLS/TUM)
%*  Stadlberger, Korbinian  * 2024-NOV-04 *     Adaption for THI lecture
%*  Stadlberger, Korbinian  * 2025-FEB-07 *     Adaption for THI lecture
%****************************************************************************************************************

classdef MultiPartitionWing < Wing
   
    properties (SetAccess = private)
        
         % The following properties are inherited from abstract class Wing
         %
         % aspectRatio            % Aspect ratio                                double           (1x1)       [-]           
         % taperRatio             % Taper ratio / meaning depends	double           (1x1)       [-] 
         % 						  % on the class inherited from this class                
         % nPartitions            % Number of partitions used for the halfspan  int              (1x1)       [-]                                                            
         % relSpanPositions       % Span positions relative to the halfspan     double           (1xn+1)     [-]             
         % chordLengths           % Chord lengths at relative span positions    double           (1xn+1)     [m]              
         % twistAngles            % Twist angles for each partition	            double           (1x1)       [rad]       
         % dihedralAngles         % Dihedral angles for each span position      double           (1x1)       [rad]       
         % sweepAngles            % Sweep angles for each span position         double           (1x1)       [rad]       
         % incidenceAngle = 0     % Incidence angle of entire wing              double           (1x1)       [rad]       
         % airfoilDB              % Airfoil database                            Airfoil          (1x1)       [-]       
         % airfoilsUsed           % Airfoil IDs for each span position          double           (1x1)       [-]       
         % flapChordLengths       % Relative flap chord length for each         double           (1x1)       [-]       
         % 						  % partition                                                    		       
         % flapDeflections        % Flap deflection for each partition          double           (1x1)       [rad]       
         % deflectionModes 		  % Deflection Modes for each partition         int              (1x1)       [-]   
         % 						  % (1)symmetric (-1)antisymmetric                               		       
         % flapControlIDs		  % Flap Control IDs for each partition         int              (1x1)       [-]   
         % dihedralKinks          % Collection of dihedral kinks                struct           (1x1)       [-]      
         % sweepKinks             % Collection of sweep kinks                   struct           (1x1)       [-]      
         % ms_carbon              % Specific mass of used carbon fibres         double           (1x1)       [kg/m³]      
         % b_spar_i               % Width of spar                               double           (1x1)       [mm]         
        
    end
    
    methods 
        
        function this = MultiPartitionWing(chordTable, airfoilDB_wing)
            % Constructor of class MultiPartitionWing.
            
%             Input variables:
%             	chordTable:			table of chord lengths and span position, double (2xm)
%             	airfoilDB_wing:     Airfoil used for entire wing, object of type 'Airfoil'
%             	nPartitions:		Number of partitions used to model half wing
%             	
%             Output variables:
%             	this: 				object of type 'MultiPartitionWing'
%             
%             Called functions:
%             	Wing() of abstract heritage class 'Wing'
%             	calculateRelSpanPositions() inherited from class 'Wing'
%             	calculateChordLengths() of this class and inherited from abstract class 'Wing'
%             	initialiseTwistAngles() inherited from class 'Wing'
%             	initialiseDihedralAngles() inherited from class 'Wing'
% 					initialiseSweepAngles() inherited from class 'Wing'    
% 					initialiseAirfoils() inherited from class 'Wing'       
% 					initialiseFlapChordLengths() inherited from class 'Wing'
% 					initialiseFlapDeflections() inherited from class 'Wing'
% 					initialiseDeflectionModes() inherited from class 'Wing'
% 					initialiseFlapControlID() inherited from class 'Wing'  
% 					refreshLref() inherited from class 'AerodynamicGeometry'

            this = this@Wing();
            
            this.name = 'multi_partition_wing';
            this.airfoilDB = airfoilDB_wing;
            
            chordTable = sortrows(chordTable', 1)';
            [~, tempm, ~] = unique(chordTable(1,:));
            if length(tempm) < size(chordTable, 2)
                fprintf('WARNING: Span positions are not unique! Doubled entries deleted -> Information lost.\n');
            end
            chordTable = chordTable(:,tempm);
            
            for i = 1:size(chordTable, 2)
                this.chordKinks(i).relYPos = chordTable(1,i) / chordTable(1,end);
                this.chordKinks(i).chordLength = chordTable(2,i);
            end
            
            % set wing parameters
            this.Sref = 0;
            for i = 1:size(chordTable,2)-1
                this.Sref = this.Sref + sum(chordTable(2,[i,i+1])) * (chordTable(1,i+1) - chordTable(1,i));
            end
            
            this.aspectRatio = (2 * chordTable(1,end))^2 / this.Sref;
            this.taperRatio = chordTable(2,end) / chordTable(2,1);

            
            this.calculateRelSpanPositions();
            this.calculateChordLengths();
            
            this.initialiseTwistAngles();
            this.initialiseDihedralAngles();
            this.initialiseSweepAngles();
            this.initialiseAirfoils();
            this.initialiseFlapChordLengths();
            this.initialiseFlapDeflections();
            this.initialiseDeflectionModes();
            this.initialiseFlapControlID();
            
            this.refreshLref();
            
        end
        
        function calculateChordLengths(this)
            % Calculates the chord lengths of each section along the halfspan. This function 
            % inherited from abstract class 'Wing'. This function implements the specific type of 
            % planform this class represents to obtain a certain shape of the wing. In this case the 
            % planform has an abritrary shape of combined tapered wing partitions.
            % 
            	
            this.chordLengths = interp1([this.chordKinks.relYPos], [this.chordKinks.chordLength], this.relSpanPositions);
        end
        
       
%         function refreshWing(this)
%             % Refreshes all important properties, which define the wing shape.
%             % 
%             % Called functions:
% 	        %    calculateRelSpanPositions() inherited from class 'Wing'
%             % 	calculateChordLengths() of this class
%             % 	setSweepAngles() inherited from class 'Wing'
%             % 	setDihedralAngles() inherited from class 'Wing'   
%             % 	refreshLref()  inherited from class 'AerodynamicGeometry'            
%             
%             this.calculateRelSpanPositions();
%             this.calculateChordLengths();
%             this.setSweepAngles();
%             this.setDihedralAngles();
%             this.refreshLref();
%         end
             
        function refreshRelCGEstimated(this)
            
            relCG_Position = 0.3; % 30% of MAC
            
            this.refreshLref;
            ac = this.getACPosition();
            cg = ac(1) + (relCG_Position - 0.25) * this.Lref;
            this.relCG = [cg;0;0];
        end
        
    end %methods public
    
    methods (Static = true)
        
    end
end
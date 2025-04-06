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
%*  Stadlberger, Korbinian  * 2024-NOV-04 *     Adaption for THI lecture
%****************************************************************************************************************

classdef TaperedWing < Wing
   
    properties (SetAccess = private)
        
         % The following properties are inherited from abstract class Wing
         %
         % aspectRatio            % Aspect ratio                                double           (1x1)       [-]           
         % taperRatio             % Taper ratio / meaning depends	double           (1x1)       [-] 
         % 						  % on the class inherited from this class                
         % nPartitions            % Number of partitions used for the halfspan  int              (1x1)       [-]  
         % nPanels                % Global number of panels the partitions are  int              (1x1)       [-]       
         % 						  % divided into                                                             
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
        
        function this = TaperedWing(Sref, airfoilDB_wing, AR, TR)
            % Constructor of class OverellipticalWing.
            
%             Input variables:
%             	Sref:					Reference area in m², double (1x1)
%             	airfoilDB_wing:	Airfoil used for entire wing, object of type 'Airfoil'
%             	AR:					Wing aspect ratio, double (1x1)
%             	TR:					Wing taper ratio, double (1x1)
%             	nPartitions:		Number of partitions used to model half wing
%             	
%             Output variables:
%             	this: 				object of type 'OverellipticalWing'
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
            
            this.airfoilDB = airfoilDB_wing;
            this.aspectRatio = AR;
            this.taperRatio = TR;
            this.Sref = Sref;
            
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
            % planform has a trapezidoal shape and represents a tapered wing.
            % 
            % Called functions:
            % 	getSpan() inherited from class 'Wing'
            	
            this.chordLengths = [];
            %calculation of root chord length
            span = this.getSpan();
            TR = this.taperRatio;
            rootChord = 2 * this.Sref / (1+TR) / span ;
            
            %calculation of chord length at every span position
            for i = 1:1:length(this.relSpanPositions)
                this.chordLengths(i) = rootChord * (1 - this.relSpanPositions(i) * (1 - TR));
            end
        end

        function setWingSize(this, newSref)
            % Setter for new wing reference area (property Sref).
            %
            % Input variables:
            %	newSref: Reference area in m²
            %
            % Called functions:
            %	refreshWing (abstract) of any class inherited from this class
            
            this.Sref = newSref;
            this.calculateRelSpanPositions();
            this.calculateChordLengths();
            this.refreshLref();
        end
        
        function angle = getSweepAngleOfPercLine(this, eta)
           
            phi25 = this.getEquivalentSweepAngle();
            phi0 = 0.5 * pi - atan((1/tan(0.5 * pi - phi25) - 0.25 * 4 * (this.taperRatio - 1) / (this.aspectRatio * (this.taperRatio + 1)))^-1);
            angle = atan((1/tan(0.5 * pi - phi0) + eta * 4 * (this.taperRatio - 1) / (this.aspectRatio * (this.taperRatio + 1)))^-1);
            angle = sign(angle) * 0.5 * pi - angle;
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
%             
%             this.setSweepAngles();
%             this.setDihedralAngles();
%             this.refreshLref();
%         end
        
        function refreshMassEstimated(this, mtom, m_surfacectrls, n_max, V_D_max)
           
            % Unit conversion factors
            kg_lb = 2.204624;
            lb_kg = 0.453592;
            ft_m = 0.3048;
            m_ft = 1 / ft_m;
            sqm_sqft = m_ft^2;
            mps_kts = 1.943844492;
            
            if this.aspectRatio >= 3
                % Nicolai Light Utility Aircraft
                massWing = 96.948 * ((mtom * kg_lb * n_max / 10^5)^0.65...
                                  * (this.aspectRatio / cos(this.getEquivalentSweepAngle))^0.57...
                                  * (this.Sref * sqm_sqft / 100)^0.61...
                                  * ((1 + this.taperRatio) / (2 * this.t_c_max_global))^0.36...
                                  * (1 + V_D_max * mps_kts / 500)^0.5)^0.993...
                                  * (1 + this.k_retractable-1)...
                                  * lb_kg;
            else
                % Nicolai conventional metal aircraft - moderate subsonic to supersonic performance
                massWing = 3.08   * ((mtom * kg_lb * n_max / this.t_c_max_global) * 10^-6 ...
                                  * ((tan(this.getSweepAngleOfPercLine(0)) - 2 * (1 - this.taperRatio) / (this.aspectRatio * (1 + this.taperRatio)))^2 + 1))^0.593 ...
                                  * ((1 + this.taperRatio) * this.aspectRatio)^0.89 * (this.Sref * sqm_sqft)^0.741 ...
                                  * (1 + this.k_retractable-1)...
                                  * lb_kg;
            end
           
            this.mass = massWing * (1 - this.k_FRP) + m_surfacectrls;
        end
        
        function refreshRelCGEstimated(this)
            
            relCG_Position = 0.3; % 30% of MAC
            
            this.refreshLref;
            ac = this.getACPosition();
            cg = ac(1) + (relCG_Position - 0.25) * this.Lref;
            this.relCG = [cg;0;0];
        end

        function filename = writeDATCOMInputFile(this, mach2calc, alt2calc, aoa2calc)

            % mach number, number of angles of attack limited to 20 values
            caseName = 'TaperedWing';

            unit_dimensions = 'M';
            unit_derivatives = 'RAD'; % causes the static and dynamic stability derivatives to be output in radian measure. The output will be in degree measure unless this flag is set.
            loop = 2.0; % Program looping control (1 vary altitude and mach together, default, 2 vary mach at fixed altitude, 3 vary altitude at fixed mach)
            stmach = 0.6; % upper limit of mach numbers for subsonic analysis (0.6 <= stmach <= 0.99). default to 0.6 if not input
            tsmach = 1.4; % lower limit of mach numbers for supersonic analysis analysis (1.01 <= tsmach <= 1.4). default to 1.4 if not input
            transitionFlag = 0.0; % drag due to lift transition flag, for regression analysis of wing-body configurations (0 for no transition, default; 1.0 for transition strips or full scale flight)
            aoa2calc_deg = aoa2calc /pi*180;
            
            Sref = this.Sref; % reference area, value of theoretical wing area used by program if not input
            c_bar = this.Lref; % longitudinal reference length value of theoretical wing mean aerodynamic chord used by program if not input
            b_l_ref = this.getSpan(); % lateral reference length value of wing span used by program if not input
            pos_AC = this.getACPosition(); % get position of the wing aerodynamic center 
            x_cg = pos_AC(1); % longitudinal location of cg = moment ref. center
            z_cg = 0.0; % vertical location of cg relative to reference plane
            x_apex = 0.0; % longitudinal location of theoretical wing apex
            z_apex = 0.0; % vertical location of theoretical wing apex relative to reference plane
            
            c_root = this.getRootChord(); % root chord
            c_tip = this.getTipChord(); % tip chord
            semispan = b_l_ref/2; % semi-span theoretical panel from theoretical root chord
            expSemispan = semispan; % semi-span exposed panel
            sweepAngle25Perc = this.getSweepAngle(1) /pi*180; % wing panel sweep angle
            refXCsweep = 0.25; % reference chord station for panel sweep angles, fraction of chord
            twist = this.getTwistAngle(1) /pi*180; % Twist angle, negative leading edge rotated down (from exposed root to tip)
            dihedral = this.getDihedralAngle(1) /pi*180;
            wingType = 1.0; % 1: straight tapered planform, 2: double delta planform (AR <= 3), 3: cranked planform (AR > 3)
            
            airfoilType = this.airfoilDB(this.airfoilsUsed(1)).type; % 1: 1-series, 4: 4-digit, 5: 5-digit, 6: 6-series, S: supersonic
            if strcmp(this.airfoilDB(this.airfoilsUsed(1)).type, 'Supersonic')
                airfoilTypeDATCOM = 'S';
            elseif strcmp(this.airfoilDB(this.airfoilsUsed(1)).type, 'NACA')
                switch length(this.airfoilDB(this.airfoilsUsed(1)).geoCode)
                    case 4
                        airfoilTypeDATCOM = '4';
                    case 5
                        airfoilTypeDATCOM = '5';
                    case 6
                        airfoilTypeDATCOM = '6';
                    otherwise
                        error('Airfoil designation problem: %s', this.airfoilDB(this.airfoilsUsed(1)).geoCode)
                end
            else
                error('Unknown airfoil type: %s', airfoilType);
            end
            airfoilDesignation = this.airfoilDB(this.airfoilsUsed(1)).geoCode;
            
            caseID = ['CASEID ',caseName];


            disp('Writing DATCOM input file...')

            % Namelists pre-process
            varindent = '         ';    % 9 spaces
            
            % Intro
            declar1 = 'DIM %s\n';     % M or FT
            declar2 = 'DERIV %s\n\n';   % DEG or RAD
            
            % Flight Conditions
            fconopen = ' $FLTCON LOOP=%.1f,\n';
            nmach = [varindent,'NMACH=%.1f,\n'];
            mach1 = ['MACH(1)=',num2str(mach2calc,'%.1f, ')];
            mach1 = dcmArraySplit(mach1,varindent); % split char array within FORTRAN limit
            nalt = [varindent,'NALT=%.1f,\n'];
            altitude1 = ['ALT(1)=',num2str(alt2calc,'%.1f, ')];
            altitude1 = dcmArraySplit(altitude1,varindent); % split char array within FORTRAN limit
            nalpha = [varindent,'NALPHA=%.2f,\n'];
            alpha1 = ['ALSCHD(1)=',num2str(aoa2calc_deg,'%.1f, ')];
            alpha1 = dcmArraySplit(alpha1,varindent); % split char array within FORTRAN limit
            fconclose = [varindent,'STMACH=%.2f, TSMACH=%.2f, TR=%.2f$\n\n'];
            
            % Reference Values
            optins = ' $OPTINS SREF=%.2f, CBARR=%.2f, BLREF=%.2f$\n\n';
            synths = ' $SYNTHS XCG=%.2f, ZCG=%.2f, XW=%.2f, ZW=%.2f$\n\n';
            
            % Wing Planform
            wingopen = ' $WGPLNF CHRDR=%.2f, CHRDTP=%.2f,\n';
            halfspan = [varindent,'SSPN=%.2f, SSPNE=%.2f,\n'];
            sweeptwist = [varindent,'SAVSI=%.2f, CHSTAT=%.2f, TWISTA=%.1f,\n'];
            wingclose = [varindent,'DHDADI=%.2f, TYPE=%.2f$\n\n'];
            
            % Airfoil
            wingprof = 'NACA W %s %s\n\n';
            
            % Write Digital Datcom input file
            filename = caseName;
            filename = regexprep(filename, ' ', '_');   % check if there are unvalid
            filename = regexprep(filename, '\', '_');   % characters in file name
            filename = regexprep(filename, '=', '_');
            filename = regexprep(filename, '°', '');
            filename = regexprep(filename,'_+','_');    % suppress consecutive underscores
            if ~strcmp(caseName,filename)
               warning([caseName,' contains unvalid characters. ',...
                   'The name tag will be preserved, but the file has been renamed as ',...
                   filename]) 
            end
%             filename = 'for005.dat';
            
            fid = fopen([filename, '.dcm'],'w');
            fprintf(fid,'*\n');
            fprintf(fid,'*   %s\n',caseName);
            fprintf(fid,'*\n\n');
            fprintf(fid,declar1,unit_dimensions);
            fprintf(fid,declar2,unit_derivatives);
            fprintf(fid,fconopen,loop);
            
            fprintf(fid,nmach,length(mach2calc));
            if iscell(mach1)
                for idx = 1:length(mach1)
                    fprintf(fid,mach1{idx});
                end
            else
                fprintf(fid,mach1);
            end
            
            fprintf(fid,nalt,length(alt2calc));
            if iscell(altitude1)
                for idx = 1:length(altitude1)
                    fprintf(fid,altitude1{idx});
                end
            else
                fprintf(fid,altitude1);
            end
            
            fprintf(fid,nalpha,length(aoa2calc_deg));
            if iscell(alpha1)
                for idx = 1:length(alpha1)
                    fprintf(fid,alpha1{idx});
                end
            else
                fprintf(fid,alpha1);
            end
            
            fprintf(fid,fconclose,stmach,tsmach,transitionFlag);
            fprintf(fid,optins,Sref,c_bar,b_l_ref);
            fprintf(fid,synths,x_cg,z_cg,x_apex,z_apex);
            fprintf(fid,wingopen,c_root,c_tip);
            fprintf(fid,halfspan,semispan,expSemispan);
            fprintf(fid,sweeptwist,sweepAngle25Perc,refXCsweep,twist);
            fprintf(fid,wingclose,dihedral,wingType);
            fprintf(fid,wingprof,airfoilTypeDATCOM,airfoilDesignation);
            fprintf(fid,caseID);
            fclose(fid);
            
            disp(['DATCOM input file written into: ', filename,'.dcm'])
            
            
            
            
            function splitArray = dcmArraySplit(myArray,varindent)
            % Split char array within FORTRAN limit
            
            chunk = 72-length(varindent); % original indentation
            rows = ceil((length(myArray)+length(varindent))/chunk);
            if rows > 1
                c = 1;
                for i = 1:rows-1
                    chunk = 72-length(varindent); % restore original indentation
                    
                    % change index where to split array if last char is not a delimiter
                    while ~strcmp(myArray(c+chunk),',') && ~strcmp(myArray(c+chunk),' ')
                        chunk = chunk - 1;
                    end
                    
                    splitArray{i} = myArray(c:c+chunk-1);
                    c = c + chunk;
                end
                splitArray{i+1} = myArray(c:end);
                
                % Concatenate char arrays to include initial indentation
                for i = 1:rows
                    splitArray{i} = [varindent, splitArray{i},'\n'];
                end
                
            else
                splitArray = [varindent,myArray,'\n'];
            end
            
            end
        end

    end %methods public
end
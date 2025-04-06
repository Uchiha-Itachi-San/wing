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
%*                                         AerodynamicProperties                                                *
%*                                                                                                              *
%****************************************************************************************************************
%*                                                                                                              *
%* Type        :   Class Definition                                                                             *
%*                                                                                                              *
%* Circulation :                                                                    							*
%*                                                                                                              *
%* Purpose     :   Definition of the class AerodynamicProperties                                                *
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


classdef AerodynamicProperties < handle
   
    %************************************************************************************************************
%%  %*                                         Properties                                                       *
    %************************************************************************************************************
    %*  Stadlberger, Korbinian           2010-JUN-22                                                            *
    %************************************************************************************************************
    %*                                                                                                          *
    %* Property                		Description                        		Type       					Dim      Unit* 
   properties (SetAccess = private)
       
       CD 									% Drag coefficient							object of class MDT		(1x1)    [-]
       CDi 									% Induced drag coefficient					object of class MDT		(1x1)    [-]
       CC 									% Sideslip coefficient						object of class MDT     (1x1)    [-]
       CL 									% Lift coefficient							object of class MDT     (1x1)    [-]
       Cl 									% Roll moment coefficient					object of class MDT     (1x1)    [-]
       Cm 									% Pitching moment coefficient				object of class MDT     (1x1)    [-]
       Cn 									% Yawing moment coefficient				object of class MDT     (1x1)    [-]
       panelData                            % aerodynamic data of single panels         struct              (nxm)   [-]
       flightstate
       derivatives 						% Derivatives of first order 				struct                  (1x1)    [-]
       CL_max 								% Maximum Lift coefficient					double                  (1x1)    [-]
       alpha_zero_lift 					% Angle of attack at zero lift			double                  (1x1)    [rad]
       trimmedPolar 						% Aerodynamic data for trimmed states	struct                  (1x1)    [-]
       customData                           % custom data from any source
   end
   
       %************************************************************************************************************
%%  %*                                          Methods                                                         *
    %************************************************************************************************************
    %*  Stadlberger, Korbinian           2010-JUN-21                                                            *
    %************************************************************************************************************
    %*                                                                                                          *
   methods
       
       function this = AerodynamicProperties()
       		% Constructor of abstract class AerodynamicProperties. No initialisation performed! 
       end
       
       function setCD(this, newCDData)
           % Setter for property CD
           %
           % Input variables:
           % newCDData: double or object of type 'MDT'
           
           this.CD = newCDData;
       end
       
       function setCDi(this, newCDiData)
           % Setter for property CDi
           %
           % Input variables:
           % newCDiData: double or object of type 'MDT'
           
           this.CDi = newCDiData;
       end
       
       function setCC(this, newCCData)
           % Setter for property CC
           %
           % Input variables:
           % newCCData: double or object of type 'MDT'
           
           this.CC = newCCData;
              end
       
       function setCL(this, newCLData)
           % Setter for property CL
           %
           % Input variables:
           % newCLData: double or object of type 'MDT'
           
           this.CL = newCLData;
       end
       
       function setCl(this, newClData)
           % Setter for property Cl
           %
           % Input variables:
           % newClData: double or object of type 'MDT'
           
           this.Cl = newClData;
       end
       
       function setCm(this, newCmData)
           % Setter for property CC
           %
           % Input variables:
           % newCmData: double or object of type 'MDT'
           
           this.Cm = newCmData;
       end       
       
       function setCn(this, newCnData)
           % Setter for property Cn
           %
           % Input variables:
           % newCnData: double or object of type 'MDT'
           
           this.Cn = newCnData;
       end
       
       function setPanelData(this, newPanelData)
           % Setter for property panelData
           
           this.panelData = newPanelData;
       end
       
       function setDerivatives(this, newDerivatives)
           % Setter for property derivatives
           %
           % Input variables:
           % newDerivatives: struct with fields for various derivatives
           
           this.derivatives = newDerivatives;
       end
       
       function setFlightState(this, newFlightState)
           
           this.flightstate = newFlightState;
       end
       
       function setCustomData(this, newData)
           
           this.customData = newData;
       end
       
       function refreshCL_max(this)
           % Refreshes the property CL_max by searching the maximum value stored in the property CL.
           
           [CLmax, index] = max(this.CL.OutputData);
           
           this.CL_max = CLmax;
       end
       
       function setCLMax(this, newCL_max)
           % Setter for property CL_max
           %
           % Input variables:
           % newCL_max: double
           
           this.CL_max = newCL_max;
       end
       
       function alpha = getAlphaFromCL(this, CL_given)
           % Returns the angle of attack corresponding to a given CL value.
           % 
           % Input variables: 
           % 		CL_given: Given CL value
           % 		
           % Output variables:
           % 		alpha: Angle of attack in radiants
           %
           % Called functions:
           % 	getCLValue(alpha) of this class
           
           function difference = myfun(x)
              	% Calculates the differece between the given and an arbitrary CL value. This function 
              	%will be "optimised" to zero.
              	
               difference = this.getCLValue(x) - CL_given;
           end
           
           options = struct(...
                    'Display', 'off',...
                    'TolX', 1e-6);
           % Zero search function to find zero of function above (MATLAB implemented)
           alpha = fzero(@myfun, 0, options);
       end
       
       function alpha = getAlphaAtCLMax(this)
           % Returns the angle of attack at maximum CL.
           %
           % Output variables:
           %	alpha: angle of attack in radiants
           % 
           % Called functions:
           % 	getCLValue(alpha) of this class
           
           alpha = fminbnd(@(x)(-1*this.getCLValue(x)), this.CL.InputMin, this.CL.InputMax);
       end
       
       function refreshAlphaZeroLift(this)
           % Refreshes the property alpha_zero_lift. If the calculated data points do not intersect 
           % the x-axis an error message is returned (linear extrapolation could avoid the error 
           % message but would not be more robust either).
           
           alpha0 = fzero(@(x)this.getCLValue(x), -0.175);
           
           if (isempty(alpha0))
               error('Too few alpha angles on beta=0 calculated to find zero lift!');
           else
               this.alpha_zero_lift = alpha0;
           end

       end
       
       function Cm0 = getCm0(this)
           % Returns the pitching moment coefficient at zero lift.
           % 
           % Output variables:
           % 	Cm0: pitching moment coefficient at zero lift
           % 
           % Called functions:
           % 	refreshAlphaZeroLift() of this class
           % 	getCmValue(alpha) of this class
           	
           this.refreshAlphaZeroLift();
           
           Cm0 = this.getCmValue(this.alpha_zero_lift);
       end
       
       function CL = getCLValue(this, alpha, beta)
           % Returns the corresponding lift coefficient value to a given angle of attack (and yaw 
           % angle).
           % 
           % Input variables:
           % 	alpha: angle of attack in radiants
           % 	beta (optional): yaw angle in radiants (is set to zero if not given)
           % 	
           % Output variables:
           % 	CL: CL value interpolated from stored aerodynamic data
           %
           % Called functions:
           %	GetIPValue(value) of class MDT
           
           if nargin == 2
               if size(this.CL.InputBreakpoints,2) == 2
                    [CL Clipped] = this.CL.GetIPValue([alpha 0]);
               else
                    [CL Clipped] = this.CL.GetIPValue(alpha);
                    
               end
           elseif nargin == 3
               
               [CL Clipped] = this.CL.GetIPValue([alpha beta]);
           end
       end
       
       function CD = getCDValue(this, alpha, beta)
           % Returns the corresponding drag coefficient value to a given angle of attack (and yaw 
           % angle).
           % 
           % Input variables:
           % 	alpha: angle of attack in radiants
           % 	beta (optional): yaw angle in radiants (is set to zero if not given)
           % 	
           % Output variables:
           % 	CD: CD value interpolated from stored aerodynamic data
           %
           % Called functions:
           %	GetIPValue(value) of class MDT
           
           if nargin == 2
               if size(this.CD.InputBreakpoints,2) == 2
                    [CD Clipped] = this.CD.GetIPValue([alpha 0]);
               else
                    [CD Clipped] = this.CD.GetIPValue(alpha);
                    
               end
           elseif nargin == 3
               
               [CD Clipped] = this.CD.GetIPValue([alpha beta]);
               
           end
       end
       
       function Cm = getCmValue(this, alpha, beta)
           % Returns the corresponding drag coefficient value to a given angle of attack (and yaw 
           % angle).
           % 
           % Input variables:
           % 	alpha: angle of attack in radiants
           % 	beta (optional): yaw angle in radiants (is set to zero if not given)
           % 	
           % Output variables:
           % 	CD: CD value interpolated from stored aerodynamic data
           %
           % Called functions:
           %	GetIPValue(value) of class MDT
           
           if nargin == 2
               if size(this.CL.InputBreakpoints,2) == 2
                    [Cm Clipped] = this.Cm.GetIPValue([alpha 0]);
               else
                    [Cm Clipped] = this.Cm.GetIPValue(alpha);
               end
           elseif nargin == 3
               
               [Cm Clipped] = this.Cm.GetIPValue([alpha beta]);

           end
       end
       
       function Cn = getCnValue(this, alpha, beta)
           % Returns the corresponding yawing moment coefficient value to a given angle of attack (and yaw 
           % angle).
           % 
           % Input variables:
           % 	alpha: angle of attack in radiants
           % 	beta (optional): yaw angle in radiants (is set to zero if not given)
           % 	
           % Output variables:
           % 	Cn: Cn value interpolated from stored aerodynamic data
           %
           % Called functions:
           %	GetIPValue(value) of class MDT
           
           if nargin == 2
               if size(this.CL.InputBreakpoints,2) == 2
                    [Cn Clipped] = this.Cn.GetIPValue([alpha 0]);
               else
                    [Cn Clipped] = this.Cn.GetIPValue(alpha);
               end
           elseif nargin == 3
               
               [Cn Clipped] = this.Cn.GetIPValue([alpha beta]);

           end
       end
       
       function addTrimmedPolar(this, polar)
           % Setter for the property trimmedPolar.
           % 
           % Input variables:
           %	polar: struct containing the aerodynamic data
           
           if ~isempty(this.trimmedPolar)
               index2store = this.findTrimmedPolar(polar.ID);
               if isempty(index2store)
                    this.trimmedPolar(size(this.trimmedPolar,2)+1) = polar;
               else
                    this.trimmedPolar(index2store) = polar;
               end
           else
               this.trimmedPolar = polar;
           end
       end
       
       function deleteTrimmedPolar(this, index)
           
           this.trimmedPolar(index) = [];
       end
       
       function CD = getCDfromCLTrimmed(this, CL, polarID)
           
           if nargin < 3
               index = this.findTrimmedPolar('default');
               if isempty(index)
                   index = 1;
               end
           else
               if isnumeric(polarID)
                   index = polarID;
               elseif ischar(polarID)
                    index = this.findTrimmedPolar(polarID);
               else
                   index = 1;
               end
           end
           [~, I] = max(this.trimmedPolar(index).CL);
           CD = interp1q(this.trimmedPolar(index).CL(1:I)', this.trimmedPolar(index).CD(1:I)', CL);
       end
       
       function CL = getCLfromAlphaTrimmed(this, alpha, polarID)
           
           if nargin < 3
               index = this.findTrimmedPolar('default');
               if isempty(index)
                   index = 1;
               end
           else
               index = this.findTrimmedPolar(polarID);
           end
           CL = interp1q(this.trimmedPolar(index).alpha', this.trimmedPolar(index).CL', alpha);
       end
       
       function CLMax = getCLMaxTrimmed(this)
          
           if isempty(this.trimmedPolar)
               disp('No polar data available.');
           else
               if nargin < 3
                   index = this.findTrimmedPolar('default');
                   if isempty(index)
                       index = 1;
                   end
               else
                   index = this.findTrimmedPolar(polarID);
               end
               CLMax = max(this.trimmedPolar(index).CL);
           end
       end
       
       function alpha = getAlphaFromCLTrimmed(this, CL_given, polarID)
           % Returns the angle of attack corresponding to a given CL value of the trimmed configuration.
           % 
           % Input variables: 
           % 		CL_given: Given CL value
           % 		
           % Output variables:
           % 		alpha: Angle of attack in radiants
           if nargin < 3
               index = this.findTrimmedPolar('default');
               if isempty(index)
                   index = 1;
               end
           else
               if isnumeric(polarID)
                   index = polarID;
               elseif ischar(polarID)
                    index = this.findTrimmedPolar(polarID);
               else
                   fprintf('No matching polar found! First polar used.\n');
                   index = 1;
               end
           end
           
           [~, I] = max(this.trimmedPolar(index).CL);
           alpha = interp1q(this.trimmedPolar(index).CL(1:I)', this.trimmedPolar(index).alpha(1:I)', CL_given);
       end
       
       function index = findTrimmedPolar(this, ID)
          
           index = [];
           for i = 1:size(this.trimmedPolar, 2)
               if strcmp(this.trimmedPolar(i).ID, ID)
                   index = i;
                   break;
               end
           end
       end
       
   end %methods public
   
   methods (Access = private)
   
       
   end
    %                                   																								 *
    %*                                                                                                          *
    %************************************************************************************************************
end
%=============================================================================================================EOF
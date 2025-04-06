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
%*                                            AerodynamicGeometry		                                        *
%*                                                                                                              *
%****************************************************************************************************************
%*                                                                                                              *
%* Type        :   Abstract Class Definition                                                                    *
%*                                                                                                              *
%* Circulation :                                                                    							*
%*                                                                                                              *
%* Purpose     :   Definition of abstract class AerodynamicGeometry                                             *
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


classdef AerodynamicGeometry < Geometry
    %************************************************************************************************************
%%  %*                                         Properties                                                       *
    %************************************************************************************************************
    %*  Stadlberger, Korbinian           2010-JUN-22                                                            *
    %************************************************************************************************************
    %*                                                                      %                                   *
    %*  Property                		Description                        Type       					Dim         Unit* 
    properties (SetAccess = protected)%                                                            				 *
        Sref                        % Reference area                   double     					(1x1)       [m²]*
        Lref                        % Refernece Length                 double     					(1x1)       [m] *
        aeroProp                    % Aerodynamic Properties           AerodynamicProperties        (1x1)       [-] *
        refPoint_mom                % Reference point for moments      double                       (3x1)       [m]
    end
    %*                                                                      %                                   *
    %************************************************************************************************************



    %************************************************************************************************************
%%  %*                                          Methods                                                         *
    %************************************************************************************************************
    %*  Stadlberger, Korbinian           2010-JUN-21                                                            *
    %************************************************************************************************************
    %*                                                                                                          *
    methods (Abstract)
        
        refreshSref(this);
        % This function is meant to refresh the reference area of the aerodynamic geometry / to 
        % update the property 'Sref' which is dependent of other properties changing during runtime.
        
        refreshLref(this);
        % This function is meant to refresh the reference length of the aerodynamic geometry / to 
        % update the property 'Lref' which is dependent of other properties changing during runtime.
        
    end % methods abstract   
    
    methods 
        
        function this = AerodynamicGeometry()
        		% Constructor of abstract class AerodynamicGeometry which initialises its properties.
        		%
        		% Called functions:
        		%	AerodynamicProperties() of class 'AerodynamicProperties'
        		
            this = this@Geometry();
            this.Sref = 0;
            this.Lref = 0;
            this.aeroProp = []; % an AerodynamicProperties object might be used
            this.refPoint_mom = this.refPoint;
        end
        
        function this = setAeroProp(this, aeroProp) 
            % Setter for the property aeroProp which contains aerodynamic data.
            
            this.aeroProp = aeroProp;
        end 
        
        function setRefPointMom(this, newRefPoint)
            
            this.refPoint_mom = newRefPoint;
        end
    end % methods public
    %                                   *
    %*                                                                      %                                   *
    %************************************************************************************************************
end
%=============================================================================================================EOF
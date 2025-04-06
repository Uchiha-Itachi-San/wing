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
%*                                               Geometry		                                                *
%*                                                                                                              *
%****************************************************************************************************************
%*                                                                                                              *
%* Type        :   Abstract Class Definition                                                                    *
%*                                                                                                              *
%* Circulation :                                                                    							*
%*                                                                                                              *
%* Purpose     :   Definition of abstract class Geometry                                                        *
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


classdef Geometry < handle
    %************************************************************************************************************
%%  %*                                         Properties                                                       *
    %************************************************************************************************************
    %*  Stadlberger, Korbinian           2010-JUN-22                                                            *
    %************************************************************************************************************
    %*                                                                                                          *
    %*      Property                Description                                       Type       Dim        Unit* 
    properties (SetAccess = protected)%                                                                         *
        name                        % Name of the geometry                            string     (1x1)      [-]
        description                 % Description of the geometry                     string     (1x1)      [-]
        refPoint                    % ReferencePoint                                  double     (3x1)      [m]
        mass                        % Geometry mass                                   double     (1x1)      [kg]
        relCG                       % Relative position of c.g. according             double     (3x1)      [m]                      
        							% to reference point   									  
        size                        % Vector defining the size [L, W, H]				  double     (1x3)      [m]
        type                        % String defining the geometry type					  string		 (1x1)		[-]
        mcg_display                 % Struct containing information for list generation struct (1x1) [-]
        
    end
    %*                                                                                                          *
    %************************************************************************************************************



    %************************************************************************************************************
%%  %*                                          Methods                                                         *
    %************************************************************************************************************
    %*  Stadlberger, Korbinian           2010-JUN-21                                                            *
    %************************************************************************************************************
    %*                                                                                                          *
    methods (Abstract)
        
        refreshMass(this);
        % This function is meant to refresh the geometry's mass / to update the property 'mass' 
        % which is dependent on other properties changing during runtime.
        
        refreshRelCG(this);
        % This function is meant to refresh the geometry's relative center of gravity position /
        % to update the property 'relCG' which is dependent of other properties changing during runtime.
        
        refreshSize(this);
        % This function is meant to refresh the geometry's size / to update the property 'size' 
        % which is dependent on other properties changing during runtime.
        
        plotGeometry(this, shiftVector);
        % This function is meant to plot the geometry. If no shiftVector is used the geometry is 
        % plotted with its reference point at [0;0;0]. Otherwise the geometry is plotted with its 
        % reference point shown at the shiftVector position.
        
    end % methods abstract
    
    
    methods
        
        function this = Geometry()
            % Constructor of abstract class 'Geometry'. Can only be called by inherited classes.
            
            this.refPoint = [0;0;0];
            this.mass = 0;
            this.relCG = [0;0;0];
        end
        
        function setName(this, newName)
            % Setter of the property name.
            %
            % Input variables:
            %	newName: string to be assigned as new name
            this.name = newName;
        end
        
        function setDescription(this, newDescription)
            
            this.description = newDescription;
        end
        
        function setRefPoint(this, newRefPoint)
            % Setter of the property refPoint.
            %
            % Input variables:
            %	newRefPoint: column vector (3x1) defining the cartesian coordinates of the geometry's  
            %  reference point in m
            
            if size(newRefPoint, 1) == 3 && size(newRefPoint, 2) == 1
                this.refPoint = newRefPoint;
            else
                error('Dimension mismatch: new reference point must have dimension 3x1');
            end
        end
        
        function setMass(this, newMass)
            % Setter of the property mass.
            %
            % Input variables:
            %	newMass: mass of geometry in kg
            
            this.mass = newMass;
        end
        
        function setRelCG(this, newRelCG)
            % Setter of the property relCG.
            %
            % Input variables:
            %	newRelCG: column vector (3x1) defining the cartesian coordinates of the geometry's  
            %  center of gravity according to its reference point in m
            
            this.relCG = newRelCG;
        end
        
        function setSize(this, newSize)
            % Setter of the property size.
            %
            % Input variables:
            %	newSize: row vector (1x3) defining the size in its three dimensions in m  

            this.size = newSize;
        end
        
        function setHeight(this, height)
            % Setter of the geometry's heigth stored in the property size.
            %
            % Input variables:
            %	height: heigth in m
            
           this.size(3) = height; 
       end
       
       function setWidth(this, width)
            % Setter of the geometry's width stored in the property size.
            %
            % Input variables:
            %	width: width in m
            
           this.size(2) = width;
       end
       
       function setLength(this, length)
            % Setter of the geometry's length stored in the property size.
            %
            % Input variables:
            %	length: length in m
            
           this.size(1) = length;
       end
       
       function setDiameter(this, diameter)
           
           this.size(2) = diameter;
           this.size(3) = diameter;
       end
       
       function setType(this, newType)
           
           this.type = newType;
       end
       
       function info = getInfoMassesAndCG(this, option)
            
            J = this.getInertia();
           
            data =   {this.name,...
                      this.mass,...
                      this.relCG(1),...
                      this.relCG(2),...
                      this.relCG(3),...
                      this.refPoint(1),...
                      this.refPoint(2),...
                      this.refPoint(3),...
                      this.relCG(1),...
                      this.relCG(2),...
                      this.relCG(3),...
                      J(1),...
                      J(2),...
                      J(3),...
                      []};
                  
            if ~isempty(this.mcg_display)
               
                nObjects = size(this.mcg_display.relevantObjects, 1);
                %tempinfo = cell(nObjects, size(data, 2));
                subitems = [];
                for i = 1:nObjects
                   
                    % Get data line
                    tempinfo = this.mcg_display.relevantObjects{i}.getInfoMassesAndCG(1);
                    for j = 1:size(tempinfo, 1)
                        % Indent subitems
                        tempinfo{j,1} = strcat([this.mcg_display.indentString, tempinfo{j,1}]);
                        % Add reference point to absolute c.g. position
                        tempinfo{j,3} = tempinfo{j,3} + this.mcg_display.relevantObjects{i}.refPoint(1);
                        tempinfo{j,4} = tempinfo{j,4} + this.mcg_display.relevantObjects{i}.refPoint(2);
                        tempinfo{j,5} = tempinfo{j,5} + this.mcg_display.relevantObjects{i}.refPoint(3);
                    end
                    subitems = vertcat(subitems, tempinfo);
                end
                data = vertcat(data, subitems);
            end
                      
            if nargin == 1 || option == 0
                
                header = {'Component',...
                          'Mass [kg]',...
                          'X position c.g. (absolute) [m]',...
                          'Y position c.g. (absolute) [m]',...
                          'Z position c.g. (absolute) [m]',...
                          'X position reference point [m]',...
                          'Y position reference point [m]',...
                          'Z position reference point [m]',...
                          'X position c.g. (relative) [m]',...
                          'Y position c.g. (relative) [m]',...
                          'Z position c.g. (relative) [m]',...
                          'Inertia around x-axis [kgm²]',...
                          'Inertia around y-axis [kgm²]',...
                          'Inertia around z-axis [kgm²]',...
                          'Comments'};
                      
                info = vertcat(header, data);
                dispnxncellarray(info, 8);
            else
                info = data;
            end
       end
        
       function J = getInertia(this, point)
           
            % This function calculates and returns the moments of inertia around the mid point of
            % the geometry defined by the property 'size'. For a homogeniously distributed density
            % this point is equal to the centre of gravity. 
            %
            % Input (optional):
                %   point: vector (3x1) to specify the rotation point
                
            % Called functions:
				%	refreshMass: refreshing the mass of this object
                %	refreshSize: refreshing the size of this object
            %    
            % Output:
                %   J: vector (3x1) containing the moments of inertia around the three axis
                
            %this.refreshMass;
            this.refreshSize;
            
            J = zeros(3,1);
            
            if strcmp(this.type, 'cubic')
                J(1) = (this.size(2)^2 + this.size(3)^2) * this.mass / 12;
                J(2) = (this.size(1)^2 + this.size(3)^2) * this.mass / 12;
                J(3) = (this.size(1)^2 + this.size(2)^2) * this.mass / 12;
            elseif strcmp(this.type, 'cylindric')
                J(1) = (0.5 * this.size(2))^2 * this.mass * 0.5;
                J(2) = this.mass * ((0.25 * this.size(2))^2 * 0.25 + this.size(1)^2 / 12);
                J(3) = this.mass * ((0.25 * this.size(3))^2 * 0.25 + this.size(1)^2 / 12);
            else
                fprintf('%s: Unknown geometry type "%s" to calculate inertia!', this.name, this.type);
            end
            
            if nargin > 1
                delta = point - this.relCG;
                dist = [sqrt(delta(2)^2 + delta(3)^2);
                        sqrt(delta(1)^2 + delta(3)^2);
                        sqrt(delta(1)^2 + delta(2)^2)];
                J = J + this.mass * dist.^2;
            end
       end

    end 

    methods (Static)
        
       function varargout = axelRot(varargin)
        %Generate roto-translation matrix for the rotation around an arbitrary line in 3D.
        %The line need not pass through the origin. Optionally, also, apply this
        %transformation to a list of 3D coordinates.
        %
        %SYNTAX 1:
        %
        %    M=AxelRot(deg,u,x0)
        %
        %
        %in:
        %
        %  u, x0: 3D vectors specifying the line in parametric form x(t)=x0+t*u 
        %         Default for x0 is [0,0,0] corresponding to pure rotation (no shift).
        %         If x0=[] is passed as input, this is also equivalent to passing
        %         x0=[0,0,0].
        %
        %  deg: The counter-clockwise rotation about the line in degrees. 
        %       Counter-clockwise is defined using the right hand rule in reference
        %       to the direction of u.
        %
        %
        %out:
        %
        % M: A 4x4 affine transformation matrix representing
        %    the roto-translation. Namely, M will have the form
        % 
        %                 M=[R,t;0 0 0 1] 
        %   
        %    where R is a 3x3 rotation and t is a 3x1 translation vector. 
        % 
        %
        %
        %SYNTAX 2:
        %
        %       [R,t]=AxelRot(deg,u,x0)
        %
        % Same as Syntax 1 except that R and t are returned as separate arguments.
        % 
        %
        %
        %SYNTAX 3: 
        %
        % This syntax requires 4 input arguments be specified, 
        %
        %   [XYZnew, R, t] = AxelRot(XYZold, deg, u, x0)
        % 
        % where the columns of the 3xN matrix XYZold specify a set of N point
        % coordinates in 3D space. The output XYZnew is the transformation of the
        % columns of XYZold by the specified rototranslation about the axis. All 
        % other input/output arguments are as before.
        %
        %   by Matt Jacobson
        %  
        %   Copyright, Xoran Technologies, Inc. 2011
        
        
        if nargin>3
            
           XYZold=varargin{1};
           varargin(1)=[];
            
           [R,t]=this.axelRot(varargin{:});
            
           XYZnew=bsxfun(@plus,R*XYZold,t);
           
           varargout={XYZnew, R,t};
           
           return; 
           
        end
        
            [deg,u]=deal(varargin{1:2});
            
            if nargin>2, x0=varargin{3}; end
        
            R3x3 = nargin>2 && isequal(x0,'R');
        
            if nargin<3 || R3x3 || isempty(x0), 
                x0=[0;0;0]; 
            end
        
            x0=x0(:); u=u(:)/norm(u);
        
            AxisShift=x0-(x0.'*u).*u;
        
        
        
        
            Mshift=mkaff(eye(3),-AxisShift);
        
            Mroto=mkaff(R3d(deg,u));
        
            M=inv(Mshift)*Mroto*Mshift;
        
            varargout(1:2)={M,[]};
            
            if R3x3 || nargout>1 
              varargout{1}=M(1:3,1:3);
            end
            
            if nargout>1,
              varargout{2}=M(1:3,4);  
            end

        function R=R3d(deg,u)
        %R3D - 3D Rotation matrix counter-clockwise about an axis.
        %
        %R=R3d(deg,axis)
        %
        % deg: The counter-clockwise rotation about the axis in degrees.
        % axis: A 3-vector specifying the axis direction. Must be non-zero
        
            R=eye(3);
            u=u(:)/norm(u);
            x=deg; %abbreviation
        
            for ii=1:3
        
                v=R(:,ii);
        
                R(:,ii)=v*cosd(x) + cross(u,v)*sind(x) + (u.'*v)*(1-cosd(x))*u;
                  %Rodrigues' formula

            end
        end

        function M=mkaff(varargin)
        % M=mkaff(R,t)
        % M=mkaff(R)
        % M=mkaff(t)
        %
        %Makes an affine transformation matrix, either in 2D or 3D. 
        %For 3D transformations, this has the form
        %
        % M=[R,t;[0 0 0 1]]
        %
        %where R is a square matrix and t is a translation vector (column or row)
        %
        %When multiplied with vectors [x;y;z;1] it gives [x';y';z;1] which accomplishes the
        %the corresponding affine transformation
        %
        % [x';y';z']=R*[x;y;z]+t
        %
        
        
        
            if nargin==1
        
               switch numel(varargin{1}) 
        
                   case {4,9} %Only rotation provided, 2D or 3D
        
                     R=varargin{1}; 
                     nn=size(R,1);
                     t=zeros(nn,1);
        
                   case {2,3}
        
                     t=varargin{1};
                     nn=length(t);
                     R=eye(nn); 
        
               end
            else
        
                [R,t]=deal(varargin{1:2});
                nn=size(R,1);
            end
        
            t=t(:); 
        
            M=eye(nn+1);
        
            M(1:end-1,1:end-1)=R;
            M(1:end-1,end)=t(:); 
        end
       end
    end
    %                                   																								 *
    %*******************************************************************************************************
end
%=============================================================================================================EOF
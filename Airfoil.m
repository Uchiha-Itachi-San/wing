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
%*                                                Airfoil	                                                    *
%*                                                                                                              *
%****************************************************************************************************************
%*                                                                                                              *
%* Type        :   Abstract Class Definition                                                                    *
%*                                                                                                              *
%* Circulation :                                                                    							*
%*                                                                                                              *
%* Purpose     :   Definition of abstract class Airfoil                                                         *
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
%*  Stadlberger, Korbinian  * 2024-NOV-04 *     Initial Release (THI)
%*  Stadlberger, Korbinian  * 2025-FEB-07 *     Functions added for THI lecture
%****************************************************************************************************************


classdef Airfoil < handle
    %AIRFOIL Simple respresentation of an airfoil
    %   By now only contains the name of the airfoil
    
    properties
        name            % generated airfoil name
        type            % airfoil type: 'NACA', 'Supersonic', 'Custom'
        geoCode         % geometry code
        maxThickness
        posMaxThickness
        coordinates
    end
    
    methods
        function this = Airfoil(type, geoCode)
            %AIRFOIL Construct an instance of this class
            %   
            if any(strcmp(type, {'NACA', 'Supersonic', 'Custom'}))
                this.type = type;
            else
                error('Airfoil type must be "NACA", "Supersonic", or "Custom"')
            end
            
            if ~strcmp(this.type, 'Custom')
                this.geoCode = geoCode;
                this.name = [this.type, '-', this.geoCode];
                if strcmp(this.type, 'NACA') && (length(this.geoCode) == 4 || length(this.geoCode) == 5)
                    
                    [Upper, Lower,~,~] = this.generateNACA45Coordinates(['NACA', strrep(this.geoCode, 'NACA', '')], 1, 200, 0);
                    this.coordinates.upperPointArray = Upper;
                    this.coordinates.lowerPointArray = Lower;
                    this.coordinates.contourPointArray = [flipud(Upper);Lower(2:end,:)];
                    this.coordinates.camberLinePointArray = this.getCamberLine();
                    [this.maxThickness, this.posMaxThickness] = this.getMaximumThickness();
                end
            else
                this.name = 'Custom data point airfoil';
            end
        end

        function plotAirfoil(this)

                figureA = ['Plot of airfoil contour ', this.name];
                handleA = findobj('type', 'figure', 'Name', figureA);
                if isempty(handleA)
                    handleA = figure('Name', figureA);
                end
                figure(handleA);
                clf
                hold off
    
                plot(this.coordinates.contourPointArray(:,1), this.coordinates.contourPointArray(:,2), 'Marker','o');
                grid on;
                xlabel('x/c [-]');
                ylabel('y/c [-]');
    
                hold on
                plot(this.coordinates.camberLinePointArray(:,1), this.coordinates.camberLinePointArray(:,2), 'Marker','o');
                axis equal
        end

        function importAirfoilCoordinates(this, filepath)
           % Import airfoil coordinates from text file with Selig format
           if nargin == 1
               filepath = [];
           end
           
           % if no path is given or file cannot be found a input window will
           % be opened
           if isempty(filepath) || exist(filepath, 'file') == 0
               % Input window
               [fileName, pathName, ~] = uigetfile('*.dat','Select the text file containing the normalised airfoil coordinates');
               filepath = [pathName, fileName];
           end
           
           fileID = fopen(filepath);
           
           tempLine = fgetl(fileID);
           % if there is no name string in the first line
           if isempty(str2num(tempLine))
               this.name = strtrim(tempLine);
               tempLine = fgetl(fileID);
           end
           
           % read coordinates
           i = 1;
           while tempLine ~= -1
                data(i,:) = str2num(tempLine);
                i = i + 1;
                tempLine = fgetl(fileID);
           end
           
            dataStruct.contourPointArray = data;

           % find point at nose
           indexNose = find(data(:,1) == 0);
           nPoints = size(data, 1);
           
           
           if data(indexNose+1, 2) > 0

               dataStruct.upperPointArray = [data(indexNose:nPoints,1), data(indexNose:nPoints,2)];
               dataStruct.lowerPointArray = [data(indexNose:-1:1,1), data(indexNose:-1:1,2)];
           else
               dataStruct.upperPointArray = [data(indexNose:-1:1,1), data(indexNose:-1:1,2)];
               dataStruct.lowerPointArray = [data(indexNose:nPoints,1), data(indexNose:nPoints,2)];
           end
           
           this.coordinates = dataStruct;
           this.coordinates.camberLinePointArray = this.getCamberLine();
           [this.maxThickness, this.posMaxThickness] = this.getMaximumThickness();
           
           fclose(fileID);
        end

        function pointArray = getCamberLine(this, nPoints)

            if nargin < 2
                nPoints = 100;
            end

            % Interpolate the mean camber line
            x_camberLine = linspace(min(this.coordinates.contourPointArray(:,1)), max(this.coordinates.contourPointArray(:,1)), nPoints)';
            y_camberLine = 0.5 * (interp1(this.coordinates.upperPointArray(:,1), this.coordinates.upperPointArray(:,2), x_camberLine, 'pchip','extrap') ...
                                + interp1(this.coordinates.lowerPointArray(:,1), this.coordinates.lowerPointArray(:,2), x_camberLine, 'pchip','extrap'));

            pointArray = [x_camberLine, y_camberLine];
        end

        function [t_c_max, x_t_max] = getMaximumThickness(this)

            if ~isempty(this.coordinates)
                nPoints = 100;
                x_thickness = linspace(min(this.coordinates.contourPointArray(:,1)), max(this.coordinates.contourPointArray(:,1)), nPoints)';
                t_c = (interp1(this.coordinates.upperPointArray(:,1), this.coordinates.upperPointArray(:,2), x_thickness, 'pchip','extrap') ...
                    - interp1(this.coordinates.lowerPointArray(:,1), this.coordinates.lowerPointArray(:,2), x_thickness, 'pchip','extrap')) ...
                    / max(this.coordinates.contourPointArray(:,1));
            
                [t_c_max, indexMax] = max(t_c);
                x_t_max = x_thickness(indexMax);
            else
                error('No coordinate data available.');
            end
        end

    end

    methods (Static)

        function [Upper, Lower,dx,dy] = generateNACA45Coordinates(NACA, c, N, varargin )
            %NACA_Airfoil: Obtains the characteristic for a NACA 4-digit and 5-digit 
            %               series and plot the shape with the cord length.
            %               
            %           Ex. [Upper, Lower,dx,dy] = airfoil_generator( NACA,c,N,option)
            %               
            %               Where the first inpute is the NACA followed by a 4 or 5 
            %               number digit as a string, NACAXXXX. The second input is the
            %               cord length of the airfoil, and the third is the number of 
            %               points to plot the airfoil. The options are optional 
            % 
            %               first option  | plot          | 1 yes (default), 0 no
            %               second option | trailing edge | 'close' (default), 'open'
            % 
            % Ex:
            %   NACA_Airfoil( 'NACA0012',2,500 )
            % is a NACA 0012 of cord length 2 meters with 500 equdistance points.
            % 
            % ------------------------------------------------------------------------
            % Author: Cesar Galan
            % Date created: 09/30/2015
            % Last date modified: 09/30/2015
            % 
            %% Error
            L = length(NACA);
            %% Extract data
            if L == 8
            % 4-digit NACA airfoil
            m = str2double(NACA(5))/100; % maximun chamber
            p = str2double(NACA(6))/10; % maximun chamber along the cord
            t = str2double(NACA(7:8))/100; % max thickness
            elseif L == 9
            % 5-digit NACA airfoil
            cl = str2double(NACA(5))*3/20; % design coefficeint of lift
            p = str2double(NACA(6))/20; % position of maximun chamber
            q = str2double(NACA(7)); % check if it is a reflex camber
            t = str2double(NACA(8:9))/100; % max thickness
            end
            % x = linspace(0,1,N);
            % Modification by K. Stadlberger: Refinement at the leading edge
            x = 1 - cos(linspace(0,0.5*pi,N));
            % Checks for trailing edge close or open and if it wants to be ploted
            if ~isempty(varargin)
                PL = varargin{1};
                if length(varargin) == 2
                    TE = varargin{2};
                elseif length(varargin) > 2
                    error('Too manny input argumets')
                else
                    TE = '';
                end
            else
                PL = 1;
                TE = '';
            end
            %% Shape of mean camber
            a0 = 0.2969;
            a1 = -0.1260;
            a2 = -0.3516;
            a3 = 0.2843;
            if strcmp(TE,'') || strcmp(TE,'close')
                a4 = -0.1036; % if close
            elseif strcmp(TE,'open')
                a4 = -.1015; % if open
            end
            yt = (t/0.2).*(a0.*sqrt(x) + a1.*(x) + a2.*(x).^2 +...
                a3.*(x).^3 + a4.*(x).^4);
            %% surface of the camber
            if L == 8
                [yc,zeta] = camber('4digit',p,x,N,m);
            elseif L == 9
                [yc,zeta] = camber('5digit',p,x,N,q,cl);
            end
            x_U = c.*(x - yt.*sin(zeta));
            x_L = c.*(x + yt.*sin(zeta));
            y_U = c.*(yc + yt.*cos(zeta));
            y_L = c.*(yc - yt.*cos(zeta));
            %% output
            dy = diff(yc);
            dx = diff(x);
            Upper = [x_U; y_U];
            Lower = [x_L; y_L];
            Upper = Upper';
            Lower = Lower';

            %% Plot
            if PL == 1
            figureA = ['Plot of Airfoil Geometry NACA', NACA];
            handleFigure = findobj('type', 'figure', 'Name', figureA);
            if isempty(handleFigure)
                handleFigure = figure('Name', figureA);
            end
            figure(handleFigure);
            clf
            plot(x_U,y_U,'b'); hold on
            plot(x_L,y_L,'b');
            title([NACA(1:4) ' ' NACA(5:end)])
            % domain
            % xmin = -.1;
            % xmax = c*.1+c;
            % ymin = -2*max(Upper(:,2));
            % ymax = 2*max(Upper(:,2));
            % axis([xmin xmax ymin ymax])
            axis equal
            end
            
            %% Compute the camber for different airfoils
            function [yc,zeta] = camber(naca,p,x,N,varargin)
            I = find(x > p,1);
            yc = zeros(1,N);
            dyc_dx = zeros(1,N);
            if strcmp(naca,'4digit') % 4-digit equation
                m = varargin{1};
                
                if m ~= 0
                    yc(1:I) = (m/p^2)*(2*p.*x(1:I) - x(1:I).^2);
                    yc(I+1:end) = (m/(1-p)^2)*(1 - 2*p + 2*p.*x(I+1:end) - x(I+1:end).^2);
                    
                    dyc_dx(1:I) = (2*m/p^2)*(p-x(1:I));
                    dyc_dx(I+1:end) = (2*m/(1-p)^2)*(p-x(I+1:end));
                end
                
            elseif strcmp(naca,'5digit')  % 5-digit equation 
                q = varargin{1};
                cl = varargin{2};
                
                if q == 0 % Not reflex
                    [r,k1] = constants(p,q,cl);
                    yc(1:I) = (k1/6)*(x(1:I).^3 - 3*r.*x(1:I).^2 + x(1:I).*r^2*(3-r));
                    yc(I+1:end) = (k1*r.^3)./6 .* (1 - x(I+1:end));
                    
                    dyc_dx(1:I) = (k1/6).* (3*x(1:I).^2 - 6*r.*x(1:I) + r^2 *(3-r));
                    dyc_dx(I+1:end) = -k1*r.^3 ./ 6;
                    
                elseif q == 1 % reflected camber
                    [r,k1,k2] = constants(p,q,cl);
                    yc(1:I) = (k1/6)*((x(1:I)-r ).^3 - x(1:I)*k2*(1-r)^3 - ...
                        x(1:I)*r^3 + r^3);
                    yc(I+1:end) = (k1/6)*(k2*(x(I+1:end) - r).^3 - ...
                        x(I+1:end).*k2*(1-r)^3 - x(I+1:end).*r^3 + r^3);
                    
                    dyc_dx(1:I) = (k1/6).* (3.*(x(1:I)-r).^2 - k2*(1-r)^3 - r^3);
                    dyc_dx(I+1:end) = (k1/6).*(3.*k2.*(x(I+1:end) - r).^2 - k2.*(1 - r)^3 - r^3);
                else
                    error('not a valid NACA 5 digit airfoil');
                end
                
                
            end
            if dyc_dx ~= zeros(1,N);
                zeta = atan(dyc_dx);
            else
                zeta = 0;
            end
            end
            %% get the constanst for different airfoil
            function [m,k1,k2] = constants(p,q,cl)
            if q == 0 && nargout == 2
                A = [.05 0.0580 361.400 ;
                     .10 .1260  51.640 ;
                     .15 .2025  15.957 ;
                     .20 .2900  6.643  ;
                     .25 .3910  3.230] ;
                 
                 I = find(p == A(:,1),1);
                if isempty(I)
                    i = round(I/5);
                    A1 = A(i,:);
                    m = (p/A1(1))*A1(2);
                    k1 = (p/A1(1))*A1(3);
                else
                    m = A(I,2);
                    k1 = A(I,3);
                end
            elseif q == 1 && nargout == 3
                A = [.10 .1300 51.990 .000764;
                     .15 .2170 15.793 .00677;
                     .20 .3180 6.520  .0303;
                     .25 .4410 3.191  .1355];
                 
                 I = find(p == A(:,1),1);
                if isempty(I)
                    i = round(I/5);
                    A1 = A(i,:);
                    m = (p/A1(1))*A1(2);
                    k1 = (p/A1(1))*A1(3);
                else
                    m = A(I,2);
                    k1 = A(I,3);
                    k2 = A(I,4);
                end
            end
            % scaling 
            if cl ~= .3 && nargout == 2
                m = m*cl/.3;
                k1 = k1*cl/.3;
            elseif cl ~= .3 && nargout == 3
                m = m*cl/.3;
                k1 = k1*cl/.3;
                k2 = k2*cl/.3;
            end
            end
        end
    end
end


classdef apex
    %APEX Derive apex localization
    %   Apices are derived from box and probe coordinates.
    
    properties
        limit
        xI, zI;
        rI, thetaI;
    end
    
    methods
        function obj = apex(box, probe, thetaS)
            % Initialize apex object derived from box and probe coordinates
            % plus a steering direction.
            
            xE = probe.xE; zE = probe.zE;
            %rE = probe.rE; thetaE = probe.thetaE;
            xF = probe.xF; zF = probe.zF;
            %rF = probe.rF; thetaF = probe.thetaF;

            % Compute limiting directions
            thetaEC = atan( (box.xC-xE) / (box.zC-zE) );
            thetaFD = atan( (box.xD-xF) / (box.zD-zF) );

            obj.limit = 'Infinite';
            % Initialize apex in case of no intersection
            rApex  = 10e3;
            obj.xI        = box.xG*(1-rApex/box.rG*cos(thetaS))-box.zG*rApex/box.rG*sin(thetaS);
            obj.zI        = box.zG*(1-rApex/box.rG*cos(thetaS))+box.xG*rApex/box.rG*sin(thetaS);
            zILeft = obj.zI; zIRight = obj.zI; 
            xILeft = obj.xI; xIRight = obj.xI;

            if (box.thetaG+thetaS>thetaEC)
                obj.limit = 'Left';
                % Compute intersection point I coordinates
                Num_zI = -(zE*box.xC-box.zC*xE)*cos(box.thetaG+thetaS) + (box.zC-zE)*box.rG*sin(thetaS);
                Num_xI = -(zE*box.xC-box.zC*xE)*sin(box.thetaG+thetaS) + (box.xC-xE)*box.rG*sin(thetaS);
                Den_I  = (box.zC-zE)*sin(box.thetaG+thetaS) - (box.xC-xE)*cos(box.thetaG+thetaS);

                zILeft = Num_zI / Den_I; xILeft = Num_xI / Den_I;
            end

            if (box.thetaG+thetaS<thetaFD)
                obj.limit = 'Right';
                % Compute intersection point I coordinates
                Num_zI = -(zF*box.xD-box.zD*xF)*cos(box.thetaG+thetaS) + (box.zD-zF)*box.rG*sin(thetaS);
                Num_xI = -(zF*box.xD-box.zD*xF)*sin(box.thetaG+thetaS) + (box.xD-xF)*box.rG*sin(thetaS);
                Den_I  = (box.zD-zF)*sin(box.thetaG+thetaS) - (box.xD-xF)*cos(box.thetaG+thetaS);

                zIRight = Num_zI / Den_I; xIRight = Num_xI / Den_I;
            end


            ILeftGVect = [zILeft; xILeft]; %- [box.zG; box.xG];
            IRightGVect = [zIRight; xIRight]; %- [box.zG; box.xG];
            if  norm(ILeftGVect) < norm(IRightGVect)
                obj.limit = 'Left';
                obj.zI = zILeft; obj.xI = xILeft;
            elseif  norm(ILeftGVect) > norm(IRightGVect)
                obj.limit = 'Right';
                obj.zI = zIRight; obj.xI = xIRight;
            else
                % Nothing to be done
            end

            obj.rI = sqrt( obj.xI^2 + obj.zI^2 );
            obj.thetaI = 2*atan( obj.xI / (obj.zI + sqrt( obj.xI^2 + obj.zI^2 )) );
        
        end
       
        function [center, width] = computeTXAperture(obj, box, probe)

            Rp = probe.radius;
            thetaY = 0;
            
            % Check which line is limiting on the other side
            if strcmp(obj.limit, 'Left')
                thetaIB = atan( (obj.xI-box.xB) / (obj.zI-box.zB) );
                thetaID = atan( (obj.xI-box.xD) / (obj.zI-box.zD) );
                if thetaIB < thetaID
                    pointChar = 'D';
                else
                    pointChar = 'B';
                end
            elseif strcmp(obj.limit, 'Right')
                thetaIA = atan( (obj.xI-box.xA) / (obj.zI-box.zA) );
                thetaIC = atan( (obj.xI-box.xC) / (obj.zI-box.zC) );
                if thetaIA < thetaIC
                    pointChar = 'A';
                else
                    pointChar = 'C';
                end
            else
                clearvars pointChar;
            end

            if exist('pointChar', 'var')
                disp(strcat('pointChar=', pointChar));

                % Find intersection point of curved array with that other line
                xX = eval(strcat('box.x', pointChar)); zX = eval(strcat('box.z', pointChar));
                rX = eval(strcat('box.r', pointChar)); thetaX = eval(strcat('box.theta', pointChar));
                thetaIX = eval(strcat('thetaI', pointChar));
                rho = - (obj.xI*sin(thetaIX)+obj.zI*cos(thetaIX)) ...
                      + sqrt((obj.xI*sin(thetaIX)+obj.zI*cos(thetaIX))^2-(obj.rI^2-Rp^2));
                disp(strcat('imag(rho)>0 == ', num2str(imag(rho)>0)));

                zY = obj.zI + rho*cos(thetaIX); xY = obj.xI + rho*sin(thetaIX);
                rY = Rp; thetaY = 2*atan( xY / (rY + zY) );

                idY = floor( Rp*(probe.deltaTheta+thetaY)/system.probe.Pitch );
                idF = ceil( Rp*(probe.deltaTheta+probe.thetaF)/system.probe.Pitch );
            end

            center = Rp*(probe.deltaTheta+(thetaY+probe.thetaF)/2);
            width = Rp*abs(probe.thetaF-thetaY);
        end
        
        function [] = drawLines(obj, box, probe, thetaS, figure_id)

            Lines(1) = tools.line(probe.rE, probe.thetaE, box.rC, box.thetaC);
            Lines(2) = tools.line(probe.rF, probe.thetaF, box.rD, box.thetaD);
            Lines(3) = tools.line(box.rG, box.thetaG, nan, thetaS);

            % Check which line is limiting on the other side
            if strcmp(obj.limit, 'Left')
                thetaIB = atan( (obj.xI-box.xB) / (obj.zI-box.zB) );
                thetaID = atan( (obj.xI-box.xD) / (obj.zI-box.zD) );
                if thetaIB < thetaID
                    pointChar = 'D';
                    Lines(4) = tools.line(obj.rI, obj.thetaI, rD, thetaD);
                else
                    pointChar = 'B';
                    Lines(4) = tools.line(obj.rI, obj.thetaI, rB, thetaB);
                end
            elseif strcmp(obj.limit, 'Right')
                thetaIA = atan( (obj.xI-box.xA) / (obj.zI-box.zA) );
                thetaIC = atan( (obj.xI-box.xC) / (obj.zI-box.zC) );
                if thetaIA < thetaIC
                    pointChar = 'A';
                    Lines(4) = tools.line(obj.rI, obj.thetaI, box.rA, box.thetaA);
                else
                    pointChar = 'C';
                    Lines(4) = tools.line(obj.rI, obj.thetaI, box.rC, box.thetaC);
                end
            else
                clearvars pointChar;
                if length(Lines)>3
                    Lines(4) = [];
                end
            end

            if exist('pointChar', 'var')
                disp(strcat('pointChar=', pointChar));

                % Find intersection point of curved array with that other line
%                 xX = eval(strcat('x', pointChar)); zX = eval(strcat('z', pointChar));
%                 rX = eval(strcat('r', pointChar)); thetaX = eval(strcat('theta', pointChar));
                thetaIX = eval(strcat('thetaI', pointChar)); Rp = system.probe.Radius;
                rho = - (obj.xI*sin(thetaIX)+obj.zI*cos(thetaIX)) ...
                      + sqrt((obj.xI*sin(thetaIX)+obj.zI*cos(thetaIX))^2-(obj.rI^2-Rp^2));
                disp(strcat('imag(rho)>0 == ', num2str(imag(rho)>0)));

                zY = obj.zI + rho*cos(thetaIX); xY = obj.xI + rho*sin(thetaIX);
                rY = Rp; thetaY = 2*atan( xY / (rY + zY) );

%                 idY = floor( Rp*(deltaThetaProbe+thetaY)/system.probe.Pitch );
%                 idF = ceil( Rp*(deltaThetaProbe+thetaF)/system.probe.Pitch );

                [ZZ, XX] = pol2cart(linspace(thetaY, probe.thetaF, 100), repmat(rY, [1, 100]));
                figure(figure_id); plot(XX*1e-3, ZZ*1e-3, 'b');
                figure(figure_id); plot(xY*1e-3, zY*1e-3, 'b.', 'MarkerSize', 26);
            end

            for j=1:1:length(Lines)
                Lines(j).draw( figure_id );
            end
            
        end
        
        function [] = drawPoint(obj, figure_id)
            figure(figure_id); plot(obj.xI*1e-3, obj.zI*1e-3, 'r.', 'MarkerSize', 10);
        end
       
    end
    
end


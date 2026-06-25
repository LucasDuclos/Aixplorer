classdef polarbox
    %POLARBOX Polar coordinates obj (e.g. color or elasto obj)
    %   Curved arrays boxes need polar kind of boxes.
    
    properties
        r0
        theta0
        delta_r
        delta_theta
        radii
        angulae
        mesh
        Z
        X
        top
        bot
        middle
        xA, zA, rA, thetaA;
        xB, zB, rB, thetaB;
        xC, zC, rC, thetaC;
        xD, zD, rD, thetaD;
        xG, zG, rG, thetaG;
    end
    
    methods
        function obj = polarbox(r0, theta0, delta_r, delta_theta)
            % Prototype: polarbox(r0, theta0, delta_r, delta_theta)
            % Description: compute position arrays of a color obj
            obj.r0 = r0;
            obj.theta0 = theta0;
            obj.delta_r = delta_r;
            obj.delta_theta = delta_theta;
            
            % Radii array
            obj.radii = linspace(r0, r0+delta_r);
            % Anguli array
            obj.angulae = linspace( -delta_theta/2, delta_theta/2 ) + theta0;

            % Register edge/corner positions
            [obj.top.left.Theta, obj.top.left.R] = deal(min(obj.angulae), min(obj.radii));
            [obj.top.left.Z, obj.top.left.X] = pol2cart(min(obj.angulae), min(obj.radii));
            
            [obj.top.right.Theta, obj.top.right.R] = deal(max(obj.angulae), min(obj.radii));
            [obj.top.right.Z, obj.top.right.X] = pol2cart(max(obj.angulae), min(obj.radii));
            
            [obj.bot.left.Theta, obj.bot.left.R] = deal(min(obj.angulae), max(obj.radii));
            [obj.bot.left.Z, obj.bot.left.X] = pol2cart(min(obj.angulae), max(obj.radii));
            
            [obj.bot.right.Theta, obj.bot.right.R] = deal(max(obj.angulae), max(obj.radii));
            [obj.bot.right.Z, obj.bot.right.X] = pol2cart(max(obj.angulae), max(obj.radii));
           
            [obj.middle.Theta, obj.middle.R] = deal(mean(obj.angulae), mean(obj.radii));
            [obj.middle.Z, obj.middle.X] = pol2cart(mean(obj.angulae), mean(obj.radii));
            
            % Colorbox (shortcut) points for side coordinates
            % A is top left, B is top right
            % C is bot left, D is bot right
            obj.rA = obj.top.left.R;  obj.thetaA = obj.top.left.Theta;
            obj.zA = obj.top.left.Z;  obj.xA = obj.top.left.X;
            obj.rB = obj.top.right.R; obj.thetaB = obj.top.right.Theta;
            obj.zB = obj.top.right.Z; obj.xB = obj.top.right.X;
            obj.rC = obj.bot.left.R;  obj.thetaC = obj.bot.left.Theta;
            obj.zC = obj.bot.left.Z;  obj.xC = obj.bot.left.X;
            obj.rD = obj.bot.right.R; obj.thetaD = obj.bot.right.Theta;
            obj.zD = obj.bot.right.Z; obj.xD = obj.bot.right.X;
            % G is the colorbox center
            obj.rG = obj.middle.R; obj.thetaG = obj.middle.Theta;
            obj.zG = obj.rG*cos(obj.thetaG); obj.xG = obj.rG*sin(obj.thetaG);

            % Full obj mesh
            obj.mesh = cat(0, [obj.radii; obj.angulae(1)*ones(size(obj.radii))]', ...
                              [obj.radii(end)*ones(size(obj.angulae)); obj.angulae]', ...
                              [obj.radii(end:-1:1); obj.angulae(end)*ones(size(obj.radii))]', ...
                              [obj.radii(1)*ones(size(obj.angulae)); obj.angulae(end:-1:1)]');

            % Convert to cartesian coordinates
            [obj.Z, obj.X] = pol2cart(obj.mesh(:, 2), obj.mesh(:, 1));
        end
        
        function [hdl] = draw(obj, varargin)
            % Prototype: draw(obj, figure_id, figure_hdl)
            % Description: draw the color obj on the specified figure
            p = inputParser;

            %addRequired(p,'obj');
            addOptional(p, 'figure_id', 4);
            addOptional(p, 'figure_hdl', nan);

            parse(p, varargin{:});

            hdl = p.Results.figure_hdl;
            id = p.Results.figure_id;

            figure(id);
            if not(isnan(hdl))
                hold on;
                set(hdl, 'XData', obj.X*1e-3, 'YData', obj.Z*1e-3);
            else
                hdl = plot(obj.X*1e-3, obj.Z*1e-3, 'w');
            end
            set(gca,'color','black')

        end
        
        function [] = drawPoints(obj, figure_id)
            figure(figure_id); plot(obj.xA*1e-3, obj.zA*1e-3, 'r.', 'MarkerSize', 26);
            figure(figure_id); plot(obj.xB*1e-3, obj.zB*1e-3, 'g.', 'MarkerSize', 26);
            figure(figure_id); plot(obj.xC*1e-3, obj.zC*1e-3, 'r.', 'MarkerSize', 26);
            figure(figure_id); plot(obj.xD*1e-3, obj.zD*1e-3, 'g.', 'MarkerSize', 26);
            figure(figure_id); plot(obj.xG*1e-3, obj.zG*1e-3, 'w.', 'MarkerSize', 26);
        end
        
    end
    
end


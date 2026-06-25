classdef line
    %LINE Polar coordinates line
    %   Definition with a point and a direction or two points
    
    properties
        Z
        X
    end
    
    methods
        function obj = line(r0, theta0, varargin)
            % Prototype: line(r0, theta0, r1, theta1, 'Type', line_type)
            % Description: compute position arrays of a line 
            %              ('Type', 'Two-points') going through 2 points
            %              ('Type', 'Steering') going through 2 points

                % Input parsing
                p = inputParser;

                addRequired(p,'r0');
                addRequired(p,'theta0');
                addOptional(p, 'r1', nan);
                addOptional(p, 'theta1', nan);

                parse(p, r0, theta0, varargin{:});
                % Parameter parsing
                [r0, theta0] = deal(p.Results.r0, p.Results.theta0);
                [r1, theta1] = deal(p.Results.r1, p.Results.theta1);

                if not(isnan(r1))

                    % Function core
                    weights = linspace(-20, 20);

                    [X1, Y1] = pol2cart(theta0, r0);
                    [X2, Y2] = pol2cart(theta1, r1);

                    rawobj = (weights'*[X1, Y1] + (1-weights')*[X2, Y2]);
                    [obj.Z, obj.X] =  deal(rawobj(:, 1), rawobj(:, 2));

                else

                    % Function core
                    weights = linspace(-20, 20);

                    [X1, Y1] = pol2cart(theta0, r0);
                    [X2, Y2] = pol2cart(theta0+theta1, weights*(r0+1));

                    rawobj = (ones(1, length(weights))'*[X1, Y1] + [X2; Y2]');
                    [obj.Z, obj.X] =  deal(rawobj(:, 1), rawobj(:, 2));

                end

        end
        
        function [hdl] = draw(obj, varargin)
            % Prototype: draw(obj, figure_id, figure_hdl)
            % Description: draw a line on the specified figure
            p = inputParser;

            %addRequired(p, 'obj');
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
                hdl = plot(obj.X*1e-3, obj.Z*1e-3, 'w-.');
            end
            set(gca,'color','black')
            
        end
        
    end
    
end


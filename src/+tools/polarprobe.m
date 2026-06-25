classdef polarprobe
    %POLARPROBE Container for curved probe coordinates.
    
    properties
        radius
        deltaTheta
        rE, rF, thetaE, thetaF;
        xE, xF, zE, zF;
    end
    
    methods
        function obj = polarprobe(probeRadius, probeDeltaTheta)
            
            obj.radius = probeRadius;
            obj.deltaTheta = probeDeltaTheta;
            
            % Compute first and last piezo coordinates
            % E is the position of the first piezo
            obj.rE = probeRadius; obj.thetaE = -probeDeltaTheta;
            obj.zE = obj.rE*cos(obj.thetaE); obj.xE = obj.rE*sin(obj.thetaE);
            % F is the position of the last piezo
            obj.rF = probeRadius; obj.thetaF = probeDeltaTheta;
            obj.zF = obj.rF*cos(obj.thetaF); obj.xF = obj.rF*sin(obj.thetaF);
        end
        
        function [] = drawPoints(obj, figure_id)
            figure(figure_id); plot(obj.xE*1e-3, obj.zE*1e-3, 'b.', 'MarkerSize', 26);
            figure(figure_id); plot(obj.xE*1e-3, obj.zE*1e-3, 'r.', 'MarkerSize', 10);
            figure(figure_id); plot(obj.xF*1e-3, obj.zF*1e-3, 'b.', 'MarkerSize', 26);
            figure(figure_id); plot(obj.xF*1e-3, obj.zF*1e-3, 'g.', 'MarkerSize', 10);
        end
    end
    
end


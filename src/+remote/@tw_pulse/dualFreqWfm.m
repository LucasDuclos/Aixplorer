function Waveform_Digital = dualFreqWfm( obj, Nsamples, NbHcycle, NDCycle)
    %DIFFTHRESHOLD Summary of this function goes here
    %   Detailed explanation goes here

    cter = 0;
    step = 0.01;
    stop = 100;
    
    Waveform_Sin1 = sin((1 : Nsamples * NbHcycle) / Nsamples * pi);
    Waveform_Cos2 = cos(2*(1 : Nsamples * NbHcycle) / Nsamples * pi);

    Waveform_Analog = 0.5*Waveform_Sin1 + 0.5*Waveform_Cos2;

    threshold = 1 - NDCycle/Nsamples;
    minus_thresh = threshold;
    plus_thresh = threshold;
    
    Waveform_Digital = applyThresholds( Waveform_Analog, Nsamples, NbHcycle, minus_thresh, plus_thresh);

    plus_occur = length(find(Waveform_Digital==1));
    minus_occur = length(find(Waveform_Digital==-1));
    NdutyC = (plus_occur + minus_occur)/NbHcycle;
    
    while (abs(plus_occur-minus_occur)>1 || abs(NDCycle - NdutyC)/(Nsamples) > 0.02) && cter<stop
        
        if plus_occur < minus_occur
            if NdutyC > NDCycle
                minus_thresh = minus_thresh + step;
            elseif NdutyC < NDCycle
                plus_thresh = plus_thresh - step;
            else
                minus_thresh = minus_thresh + step/2;
                plus_thresh = plus_thresh - step/2;
            end
        elseif plus_occur > minus_occur
            if NdutyC > NDCycle
                plus_thresh = plus_thresh + step;
            elseif NdutyC < NDCycle
                minus_thresh = minus_thresh - step;
            else
                minus_thresh = minus_thresh - step/2;
                plus_thresh = plus_thresh + step/2;
            end
        end
        
        Waveform_Digital = applyThresholds( Waveform_Analog, Nsamples, NbHcycle, minus_thresh, plus_thresh);

        plus_occur = length(find(Waveform_Digital==1));
        minus_occur = length(find(Waveform_Digital==-1));
        NdutyC = plus_occur + minus_occur;
        
        cter = cter + 1;
        
    end
    
end


function Waveform_Digital = applyThresholds( Waveform_Analog, Nsamples, NbHcycle, minus_thresh, plus_thresh)
    %DIFFTHRESHOLD Summary of this function goes here
    %   Detailed explanation goes here
    
    Waveform_Digital = zeros(1, Nsamples * NbHcycle);
    
        % Discretize analog waveform into digital waveform
    for i=1:length(Waveform_Analog)
        if Waveform_Analog(i) > plus_thresh
            Waveform_Digital(i) = 1;
        elseif Waveform_Analog(i) < -minus_thresh
            Waveform_Digital(i) = -1;
        else
            Waveform_Digital(i) = 0;

        end
    end
    
    % Set beginning at the right phase
    if Waveform_Digital(end)==1
        i=1;
        while Waveform_Digital(end-i) == 1
            i = i+1;
        end
        Waveform_Digital = [Waveform_Digital(end-i:end) Waveform_Digital(1:end-i)];
    end
    
end

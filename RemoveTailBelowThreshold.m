%% Removes trailing values below threshold
function output_size = RemoveTailBelowThreshold(IR, input_size, threshold_dB)
    output_size = input_size;
    threshold = threshold_dB;
    threshold_mag = db2mag(threshold);
    for i = input_size:-1:1
       if abs(IR(i)) > threshold_mag
          output_size = i;
          break 
       end
    end


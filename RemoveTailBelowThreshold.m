%% Removes trailing values below threshold
function IR = RemoveTailBelowThreshold(IR, threshold_dB)
    threshold = threshold_dB;
    threshold_mag = db2mag(threshold);
    for i = size(IR,1):-1:1
       if abs(mean(IR(i,1,:))) > threshold_mag
          IR = IR(1:i,:);
          break 
       end
    end


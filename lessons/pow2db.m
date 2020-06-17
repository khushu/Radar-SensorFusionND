function ydb = pow2db(y)
%pow2db: This is quivalent signal processing function for non toolbox users
%Convert power to decibels, returns expresses in decibels (dB) the power 
% measurements specified in y. The relationship between power and decibels 
%is ydb = 10 log10(y).
    ydb = 10*log10(y);
end


function y = db2pow(ydb)
%db2pow: This is quivalent signal processing function for non toolbox users
% Convert decibels to power, returns the power measurements, y, that correspond to the decibel (dB)
% values specified in ydb. The relationship between power and decibels is
% ydb = 10 log10(y)
    y = 10.^(ydb/10);
end


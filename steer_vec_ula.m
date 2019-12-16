function s = steer_vec_ula(ula,lambda,angles)
%STREER_VEC Summary of this function goes here
%   Detailed explanation goes here
k = 2*pi/lambda;
txEl = ula.getNumElements;
d = ula.ElementSpacing;
s = zeros(txEl,length(angles));
for i = 1: length(angles)
    for n = 1:txEl
        s(txEl - n+1,i) = exp(-1i*k*(n-1) * d * sind(angles(i)) );
    end
end
end


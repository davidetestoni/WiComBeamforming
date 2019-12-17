function s = steer_vec_ula(ula,lambda,angles)
%STREER_VEC Summary of this function goes here
%   Detailed explanation goes here
k = 2*pi/lambda;
txEl = ula.getNumElements;
d = ula.ElementSpacing;
s = zeros(txEl,length(angles));
for j = 1: length(angles)
    for ii = 1:txEl
        n = (ii-1) - ((txEl-1)/2);
        s(txEl - (ii-1),j) = exp(-1i*k*n * d * sind(angles(j)) );
    end
end
end


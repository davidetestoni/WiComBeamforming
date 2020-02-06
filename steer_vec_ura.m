function [s] = steer_vec_ura(ura,lambda,angles)
%STREER_VEC Summary of this function goes here
%   Detailed explanation goes here

k = 2*pi/lambda;

d = ura.ElementSpacing(1);
%b_angles = az2broadside(angles(1,:),angles(2,:)); 
N_tx = ura.getNumElements;
s = zeros(N_tx,size(angles,2));

for j = 1:size(angles,2)
        n1 = (ura.Size(1) - 1)/2 : -1 : -(ura.Size(1) - 1)/2;
        
        s1 = exp(-1i * d  * k * sind(az2broadside(angles(1,j),angles(2,j))) .* n1);
        
        n2 = (ura.Size(2) - 1)/2 : -1 :-(ura.Size(1) - 1)/2;
        s2 = exp(-1i *  d * k * sind(angles(2,j)) .* n2 );        

        
        s0 = (s1' * s2)';
        s(:,j) = s0(:);
end


end


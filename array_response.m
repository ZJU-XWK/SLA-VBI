function [aR] = array_response(theta,N)
   n = 0:1:N-1;
   aR = exp(1j*pi*sin(theta)*n.')/sqrt(N);
end


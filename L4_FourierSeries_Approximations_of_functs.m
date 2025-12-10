%The following will be an attempt to represent various
%Functions using the Fourier series

clc, clf;
clear all, close all;

dx = 0.005;
L = 1;
x = dx:dx:L; %starts at approx 0 but not
%Another reason is because length is now 200 not 201
%And there would be double counting
f = ones(size(x));
%turn f from 1 from 0to L to be 0 from 0 to L/2 then 1 from L/2 to L
f(1:length(f)/2) = 0*f(1:length(f)/2);
plot(x,f)
%note could comment away lines 14 and 15 and replace with:
%f = sin(@*pi*x)
axis([-0.5, 1.5, -0.5, 1.5])
hold on

%Now for FS approximations:
fFS = zeros(size(x));
%DC offset
A0 = (2/L)*sum(f.*ones(size(x)))*dx;
%The 'm' stuff is for visual of seeing how u increase n the approx gets
%better
for m = 1:100 
    hold off
    plot(x,f)
    axis([-0.5, 1.5, -0.5, 1.5])
    
    %Initially starting at n = 0:
    fFS = A0/2;
    
    %Now for the series rep:
    
    for n = 1:m
        An = (2/L)*sum(f.*cos(2*pi*n*x/L))*dx;
        Bn = (2/L)*sum(f.*sin(2*pi*n*x/L))*dx;
        
        fFS = fFS + An*cos(2*pi*n*x/L) + Bn*sin(2*pi*n*x/L);    
    end
    hold on
    plot(x,fFS)
    legend("Original unit step", "FS approximate for n = " + m)

    drawnow
    pause(0.01)
end

%Note at m(n) = 100, we get an exact replica. This is due to aliasing due
%to the specific dx = 0.005 and L = 1 (ie 200 pts)

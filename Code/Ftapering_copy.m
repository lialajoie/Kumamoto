function Yf2 = Ftapering(Yf1,percent)
% Effectue un tapering du signal vers 0 à ses 2 extrémités, sur
% une distance de x 'percent' de la longueur totale du signal

% Perform signal tapering to 0 at both ends over a distance of x 'percent'
% of the total length of the signal

N       = length(Yf1);
M1      = ceil(N*percent);
xtap    = [1:1:N];
ANG     = pi/M1;
cs(1:M1) = (1.-cos(xtap(1:M1)*ANG))/2;
M2           = N-M1;
ANG          = pi/M1;
cs((M2+1):N) = fliplr((1.-cos(xtap(1:M1)*ANG))/2);
cs(M1+1:M2)  = 1;
cs2          = cs.*length(cs)./sum(cs);
[S1 S2]      = size(Yf1);
if(S1>S2)
    cs=cs';
end
cs           = cs.*length(cs)./sum(cs);
Yf2          = Yf1.*cs;

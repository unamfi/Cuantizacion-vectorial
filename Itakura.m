%% ITAKURA

function d = Itakura(y,Ry,x)
    d = log((x*Ry*transpose(x))/(y*Ry*transpose(y)));
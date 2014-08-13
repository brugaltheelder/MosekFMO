obj = 0;

diff = (dose-thresh);
over = diff;
under = - diff;
over(over<0) = 0;
under(under<0) = 0;

obj = sum(aOver.*(over.^2)+aUnder.*(under.^2))/2;
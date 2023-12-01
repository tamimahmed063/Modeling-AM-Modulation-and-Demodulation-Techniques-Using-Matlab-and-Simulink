function y = triangl(t)
y=(1-abs(t)).*(t>=-1).*(t<1); % i.e. setting y to 1 -|t|  if  |t|<1 and to 0 if not
end

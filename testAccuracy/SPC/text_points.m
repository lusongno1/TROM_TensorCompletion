function text_points(x,y,R)
sz = length(x);
if(size(x,2)>1)
    x = x.';
end
if(size(y,2)>1)
    y = y.';
end
if(size(R,2)>1)
    R = R.';
end
strs = [num2str(R)];
text(x,y,strs);
end
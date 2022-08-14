function text_point(x,y)
plot(x,y,'o');
sz = length(x);
if(size(x,2)>1)
    x = x.';
end
if(size(y,2)>1)
    y = y.';
end
strs = [repmat('(',sz,1) num2str(x) repmat(',',sz,1) num2str(y) repmat(')',sz,1)];
%strs = [num2str(y)];
text(x,y,strs);
end
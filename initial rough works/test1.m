data=[1,2,3,1;
      2,3,4,.5
      4,5,6,.3];
h=axes;hold on;
for k=1:size(data,1);        [x,y,z]=ellipsoid(data(k,1),data(k,2),data(k,3),data(k,4),data(k,4),data(k,4));
surf(h,x,y,z);
end
view(3);grid on;
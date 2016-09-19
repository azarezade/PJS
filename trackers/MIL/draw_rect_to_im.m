function im=draw_rect_to_im(im, state)
%   ---------------------------------> x
%   |
%   |
%   |
%   |
%   |
%   |
%   |
%   \/
%   y
rect_color=[25 255 255]; 

cx=state(1,1); 
cy=state(1,2); 
hw=state(1,3); 
hh=state(1,4);

%top left corner
xx1= cx;
yy1=cy;
 
%bottom right corner
xx2=cx+hw;  
yy2=cy+hh;

im(yy1, xx1:xx2,1)=rect_color(1);im(yy1, xx1:xx2,2)=rect_color(2);im(yy1, xx1:xx2,3)=rect_color(3);
im(yy2, xx1:xx2,1)=rect_color(1);im(yy2, xx1:xx2,2)=rect_color(2);im(yy2, xx1:xx2,3)=rect_color(3);
im(yy1:yy2,xx1,1)=rect_color(1);im(yy1:yy2,xx1,2)=rect_color(2);im(yy1:yy2,xx1,3)=rect_color(3);
im(yy1:yy2,xx2,1)=rect_color(1);im(yy1:yy2,xx2,2)=rect_color(2);im(yy1:yy2,xx2,3)=rect_color(3);
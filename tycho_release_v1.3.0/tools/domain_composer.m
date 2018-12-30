

%%%%%%   TYCHO DOMAIN COMPOSER 2D %%%%%%%%%%%% 
%This octave script allows you to define obstacles in a                   %
%computational domain with resolution defined at the beginning  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;

disp "========================"
disp "This is the 2D TYCHO DOMAIN COMPOSER"
disp "========================"


resolution_x = input("The Resolution in X direction: ", "s");
res_x = str2num(resolution_x);
resolution_y = input("The Resolution in Y direction: ", "s");
res_y = str2num(resolution_y);

form1 = input("[o]bstacle or [w]indemitter or [s]top", "s");

if (form1=="w")
   form2= input("In which direction does the wind blow? [2]..+x  [3]..-x [4]..+y [5]..-y [6]..+z [7]..-z"); 
endif

form = input("[r]ectangle, [l]ine, [c]ircle, [h]orizontal line, [v]ertical line", "s");

img=zeros(res_y,res_x);
ih=imagesc(img);
axis image;
set(gca,'YDir','normal')
caxis([0,1]);
colormap(jet(2));

while (form1 ~= "s")
	
if (form1 == "o")
    constant = 1;
elseif (form1 == "w")
    constant = form2;
endif 

[x,y,button] = ginput(2);

for i=1:length(x)
    x_round(i) = round(x(i));
    if (x_round(i)<0)
        x_round(i) = 0.0;
    endif
     if (x_round(i)>res_x)
        x_round(i) = res_x;
    endif
    y_round(i) = round(y(i));
     if (y_round(i)<0)
        y_round(i) = 0.0;
    endif
     if (y_round(i)>res_y)
        y_round(i) = res_y;
    endif
end

 

delta_x = x_round(2)-x_round(1)
delta_y = y_round(2)-y_round(1)

abs_delta_x = abs(delta_x);
abs_delta_y = abs(delta_y);

%%%%RECTANGLE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (form == "r")
for i=1:abs_delta_x
    for j=1:abs_delta_y
        if ((delta_x>0)&&(delta_y>0))
            img((y_round(1)+j),(x_round(1)+i))=constant;
        endif
        if ((delta_x<0)&&(delta_y>0))
            img((y_round(2)+j),(x_round(1)+i))=constant;
        endif
        if ((delta_x>0)&&(delta_y<0))
            img((y_round(2)+j),(x_round(1)+i))=constant;
        endif
        if ((delta_x<0)&&(delta_y<0))
            img((y_round(2)+j),(x_round(2)+i))=constant;
        endif
    endfor
endfor
endif
%%%%RECTANGLE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%CIRCLE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (form == "c")
for i=-res_x:2*res_x           
    for j=-res_y:2*res_y
         if (i^2+j^2  <= (delta_x^2+delta_y^2))
         	if (((y_round+j)<0)&&((x_round+i)>0))
                        if (((y_round(1)+abs(j))<=res_y)&&((x_round(1)+i)<=res_x))
	    img((y_round(1)+abs(j)),(x_round(1)+i))=constant;
	    endif
	endif
         	if (((y_round+j)>0)&&((x_round+i)<0))
                        if (((y_round(1)+j)<=res_y)&&((x_round(1)+abs(i))<=res_x))
	    img((y_round(1)+j),(x_round(1)+abs(i)))=constant;
	    endif
	endif
         	if (((y_round+j)<0)&&((x_round+i)<0))
                       if (((y_round(1)+abs(j))<=res_y)&&((x_round(1)+abs(i))<=res_x))
	    img((y_round(1)+abs(j)),(x_round(1)+abs(i)))=constant;
	    endif
	endif
	if (((y_round+j)>0)&&((x_round+i)>0))
                       if (((y_round(1)+j)<=res_y)&&((x_round(1)+i)<=res_x))
	    img((y_round(1)+j),(x_round(1)+i))=constant;
	    endif
	endif
        endif
    endfor
endfor
endif
%%%%CIRCLE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%LINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (form == "l")
	
if ((abs_delta_x > abs_delta_y)&&(abs_delta_y>0)&&(abs_delta_x>0))

counter = round(abs_delta_x/abs_delta_y);
y_coord = y_round(1);		
	
for i=1:abs_delta_x
      if ((delta_x>0)&&(delta_y>0))
            if (mod(i,counter)==0)	
                y_coord++;
           endif
            img((y_coord),(x_round(1)+i))=constant;
      endif
        if ((delta_x<0)&&(delta_y>0))
            if (mod(i,counter)==0)	
                y_coord++;
           endif        	
            img(y_coord,(x_round(1)-i))=constant;
        endif
        if ((delta_x>0)&&(delta_y<0))
            if (mod(i,counter)==0)	
                y_coord--;
           endif        	        	
            img((y_coord),(x_round(1)+i))=constant;
        endif
        if ((delta_x<0)&&(delta_y<0))
            if (mod(i,counter)==0)	
                y_coord--;
           endif        	        	        	
            img(y_coord,(x_round(1)-i))=constant;
        endif         
endfor
endif

if ((abs_delta_x < abs_delta_y)&&(abs_delta_y>0)&&(abs_delta_x>0))
	
counter = round(abs_delta_y/abs_delta_x);
x_coord = x_round(1);
	
for j=1:abs_delta_y
      if ((delta_x>0)&&(delta_y>0))
            if (mod(j,counter)==0)	
                x_coord++;
           endif
            img((y_round(1)+j),x_coord)=constant;
      endif
        if ((delta_x<0)&&(delta_y>0))
            if (mod(j,counter)==0)	
                x_coord--;
           endif        	
            img((y_round(1)+j),x_coord)=constant;
        endif
        if ((delta_x>0)&&(delta_y<0))
            if (mod(j,counter)==0)	
                x_coord++;
           endif        	        	
            img((y_round(1)-j),x_coord)=constant;
        endif
        if ((delta_x<0)&&(delta_y<0))
            if (mod(j,counter)==0)	
                x_coord--;
           endif        	        	        	
            img((y_round(1)-j),x_coord)=constant;
        endif         
endfor
endif

if ((delta_y==0)&&(delta_x>0))
for i=1:abs_delta_x
      img(y_round(1),x_round(1)+i)=constant;
endfor
endif	

if ((delta_y==0)&&(delta_x<0))
for i=1:abs_delta_x
      img(y_round(1),x_round(1)-i)=constant;
endfor
endif

if ((delta_y>0)&&(delta_x==0))
for j=1:abs_delta_y
      img(y_round(1) +j,x_round(1))=constant;
endfor
endif

if ((delta_y<0)&&(delta_x==0))
for j=1:abs_delta_y
      img(y_round(1) -j,x_round(1))=constant;
endfor
endif

endif
%%%%LINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%H-LINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (form == "h")
	
if (delta_x>0)
for i=1:abs_delta_x
      img(y_round(1),x_round(1)+i)=constant;
endfor
endif	

if (delta_x<0)
for i=1:abs_delta_x
      img(y_round(1),x_round(1)-i)=constant;
endfor
endif
endif
%%%%H-LINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%V-LINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (form == "v")

if (delta_y>0)
for j=1:abs_delta_y
      img(y_round(1) +j,x_round(1))=constant;
endfor
endif

if (delta_y<0)
for j=1:abs_delta_y
      img(y_round(1) -j,x_round(1))=constant;
endfor
endif

endif
%%%%V-LINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear x,y;

clear x_round;
clear y_round;
clear delta_x;
clear delta_y;

ih=imagesc(img);
axis image;
set(gca,'YDir','normal')
caxis([0,2]);
colormap(jet(3));

clear form;
clear form1;

form1 = input("[o]bstacle or [w]indemitter or [s]top", "s");

if (form1=="w")
   form2= input("In which direction does the wind blow? [2]..+x  [3]..-x [4]..+y [5]..-y [6]..+z [7]..-z"); 
endif

if (form1 ~= "s")
form = form = input("[r]ectangle, [l]ine, [c]ircle, [h]orizontal line, [v]ertical line", "s");
else
disp("Done");	
endif
	
endwhile

ih=imagesc(img);
axis image;
set(gca,'YDir','normal');
caxis([0,2]);
colormap(gray(3));
colorbar;200


%now write the TYCHO Domain File:
fd=fopen('domain_ic.tyc','wb');

for i=1:res_y
    for j=1:res_x
        if (img(i,j)==1)
           fwrite(fd,img(i,j),'int32');
        else
           fwrite(fd,0,'int32');          
       endif
    end
end
fclose(fd);

disp('The domain_ic.tyc file is written.');

%now write the TYCHO Wind File:
fd=fopen('wind_ic.tyc','wb');

for i=1:res_y
    for j=1:res_x
        if (img(i,j)~=1)
           fwrite(fd,img(i,j),'int32');
        else
           fwrite(fd,0,'int32');          
       endif
    end
end
fclose(fd);

disp('The wind_ic.tyc file is written.');
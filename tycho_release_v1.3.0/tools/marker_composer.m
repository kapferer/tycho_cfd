%%%%%%   TYCHO Marker COMPOSER 2D %%%%%%%%%%%%% 
%This octave script allows you to define a Marker field in the         %
%computational domain with resolution defined at the beginning  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;

disp "==========================="
disp "This is the 2D TYCHO MARKER COMPOSER"
disp "==========================="
disp "Do you want to read in a obstacle distribution?"
disp "==========================="

read_in=input("Read in an obstacle distribution [y]es - [n]o: ", "s");

resolution_x = input("The Resolution in X direction: ", "s");
res_x = str2num(resolution_x);
resolution_y = input("The Resolution in Y direction: ", "s");
res_y = str2num(resolution_y);

if (read_in=="y")
%now write the TYCHO Domain File:
fd=fopen('domain_ic.tyc','r');
for i=1:res_y
    for j=1:res_x
      img(i,j)= fread(fd,1,'int32');
    end
end
fclose(fd);
endif

if (read_in=="n")
    img=zeros(res_x, res_y);
endif

disp "Please choose beteen different forms"
disp "------------------------------------------------"
disp "r......rectangle"
disp "l......line"

form = input("[r]ectangle, [l]ine, [c]ircle, [h]orizontal line, [v]ertical line or [s]top", "s");

ih=imagesc(img);
axis image;
set(gca,'YDir','normal')
caxis([0,1]);
colormap(jet(2));

while (form ~= "s")

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
            img((y_round(1)+j),(x_round(1)+i))=2;
        endif
        if ((delta_x<0)&&(delta_y>0))
            img((y_round(2)+j),(x_round(1)+i))=2;
        endif
        if ((delta_x>0)&&(delta_y<0))
            img((y_round(2)+j),(x_round(1)+i))=2;
        endif
        if ((delta_x<0)&&(delta_y<0))
            img((y_round(2)+j),(x_round(2)+i))=2;
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
	    img((y_round(1)+abs(j)),(x_round(1)+i))=2;
	    endif
	endif
         	if (((y_round+j)>0)&&((x_round+i)<0))
                        if (((y_round(1)+j)<=res_y)&&((x_round(1)+abs(i))<=res_x))
	    img((y_round(1)+j),(x_round(1)+abs(i)))=2;
	    endif
	endif
         	if (((y_round+j)<0)&&((x_round+i)<0))
                       if (((y_round(1)+abs(j))<=res_y)&&((x_round(1)+abs(i))<=res_x))
	    img((y_round(1)+abs(j)),(x_round(1)+abs(i)))=2;
	    endif
	endif
	if (((y_round+j)>0)&&((x_round+i)>0))
                       if (((y_round(1)+j)<=res_y)&&((x_round(1)+i)<=res_x))
	    img((y_round(1)+j),(x_round(1)+i))=2;
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
            img((y_coord),(x_round(1)+i))=2;
      endif
        if ((delta_x<0)&&(delta_y>0))
            if (mod(i,counter)==0)	
                y_coord++;
           endif        	
            img(y_coord,(x_round(1)-i))=2;
        endif
        if ((delta_x>0)&&(delta_y<0))
            if (mod(i,counter)==0)	
                y_coord--;
           endif        	        	
            img((y_coord),(x_round(1)+i))=2;
        endif
        if ((delta_x<0)&&(delta_y<0))
            if (mod(i,counter)==0)	
                y_coord--;
           endif        	        	        	
            img(y_coord,(x_round(1)-i))=2;
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
            img((y_round(1)+j),x_coord)=2;
      endif
        if ((delta_x<0)&&(delta_y>0))
            if (mod(j,counter)==0)	
                x_coord--;
           endif        	
            img((y_round(1)+j),x_coord)=2;
        endif
        if ((delta_x>0)&&(delta_y<0))
            if (mod(j,counter)==0)	
                x_coord++;
           endif        	        	
            img((y_round(1)-j),x_coord)=2;
        endif
        if ((delta_x<0)&&(delta_y<0))
            if (mod(j,counter)==0)	
                x_coord--;
           endif        	        	        	
            img((y_round(1)-j),x_coord)=2;
        endif         
endfor
endif

if ((delta_y==0)&&(delta_x>0))
for i=1:abs_delta_x
      img(y_round(1),x_round(1)+i)=2;
endfor
endif	

if ((delta_y==0)&&(delta_x<0))
for i=1:abs_delta_x
      img(y_round(1),x_round(1)-i)=2;
endfor
endif

if ((delta_y>0)&&(delta_x==0))
for j=1:abs_delta_y
      img(y_round(1) +j,x_round(1))=2;
endfor
endif

if ((delta_y<0)&&(delta_x==0))
for j=1:abs_delta_y
      img(y_round(1) -j,x_round(1))=2;
endfor
endif

endif
%%%%LINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%H-LINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (form == "h")
	
if (delta_x>0)
for i=1:abs_delta_x
      img(y_round(1),x_round(1)+i)=2;
endfor
endif	

if (delta_x<0)
for i=1:abs_delta_x
      img(y_round(1),x_round(1)-i)=2;
endfor
endif
endif
%%%%H-LINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%V-LINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (form == "v")

if (delta_y>0)
for j=1:abs_delta_y
      img(y_round(1) +j,x_round(1))=2;
endfor
endif

if (delta_y<0)
for j=1:abs_delta_y
      img(y_round(1) -j,x_round(1))=2;
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
caxis([0,1]);
colormap(jet(2));

clear form;
form = input("[r]ectangle, [l]ine, [c]ircle or [s]top", "s");

endwhile

ih=imagesc(img);
axis image;
set(gca,'YDir','normal');
caxis([0,2]);
colormap(gray(3));
colorbar;200


%now write the TYCHO Marker File:
fd=fopen('marker_ic.tyc','wb');
data = 0;

for i=1:res_y
    for j=1:res_x
    	if (img(i,j)==2)
    	data=1;
    	else
	data=0;
    	endif	
        fwrite(fd,data,'int32');
    end
end
fclose(fd);

disp('The marker_ic.tyc file is written.');

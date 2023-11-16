function alldata = All_Data() 

close all; clear all;

load('Data_Wen-Chuang.mat')

Nsubj     = length(rec_all);
N_trials  = length(rec_all(1).set);

for k = 1: length(rec_all)%load(files_all(k).name)
    
    rec = rec_all(k);
    
    field1 = 'set_size';
    value1 = [rec.set];
    
    field2 = 'stim_order';
    value2 = [rec.ord];
    
    field3 = 'spatial_dist';
    value3 = [rec.dis];
    
    field4 = 'col_val';
    value4 = [rec.c];
    
    field5 = 'response';
    value5 = [rec.res];
    
    field6 = 'reaction_time';
    value6 = [rec.rt];
    
    for j = 1:N_trials
        
        
      
        %rec.c is between 0 and 180
        color_target(j)    = 2*(rec.c(rec.ord(:,j) == 1,j)/180*pi-pi/2);% because color values are 2 degrees apart
        color_nontarget(j) = 2*(rec.c(rec.ord(:,j) == 2,j)/180*pi-pi/2);
        % both unif between (-pi,pi)
    end
    
    field7 = 'col_dist';
    value7 = abs(color_target - color_nontarget);
    value7(value7>pi) = 2*pi - value7(value7>pi); % uniform between 0 and pi
    
    data = struct(field1, value1, field2, value2, field3, value3, field4, value4, ...
        field5, value5, field6, value6, field7, value7);
    
    alldata(k).data = data;
end


end



function [out,GPS_wk,GPS_sec_wk] = read_pos_auto_with_velocity()

if nargin == 0
    [file, path] = uigetfile('*.*');
    filename = fullfile(path,file);  
end

    f = fopen(filename);
    k = 0;
    
    for i = 1:14
        s=fgetl(f);
    end
    date = [];
    while feof(f)==0 
        
        s=fgetl(f);
        if s(1) ~= '%'
            S = split(s);
            
            k = k + 1;
            [GPS_wk(k), GPS_sec_wk(k)] = GPSweek(str2num(S{1,1}(1:4)),str2num(S{1,1}(6:7)),str2num(S{1,1}(9:10)),str2num(S{2,1}(1:2)),str2num(S{2,1}(4:5)),str2num(S{2,1}(7:end)));
%             date = [date; str2num(S{2,1}(1:2)) * 3600 + str2num(S{2,1}(4:5)) * 60 + str2num(S{2,1}(7:end))];
            t = str2num(S{2,1}(1:2)) * 3600 + str2num(S{2,1}(4:5)) * 60 + str2num(S{2,1}(7:end));
            out(:,k) = [t; str2num(S{3,1}); str2num(S{4,1}); str2num(S{5,1});str2num(S{16,1});str2num(S{17,1});str2num(S{18,1})];
%             out(:,k) = [t; str2num(S{3,1}); str2num(S{4,1}); str2num(S{5,1})];
        end
           
    end
    
    
fclose(f);

end




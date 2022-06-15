function [out_new] = read_and_parse_mdek_raw_log(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0
    [file, path] = uigetfile('*.*');
    filename = fullfile(path,file);  
end
fid = fopen(filename);%create file variable and open file fname
k = 0;

while feof(fid) == 0% a~=-1
    current_string = fgetl(fid); %read new line
    k = k + 1;
    out(k).raw_data = current_string;
end

%%
% КАРОЧЕ НАДО СДЕЛАТЬ СНАЧАЛА ПРОВЕРКУ ПО ВСЕМ НАБЛЮДЕН�?ЯМ, НАЙТ�? МАКС�?МУМ
% Ч�?СЛА Д�?СТОВ, ТАМ СДЕЛАТЬ ПАРС ВСЕХ ЯКОРЕЙ �? ЗАП�?САТЬ �?Х ПОСЫ �? АЙД�?ШН�?К�?
% ПО СТРУКТУРЕ �? ДАЛЬШЕ В СВ�?ТЧ-КЕЙСЕ ПЕРЕБ�?РАТЬ �?Х В Ц�?КЛЕ ФОР ПО ВСЕМ
% ACN_IND
dist_count = nan(1,length(out));
out_struct = struct('PC_time',[],'anc_id',[],'anc_pos',[],'anc_dist',[],'pos',[]);
for k = 1:length(out)
    data_mass = strsplit(out(k).raw_data,',');
    if any( (contains(data_mass,'DIST')))
        dist_count(k) = str2num(data_mass{3});
    end
end
[val,ind] = max(dist_count);
data_mass = strsplit(out(ind).raw_data,',');
anc_ind_mass = find(contains(data_mass, 'AN'));
[anc_ids_mass,ind] = sort(data_mass(anc_ind_mass+1));
%%
out_struct.anc_id = anc_ids_mass;
out_new = repmat(out_struct,1,length(out));
%%
for k = 1:length(out)
    data_mass = strsplit(out(k).raw_data,',');
    out_new(k).PC_time = str2num(data_mass{1});
    anc_ind = find(contains(data_mass, 'AN'));
%     [anc_id,ind] = sort(data_mass(anc_ind+1));
    anc_id = data_mass(anc_ind+1);
    out_new(k).anc_pos = zeros(4,3);
    out_new(k).anc_dist = zeros(4,1);
    for i = 1:length(anc_ind)
        switch anc_id{i}
            case out_new(k).anc_id{1}
                out_new(k).anc_pos(1,:) = [str2num(data_mass{anc_ind(i)+2}) str2num(data_mass{anc_ind(i)+3}) str2num(data_mass{anc_ind(i)+4})];
                out_new(k).anc_dist(1,1) = str2num(data_mass{anc_ind(i)+5});
            case out_new(k).anc_id{2}
                out_new(k).anc_pos(2,:) = [str2num(data_mass{anc_ind(i)+2}) str2num(data_mass{anc_ind(i)+3}) str2num(data_mass{anc_ind(i)+4})];
                out_new(k).anc_dist(2,1) = str2num(data_mass{anc_ind(i)+5});
            case out_new(k).anc_id{3}
                out_new(k).anc_pos(3,:) = [str2num(data_mass{anc_ind(i)+2}) str2num(data_mass{anc_ind(i)+3}) str2num(data_mass{anc_ind(i)+4})];
                out_new(k).anc_dist(3,1) = str2num(data_mass{anc_ind(i)+5});
            case out_new(k).anc_id{4}
                out_new(k).anc_pos(4,:) = [str2num(data_mass{anc_ind(i)+2}) str2num(data_mass{anc_ind(i)+3}) str2num(data_mass{anc_ind(i)+4})];
                out_new(k).anc_dist(4,1) = str2num(data_mass{anc_ind(i)+5});
        end
%         out(k).anc_id(i,:) = anc_id{i};
%         out(k).anc_pos(i,:) = [str2num(data_mass{anc_ind(ind(i))+2}) str2num(data_mass{anc_ind(ind(i))+3}) str2num(data_mass{anc_ind(ind(i))+4})];
%         out(k).anc_dist(i,1) = str2num(data_mass{anc_ind(ind(i))+5});

%         anc_ids = data_mass(anc_ind+1);
%         anc_pos = [];
%         anc_dist = data_mass{anc_ind+5};
    end
    if any(contains(data_mass,'POS'))
        pos_ind = find(strcmp(data_mass, 'POS'));
        out_new(k).pos = [str2num(data_mass{pos_ind+1});str2num(data_mass{pos_ind+2});str2num(data_mass{pos_ind+3})];
    end
end
% out = rmfield(out,'raw_data');
end

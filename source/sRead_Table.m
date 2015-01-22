% This reads a .CSV or a tab-separated table.
% The first line is expected to be the header
%
% Written by Pejman Mohammadi
% pejman.mohammadi@bsse.ethz.ch
%
%

function Data = sRead_Table(Path2File)
DLM = '[,\t]'; % list of potential delimiters
Fin = fopen(Path2File, 'r');
% Identify Names
RawHeader = fgetl(Fin);
dlm    = regexp(RawHeader, DLM, 'match', 'once');
Header = regexp(RawHeader, dlm, 'split');
Header = Remove_nonAlphanumerics(Header);
noHeaderLines = 1; % Number of header lines
%% Identify File format
Firstrow = regexp(fgetl(Fin), dlm, 'split');
fclose(Fin);
if Firstrow{1}(1)=='#'
    % This is a header row with values to guide correct format specification.
    % For example if the fisrt row of a numerical column by chance is "NaN",
    % auto formating would recognize it as "Text", this extra line that startes
    % with "#' holds values that ensure correct format recognition, without
    % being included in the final data.
    Firstrow{1}(1)=[];
    noHeaderLines = 2;
end
FormatS = '';
Report  = ''; 
if length(Header)<length(Firstrow)
    warning('The header line does not match the data')
    for c = length(Header)+1:length(Firstrow)
        Header{c} = ['UnLabeled_C' int2str(c)];
    end
end
for i = 1:length(Firstrow)
    if ~isempty(regexp(Firstrow{i}, '^[+-]?\d+$'))
        % It's an integer
        FormatS = [FormatS '%f' ];
        Report  = [Report  Header{i} '(int)\t'];
    elseif ~isempty(regexp(Firstrow{i}, '^[+-]?\d+(\.)?\d*[eE]?[+-]?\d+$'))
        %read it as Float
        FormatS = [FormatS '%f' ];
        Report  = [Report  Header{i} '(float)\t'];
    else
        %read it as string
        FormatS = [FormatS '%s' ];
        Report  = [Report  Header{i} '(str)\t'];
    end
end
% fprintf(['Input file: %s \nRecognized format:\n' Report '\n'], Path2File);

%% Read file
Fin = fopen(Path2File, 'r');
DES_Res = textscan(Fin, FormatS, 'HeaderLines',noHeaderLines,'Delimiter', dlm);%, 'bufsize', 1E+5); 
for i = 1:length(DES_Res)
    Data.(Header{i}) = DES_Res{i}; DES_Res{i} = [];
end  

fclose(Fin);
end


function S = Remove_nonAlphanumerics(S)
% this removed all but alphanumerics and substitutes them with '_'. Extra
% underscores at the begining and the end are removed.
for i = 1:length(S)
    tmpS = S{i};
    Valid = ...
        (tmpS >='0' & tmpS <= '9')  | ...
        (tmpS >='a' & tmpS <= 'z')  | ...
        (tmpS >='A' & tmpS <= 'Z')  | ...
        (tmpS == '_');
    tmpS(~Valid)= '_';
    while tmpS(1)  =='_'; tmpS(1)  =[];end
    while tmpS(end)=='_'; tmpS(end)=[];end
    S{i} = tmpS;
end
end

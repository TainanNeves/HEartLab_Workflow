function t = CalcTime(path)
%Function to calculate the initial recording time, from sync file from Open
%Ephys

filename = 'sync_messages.txt';
fullfilename = fullfile(path,filename);
f = fopen(fullfilename,"r");
f = fscanf(f, '%s');
s = textscan(f, 'Software Time (milliseconds since midnight Jan 1st 1970 UTC): %d64');
mil = milliseconds(s{1});

d = "01-01-1970 00:00:00.000";
t = datetime(d, "Format", "MM-dd-yyyy HH:mm:ss.SSS", 'TimeZone', 'UTC');
t = t + mil;
t.TimeZone = 'America/Sao_Paulo';




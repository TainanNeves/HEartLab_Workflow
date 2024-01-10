% rewrites DATA file into Evolve camera file format. 
% syntax  opwrite(file, DATA, filenew, test)
% file - old filename, header from this file is attached to filenew
% DATA - matlab data variable to be saved
% filenew - new filename in Evolve camera format
% test, bonary flag, when 1 after saving, saved file is readout and
% compared with DATA variable
% this script uses opread function to tests
function [v,h,times] = opread(file, DATA, filenew, test) % this is syntax of the function

h = opheader(file);
disp('File header:');
disp(h);



[fid,msg] = fopen(file,'rb',h.machineformat);
if fid == -1
   error('Cannot open file %s -> %s',file,msg);
end

[fid2,msg] = fopen(filenew,'w');
if fid2 == -1
   error('Cannot open file %s -> %s',file,msg);
end

frewind(fid);
Header = uint8(fread( fid,1024, 'uint8' ));

frewind(fid2);
fwrite(fid2, Header(1:length(Header)), 'uint8' );
fseek(fid2,5,'bof'); fwrite(fid2,uint32(size(DATA,3)),'uint32');
fseek(fid2,9,'bof');
fwrite(fid2, size(DATA,2) );
fseek(fid2,13,'bof');
fwrite(fid2,size(DATA,1));
fseek(fid2,1024,'bof');

for loop = 1:size(DATA,3)
    image = DATA(:,:,loop)'; image=reshape(image, [1, size(image,1)*size(image,2)]);
    fwrite(fid2, image, 'uint16'); 
    fwrite(fid2, 1, '*uint64');
	if (mod(loop,50) == 0)
		stop = progressbar(loop/size(DATA,3));
		if (stop) 
			disp('OPREAD aborted by user.');
			break; 
		end
	end
end
ftell(fid2)
progressbar(1); % finalize the progress bar.
close(gcf); % close the progressbar figure.
fclose(fid);fclose(fid2);


if test==1
    TEST = opread(filenew);
    TEST = sum(sum(sum(TEST,1),2),3);
    DATA = sum(sum(sum(DATA,1),2),3);
    if (DATA - TEST) ==0
        display('Write test ok!')
    else
        display('Write test failed')
    end
end

return;

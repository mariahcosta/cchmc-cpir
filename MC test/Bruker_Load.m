function raw_data = Bruker_Load(file_name);

%Get the file
if nargin>0
      filename= file_name;  
else
    [filename,path]=uigetfile('*','Select file');
    cd(path)
end

% open data
fd=fopen(char(filename),'r','l');
fseek(fd,0,'bof');  
raw=fread(fd,inf,'int32');
raw_data = complex(raw(1:2:end),raw(2:2:end));
fclose(fd);

end

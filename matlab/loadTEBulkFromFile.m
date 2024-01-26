function bulkTE_in = loadTEBulkFromFile(file,path);

	data    		= importdata(strcat(path,file), ' ');
	bulk_data       = data.data;
	bulk_TE     	= data.colheaders;

	if length(data.data) ~= 28
		disp('Error, the number of provided trace elements is wrong');
	end

	bulkTE_in = table(bulk_TE',bulk_data','VariableNames',{'El','C [ppm]'});

end

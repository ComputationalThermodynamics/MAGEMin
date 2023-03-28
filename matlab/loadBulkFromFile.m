function bulk_in = loadBulkFromFile(sysunit,file,path,db);

	data    = importdata(strcat(path,file), ' ');

	bulk        = data.data;
	bulk_ox     = data.colheaders;

	ref_ox          = {'SiO2', 'Al2O3', 'CaO', 'MgO' ,'FeO', 'Fe2O3', 'K2O','Na2O', 'TiO2', 'O', 'Cr2O3', 'MnO', 'H2O'};
	ref_MolarMass   = [60.08 101.96 56.08 40.30 71.85 79.85 94.2 61.98 79.88 16.0 151.99,70.94,18.015];      %Molar mass of oxides

	if strcmp(db,'ig')
		MAGEMin_bulk    = zeros(11,1);
		MAGEMin_ox      = {'SiO2', 'Al2O3', 'CaO', 'MgO' ,'FeO', 'K2O','Na2O', 'TiO2', 'O', 'Cr2O3', 'H2O'};
	end
	if strcmp(db,'mp')
		MAGEMin_bulk    = zeros(11,1);
		MAGEMin_ox      = {'SiO2', 'Al2O3', 'CaO', 'MgO' ,'FeO', 'K2O','Na2O', 'TiO2', 'O', 'MnO', 'H2O'};
	end

	% convert to mol, if system unit = wt
	if strcmp(sysunit, 'wt')
		for i=1:length(bulk_ox)
			id = find(strcmp(ref_ox, bulk_ox(i)));
			bulk(i) = bulk(i)/ref_MolarMass(id);
		end
	end
	bulk = normalize(bulk,'norm',1); 

	for i=1:length(MAGEMin_ox)
		id = find(strcmp(bulk_ox, MAGEMin_ox(i)));
		if id
			MAGEMin_bulk(i) = bulk(id);
		end
	end
	idFe2O3 = find(strcmp(bulk_ox, 'Fe2O3'));
	if idFe2O3
		%disp('convert Fe2O3 to FeOt and O...');  
		
		idFeO = find(strcmp(MAGEMin_ox, 'FeO'));
		MAGEMin_bulk(idFeO) = MAGEMin_bulk(idFeO) + bulk(idFe2O3)*2;
		
		idO = find(strcmp(MAGEMin_ox, 'O'));
		MAGEMin_bulk(idO) = MAGEMin_bulk(idO) + bulk(idFe2O3);
	end

	MAGEMin_bulk = normalize(MAGEMin_bulk,'norm',1);

	idNonH2O = find( strcmp(MAGEMin_ox, 'H2O') == 0);
	MAGEMin_bulk(find(MAGEMin_bulk(idNonH2O) == 0)) = 1e-4;
	MAGEMin_bulk = normalize(MAGEMin_bulk,'norm',1)*100.0;
	
	bulk_in = table(MAGEMin_ox',MAGEMin_bulk,MAGEMin_bulk,'VariableNames',{'Oxide','mol %','mol2 %'});

end

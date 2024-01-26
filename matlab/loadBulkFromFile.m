function bulk_in = loadBulkFromFile(sysunit,file,path,db);

	data    		= importdata(strcat(path,file), ' ');
	bulk_txt       	= data.data;
	bulk_ox     	= data.colheaders;
	ref_ox          = {'SiO2', 'Al2O3', 'CaO', 'MgO' ,'FeO', 'Fe2O3', 'K2O','Na2O', 'TiO2', 'O', 'Cr2O3', 'MnO', 'H2O', 'CO2', 'S'};
	ref_MolarMass   = [60.08, 101.96, 56.08, 40.30, 71.85,  159.69, 94.2,  61.98,  79.88,  16.0, 151.99, 70.94, 18.015, 44.01,32.06	];      %Molar mass of oxides
	nB 				= size(bulk_txt,1);

	if nB == 1 
		bulk = bulk_txt;
	else
		bulk  = bulk_txt(1,:);
		bulk2 = bulk_txt(2,:);
	end

	if strcmp(db,'ig')
		MAGEMin_bulk    = zeros(11,1);
		MAGEMin_ox      = {'SiO2', 'Al2O3', 'CaO', 'MgO' ,'FeO', 'K2O','Na2O', 'TiO2', 'O', 'Cr2O3', 'H2O'};
		nonZero 		= [1 2 3 4 5 6 7]';
	elseif strcmp(db,'igd')
		MAGEMin_bulk    = zeros(11,1);
		MAGEMin_ox      = {'SiO2', 'Al2O3', 'CaO', 'MgO' ,'FeO', 'K2O','Na2O', 'TiO2', 'O', 'Cr2O3', 'H2O'};
		nonZero 		= [1 2 3 4 5 6 7]';
	elseif strcmp(db,'alk')
		MAGEMin_bulk    = zeros(11,1);
		MAGEMin_ox      = {'SiO2', 'Al2O3', 'CaO', 'MgO' ,'FeO', 'K2O','Na2O', 'TiO2', 'O', 'Cr2O3', 'H2O'};
		nonZero 		= [1 2 3 4 5 6 7]';
	elseif strcmp(db,'mp')
		MAGEMin_bulk    = zeros(11,1);
		MAGEMin_ox      = {'SiO2', 'Al2O3', 'CaO', 'MgO' ,'FeO', 'K2O','Na2O', 'TiO2', 'O', 'MnO', 'H2O'};
		nonZero 		= [1 2 3 4 5 6 7]';
	elseif strcmp(db,'mb')
		MAGEMin_bulk    = zeros(10,1);
		MAGEMin_ox      = {'SiO2', 'Al2O3', 'CaO', 'MgO' ,'FeO', 'K2O','Na2O', 'TiO2', 'O', 'H2O'};
		nonZero 		= [1 2 3 4 5 6 7]';
	elseif strcmp(db,'um')
		MAGEMin_bulk    = zeros(7,1);
		MAGEMin_ox      = {'SiO2', 'Al2O3', 'MgO' ,'FeO', 'O', 'H2O', 'S'};
		nonZero 		= [1 2 3 4 5 6 7]';
	else
		disp('wrong database id...')
	end

	if nB == 2 
		MAGEMin_bulk2    = zeros(size(MAGEMin_bulk));
	end

	% convert to mol, if system unit = wt
	if strcmp(sysunit, 'wt')
		for i=1:length(bulk_ox)
			id = find(strcmp(ref_ox, bulk_ox(i)));
			bulk(i) = bulk(i)/ref_MolarMass(id);
			if nB == 2 
				bulk2(i) = bulk2(i)/ref_MolarMass(id);
			end
		end
	end
	bulk = normalize(bulk,'norm',1);
	if nB == 2 
		bulk2 = normalize(bulk2,'norm',1);
	end
	for i=1:length(MAGEMin_ox)
		id = find(strcmp(bulk_ox, MAGEMin_ox(i)));
		if id
			MAGEMin_bulk(i) = bulk(id);
			if nB == 2 
				MAGEMin_bulk2(i) = bulk2(id);
			end
		end
	end

	idFe2O3 = find(strcmp(bulk_ox, 'Fe2O3'));
	if idFe2O3
		%disp('convert Fe2O3 to FeOt and O...');  
		
		idFeO = find(strcmp(MAGEMin_ox, 'FeO'));
		MAGEMin_bulk(idFeO) = MAGEMin_bulk(idFeO) + bulk(idFe2O3)*2.0;
		if nB == 2 
			MAGEMin_bulk2(idFeO) = MAGEMin_bulk2(idFeO) + bulk2(idFe2O3);
		end
		
		idO = find(strcmp(MAGEMin_ox, 'O'));
		MAGEMin_bulk(idO) = MAGEMin_bulk(idO) + bulk(idFe2O3);
		if nB == 2 
			MAGEMin_bulk2(idO) = MAGEMin_bulk2(idO) + bulk2(idFe2O3);
		end
	end

	MAGEMin_bulk = normalize(MAGEMin_bulk,'norm',1);
	if nB == 2 
		MAGEMin_bulk2 = normalize(MAGEMin_bulk2,'norm',1);
	end
	% idNonH2O = find( strcmp(MAGEMin_ox', 'H2O') == 0);
	MAGEMin_bulk(find(MAGEMin_bulk(nonZero) == 0)) = 1e-4;
	if nB == 2 
		MAGEMin_bulk2(find(MAGEMin_bulk2(nonZero) == 0)) = 1e-4;
	end
	MAGEMin_bulk = normalize(MAGEMin_bulk,'norm',1)*100.0;
	if nB == 2 
		MAGEMin_bulk2 = normalize(MAGEMin_bulk2,'norm',1)*100.0;
	end

	if nB == 1 
		bulk_in = table(MAGEMin_ox',MAGEMin_bulk,MAGEMin_bulk,'VariableNames',{'Oxide','mol %','mol2 %'});
	else
		bulk_in = table(MAGEMin_ox',MAGEMin_bulk,MAGEMin_bulk2,'VariableNames',{'Oxide','mol %','mol2 %'});
	end


end

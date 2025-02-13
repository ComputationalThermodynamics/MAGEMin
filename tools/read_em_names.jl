using JSON

# Read the JSON file
file_path   = "em_name.json"
dict_em     = JSON.parsefile(file_path)


dict_ss = Dict{String, Any}()
dict_ss["sp"] = "spinel",[ "sp" "spinel";
                            "mt" "magnetite" ] 

dict_ss["spl"] = "spinel",[ "spl" "spinel";
                            "cm" "chromite";
                            "usp" "uvospinel";
                            "mgt" "magnetite" ] 

dict_ss["fsp"] = "feldspar",[   "afs" "alkali-feldspar";
                                "pl" "plagioclase" ] 

dict_ss["mu"] = "muscovite",[   "pat" "paragonite";
                                "mu" "muscovite" ] 

dict_ss["amp"] = "amphibole",[  "gl" "glaucophane";
                                "act" "actinolite";
                                "cumm" "cummingtonite";
                                "tr" "tremolite";
                                "amp" "amphibole" ] 

dict_ss["ilm"] = "ilmenite",[   "hem" "hematite";
                                "ilm" "ilmenite" ] 

dict_ss["ilmm"] = "ilmenite",[   "hemm" "hematite";
                                "ilmm" "ilmenite" ] 

dict_ss["nph"] = "nepheline",[  "K-nph" "K-rich nepheline";
                                "nph" "nepheline" ] 

dict_ss["cpx"] = "clinopyroxene",[  "pig" "pigeonite";
                                    "Na-cpx" "Na-rich clinopyroxene";
                                    "cpx" "clinopyroxene" ] 

dict_ss["dio"] = "diopside",[   "dio" "diopside";
                                "omph" "omphacite";
                                "jd" "jadeite" ] 

dict_ss["occm"] = "carbonate",[ "sid" "siderite";
                                "ank" "ankerite";
                                "mag" "magnesite";
                                "cc" "calcite"] 

dict_ss["opx"]  = "orthopyroxene",[ "opx" "orthopyroxene"] 
dict_ss["liq"]  = "liquid",[ "liq" "liquid"] 
dict_ss["ol"]   = "olivine",[ "ol" "olivine"] 
dict_ss["ep"]   = "epidote",[ "ep" "epidote"] 
dict_ss["g"]    = "garnet",[ "g" "garnet"] 
dict_ss["chl"]  = "chlorite",[ "chl" "chlorite"] 
dict_ss["bi"]   = "biotite",[ "bi" "biotite"] 
dict_ss["aug"]  = "augite",[ "aug" "augite"] 
dict_ss["abc"]  = "abc",[ "abc" "abc"] 
dict_ss["ctd"]  = "chloritoid",[ "ctd" "chloritoid"] 
dict_ss["ma"]   = "margarite",[ "ma" "margarite"] 
dict_ss["st"]   = "staurolite",[ "st" "staurolite"] 
dict_ss["sa"]   = "sapphirine",[ "sa" "sapphirine"] 
dict_ss["cd"]   = "cordierite",[ "cd" "cordierite"] 
dict_ss["mt"]   = "magnetite",[ "mt" "magnetite"] 
dict_ss["fl"]   = "fluid",[ "fl" "fluid"] 
dict_ss["lct"]  = "leucite",[ "lct" "leucite"] 
dict_ss["kals"] = "kalsilite",[ "kals" "kalsilite"] 
dict_ss["mel"]  = "melilite",[ "mel" "melilite"] 
dict_ss["br"]   = "brucite",[ "br" "brucite"] 
dict_ss["atg"]  = "antigorite",[ "atg" "antigorite"] 
dict_ss["ta"]   = "talc",[ "ta" "talc"] 
dict_ss["anth"] = "anthophyllite",[ "anth" "anthophyllite"] 
dict_ss["po"]   = "pyrrhotite",[ "po" "pyrrhotite"] 

dict_ss["fp"]   = "ferropericlase",[ "fp" "ferropericlase"] 
dict_ss["mpv"]  = "MgSi-perovskite",[ "mpv" "MgSi-perovskite"] 
dict_ss["cpv"]  = "CaSi-perovskite",[ "cpv" "CaSi-perovskite"] 
dict_ss["crn"]  = "Mg-Al corundum",[ "crn" "Mg-Al corundum"] 
dict_ss["cf"]   = "Ca-Ferrite",[ "cf" "Ca-Ferrite"] 
dict_ss["nal"]  = "NAL-phase",[ "nal" "NAL-phase"] 
dict_ss["aki"]  = "akimotoite",[ "aki" "akimotoite"] 
dict_ss["wad"]  = "Fe-Mg wadsleyite",[ "wad" "Fe-Mg wadsleyite"] 
dict_ss["ring"] = "Fe-Mg ringwoodite",[ "ring" "Fe-Mg ringwoodite"] 
dict_ss["hpx"]  = "high-P orthopyroxene",[ "hpx" "high-P orthopyroxene"] 

# Convert the dictionary to a JSON string
json_str = JSON.json(dict_ss)

# Write the JSON string to a file
file_path = "ss_name.json"
open(file_path, "w") do file
    write(file, json_str)
end
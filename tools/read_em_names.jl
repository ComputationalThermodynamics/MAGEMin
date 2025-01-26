using JSON

# Read the JSON file
file_path   = "em_name.json"
dict_em     = JSON.parsefile(file_path)

println(sorted_dict)


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
                                "omph" "omphacite" ] 

dict_ss["occm"] = "carbonate",[ "sid" "siderite";
                                "ank" "ankerite";
                                "mag" "magnesite";
                                "cc" "calcite"] 

# Convert the dictionary to a JSON string
json_str = JSON.json(dict_ss)

# Write the JSON string to a file
file_path = "ss_name.json"
open(file_path, "w") do file
    write(file, json_str)
end
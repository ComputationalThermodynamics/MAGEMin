#=
    NR 21-07-2023
    The goal here is to copy and paste table from published experimental work into an excel spreadsheet and extract/sort all needed data
    Although excel is evil, it works decently to copy from and paste from pdf

    *xls file must be formatted correctly! (check xls files available)
    *one spreadsheet for mode
    *one spreadsheet for mineral composition
    *make sure the oxide names are identical between the list and the xls file containing the tables extracted from the publication
    *If a phase is present but its fraction is not measured set the value to NaN
    *Make sure experiments are labelled with "Run"

    SPECIFIC NOTES
    ==============
    [L04]
        *Oxygen fugacity between the magnetite–wustite and the iron–wustite buffers (fO2=10^-15.91bar)
        *-3.1 to be adjusted as function of QFM
        *check -> https://fo2.rses.anu.edu.au/fo2app/ to calculate buffers
        *import table from pictures -> https://nanonets.com/free-tools/extract-table-from-image (works quite well and saves formatting time)
        *reg expression for "everything between parenthesis" -> \(([^\)]+)\)


=#
mutable struct dataArticle{ _T  } 
    file        ::  String                  # excel file name
    wox         ::  Vector{String}          # working set of oxides
    X_unit      ::  String                  # system unit for the composition (wt% or mol%)
    P_unit      ::  String                  # pressure unit ((G)Pa or (k)bar)
    T_unit      ::  String                  # temperature unit (K or C)
end


list_articles = [
    dataArticle{Float64}(
    "W98.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "GPa",
    "°C"),

    dataArticle{Float64}(
    "M03.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "GPa",
    "°C"),

    dataArticle{Float64}(
    "L04.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "GPa",
    "°C"),

    dataArticle{Float64}(
    "W03.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "GPa",
    "°C"),

    dataArticle{Float64}(
    "F99.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "GPa",
    "°C"),

    dataArticle{Float64}(
    "B94.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "GPa",
    "°C"),

    dataArticle{Float64}(
    "T93.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "GPa",
    "°C"),

    dataArticle{Float64}(
    "J90.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Fe2O3", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "kbar",
    "°C"),

    dataArticle{Float64}(
    "F88.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Fe2O3", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "kbar",
    "°C"),

    dataArticle{Float64}(
    "P00.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "GPa",
    "°C"),

    dataArticle{Float64}(
    "K98.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "GPa",
    "°C"),

    dataArticle{Float64}(
    "S01.xlsx",
    ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"],       # chemical system you want to work with
    "wt%",
    "GPa",
    "°C"),

    ]

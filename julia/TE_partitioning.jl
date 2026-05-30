#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#
# Trace element partitioning prediction routine
# Note that only metapelite, metabasite and igneous database can be used for trace element prediction
# NR 17/06/2023


"""
    get_TE_database(tedb="OL12")

    Return the built-in trace element (TE) partitioning coefficient database.

    Parameters
    ----------
    tedb : String, optional
        Database identifier (default: "OL12"). Currently only "OL12" is supported (Laurent, 2012).

    Returns
    -------
    db : custom_KDs_database
        Trace element partitioning coefficient database containing element names, phase names, and partition coefficient expressions.
"""
function get_TE_database(tedb :: String = "OL12")
    if tedb == "OL12"
        infos               = "Laurent, O. (2012). Les changements géodynamiques à la transition Archéen-Protérozoïque : étude des granitoïdes de la marge Nord du craton du Kaapvaal (Afrique du Sud). PhD, Université Blaise Pascal, Clermont-Ferrand."

        el                 = ["Rb", "Ba", "Th", "U", "Nb", "Ta", "La", "Ce", "Pb", "Pr", "Sr", "Nd", "Zr", "Hf", "Sm", "Eu", "Gd", "Tb", "Dy", "Y", "Ho", "Er", "Tm", "Yb", "Lu", "V", "Sc"]
        ph                 = ["all", "amp", "ap", "bi", "cd", "cpx", "FeTiOx", "g", "afs", "mgt", "smt", "ol", "opx", "pl", "q", "ru", "sp", "spl", "ttn", "zrc", "ep", "and", "sill", "mu"]
        KDs                = ["if [:liq].Comp_wt[1] <= 52.0 0.0632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0632455532033676 else 0.0632455532033676 end" "if [:liq].Comp_wt[1] <= 52.0 3.46410161513776 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.46410161513776 else 3.46410161513776 end" "if [:liq].Comp_wt[1] <= 52.0 424.264068711929 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 324.037034920393 else 324.037034920393 end" "if [:liq].Comp_wt[1] <= 52.0 20.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 13.228756555323 else 13.228756555323 end" "if [:liq].Comp_wt[1] <= 52.0 0.447213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.836660026534076 else 0.836660026534076 end" "if [:liq].Comp_wt[1] <= 52.0 2.73861278752583 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.836660026534076 else 0.836660026534076 end" "if [:liq].Comp_wt[1] <= 52.0 1549.19333848297 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1414.2135623731 else 1414.2135623731 end" "if [:liq].Comp_wt[1] <= 52.0 1224.74487139159 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1174.73401244707 else 1174.73401244707 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0 else 1.0 end" "if [:liq].Comp_wt[1] <= 52.0 1005.98210719674 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 948.683298050515 else 948.683298050515 end" "if [:liq].Comp_wt[1] <= 52.0 1.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0 else 1.0 end" "if [:liq].Comp_wt[1] <= 52.0 793.725393319378 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 612.372435695795 else 612.372435695795 end" "if [:liq].Comp_wt[1] <= 52.0 0.173205080756888 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.173205080756888 else 0.173205080756888 end" "if [:liq].Comp_wt[1] <= 52.0 12.2474487139159 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 11.180339887499 else 11.180339887499 end" "if [:liq].Comp_wt[1] <= 52.0 447.213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 300.0 else 300.0 end" "if [:liq].Comp_wt[1] <= 52.0 97.9795897113271 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 34.6410161513775 else 34.6410161513775 end" "if [:liq].Comp_wt[1] <= 52.0 254.950975679639 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 116.189500386222 else 116.189500386222 end" "if [:liq].Comp_wt[1] <= 52.0 169.558249578132 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 58.309518948453 else 58.309518948453 end" "if [:liq].Comp_wt[1] <= 52.0 118.321595661992 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 31.6227766016838 else 31.6227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 95.3939201416947 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 31.6227766016838 else 31.6227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 84.8528137423857 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 18.9736659610103 else 18.9736659610103 end" "if [:liq].Comp_wt[1] <= 52.0 63.2455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 13.4164078649987 else 13.4164078649987 end" "if [:liq].Comp_wt[1] <= 52.0 43.301270189222 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 10.2469507659596 else 10.2469507659596 end" "if [:liq].Comp_wt[1] <= 52.0 26.8328157299975 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 7.74596669241484 else 7.74596669241484 end" "if [:liq].Comp_wt[1] <= 52.0 17.7482393492988 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 6.70820393249937 else 6.70820393249937 end" "if [:liq].Comp_wt[1] <= 52.0 1.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0 else 1.0 end" "if [:liq].Comp_wt[1] <= 52.0 57.0087712549569 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 57.0087712549569 else 57.0087712549569 end"; "if [:liq].Comp_wt[1] <= 52.0 0.0547722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.122474487139159 else 0.223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.189736659610103 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.223606797749979 else 0.31224989991992 end" "if [:liq].Comp_wt[1] <= 52.0 0.212132034355964 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.14142135623731 else 0.0141421356237309 end" "if [:liq].Comp_wt[1] <= 52.0 0.821583836257749 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.1 else 0.0193649167310371 end" "if [:liq].Comp_wt[1] <= 52.0 2.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.559016994374947 else 0.158113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 0.4 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.273861278752583 else 0.220794021658196 end" "if [:liq].Comp_wt[1] <= 52.0 0.559016994374947 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.335410196624968 else 0.0724568837309472 end" "if [:liq].Comp_wt[1] <= 52.0 1.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.412310562561766 else 0.144913767461894 end" "if [:liq].Comp_wt[1] <= 52.0 0.632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.223606797749979 else 0.0632455532033676 end" "if [:liq].Comp_wt[1] <= 52.0 1.59687194226713 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.574456264653803 else 0.234520787991171 end" "if [:liq].Comp_wt[1] <= 52.0 0.387298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.273861278752583 else 0.3 end" "if [:liq].Comp_wt[1] <= 52.0 2.64575131106459 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.866025403784439 else 0.346410161513775 end" "if [:liq].Comp_wt[1] <= 52.0 0.632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.4 else 0.335410196624968 end" "if [:liq].Comp_wt[1] <= 52.0 0.707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.547722557505166 else 0.391152144312159 end" "if [:liq].Comp_wt[1] <= 52.0 3.46410161513776 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.36930639376292 else 0.447213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 3.46410161513776 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0 else 0.524404424085076 end" "if [:liq].Comp_wt[1] <= 52.0 4.74341649025257 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.80277563773199 else 0.589915248150105 end" "if [:liq].Comp_wt[1] <= 52.0 5.37354631505117 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.21359436211787 else 0.634822809924155 end" "if [:liq].Comp_wt[1] <= 52.0 5.42217668469038 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.29128784747792 else 0.66932802122726 end" "if [:liq].Comp_wt[1] <= 52.0 5.74456264653803 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.12132034355964 else 0.612372435695794 end" "if [:liq].Comp_wt[1] <= 52.0 5.24404424085076 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.08566536146142 else 0.681175454637056 end" "if [:liq].Comp_wt[1] <= 52.0 4.96990945591567 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.8165902124585 else 0.648074069840786 end" "if [:liq].Comp_wt[1] <= 52.0 4.5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.51657508881031 else 0.603738353924943 end" "if [:liq].Comp_wt[1] <= 52.0 3.74165738677394 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.40712472794703 else 0.53619026473818 end" "if [:liq].Comp_wt[1] <= 52.0 2.95803989154981 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.30384048104053 else 0.458257569495584 end" "if [:liq].Comp_wt[1] <= 52.0 6.70820393249937 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 5.47722557505166 else 4.89897948556636 end" "if [:liq].Comp_wt[1] <= 52.0 14.142135623731 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 10.0 else 3.87298334620742 end"; "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.273861278752583 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0790569415042095 else 0.0790569415042095 end" "if [:liq].Comp_wt[1] <= 52.0 1.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.790569415042095 else 0.790569415042095 end" "if [:liq].Comp_wt[1] <= 52.0 0.948683298050514 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.387298334620742 else 0.387298334620742 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00353553390593274 else 0.00353553390593274 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00707106781186547 else 0.00707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 13.228756555323 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.87298334620742 else 3.87298334620742 end" "if [:liq].Comp_wt[1] <= 52.0 21.2132034355964 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 7.07106781186548 else 7.07106781186548 end" "if [:liq].Comp_wt[1] <= 52.0 0.158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.316227766016838 else 0.316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 30.5777697028413 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 11.0679718105893 else 11.0679718105893 end" "if [:liq].Comp_wt[1] <= 52.0 7.07106781186548 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.93649167310371 else 1.93649167310371 end" "if [:liq].Comp_wt[1] <= 52.0 40.3112887414928 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 15.0 else 15.0 end" "if [:liq].Comp_wt[1] <= 52.0 0.387298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.193649167310371 else 0.193649167310371 end" "if [:liq].Comp_wt[1] <= 52.0 0.387298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0707106781186547 else 0.0707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 47.4341649025257 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 17.3205080756888 else 17.3205080756888 end" "if [:liq].Comp_wt[1] <= 52.0 21.2132034355964 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 11.180339887499 else 11.180339887499 end" "if [:liq].Comp_wt[1] <= 52.0 52.9150262212919 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 18.3711730708738 else 18.3711730708738 end" "if [:liq].Comp_wt[1] <= 52.0 52.9150262212919 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 17.8885438199983 else 17.8885438199983 end" "if [:liq].Comp_wt[1] <= 52.0 45.8257569495584 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 15.6524758424985 else 15.6524758424985 end" "if [:liq].Comp_wt[1] <= 52.0 38.7298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 12.2474487139159 else 12.2474487139159 end" "if [:liq].Comp_wt[1] <= 52.0 37.0809924354783 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 13.4164078649987 else 13.4164078649987 end" "if [:liq].Comp_wt[1] <= 52.0 30.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 11.180339887499 else 11.180339887499 end" "if [:liq].Comp_wt[1] <= 52.0 22.9128784747792 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 8.94427190999916 else 8.94427190999916 end" "if [:liq].Comp_wt[1] <= 52.0 15.8113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 6.70820393249937 else 6.70820393249937 end" "if [:liq].Comp_wt[1] <= 52.0 12.5499003980111 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 4.47213595499958 else 4.47213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.223606797749979 else 0.223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.316227766016838 else 0.316227766016838 end"; "if [:liq].Comp_wt[1] <= 52.0 4.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.16227766016838 else 2.0 end" "if [:liq].Comp_wt[1] <= 52.0 6.70820393249937 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 5.47722557505166 else 4.47213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0474341649025257 else 0.0474341649025257 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0707106781186547 else 0.0707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 6.32455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.474341649025257 else 0.273861278752583 end" "if [:liq].Comp_wt[1] <= 52.0 1.89736659610103 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.316227766016838 else 0.187082869338697 end" "if [:liq].Comp_wt[1] <= 52.0 1.26491106406735 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0591607978309961 else 0.0144913767461894 end" "if [:liq].Comp_wt[1] <= 52.0 0.948683298050514 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.06 else 0.0189736659610103 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.0316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.790569415042095 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0648074069840786 else 0.0244948974278318 end" "if [:liq].Comp_wt[1] <= 52.0 0.223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.234520787991171 else 0.234520787991171 end" "if [:liq].Comp_wt[1] <= 52.0 0.790569415042095 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0670820393249937 else 0.031224989991992 end" "if [:liq].Comp_wt[1] <= 52.0 0.632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.229128784747792 else 0.0547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.316227766016838 else 0.0547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0707106781186547 else 0.0387298334620742 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.13228756555323 else 0.0418330013267038 end" "if [:liq].Comp_wt[1] <= 52.0 0.51234753829798 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0741619848709566 else 0.0447213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 0.521536192416212 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0774596669241483 else 0.0519615242270663 end" "if [:liq].Comp_wt[1] <= 52.0 0.519615242270663 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0793725393319377 else 0.0591607978309961 end" "if [:liq].Comp_wt[1] <= 52.0 0.774596669241483 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.111803398874989 else 0.066332495807108 end" "if [:liq].Comp_wt[1] <= 52.0 0.53851648071345 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0774596669241483 else 0.0692820323027551 end" "if [:liq].Comp_wt[1] <= 52.0 0.554977477020464 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0747663025700749 else 0.0747663025700749 end" "if [:liq].Comp_wt[1] <= 52.0 0.557225268630201 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0734846922834953 else 0.0793725393319377 end" "if [:liq].Comp_wt[1] <= 52.0 0.558569601750758 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0726636084983398 else 0.0819756061276768 end" "if [:liq].Comp_wt[1] <= 52.0 0.559016994374947 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0707106781186547 else 0.0866025403784439 end" "if [:liq].Comp_wt[1] <= 52.0 3.16227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.16227766016838 else 3.16227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 9.21954445729289 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 9.21954445729289 else 9.21954445729289 end"; "if [:liq].Comp_wt[1] <= 52.0 0.106066017177982 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.106066017177982 else 0.106066017177982 end" "if [:liq].Comp_wt[1] <= 52.0 0.0223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0223606797749979 else 0.0223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.193649167310371 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.193649167310371 else 0.193649167310371 end" "if [:liq].Comp_wt[1] <= 52.0 0.632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.632455532033676 else 0.632455532033676 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.0316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.0316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.0741619848709566 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0741619848709566 else 0.0741619848709566 end" "if [:liq].Comp_wt[1] <= 52.0 0.0848528137423857 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0848528137423857 else 0.0848528137423857 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.0316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.0989949493661167 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0989949493661167 else 0.0989949493661167 end" "if [:liq].Comp_wt[1] <= 52.0 0.187082869338697 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.187082869338697 else 0.187082869338697 end" "if [:liq].Comp_wt[1] <= 52.0 0.12 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.12 else 0.12 end" "if [:liq].Comp_wt[1] <= 52.0 0.0612372435695794 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0612372435695794 else 0.0612372435695794 end" "if [:liq].Comp_wt[1] <= 52.0 0.0774596669241483 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0774596669241483 else 0.0774596669241483 end" "if [:liq].Comp_wt[1] <= 52.0 0.173205080756888 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.173205080756888 else 0.173205080756888 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.0316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.403112887414927 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.403112887414927 else 0.403112887414927 end" "if [:liq].Comp_wt[1] <= 52.0 0.670820393249937 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.670820393249937 else 0.670820393249937 end" "if [:liq].Comp_wt[1] <= 52.0 1.14564392373896 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.14564392373896 else 1.14564392373896 end" "if [:liq].Comp_wt[1] <= 52.0 0.987420882906575 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.987420882906575 else 0.987420882906575 end" "if [:liq].Comp_wt[1] <= 52.0 1.58113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.58113883008419 else 1.58113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 2.29128784747792 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.29128784747792 else 2.29128784747792 end" "if [:liq].Comp_wt[1] <= 52.0 2.80624304008046 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.80624304008046 else 2.80624304008046 end" "if [:liq].Comp_wt[1] <= 52.0 3.16227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.16227766016838 else 3.16227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 3.51781181986757 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.51781181986757 else 3.51781181986757 end" "if [:liq].Comp_wt[1] <= 52.0 0.447213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.447213595499958 else 0.447213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 2.23606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.23606797749979 else 2.23606797749979 end"; "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.00948683298050514 end" "if [:liq].Comp_wt[1] <= 52.0 0.122474487139159 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0547722557505166 else 0.005 end" "if [:liq].Comp_wt[1] <= 52.0 0.14142135623731 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0707106781186547 else 0.00935414346693485 end" "if [:liq].Comp_wt[1] <= 52.0 0.154919333848297 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0447213595499958 else 0.00707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 0.109544511501033 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.111803398874989 else 0.00894427190999916 end" "if [:liq].Comp_wt[1] <= 52.0 0.154919333848297 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0707106781186547 else 0.0180277563773199 end" "if [:liq].Comp_wt[1] <= 52.0 0.632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.13228756555323 else 0.0474341649025257 end" "if [:liq].Comp_wt[1] <= 52.0 1.09544511501033 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.212132034355964 else 0.0866025403784439 end" "if [:liq].Comp_wt[1] <= 52.0 0.223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.316227766016838 else 0.0223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 1.48323969741913 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.290688837074973 else 0.14456832294801 end" "if [:liq].Comp_wt[1] <= 52.0 0.335410196624968 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.287228132326901 else 0.122474487139159 end" "if [:liq].Comp_wt[1] <= 52.0 1.93649167310371 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.401870625948202 else 0.216333076527839 end" "if [:liq].Comp_wt[1] <= 52.0 0.387298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.234520787991171 else 0.0707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 0.547722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.31224989991992 else 0.187082869338697 end" "if [:liq].Comp_wt[1] <= 52.0 2.32379000772445 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.5 else 0.301662062579967 end" "if [:liq].Comp_wt[1] <= 52.0 2.77488738510232 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.692820323027551 else 0.387298334620742 end" "if [:liq].Comp_wt[1] <= 52.0 3.06186217847897 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.725603197346869 else 0.460977222864644 end" "if [:liq].Comp_wt[1] <= 52.0 3.46410161513775 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.848528137423857 else 0.522494019104525 end" "if [:liq].Comp_wt[1] <= 52.0 3.57770876399966 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.953939201416946 else 0.551361950083609 end" "if [:liq].Comp_wt[1] <= 52.0 3.24037034920393 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.09544511501033 else 0.6 end" "if [:liq].Comp_wt[1] <= 52.0 3.35410196624969 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.06348483769163 else 0.585662018573853 end" "if [:liq].Comp_wt[1] <= 52.0 3.1859064644148 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.14105214604767 else 0.608276253029822 end" "if [:liq].Comp_wt[1] <= 52.0 3.07408522978788 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.18215904175369 else 0.608276253029822 end" "if [:liq].Comp_wt[1] <= 52.0 2.89827534923789 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.1861703081767 else 0.585662018573853 end" "if [:liq].Comp_wt[1] <= 52.0 2.77488738510232 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.18321595661992 else 0.559910707166777 end" "if [:liq].Comp_wt[1] <= 52.0 4.74341649025257 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.93649167310371 else 1.58113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 26.4575131106459 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 7.74596669241484 else 4.74341649025257 end"; "if [:liq].Comp_wt[1] <= 52.0 0.0223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.0223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00447213595499958 else 0.00447213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 0.387298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0158113883008419 else 0.0158113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 0.387298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0237170824512628 else 0.0237170824512628 end" "if [:liq].Comp_wt[1] <= 52.0 28.2842712474619 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 4.47213595499958 else 4.47213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 44.7213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 6.70820393249937 else 6.70820393249937 end" "if [:liq].Comp_wt[1] <= 52.0 3.35410196624968 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0223606797749979 else 0.0223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 2.82842712474619 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0306186217847897 else 0.0306186217847897 end" "if [:liq].Comp_wt[1] <= 52.0 0.447213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0158113883008419 else 0.0158113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 2.57875939164553 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0367423461417476 else 0.0367423461417476 end" "if [:liq].Comp_wt[1] <= 52.0 0.193649167310371 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00632455532033676 else 0.00632455532033676 end" "if [:liq].Comp_wt[1] <= 52.0 2.50998007960223 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0412310562561766 else 0.0412310562561766 end" "if [:liq].Comp_wt[1] <= 52.0 1.58113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.22474487139159 else 1.22474487139159 end" "if [:liq].Comp_wt[1] <= 52.0 1.58113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.5 else 1.5 end" "if [:liq].Comp_wt[1] <= 52.0 2.36643191323985 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0547722557505166 else 0.0547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0447213595499958 else 0.0447213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 2.36643191323985 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0790569415042095 else 0.0790569415042095 end" "if [:liq].Comp_wt[1] <= 52.0 1.89736659610103 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.09 else 0.09 end" "if [:liq].Comp_wt[1] <= 52.0 1.58113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0989949493661167 else 0.0989949493661167 end" "if [:liq].Comp_wt[1] <= 52.0 0.447213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0707106781186547 else 0.0707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 1.58113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.101980390271856 else 0.101980390271856 end" "if [:liq].Comp_wt[1] <= 52.0 1.58113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.111803398874989 else 0.111803398874989 end" "if [:liq].Comp_wt[1] <= 52.0 1.42302494707577 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.117260393995586 else 0.117260393995586 end" "if [:liq].Comp_wt[1] <= 52.0 1.26491106406735 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.117473401244707 else 0.117473401244707 end" "if [:liq].Comp_wt[1] <= 52.0 1.26491106406735 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.109544511501033 else 0.109544511501033 end" "if [:liq].Comp_wt[1] <= 52.0 31.6227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 11.180339887499 else 11.180339887499 end" "if [:liq].Comp_wt[1] <= 52.0 4.47213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.73205080756888 else 1.73205080756888 end"; "if [:liq].Comp_wt[1] <= 52.0 0.00790569415042094 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00836660026534076 else 0.000866025403784438 end" "if [:liq].Comp_wt[1] <= 52.0 0.0173205080756888 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0154919333848297 else 0.000387298334620741 end" "if [:liq].Comp_wt[1] <= 52.0 0.0774596669241483 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0474341649025257 else 0.002 end" "if [:liq].Comp_wt[1] <= 52.0 0.158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0935414346693485 else 0.00612372435695794 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0264575131106459 else 0.013228756555323 end" "if [:liq].Comp_wt[1] <= 52.0 0.0707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.05 else 0.0180277563773199 end" "if [:liq].Comp_wt[1] <= 52.0 0.1 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0264575131106459 else 0.00894427190999916 end" "if [:liq].Comp_wt[1] <= 52.0 0.237170824512628 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0519615242270663 else 0.0219089023002067 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.0122474487139159 end" "if [:liq].Comp_wt[1] <= 52.0 0.460977222864644 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.102469507659596 else 0.0565685424949238 end" "if [:liq].Comp_wt[1] <= 52.0 0.0223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0173205080756888 else 0.013228756555323 end" "if [:liq].Comp_wt[1] <= 52.0 0.866025403784439 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.205548047910945 else 0.122474487139159 end" "if [:liq].Comp_wt[1] <= 52.0 0.447213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.547722557505166 else 0.3 end" "if [:liq].Comp_wt[1] <= 52.0 0.707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.387298334620742 else 0.252487623459052 end" "if [:liq].Comp_wt[1] <= 52.0 1.93649167310371 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.433012701892219 else 0.264575131106459 end" "if [:liq].Comp_wt[1] <= 52.0 1.58113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.5 else 0.424264068711928 end" "if [:liq].Comp_wt[1] <= 52.0 6.12372435695795 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.58113883008419 else 0.821583836257749 end" "if [:liq].Comp_wt[1] <= 52.0 11.180339887499 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.16227766016838 else 1.3228756555323 end" "if [:liq].Comp_wt[1] <= 52.0 22.3606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 5.45435605731786 else 1.93649167310371 end" "if [:liq].Comp_wt[1] <= 52.0 25.4950975679639 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 6.12372435695795 else 2.32379000772445 end" "if [:liq].Comp_wt[1] <= 52.0 32.4037034920393 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 8.2915619758885 else 2.54950975679639 end" "if [:liq].Comp_wt[1] <= 52.0 41.2310562561766 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 9.89949493661167 else 3.1224989991992 end" "if [:liq].Comp_wt[1] <= 52.0 44.7213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 10.0623058987491 else 3.51781181986757 end" "if [:liq].Comp_wt[1] <= 52.0 44.7213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 9.35414346693486 else 4.02336923485777 end" "if [:liq].Comp_wt[1] <= 52.0 37.4165738677394 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 7.74596669241484 else 4.47213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 4.47213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.64575131106459 else 2.44948974278318 end" "if [:liq].Comp_wt[1] <= 52.0 22.3606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.16227766016838 else 3.74165738677394 end"; "if [:liq].Comp_wt[1] <= 52.0 0.948683298050514 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.948683298050514 else 0.948683298050514 end" "if [:liq].Comp_wt[1] <= 52.0 9.48683298050514 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 5.47722557505166 else 5.47722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.0144913767461894 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0173205080756888 else 0.0173205080756888 end" "if [:liq].Comp_wt[1] <= 52.0 0.0158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0273861278752583 else 0.0273861278752583 end" "if [:liq].Comp_wt[1] <= 52.0 0.0447213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0158113883008419 else 0.0158113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 0.00707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.14142135623731 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.111803398874989 else 0.111803398874989 end" "if [:liq].Comp_wt[1] <= 52.0 0.0866025403784439 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0670820393249937 else 0.0670820393249937 end" "if [:liq].Comp_wt[1] <= 52.0 0.948683298050514 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.223606797749979 else 0.223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.0547722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0447213595499958 else 0.0447213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 3.87298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.23606797749979 else 2.23606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.03 else 0.03 end" "if [:liq].Comp_wt[1] <= 52.0 0.0193649167310371 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0223606797749979 else 0.0223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.0158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0223606797749979 else 0.0223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.0187082869338697 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0212132034355964 else 0.0212132034355964 end" "if [:liq].Comp_wt[1] <= 52.0 3.87298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.35410196624969 else 3.35410196624969 end" "if [:liq].Comp_wt[1] <= 52.0 0.0187082869338697 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0167332005306815 else 0.0167332005306815 end" "if [:liq].Comp_wt[1] <= 52.0 0.0167332005306815 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.015 else 0.015 end" "if [:liq].Comp_wt[1] <= 52.0 0.02 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0144913767461894 else 0.0144913767461894 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0254950975679639 else 0.0254950975679639 end" "if [:liq].Comp_wt[1] <= 52.0 0.02 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0162018517460197 else 0.0162018517460197 end" "if [:liq].Comp_wt[1] <= 52.0 0.0223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0189736659610103 else 0.0189736659610103 end" "if [:liq].Comp_wt[1] <= 52.0 0.0234520787991171 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.02 else 0.02 end" "if [:liq].Comp_wt[1] <= 52.0 0.0234520787991171 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0232379000772445 else 0.0232379000772445 end" "if [:liq].Comp_wt[1] <= 52.0 0.0244948974278318 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0273861278752583 else 0.0273861278752583 end" "if [:liq].Comp_wt[1] <= 52.0 0.223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.111803398874989 else 0.111803398874989 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.02 else 0.02 end"; "if [:liq].Comp_wt[1] <= 52.0 0.0223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0387298334620742 else 0.0158113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 0.0223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.158113883008419 else 0.0158113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 0.221359436211787 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0935414346693485 else 0.0387298334620742 end" "if [:liq].Comp_wt[1] <= 52.0 0.447213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.13228756555323 else 0.0387298334620742 end" "if [:liq].Comp_wt[1] <= 52.0 1.67332005306815 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.136930639376291 else 0.1 end" "if [:liq].Comp_wt[1] <= 52.0 2.12132034355964 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.122474487139159 else 0.1 end" "if [:liq].Comp_wt[1] <= 52.0 1.02469507659596 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.14142135623731 else 0.0126491106406735 end" "if [:liq].Comp_wt[1] <= 52.0 1.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.189736659610103 else 0.0273861278752583 end" "if [:liq].Comp_wt[1] <= 52.0 0.547722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.316227766016838 else 0.0591607978309961 end" "if [:liq].Comp_wt[1] <= 52.0 1.06066017177982 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.234520787991171 else 0.0474341649025257 end" "if [:liq].Comp_wt[1] <= 52.0 0.0264575131106459 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0223606797749979 else 0.0122474487139159 end" "if [:liq].Comp_wt[1] <= 52.0 1.09544511501033 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.267394839142419 else 0.062449979983984 end" "if [:liq].Comp_wt[1] <= 52.0 0.670820393249937 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.220794021658196 else 0.150831031289984 end" "if [:liq].Comp_wt[1] <= 52.0 0.707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.158113883008419 else 0.2 end" "if [:liq].Comp_wt[1] <= 52.0 1.09544511501033 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.3 else 0.0724568837309471 end" "if [:liq].Comp_wt[1] <= 52.0 0.948683298050514 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.254950975679639 else 0.0360555127546399 end" "if [:liq].Comp_wt[1] <= 52.0 1.18321595661992 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.324037034920393 else 0.0774596669241483 end" "if [:liq].Comp_wt[1] <= 52.0 1.18321595661992 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.31224989991992 else 0.0702851335632223 end" "if [:liq].Comp_wt[1] <= 52.0 1.10679718105893 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.301662062579967 else 0.0666333249958307 end" "if [:liq].Comp_wt[1] <= 52.0 0.948683298050514 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.335410196624968 else 0.0591607978309961 end" "if [:liq].Comp_wt[1] <= 52.0 1.02469507659596 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.268328157299975 else 0.0620483682299543 end" "if [:liq].Comp_wt[1] <= 52.0 0.866025403784439 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.245967477524977 else 0.0565685424949238 end" "if [:liq].Comp_wt[1] <= 52.0 0.866025403784439 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.234520787991171 else 0.0547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.223606797749979 else 0.05 end" "if [:liq].Comp_wt[1] <= 52.0 0.547722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.212132034355964 else 0.0447213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 36.0555127546399 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 22.3606797749979 else 6.32455532033676 end" "if [:liq].Comp_wt[1] <= 52.0 8.66025403784439 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.16227766016838 else 1.4142135623731 end"; "if [:liq].Comp_wt[1] <= 52.0 0.0223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0387298334620742 else 0.0158113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 0.0223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.158113883008419 else 0.0158113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 0.221359436211787 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0935414346693485 else 0.0387298334620742 end" "if [:liq].Comp_wt[1] <= 52.0 0.447213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.13228756555323 else 0.0387298334620742 end" "if [:liq].Comp_wt[1] <= 52.0 1.67332005306815 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.136930639376291 else 0.1 end" "if [:liq].Comp_wt[1] <= 52.0 2.12132034355964 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.122474487139159 else 0.1 end" "if [:liq].Comp_wt[1] <= 52.0 1.02469507659596 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.14142135623731 else 0.0126491106406735 end" "if [:liq].Comp_wt[1] <= 52.0 1.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.189736659610103 else 0.0273861278752583 end" "if [:liq].Comp_wt[1] <= 52.0 0.547722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.316227766016838 else 0.0591607978309961 end" "if [:liq].Comp_wt[1] <= 52.0 1.06066017177982 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.234520787991171 else 0.0474341649025257 end" "if [:liq].Comp_wt[1] <= 52.0 0.0264575131106459 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0223606797749979 else 0.0122474487139159 end" "if [:liq].Comp_wt[1] <= 52.0 1.09544511501033 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.267394839142419 else 0.062449979983984 end" "if [:liq].Comp_wt[1] <= 52.0 0.670820393249937 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.220794021658196 else 0.150831031289984 end" "if [:liq].Comp_wt[1] <= 52.0 0.707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.158113883008419 else 0.2 end" "if [:liq].Comp_wt[1] <= 52.0 1.09544511501033 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.3 else 0.0724568837309471 end" "if [:liq].Comp_wt[1] <= 52.0 0.948683298050514 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.254950975679639 else 0.0360555127546399 end" "if [:liq].Comp_wt[1] <= 52.0 1.18321595661992 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.324037034920393 else 0.0774596669241483 end" "if [:liq].Comp_wt[1] <= 52.0 1.18321595661992 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.31224989991992 else 0.0702851335632223 end" "if [:liq].Comp_wt[1] <= 52.0 1.10679718105893 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.301662062579967 else 0.0666333249958307 end" "if [:liq].Comp_wt[1] <= 52.0 0.948683298050514 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.335410196624968 else 0.0591607978309961 end" "if [:liq].Comp_wt[1] <= 52.0 1.02469507659596 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.268328157299975 else 0.0620483682299543 end" "if [:liq].Comp_wt[1] <= 52.0 0.866025403784439 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.245967477524977 else 0.0565685424949238 end" "if [:liq].Comp_wt[1] <= 52.0 0.866025403784439 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.234520787991171 else 0.0547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.223606797749979 else 0.05 end" "if [:liq].Comp_wt[1] <= 52.0 0.547722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.212132034355964 else 0.0447213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 36.0555127546399 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 22.3606797749979 else 6.32455532033676 end" "if [:liq].Comp_wt[1] <= 52.0 8.66025403784439 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.16227766016838 else 1.4142135623731 end"; "if [:liq].Comp_wt[1] <= 52.0 0.0223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.013228756555323 else 0.00223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.03 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0122474487139159 else 0.00223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.02 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0154919333848297 else 0.000707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 0.0154919333848297 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0223606797749979 else 0.000707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 0.031224989991992 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0193649167310371 else 0.00273861278752583 end" "if [:liq].Comp_wt[1] <= 52.0 0.0335410196624968 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0273861278752583 else 0.005 end" "if [:liq].Comp_wt[1] <= 52.0 0.0670820393249937 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00273861278752583 else 0.000223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.0774596669241483 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00489897948556636 else 0.000447213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 0.158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0187082869338697 else 0.000193649167310371 end" "if [:liq].Comp_wt[1] <= 52.0 0.0916515138991168 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00790569415042094 else 0.000935414346693485 end" "if [:liq].Comp_wt[1] <= 52.0 0.13228756555323 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.01 else 0.000109544511501033 end" "if [:liq].Comp_wt[1] <= 52.0 0.106066017177982 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0116189500386222 else 0.00158113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 0.025298221281347 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0387298334620742 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.0284604989415154 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0223606797749979 else 0.00547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.121963109176505 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0152970585407783 else 0.00250998007960223 end" "if [:liq].Comp_wt[1] <= 52.0 0.14142135623731 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0187616630392937 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.13228756555323 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0228035085019828 else 0.00424264068711928 end" "if [:liq].Comp_wt[1] <= 52.0 0.154919333848297 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0270554985169374 else 0.0062928530890209 end" "if [:liq].Comp_wt[1] <= 52.0 0.175499287747842 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0324037034920393 else 0.00908295106229247 end" "if [:liq].Comp_wt[1] <= 52.0 0.126491106406735 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0346410161513775 else 0.00774596669241483 end" "if [:liq].Comp_wt[1] <= 52.0 0.193649167310371 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0369864840178138 else 0.013228756555323 end" "if [:liq].Comp_wt[1] <= 52.0 0.232379000772445 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.04 else 0.0193649167310371 end" "if [:liq].Comp_wt[1] <= 52.0 0.277488738510232 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0439545219516718 else 0.0241867732448956 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0474341649025257 else 0.0273861278752583 end" "if [:liq].Comp_wt[1] <= 52.0 0.387298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0529150262212918 else 0.0282842712474619 end" "if [:liq].Comp_wt[1] <= 52.0 0.244948974278318 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.16431676725155 else 0.0774596669241483 end" "if [:liq].Comp_wt[1] <= 52.0 0.447213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.433012701892219 else 0.287228132326901 end"; "if [:liq].Comp_wt[1] <= 52.0 0.0158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0244948974278318 else 0.00547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.05 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.00173205080756888 end" "if [:liq].Comp_wt[1] <= 52.0 0.158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.1 else 0.000316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.0005 end" "if [:liq].Comp_wt[1] <= 52.0 0.25298221281347 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.05 else 0.00346410161513776 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0670820393249937 else 0.00632455532033676 end" "if [:liq].Comp_wt[1] <= 52.0 0.158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00894427190999916 else 0.000836660026534075 end" "if [:liq].Comp_wt[1] <= 52.0 0.239791576165636 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0204939015319192 else 0.00264575131106459 end" "if [:liq].Comp_wt[1] <= 52.0 0.0894427190999916 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.111803398874989 else 0.02 end" "if [:liq].Comp_wt[1] <= 52.0 0.289827534923789 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0363318042491699 else 0.00547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.0447213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0193649167310371 else 0.0141421356237309 end" "if [:liq].Comp_wt[1] <= 52.0 0.353553390593274 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0574456264653803 else 0.0109544511501033 end" "if [:liq].Comp_wt[1] <= 52.0 0.0547722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0346410161513775 else 0.0291547594742265 end" "if [:liq].Comp_wt[1] <= 52.0 0.0948683298050513 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0547722557505166 else 0.0547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.424264068711928 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0812403840463596 else 0.0180277563773199 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.114017542509914 else 0.0223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.58309518948453 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.157480157480236 else 0.0324037034920393 end" "if [:liq].Comp_wt[1] <= 52.0 0.724568837309472 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.207846096908265 else 0.0412310562561766 end" "if [:liq].Comp_wt[1] <= 52.0 0.848528137423857 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.256124969497314 else 0.05 end" "if [:liq].Comp_wt[1] <= 52.0 0.774596669241483 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.282842712474619 else 0.0670820393249937 end" "if [:liq].Comp_wt[1] <= 52.0 0.924662100445346 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.311126983722081 else 0.0619677335393187 end" "if [:liq].Comp_wt[1] <= 52.0 1.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.351425667816112 else 0.0734846922834953 end" "if [:liq].Comp_wt[1] <= 52.0 1.16189500386223 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.4 else 0.0916515138991168 end" "if [:liq].Comp_wt[1] <= 52.0 1.2747548783982 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.458257569495584 else 0.107238052947636 end" "if [:liq].Comp_wt[1] <= 52.0 1.54919333848297 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.5 else 0.120415945787923 end" "if [:liq].Comp_wt[1] <= 52.0 3.16227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.5 else 1.0 end" "if [:liq].Comp_wt[1] <= 52.0 22.3606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 5.47722557505166 else 1.22474487139159 end"; "if [:liq].Comp_wt[1] <= 52.0 0.126491106406735 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0447213595499958 else 0.0547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.387298334620742 else 0.335410196624968 end" "if [:liq].Comp_wt[1] <= 52.0 0.0547722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.1 end" "if [:liq].Comp_wt[1] <= 52.0 0.1 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0547722557505166 else 0.0632455532033676 end" "if [:liq].Comp_wt[1] <= 52.0 0.1 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0612372435695794 else 0.0447213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 0.0670820393249937 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0447213595499958 else 0.0447213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 0.387298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.2 else 0.111803398874989 end" "if [:liq].Comp_wt[1] <= 52.0 0.284604989415154 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.154919333848297 else 0.0894427190999916 end" "if [:liq].Comp_wt[1] <= 52.0 0.670820393249937 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.5 else 0.591607978309962 end" "if [:liq].Comp_wt[1] <= 52.0 0.221359436211786 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.122474487139159 else 0.0714142842854285 end" "if [:liq].Comp_wt[1] <= 52.0 4.47213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.44948974278318 else 1.6583123951777 end" "if [:liq].Comp_wt[1] <= 52.0 0.158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.106066017177982 else 0.0591607978309961 end" "if [:liq].Comp_wt[1] <= 52.0 0.0894427190999916 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0223606797749979 else 0.01 end" "if [:liq].Comp_wt[1] <= 52.0 0.0670820393249937 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0324037034920393 else 0.0154919333848297 end" "if [:liq].Comp_wt[1] <= 52.0 0.126491106406735 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0774596669241483 else 0.0489897948556635 end" "if [:liq].Comp_wt[1] <= 52.0 2.73861278752583 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.1180339887499 else 0.707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 0.116189500386222 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0447213595499958 else 0.0387298334620742 end" "if [:liq].Comp_wt[1] <= 52.0 0.102469507659596 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0367423461417476 else 0.0338230690505755 end" "if [:liq].Comp_wt[1] <= 52.0 0.0948683298050513 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.03 else 0.03 end" "if [:liq].Comp_wt[1] <= 52.0 0.1 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0292403830344269 else 0.0357071421427143 end" "if [:liq].Comp_wt[1] <= 52.0 0.0774596669241483 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0256904651573302 else 0.0267394839142419 end" "if [:liq].Comp_wt[1] <= 52.0 0.0692820323027551 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0212132034355964 else 0.0240831891575846 end" "if [:liq].Comp_wt[1] <= 52.0 0.0574456264653803 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0182345825288105 else 0.0222261107708929 end" "if [:liq].Comp_wt[1] <= 52.0 0.0524404424085076 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0169705627484771 else 0.0203469899493758 end" "if [:liq].Comp_wt[1] <= 52.0 0.0447213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0154919333848297 else 0.0178885438199983 end" "if [:liq].Comp_wt[1] <= 52.0 0.158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.158113883008419 else 0.0670820393249937 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.0316227766016838 end"; "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.00316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00316227766016838 else 0.00316227766016838 end"; "if [:liq].Comp_wt[1] <= 52.0 0.0118321595661992 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.01 else 0.01 end" "if [:liq].Comp_wt[1] <= 52.0 0.01 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00866025403784438 else 0.00866025403784438 end" "if [:liq].Comp_wt[1] <= 52.0 0.335410196624968 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.162018517460196 else 0.162018517460196 end" "if [:liq].Comp_wt[1] <= 52.0 0.335410196624968 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.2 else 0.2 end" "if [:liq].Comp_wt[1] <= 52.0 61.2372435695795 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 33.5410196624969 else 33.5410196624969 end" "if [:liq].Comp_wt[1] <= 52.0 27.3861278752583 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 22.3606797749979 else 22.3606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.0353553390593274 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0111803398874989 else 0.0111803398874989 end" "if [:liq].Comp_wt[1] <= 52.0 0.0458257569495584 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0141421356237309 else 0.0141421356237309 end" "if [:liq].Comp_wt[1] <= 52.0 0.0316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0273861278752583 else 0.0273861278752583 end" "if [:liq].Comp_wt[1] <= 52.0 0.0494974746830583 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0156524758424985 else 0.0156524758424985 end" "if [:liq].Comp_wt[1] <= 52.0 0.0707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0316227766016838 else 0.0316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.0565685424949238 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0172481883106603 else 0.0172481883106603 end" "if [:liq].Comp_wt[1] <= 52.0 3.46410161513776 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.23606797749979 else 2.23606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 5.29150262212918 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.74165738677394 else 3.74165738677394 end" "if [:liq].Comp_wt[1] <= 52.0 0.0565685424949238 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0184932420089069 else 0.0184932420089069 end" "if [:liq].Comp_wt[1] <= 52.0 0.00111803398874989 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.000632455532033676 else 0.000632455532033676 end" "if [:liq].Comp_wt[1] <= 52.0 0.0354964786985977 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0194935886896179 else 0.0194935886896179 end" "if [:liq].Comp_wt[1] <= 52.0 0.0278208554864871 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0204939015319192 else 0.0204939015319192 end" "if [:liq].Comp_wt[1] <= 52.0 0.0256904651573302 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.022248595461287 else 0.022248595461287 end" "if [:liq].Comp_wt[1] <= 52.0 0.0707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0273861278752583 else 0.0273861278752583 end" "if [:liq].Comp_wt[1] <= 52.0 0.0241867732448956 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.024 else 0.024 end" "if [:liq].Comp_wt[1] <= 52.0 0.0217944947177034 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0264575131106459 else 0.0264575131106459 end" "if [:liq].Comp_wt[1] <= 52.0 0.0185202591774521 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0291204395571221 else 0.0291204395571221 end" "if [:liq].Comp_wt[1] <= 52.0 0.0173205080756888 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0314642654451045 else 0.0314642654451045 end" "if [:liq].Comp_wt[1] <= 52.0 0.0158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.033466401061363 else 0.033466401061363 end" "if [:liq].Comp_wt[1] <= 52.0 47.4341649025257 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.16227766016838 else 3.16227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.273861278752583 else 0.273861278752583 end"; "if [:liq].Comp_wt[1] <= 52.0 0.00547722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00547722557505166 else 0.00547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.00223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00223606797749979 else 0.00223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.00790569415042094 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00790569415042094 else 0.00790569415042094 end" "if [:liq].Comp_wt[1] <= 52.0 0.013228756555323 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.013228756555323 else 0.013228756555323 end" "if [:liq].Comp_wt[1] <= 52.0 0.0707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0707106781186547 else 0.0707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 0.0707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0707106781186547 else 0.0707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 0.00244948974278318 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00244948974278318 else 0.00244948974278318 end" "if [:liq].Comp_wt[1] <= 52.0 0.00387298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00387298334620742 else 0.00387298334620742 end" "if [:liq].Comp_wt[1] <= 52.0 0.000316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.000316227766016838 else 0.000316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.006 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.006 else 0.006 end" "if [:liq].Comp_wt[1] <= 52.0 0.00346410161513776 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00346410161513776 else 0.00346410161513776 end" "if [:liq].Comp_wt[1] <= 52.0 0.00812403840463596 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00812403840463596 else 0.00812403840463596 end" "if [:liq].Comp_wt[1] <= 52.0 0.0790569415042095 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0790569415042095 else 0.0790569415042095 end" "if [:liq].Comp_wt[1] <= 52.0 0.1 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.1 else 0.1 end" "if [:liq].Comp_wt[1] <= 52.0 0.0116189500386222 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0116189500386222 else 0.0116189500386222 end" "if [:liq].Comp_wt[1] <= 52.0 0.0106066017177982 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0106066017177982 else 0.0106066017177982 end" "if [:liq].Comp_wt[1] <= 52.0 0.0144913767461894 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0144913767461894 else 0.0144913767461894 end" "if [:liq].Comp_wt[1] <= 52.0 0.0149666295470958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0149666295470958 else 0.0149666295470958 end" "if [:liq].Comp_wt[1] <= 52.0 0.015 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.015 else 0.015 end" "if [:liq].Comp_wt[1] <= 52.0 0.0158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0158113883008419 else 0.0158113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 0.0141421356237309 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0141421356237309 else 0.0141421356237309 end" "if [:liq].Comp_wt[1] <= 52.0 0.0132664991614216 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0132664991614216 else 0.0132664991614216 end" "if [:liq].Comp_wt[1] <= 52.0 0.0130384048104053 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0130384048104053 else 0.0130384048104053 end" "if [:liq].Comp_wt[1] <= 52.0 0.012369316876853 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.012369316876853 else 0.012369316876853 end" "if [:liq].Comp_wt[1] <= 52.0 0.0116189500386222 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0116189500386222 else 0.0116189500386222 end" "if [:liq].Comp_wt[1] <= 52.0 8.36660026534076 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 8.36660026534076 else 8.36660026534076 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.316227766016838 else 0.316227766016838 end"; "if [:liq].Comp_wt[1] <= 52.0 0.00547722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00547722557505166 else 0.00547722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 0.00223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00223606797749979 else 0.00223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.00790569415042094 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00790569415042094 else 0.00790569415042094 end" "if [:liq].Comp_wt[1] <= 52.0 0.013228756555323 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.013228756555323 else 0.013228756555323 end" "if [:liq].Comp_wt[1] <= 52.0 0.0707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0707106781186547 else 0.0707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 0.0707106781186547 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0707106781186547 else 0.0707106781186547 end" "if [:liq].Comp_wt[1] <= 52.0 0.00244948974278318 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00244948974278318 else 0.00244948974278318 end" "if [:liq].Comp_wt[1] <= 52.0 0.00387298334620742 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00387298334620742 else 0.00387298334620742 end" "if [:liq].Comp_wt[1] <= 52.0 0.000316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.000316227766016838 else 0.000316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 0.006 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.006 else 0.006 end" "if [:liq].Comp_wt[1] <= 52.0 0.00346410161513776 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00346410161513776 else 0.00346410161513776 end" "if [:liq].Comp_wt[1] <= 52.0 0.00812403840463596 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.00812403840463596 else 0.00812403840463596 end" "if [:liq].Comp_wt[1] <= 52.0 0.0790569415042095 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0790569415042095 else 0.0790569415042095 end" "if [:liq].Comp_wt[1] <= 52.0 0.1 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.1 else 0.1 end" "if [:liq].Comp_wt[1] <= 52.0 0.0116189500386222 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0116189500386222 else 0.0116189500386222 end" "if [:liq].Comp_wt[1] <= 52.0 0.0106066017177982 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0106066017177982 else 0.0106066017177982 end" "if [:liq].Comp_wt[1] <= 52.0 0.0144913767461894 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0144913767461894 else 0.0144913767461894 end" "if [:liq].Comp_wt[1] <= 52.0 0.0149666295470958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0149666295470958 else 0.0149666295470958 end" "if [:liq].Comp_wt[1] <= 52.0 0.015 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.015 else 0.015 end" "if [:liq].Comp_wt[1] <= 52.0 0.0158113883008419 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0158113883008419 else 0.0158113883008419 end" "if [:liq].Comp_wt[1] <= 52.0 0.0141421356237309 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0141421356237309 else 0.0141421356237309 end" "if [:liq].Comp_wt[1] <= 52.0 0.0132664991614216 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0132664991614216 else 0.0132664991614216 end" "if [:liq].Comp_wt[1] <= 52.0 0.0130384048104053 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0130384048104053 else 0.0130384048104053 end" "if [:liq].Comp_wt[1] <= 52.0 0.012369316876853 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.012369316876853 else 0.012369316876853 end" "if [:liq].Comp_wt[1] <= 52.0 0.0116189500386222 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.0116189500386222 else 0.0116189500386222 end" "if [:liq].Comp_wt[1] <= 52.0 8.36660026534076 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 8.36660026534076 else 8.36660026534076 end" "if [:liq].Comp_wt[1] <= 52.0 0.316227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.316227766016838 else 0.316227766016838 end"; "if [:liq].Comp_wt[1] <= 52.0 0.4 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.4 else 0.4 end" "if [:liq].Comp_wt[1] <= 52.0 1.36930639376292 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.36930639376292 else 1.36930639376292 end" "if [:liq].Comp_wt[1] <= 52.0 0.223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.223606797749979 else 0.223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.2 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.2 else 0.2 end" "if [:liq].Comp_wt[1] <= 52.0 3.46410161513776 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.46410161513776 else 3.46410161513776 end" "if [:liq].Comp_wt[1] <= 52.0 8.66025403784439 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 8.66025403784439 else 8.66025403784439 end" "if [:liq].Comp_wt[1] <= 52.0 6.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 6.0 else 6.0 end" "if [:liq].Comp_wt[1] <= 52.0 7.41619848709567 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 7.41619848709567 else 7.41619848709567 end" "if [:liq].Comp_wt[1] <= 52.0 0.223606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.223606797749979 else 0.223606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 8.83176086632785 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 8.83176086632785 else 8.83176086632785 end" "if [:liq].Comp_wt[1] <= 52.0 2.23606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.23606797749979 else 2.23606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 10.2469507659596 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 10.2469507659596 else 10.2469507659596 end" "if [:liq].Comp_wt[1] <= 52.0 1.34164078649987 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.34164078649987 else 1.34164078649987 end" "if [:liq].Comp_wt[1] <= 52.0 1.73205080756888 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.73205080756888 else 1.73205080756888 end" "if [:liq].Comp_wt[1] <= 52.0 10.9544511501033 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 10.9544511501033 else 10.9544511501033 end" "if [:liq].Comp_wt[1] <= 52.0 8.06225774829855 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 8.06225774829855 else 8.06225774829855 end" "if [:liq].Comp_wt[1] <= 52.0 10.2469507659596 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 10.2469507659596 else 10.2469507659596 end" "if [:liq].Comp_wt[1] <= 52.0 9.2951600308978 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 9.2951600308978 else 9.2951600308978 end" "if [:liq].Comp_wt[1] <= 52.0 7.88035532193822 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 7.88035532193822 else 7.88035532193822 end" "if [:liq].Comp_wt[1] <= 52.0 6.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 6.0 else 6.0 end" "if [:liq].Comp_wt[1] <= 52.0 6.557438524302 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 6.557438524302 else 6.557438524302 end" "if [:liq].Comp_wt[1] <= 52.0 5.17010638188423 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 5.17010638188423 else 5.17010638188423 end" "if [:liq].Comp_wt[1] <= 52.0 3.73496987939662 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.73496987939662 else 3.73496987939662 end" "if [:liq].Comp_wt[1] <= 52.0 2.56904651573303 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.56904651573303 else 2.56904651573303 end" "if [:liq].Comp_wt[1] <= 52.0 1.67332005306815 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.67332005306815 else 1.67332005306815 end" "if [:liq].Comp_wt[1] <= 52.0 5.47722557505166 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 5.47722557505166 else 5.47722557505166 end" "if [:liq].Comp_wt[1] <= 52.0 2.23606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.23606797749979 else 2.23606797749979 end"; "if [:liq].Comp_wt[1] <= 52.0 0.632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.632455532033676 else 0.632455532033676 end" "if [:liq].Comp_wt[1] <= 52.0 0.632455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.632455532033676 else 0.632455532033676 end" "if [:liq].Comp_wt[1] <= 52.0 18.02775637732 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 18.02775637732 else 18.02775637732 end" "if [:liq].Comp_wt[1] <= 52.0 31.6227766016838 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 31.6227766016838 else 31.6227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 22.3606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 22.3606797749979 else 22.3606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 22.3606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 22.3606797749979 else 22.3606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 1.26491106406735 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.26491106406735 else 1.26491106406735 end" "if [:liq].Comp_wt[1] <= 52.0 3.46410161513776 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 3.46410161513776 else 0.1 end" "if [:liq].Comp_wt[1] <= 52.0 0.1 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.1 else 0.316227766016838 end" "if [:liq].Comp_wt[1] <= 52.0 2.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.0 else 2.0 end" "if [:liq].Comp_wt[1] <= 52.0 4.47213595499958 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 4.47213595499958 else 4.47213595499958 end" "if [:liq].Comp_wt[1] <= 52.0 2.82842712474619 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.82842712474619 else 2.82842712474619 end" "if [:liq].Comp_wt[1] <= 52.0 948.683298050515 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 948.683298050515 else 948.683298050515 end" "if [:liq].Comp_wt[1] <= 52.0 948.683298050515 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 948.683298050515 else 948.683298050515 end" "if [:liq].Comp_wt[1] <= 52.0 7.74596669241484 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 7.74596669241484 else 7.74596669241484 end" "if [:liq].Comp_wt[1] <= 52.0 2.23606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 2.23606797749979 else 2.23606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 20.4939015319192 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 20.4939015319192 else 20.4939015319192 end" "if [:liq].Comp_wt[1] <= 52.0 28.4604989415154 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 28.4604989415154 else 28.4604989415154 end" "if [:liq].Comp_wt[1] <= 52.0 44.1588043316392 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 44.1588043316392 else 44.1588043316392 end" "if [:liq].Comp_wt[1] <= 52.0 67.0820393249938 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 67.0820393249938 else 67.0820393249938 end" "if [:liq].Comp_wt[1] <= 52.0 73.4846922834954 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 73.4846922834954 else 73.4846922834954 end" "if [:liq].Comp_wt[1] <= 52.0 110.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 110.0 else 110.0 end" "if [:liq].Comp_wt[1] <= 52.0 139.64240043769 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 139.64240043769 else 139.64240043769 end" "if [:liq].Comp_wt[1] <= 52.0 173.205080756888 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 173.205080756888 else 173.205080756888 end" "if [:liq].Comp_wt[1] <= 52.0 223.606797749979 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 223.606797749979 else 223.606797749979 end" "if [:liq].Comp_wt[1] <= 52.0 0.1 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 0.1 else 0.1 end" "if [:liq].Comp_wt[1] <= 52.0 63.2455532033676 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 63.2455532033676 else 63.2455532033676 end"; "if [:liq].Comp_wt[1] <= 52.0 0.0045 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 0.0045 end" "if [:liq].Comp_wt[1] <= 52.0 0.408 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 0.408 end" "if [:liq].Comp_wt[1] <= 52.0 156.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 156.0 end" "if [:liq].Comp_wt[1] <= 52.0 1.29 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.29 end" "if [:liq].Comp_wt[1] <= 52.0 0.226 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 0.226 end" "if [:liq].Comp_wt[1] <= 52.0 0.226 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 0.226 end" "if [:liq].Comp_wt[1] <= 52.0 2.05 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 2.05 end" "if [:liq].Comp_wt[1] <= 52.0 2.44 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 2.44 end" "if [:liq].Comp_wt[1] <= 52.0 0.5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 0.5 end" "if [:liq].Comp_wt[1] <= 52.0 2.85475042692001 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 2.89 end" "if [:liq].Comp_wt[1] <= 52.0 2.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 2.0 end" "if [:liq].Comp_wt[1] <= 52.0 3.34 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 3.78 end" "if [:liq].Comp_wt[1] <= 52.0 0.1 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 0.1 end" "if [:liq].Comp_wt[1] <= 52.0 10.0 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 10.0 end" "if [:liq].Comp_wt[1] <= 52.0 4.22 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 4.22 end" "if [:liq].Comp_wt[1] <= 52.0 3.78 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 3.78 end" "if [:liq].Comp_wt[1] <= 52.0 4.67 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 4.67 end" "if [:liq].Comp_wt[1] <= 52.0 4.58421203698084 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 4.59 end" "if [:liq].Comp_wt[1] <= 52.0 4.5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 4.5 end" "if [:liq].Comp_wt[1] <= 52.0 4.3 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 4.3 end" "if [:liq].Comp_wt[1] <= 52.0 4.12431812546026 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 4.14 end" "if [:liq].Comp_wt[1] <= 52.0 3.78 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 3.78 end" "if [:liq].Comp_wt[1] <= 52.0 3.34496636754392 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 3.37 end" "if [:liq].Comp_wt[1] <= 52.0 2.96 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 2.96 end" "if [:liq].Comp_wt[1] <= 52.0 2.61933874283863 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 2.55 end" "if [:liq].Comp_wt[1] <= 52.0 0.1 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 0.1 end" "if [:liq].Comp_wt[1] <= 52.0 0.0001 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 0.0001 end"; "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end"; "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end"; "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end" "if [:liq].Comp_wt[1] <= 52.0 1.0e-5 elseif [:liq].Comp_wt[1] <= 56.0 && [:liq].Comp_wt[1] > 52.0 1.0e-5 else 1.0e-5 end"]
        return (infos,el,ph,KDs)
    elseif tedb == "CO"
        return get_CO_KDs_database()
    end
end

"""
    mineral_classification(out, dtb)

    Classify the stable phases from a MAGEMin minimization result into mineralogical names compatible with the trace element partitioning coefficient database.

    Solution phases that straddle a solvus (e.g., feldspar, spinel, ilmenite) are disambiguated using their compositional variables. Duplicate phases resulting from solvus splitting are merged by summing their weight fractions.

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output structure.
    dtb : String
        Database identifier (e.g., "ig", "igad", "mp", "mpe", "mb", "ume", "mbe").

    Returns
    -------
    ph : Vector{String}
        Classified phase names compatible with the TE partitioning database.
    ph_wt : Vector{Float64}
        Weight fractions corresponding to each phase.
"""
function mineral_classification(    out             :: MAGEMin_C.gmin_struct{Float64, Int64},
                                    dtb             :: String  )

    # @warn "Breaking changes in v2.1.3 by disambiguation of solvus names: 'spl': 'sp' > 'spl'; 'sp': 'mt' > 'smt'."

    ph      = Array{String}(undef, out.n_SS + out.n_PP) 
    ph_wt   = Array{Float64}(undef, out.n_SS + out.n_PP) 

    # add solution phase and classify some solution phases (spl, fsp, ilm)                             
    for i = 1:out.n_SS                             
        ss              = out.ph[i]
        mineral_name    = ss
        ph_wt[i]        = out.ph_frac_wt[i]

        x               = out.SS_vec[i].compVariables
        if dtb == "ig" || dtb == "igad"
            if ss == "spl"
                if x[3] - 0.5 > 0.0;        mineral_name = "cm";
                elseif x[4] - 0.5 > 0.0;    mineral_name = "sp";
                elseif x[2] - 0.5 > 0.0;    mineral_name = "mgt";
                else                        mineral_name = "spl";    end
            elseif ss == "pig" || ss == "Na-cpx"
                mineral_name = "cpx";
            elseif ss == "gl" || ss == "act" || ss == "amp" || ss == "cumm" || ss == "tr"
                mineral_name = "amp";
            elseif ss == "fsp"
                if x[2] - 0.5 > 0.0;        mineral_name = "afs";
                else                        mineral_name = "pl";    end
            elseif ss == "mu"
                if x[4] - 0.5 > 0.0;        mineral_name = "pat";
                else                        mineral_name = "mu";    end
            elseif ss == "ilm"
                if -x[1] + 0.5 > 0.0;       mineral_name = "hem";
                else                        mineral_name = "FeTiOx";   end 
            end
    
        elseif dtb == "mp" || dtb == "mpe" || dtb == "mb" || dtb == "ume" || db == "mbe"
            if ss == "sp"
                if x[2] - 0.5 > 0.0;        mineral_name = "sp";
                else                        mineral_name = "smt";    end
            elseif ss == "gl" || ss == "act" || ss == "amp" || ss == "cumm" || ss == "tr"
                mineral_name = "amp";
            elseif ss == "omph" || ss == "dio" || ss == "jd"
                mineral_name = "cpx";
            elseif ss == "spl"
                if x[3] - 0.5 > 0.0;        mineral_name = "cm";
                elseif x[2] - 0.5 > 0.0;    mineral_name = "mgt";
                else                        mineral_name = "spl";    end
            elseif ss == "fsp"
                if x[2] - 0.5 > 0.0;        mineral_name = "afs";
                else                        mineral_name = "pl";    end
            elseif ss == "mu"
                if x[4] - 0.5 > 0.0;        mineral_name = "pat";
                else                        mineral_name = "mu";    end
            elseif ss == "ilmm"
                if x[1] - 0.5 > 0.0;        mineral_name = "FeTiOx";
                else                        mineral_name = "hem";   end 
            elseif ss == "ilm"
                if 1.0 - x[1] > 0.5;        mineral_name = "hem";
                else                        mineral_name = "FeTiOx";   end 
            elseif ss == "occm"
                if x[2] > 0.5;              mineral_name = "sid";
                elseif x[3] > 0.5;          mineral_name = "ank";  
                elseif x[1] > 0.25 && x[3] < 0.01;         mineral_name = "mag";  
                else                        mineral_name = "cc";   end 
            end
        end

        # provide the right day
        ph[i]   = mineral_name
    end

    # add pure phases
    for i=1:out.n_PP
        ph[i+out.n_SS]      = out.ph[i+out.n_SS]
        ph_wt[i+out.n_SS]   = out.ph_frac_wt[i+out.n_SS]
    end

    unique_ph   = unique(ph)
    n_ph        = length(unique_ph)
    if n_ph    != length(ph)
        ph_     = Array{String}(undef, length(unique_ph)) 
        ph_wt_  = Array{Float64}(undef, length(unique_ph)) 

        id_ph   = [findall(x->x==ph[i],unique_ph)[1] for i=1:length(ph)]

        for i=1:n_ph
            ph_[i]      = unique_ph[i]
            ph_wt_[i]   = sum(ph_wt[id_ph .== i])
        end   

        return ph_, ph_wt_
    else
        return ph, ph_wt
    end

end

"""
    SaturationConfig

    Configuration struct for saturation models used in `TE_prediction` and `solve_with_saturation`.
    All fields default to `"none"` (disabled).

    Fields
    ------
    Zr   : String — zirconium saturation model. Options: "none", "CB", "WH", "B".
    S    : String — sulfur saturation model. Options: "none", "Liu07", "Oneill21", "<N>ppm".
    P2O5 : String — phosphate saturation model. Options: "none", "Klein26", "HWBea92", "Tollari06".
    CO2  : String — CO₂ saturation model. Options: "none", "SY26".
"""
Base.@kwdef struct SaturationConfig
    Zr   :: String = "none"
    S    :: String = "none"
    P2O5 :: String = "none"
    CO2  :: String = "none"
    # Optional per-element phase overrides.  Keys are element names (e.g. "Zr");
    # values are (phase, adjust) NamedTuples where:
    #   phase  :: String   — phase label used in the KDs database
    #   adjust :: Function — bulk-correction function with signature
    #                        (out, Cliq_el, Sat_el, liq_wt) -> (phase_wt, bulk_delta)
    # Elements absent from this dict use the _SAT_PHASE default and the built-in
    # adjust_bulk_4_* functions.
    overrides :: Dict{String, @NamedTuple{phase::String, adjust::Function}} =
                 Dict{String, @NamedTuple{phase::String, adjust::Function}}()
end

"""
    out_tepm

    Structure holding the output of the trace element (TE) partitioning routine.

    Fields
    ------
    elements : Vector{String}
        Names of the trace elements.
    C0 : Vector{Float64}
        Initial bulk trace element composition [ppm].
    Cliq : Vector{Float64}
        Trace element concentrations in the melt phase [ppm].
    Csol : Vector{Float64}
        Trace element concentrations in the bulk solid [ppm].
    Cmin : Matrix{Float64}
        Trace element concentrations in each individual mineral phase [ppm] (phases × elements).
    ph_TE : Vector{String}
        Names of the phases included in the partitioning calculation.
    ph_wt_norm : Vector{Float64}
        Normalized weight fractions of the solid phases.
    liq_wt_norm : Float64
        Normalized weight fraction of the melt.
    bulk_D : Float64
        Bulk partition coefficient.
    bulk_cor_wt : Vector{Float64}
        Corrected bulk oxide weight fractions (accounting for saturation phases).
    bulk_cor_mol : Vector{Float64}
        Corrected bulk oxide molar fractions.
    Sat_Zr_liq : Float64
        Zirconium saturation concentration in the melt [ppm] (NaN if not computed).
    zrc_wt : Float64
        Weight fraction of zircon precipitated (NaN if not computed).
    Sat_S_liq : Float64
        Sulfur saturation concentration in the melt [ppm] (NaN if not computed).
    sulf_wt : Float64
        Weight fraction of sulfide precipitated (NaN if not computed).
    Sat_P2O5_liq : Float64
        P₂O₅ saturation concentration in the melt [ppm] (NaN if not computed).
    fapt_wt : Float64
        Weight fraction of fluorapatite precipitated (NaN if not computed).
    Sat_CO2_liq : Float64
        CO₂ saturation concentration in the melt [ppm] (NaN if not computed).
    fl_CO2_wt : Float64
        Weight fraction of CO₂ fluid formed (NaN if not computed).
"""
struct out_tepm
    elements        :: Union{Float64, Vector{String}}
    C0              :: Union{Float64, Vector{Float64}}
    Cliq            :: Union{Float64, Vector{Float64}}
    Csol            :: Union{Float64, Vector{Float64}}
    Cmin            :: Union{Float64, Matrix{Float64}}
    ph_TE           :: Union{Nothing, Vector{String}}
    ph_wt_norm      :: Union{Float64, Vector{Float64}}
    liq_wt_norm     :: Union{Float64, Float64}
    bulk_D          :: Union{Float64, Float64}

    bulk_cor_wt     :: Union{Float64, Vector{Float64}}
    bulk_cor_mol    :: Union{Float64, Vector{Float64}}

    Sat_Zr_liq      :: Union{Float64, Float64}
    zrc_wt          :: Union{Float64, Float64}

    Sat_S_liq       :: Union{Float64, Float64}
    sulf_wt         :: Union{Float64, Float64}

    Sat_P2O5_liq    :: Union{Float64, Float64}
    fapt_wt         :: Union{Float64, Float64}

    Sat_CO2_liq     :: Union{Float64, Float64}
    fl_CO2_wt       :: Union{Float64, Float64}
end


"""
    custom_KDs_database

    Structure holding a trace element partitioning coefficient (KD) database.

    Fields
    ------
    infos : String
        Description or citation for the database.
    element_name : Vector{String}
        Names of the trace elements.
    phase_name : Vector{String}
        Names of the mineral phases for which KDs are defined.
    KDs_expr : Matrix{Function}
        Matrix of compiled KD functions (phases × elements). Each function takes a `gmin_struct` as input and returns a Float64.
"""
struct custom_KDs_database
    infos           #:: String
    element_name    #:: Vector{String}
    phase_name      #:: Vector{String}, String
    KDs_expr        #:: Matrix{Expr}, Vector{Expr}
end

"""
    retrieve_eval_rules_TE()

    Return the token substitution rules used when compiling KD expression strings.

    Plain tokens such as `T_C` and `P_kbar` are replaced with their fully-qualified counterparts on a `gmin_struct` (e.g., `out.T_C`), so that compiled functions can be called with a single `out` argument.

    Returns
    -------
    in_eval_TE : Vector{String}
        Tokens to search for in the expression string.
    out_eval_TE : Vector{String}
        Replacement strings referencing the `out` argument.
"""
function retrieve_eval_rules_TE()

    in_eval_TE     = ["T_C","P_kbar","oxides"]
    out_eval_TE    = ["out.T_C","out.P_kbar","out.oxides"]

    return in_eval_TE, out_eval_TE
end


"""
    convert_SS_eval_TE(str)

    Rewrite `[:phase]` subscript tokens in a KD expression string to fully-qualified `gmin_struct` accessor calls.

    The pattern `[:name]` is replaced with `out.SS_vec[out.SS_syms[:name]]`, allowing KD expressions to reference solution phase data by short name (e.g., `[:liq].compVariables[1]` → `out.SS_vec[out.SS_syms[:liq]].compVariables[1]`).

    Parameters
    ----------
    str : String
        KD expression string potentially containing `[:phase]` tokens.

    Returns
    -------
    str : String
        Expression string with all `[:phase]` tokens replaced.
"""
function convert_SS_eval_TE(str)
    pattern = r"\[:([A-Za-z_][A-Za-z0-9_]*)\]"
    matches = collect(eachmatch(pattern, str))

    seen = Set{String}()
    for m in matches
        raw = m.match
        raw in seen && continue
        push!(seen, raw)
        name = m.captures[1]
        str = replace(str, raw => "out.SS_vec[out.SS_syms[:$name]]")
    end

    return str
end

"""
    adjust_chemical_system(KDs_dtb, bulk_TE, elem_TE)

    Reorder and subset a bulk trace element composition vector to match the element order defined in a KD database.

    Elements present in `KDs_dtb` but absent from `elem_TE` are set to zero.

    Parameters
    ----------
    KDs_dtb : custom_KDs_database
        KD database whose `element_name` field defines the target element order.
    bulk_TE : Vector{Float64}
        Input bulk trace element concentrations [ppm].
    elem_TE : Vector{String}
        Element names corresponding to `bulk_TE`.

    Returns
    -------
    C0_TE : Vector{Float64}
        Bulk trace element concentrations reordered to match `KDs_dtb.element_name` [ppm].
"""
function adjust_chemical_system(    KDs_dtb     :: custom_KDs_database,
                                    bulk_TE     :: Vector{Float64},
                                    elem_TE     :: Vector{String}       )

    n_el        = length(KDs_dtb.element_name)
    C0_TE       = zeros(Float64,n_el)

    for i=1:n_el
        id = findfirst(KDs_dtb.element_name[i] .== elem_TE)
        if !isnothing(id)
            C0_TE[i] = bulk_TE[id]
        end
    end

    return C0_TE
end

"""
    create_custom_KDs_database(el_name; info)

    Create an elements-only KDs database with no phases. Phases (and zero KDs) are added
    automatically by `_augment_KDs_for_saturation` when a `SaturationConfig` is provided
    to `TE_prediction` or `solve_with_saturation`.
"""
function create_custom_KDs_database(el_name :: Vector{String};
                                    info    :: String = "Custom KDs database")
    return create_custom_KDs_database(el_name, String[],
                                      Matrix{String}(undef, 0, length(el_name));
                                      info = info)
end

"""
    create_custom_KDs_database(el_name, phase_name; info)

    Create a KDs database with the given elements and phases, all KDs set to zero.
    Useful when phases are known but partition coefficients are saturation-controlled.
"""
function create_custom_KDs_database(el_name    :: Vector{String},
                                    phase_name :: Vector{String};
                                    info       :: String = "Custom KDs database")
    n_ph = length(phase_name)
    n_el = length(el_name)
    return create_custom_KDs_database(el_name, phase_name,
                                      fill("0.0", n_ph, n_el);
                                      info = info)
end

"""
    create_custom_KDs_database(el_name, phase_name, KDs_expr_str; info="Custom KDs database")

    Create a custom trace element partitioning coefficient (KD) database from string expressions.

    Each expression in `KDs_expr_str` may reference `T_C`, `P_kbar`, `oxides`, and solution phase compositional variables using the syntax `[:phase_name]` (e.g., `[:liq].compVariables[1]`). These are resolved against the `gmin_struct` output at evaluation time.

    Parameters
    ----------
    el_name : Vector{String}
        Names of the trace elements.
    phase_name : Vector{String}
        Names of the mineral phases.
    KDs_expr_str : Union{Matrix{String}, Vector{String}}
        Matrix (phases × elements) or vector of KD expressions as strings.
    info : String, optional
        Description or citation for the database (default: "Custom KDs database").

    Returns
    -------
    db : custom_KDs_database
        Compiled KD database ready for use in `TE_prediction`.
"""
function create_custom_KDs_database(el_name         :: Vector{String},
                                    phase_name      :: Vector{String},
                                    KDs_expr_str    :: Union{Matrix{String},Vector{String}} ;
                                    info            :: String = "Custom KDs database")

    n_el    = length(el_name)
    n_ph    = length(phase_name)

    KDs_expr = Matrix{Function}(undef, n_ph, n_el)
    in_eval_TE, out_eval_TE = retrieve_eval_rules_TE()

    for i = 1:n_ph
        for j = 1:n_el

            expr_str      = KDs_expr_str[i,j]
            for i=1:length(in_eval_TE)
                expr_str = replace(expr_str, in_eval_TE[i] => out_eval_TE[i])
            end

            expr_str      = convert_SS_eval_TE(expr_str)
            KDs_expr[i,j] = eval(Meta.parse("out -> ($(expr_str) +0)"))
        end
    end

    return custom_KDs_database(info, el_name, phase_name, KDs_expr)
end

"""
    partition_TE(KDs_database, out, C0, ph, ph_wt, liq_wt; norm_TE=true)

    Apply the batch melting equation to partition trace elements between melt and solid phases for the subset of phases present in both the MAGEMin output and the KD database.

    Parameters
    ----------
    KDs_database : custom_KDs_database
        Compiled KD database.
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output used to evaluate P-T-composition-dependent KDs.
    C0 : Vector{Float64}
        Initial bulk trace element composition [ppm].
    ph : Vector{String}
        Phase names from `mineral_classification`.
    ph_wt : Vector{Float64}
        Weight fractions of each phase.
    liq_wt : Float64
        Melt weight fraction.
    norm_TE : Bool, optional
        Normalize phase fractions to the phases present in the KD database (default: true).

    Returns
    -------
    Cliq : Vector{Float64}
        Trace element concentrations in the melt [ppm].
    Cmin : Matrix{Float64}
        Trace element concentrations in each mineral phase [ppm] (phases × elements).
    Csol : Vector{Float64}
        Trace element concentrations in the bulk solid [ppm].
    ph_TE : Vector{String}
        Names of phases included in the calculation.
    ph_wt_norm : Vector{Float64}
        Normalized solid phase weight fractions.
    liq_wt_norm : Float64
        Normalized melt weight fraction.
    bulk_D : Float64
        Bulk partition coefficient.
"""
function partition_TE(  KDs_database:: custom_KDs_database,
                        out         :: MAGEMin_C.gmin_struct{Float64, Int64},   
                        C0          :: Vector{Float64},
                        ph          :: Vector{String},
                        ph_wt       :: Vector{Float64}, 
                        liq_wt      :: Float64;
                        norm_TE     :: Bool = true)


    phase_name  = KDs_database.phase_name
    TE_ph       = intersect(ph, phase_name);

    # get indexes of the phase with respect to MAGEMin output and TE_database
    MM_ph_idx   = [findfirst(isequal(x), ph) for x in TE_ph]
    TE_ph_idx   = [findfirst(isequal(x), phase_name) for x in TE_ph]

    # normalize phase fractions
    sum_ph_frac = sum(ph_wt[MM_ph_idx]);
    if norm_TE == true
        liq_wt_norm = liq_wt/(sum_ph_frac+liq_wt);
        ph_wt_norm  = ph_wt[MM_ph_idx]./sum_ph_frac;
    else
        liq_wt_norm = liq_wt
        ph_wt_norm  = ph_wt[MM_ph_idx]
    end
    ph_TE       = ph[MM_ph_idx];

    n_ph        = length(ph_wt_norm)
    n_el        = length(C0)
    KDs         = zeros(Float64, n_ph, n_el) # partitioning coefficients

    for i=1:n_ph
        for j=1:n_el
            expr        = KDs_database.KDs_expr[TE_ph_idx[i],j]
            KDs[i,j]    = Base.invokelatest(expr, out)
        end
    end

    D           = KDs'*ph_wt_norm;
    bulk_D      = sum(D);
    Cliq        = C0 ./ (D .+ liq_wt_norm.*(1.0 .- D));
    Csol        = (C0 .- Cliq .*  liq_wt_norm) ./ (1.0 .- liq_wt_norm)
    Cmin        = similar(KDs); 

    for i = 1:length(ph_wt_norm)
        Cmin[i,:] = KDs[i,:] .* Cliq;
    end

    return Cliq, Cmin, Csol, ph_TE, ph_wt_norm, liq_wt_norm, bulk_D
end


"""
    health_check_TE(C0, KDs_database)

    Validate that the initial composition vector and KD database are consistent before running TE partitioning.

    Parameters
    ----------
    C0 : Vector{Float64}
        Initial bulk trace element composition [ppm].
    KDs_database : custom_KDs_database
        KD database to validate against `C0`.

    Returns
    -------
    status : Int64
        1 if all checks pass, 0 if any check fails (errors are also thrown).
"""
function health_check_TE(C0, KDs_database)

    status = 1
    if length(C0) != length(KDs_database.element_name)
        error("C0 and KDs_database.element_name must have the same length")
        status = 0
    end
    if length(KDs_database.phase_name) == 0
        error("KDs_database.phase_name must not be empty")
        status = 0
    end
    if length(KDs_database.KDs_expr) == 0
        error("KDs_database.KDs_expr must not be empty")
        status = 0
    end
    if size(KDs_database.KDs_expr,1) != length(KDs_database.phase_name)
        error("KDs_database.KDs_expr must have the same number of rows as KDs_database.phase_name")
        status = 0
    end

    return status
end


"""
    compute_TE_partitioning(KDs_database, out, C0, ph, ph_wt, liq_wt, sol_wt; norm_TE=true)

    Partition trace elements between melt and solid phases using the supplied KD database.

    Handles three end-member cases: fully molten (`liq_wt == 1.0`), fully solid (`liq_wt == 0.0`), and mixed (`0 < liq_wt < 1`). In the mixed case the batch melting equation is applied.

    Parameters
    ----------
    KDs_database : custom_KDs_database
        Compiled trace element partitioning coefficient database.
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output used to evaluate P-T-composition-dependent KDs.
    C0 : Vector{Float64}
        Initial bulk trace element composition [ppm].
    ph : Vector{String}
        Phase names from `mineral_classification`.
    ph_wt : Vector{Float64}
        Weight fractions of each phase.
    liq_wt : Float64
        Melt weight fraction.
    sol_wt : Float64
        Solid weight fraction.
    norm_TE : Bool, optional
        Normalize phase fractions before computing KDs (default: true).

    Returns
    -------
    Cliq : Vector{Float64}
        Trace element concentrations in the melt [ppm].
    Csol : Vector{Float64}
        Trace element concentrations in the bulk solid [ppm].
    Cmin : Matrix{Float64}
        Trace element concentrations in each mineral phase [ppm].
    ph_TE : Vector{String}
        Phase names included in the calculation.
    ph_wt_norm : Vector{Float64}
        Normalized solid phase weight fractions.
    liq_wt_norm : Float64
        Normalized melt weight fraction.
    bulk_D : Float64
        Bulk partition coefficient.
"""
function compute_TE_partitioning(   KDs_database:: custom_KDs_database,
                                    out         :: MAGEMin_C.gmin_struct{Float64, Int64},
                                    C0          :: Vector{Float64},
                                    ph          :: Vector{String},
                                    ph_wt       :: Vector{Float64}, 
                                    liq_wt      :: Float64,
                                    sol_wt      :: Float64;
                                    norm_TE     :: Bool = true  )


    if liq_wt > 0.0 && liq_wt < 1.0 && sol_wt > 0.0
        
        Cliq, Cmin, Csol, ph_TE, ph_wt_norm, liq_wt_norm, bulk_D = partition_TE(    KDs_database, out, C0, 
                                                                                    ph, ph_wt, liq_wt; norm_TE=norm_TE)
    elseif liq_wt == 0.0
        Csol        = C0
        Cliq, Cmin, ph_TE, ph_wt_norm, liq_wt_norm, bulk_D = C0.*0.0, NaN, nothing, NaN, NaN, NaN

    elseif liq_wt == 1.0 || (sol_wt == 0.0 && liq_wt > 0.0) #latter means there is fluid + melt
        Cliq        = C0
        Csol, Cmin, ph_TE, ph_wt_norm, bulk_D  = C0.*0.0, NaN, nothing, NaN, NaN
        liq_wt_norm = 1.0
    else
        println("unrecognized case!")
        Cliq, Csol, Cmin, ph_TE, ph_wt_norm, liq_wt_norm, bulk_D = NaN, NaN, NaN, nothing, NaN, NaN, NaN
    end

    return Cliq, Csol, Cmin, ph_TE, ph_wt_norm, liq_wt_norm, bulk_D
end


# Default element → saturation phase mapping
const _SAT_PHASE = Dict("Zr" => "zrc", "S" => "sulf", "P2O5" => "fapt", "CO2" => "flC")

# Return the saturation phase for `el`: override → _SAT_PHASE default → nothing
_sat_phase(sat::SaturationConfig, el::String) =
    haskey(sat.overrides, el) ? sat.overrides[el].phase : get(_SAT_PHASE, el, nothing)

"""
    _augment_KDs_for_saturation(KDs_database, sat)

Internal helper. When `sat` is a `SaturationConfig`, append a column of zero-KD
functions for each active saturation phase that is not already present in
`KDs_database.phase_name`.  Returns the (possibly augmented) database; the
original is never mutated.
"""
function _augment_KDs_for_saturation(KDs_database :: custom_KDs_database,
                                     sat          :: SaturationConfig)
    extra = String[]
    for field in ("Zr", "S", "P2O5", "CO2")
        model = getfield(sat, Symbol(field))
        phase = _sat_phase(sat, field)
        if model != "none" && !isnothing(phase) && phase ∉ KDs_database.phase_name
            push!(extra, phase)
        end
    end

    isempty(extra) && return KDs_database

    n_elem   = length(KDs_database.element_name)
    zero_fn  = (_) -> 0.0
    new_KDs  = vcat(KDs_database.KDs_expr, fill(zero_fn, length(extra), n_elem))

    return custom_KDs_database(
        KDs_database.infos,
        KDs_database.element_name,
        vcat(KDs_database.phase_name, extra),
        new_KDs,
    )
end


"""
    compute_Zr_sat_n_part(out, KDs_database, Cliq, bulk_cor_wt, C0, liq_wt; ZrSat_model="CB")

    Check zircon saturation and adjust the corrected bulk composition if the melt exceeds the Zr saturation limit.

    If Zr in the melt exceeds the saturation concentration, the excess is removed from the melt and the corresponding SiO₂ and O are returned to the bulk. If there is no melt (`liq_wt == 0`), all Zr is assumed to have precipitated as zircon.

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output.
    KDs_database : custom_KDs_database
        Trace element partitioning coefficient database (must include "Zr").
    Cliq : Vector{Float64}
        Current trace element concentrations in the melt [ppm].
    bulk_cor_wt : Vector{Float64}
        Corrected bulk oxide weight fractions (modified in place).
    C0 : Vector{Float64}
        Initial bulk trace element composition [ppm].
    liq_wt : Float64
        Melt weight fraction.
    ZrSat_model : String, optional
        Zirconium saturation model (default: "CB"). Passed to `zirconium_saturation`.

    Returns
    -------
    Sat_Zr_liq : Float64
        Zr saturation concentration in the melt [ppm].
    zrc_wt : Float64
        Weight fraction of precipitated zircon.
    bulk_cor_wt : Vector{Float64}
        Updated corrected bulk oxide weight fractions.
"""
function compute_Zr_sat_n_part(     out         :: MAGEMin_C.gmin_struct{Float64, Int64},
                                    KDs_database:: custom_KDs_database,
                                    Cliq, bulk_cor_wt, C0,
                                    liq_wt      :: Float64;
                                    ZrSat_model :: String = "CB")
    Sat_Zr_liq  = NaN
    id_Zr       = findfirst(KDs_database.element_name .== "Zr")
    if liq_wt > 0.0
        Cliq_Zr     = Cliq[id_Zr]
        Sat_Zr_liq  = zirconium_saturation( out; 
                                            model = ZrSat_model)

        if Cliq_Zr > Sat_Zr_liq
            zrc_wt, SiO2_zrc_wt, O_zrc_wt       = adjust_bulk_4_zircon(Cliq_Zr, Sat_Zr_liq, liq_wt)

            SiO2_bulk_wt             = out.SS_vec[out.SS_syms[:liq]].Comp_wt[findfirst(isequal("SiO2"), out.oxides)]

            if SiO2_zrc_wt*liq_wt  > SiO2_bulk_wt
                @warn "Not enough SiO2 in the bulk composition to saturate in zircon. Increasing the Sat_Zr_liq to the available SiO2 content."
                factor                          = SiO2_bulk_wt / SiO2_zrc_wt*liq_wt 
                Sat_Zr_liq                      = Sat_Zr_liq + factor * (Cliq_Zr - Sat_Zr_liq)
                zrc_wt, SiO2_zrc_wt, O_zrc_wt   = adjust_bulk_4_zircon(Cliq_Zr, Sat_Zr_liq, liq_wt)
            end
        else
            zrc_wt, SiO2_zrc_wt, O_zrc_wt = 0.0, 0.0, 0.0
        end
    else
        C0_Zr       =  C0[id_Zr]
        zrc_wt, SiO2_zrc_wt, O_zrc_wt       = adjust_bulk_4_zircon(C0_Zr, 0.0, 1.0)
    end

    bulk_cor_wt[findfirst(out.oxides .== "SiO2")]   += SiO2_zrc_wt
    bulk_cor_wt[findfirst(out.oxides .== "O")]      += 2.0*O_zrc_wt

    return Sat_Zr_liq, zrc_wt, bulk_cor_wt
end

"""
    compute_S_sat_n_part(out, KDs_database, Cliq, bulk_cor_wt, C0, liq_wt; SSat_model="1000ppm")

    Check sulfur saturation and adjust the corrected bulk composition if the melt exceeds the S saturation limit.

    If S in the melt exceeds the saturation concentration, the excess is removed and the corresponding FeO and O are returned to the bulk. If there is no melt, all S is assumed to have precipitated as sulfide.

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output.
    KDs_database : custom_KDs_database
        Trace element partitioning coefficient database (must include "S").
    Cliq : Vector{Float64}
        Current trace element concentrations in the melt [ppm].
    bulk_cor_wt : Vector{Float64}
        Corrected bulk oxide weight fractions (modified in place).
    C0 : Vector{Float64}
        Initial bulk trace element composition [ppm].
    liq_wt : Float64
        Melt weight fraction.
    SSat_model : String, optional
        Sulfur saturation model (default: "1000ppm"). Passed to `sulfur_saturation`.

    Returns
    -------
    Sat_S_liq : Float64
        S saturation concentration in the melt [ppm].
    sulf_wt : Float64
        Weight fraction of precipitated sulfide.
    bulk_cor_wt : Vector{Float64}
        Updated corrected bulk oxide weight fractions.
"""
function compute_S_sat_n_part(      out         :: MAGEMin_C.gmin_struct{Float64, Int64},
                                    KDs_database:: custom_KDs_database,
                                    Cliq, bulk_cor_wt, C0,
                                    liq_wt      :: Float64;
                                    SSat_model  :: String = "1000ppm")
    Sat_S_liq   = NaN
    id_S        = findfirst(KDs_database.element_name .== "S")
    if liq_wt > 0.0
        Cliq_S      = Cliq[id_S]
        Sat_S_liq   = sulfur_saturation(    out; 
                                            model = SSat_model)

        if Cliq_S > Sat_S_liq
            sulf_wt, FeO_sulf_wt, O_sulf_wt     = adjust_bulk_4_sulfide(Cliq_S, Sat_S_liq, liq_wt)

            FeO_bulk_wt             = out.SS_vec[out.SS_syms[:liq]].Comp_wt[findfirst(isequal("FeO"), out.oxides)]
            # O_bulk_wt               = out.SS_vec[out.SS_syms[:liq]].Comp_wt[findfirst(isequal("O"),    out.oxides)]
            
            if FeO_sulf_wt*liq_wt > FeO_bulk_wt
                @warn "Not enough FeO in the bulk composition to saturate in sulfide. Increasing the Sat_S_liq to the available FeO content."
                factor              = FeO_bulk_wt / FeO_sulf_wt*liq_wt 
                Sat_S_liq           = Sat_S_liq + factor * (Cliq_S - Sat_S_liq)
                sulf_wt, FeO_sulf_wt, O_sulf_wt = adjust_bulk_4_sulfide(Cliq_S, Sat_S_liq, liq_wt)
            end
        else
            sulf_wt, FeO_sulf_wt, O_sulf_wt = 0.0, 0.0, 0.0
        end
    else
        C0_S       =  C0[id_S]
        sulf_wt, FeO_sulf_wt, O_sulf_wt      = adjust_bulk_4_sulfide(C0_S, 0.0, 1.0)
    end

    bulk_cor_wt[findfirst(out.oxides .== "FeO")] += FeO_sulf_wt
    bulk_cor_wt[findfirst(out.oxides .== "O")]   += O_sulf_wt

    return Sat_S_liq, sulf_wt, bulk_cor_wt
end


"""
    compute_P2O5_sat_n_part(out, KDs_database, Cliq, bulk_cor_wt, C0, liq_wt; P2O5Sat_model="Klein26")

    Check phosphate saturation and adjust the corrected bulk composition if the melt exceeds the P₂O₅ saturation limit.

    If P₂O₅ in the melt exceeds the saturation concentration, the excess is removed and the corresponding CaO is returned to the bulk. If there is no melt, all P₂O₅ is assumed to have precipitated as fluorapatite.

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output.
    KDs_database : custom_KDs_database
        Trace element partitioning coefficient database (must include "P2O5").
    Cliq : Vector{Float64}
        Current trace element concentrations in the melt [ppm].
    bulk_cor_wt : Vector{Float64}
        Corrected bulk oxide weight fractions (modified in place).
    C0 : Vector{Float64}
        Initial bulk trace element composition [ppm].
    liq_wt : Float64
        Melt weight fraction.
    P2O5Sat_model : String, optional
        Phosphate saturation model (default: "Klein26"). Passed to `phosphate_saturation`.

    Returns
    -------
    Sat_P2O5_liq : Float64
        P₂O₅ saturation concentration in the melt [ppm].
    fapt_wt : Float64
        Weight fraction of precipitated fluorapatite.
    bulk_cor_wt : Vector{Float64}
        Updated corrected bulk oxide weight fractions.
"""
function compute_P2O5_sat_n_part(   out         :: MAGEMin_C.gmin_struct{Float64, Int64},
                                    KDs_database:: custom_KDs_database,
                                    Cliq, bulk_cor_wt, C0,
                                    liq_wt      :: Float64;
                                    P2O5Sat_model   :: String = "Klein26")
    Sat_P2O5_liq   = NaN
    id_P2O5        = findfirst(KDs_database.element_name .== "P2O5")
    if liq_wt > 0.0                                    
        Cliq_P2O5     = Cliq[id_P2O5]
        Sat_P2O5_liq  = phosphate_saturation(   out; 
                                                model = P2O5Sat_model)

        if Cliq_P2O5 > Sat_P2O5_liq
            fapt, CaO_fpat_wt     = adjust_bulk_4_fapatite(Cliq_P2O5, Sat_P2O5_liq, liq_wt)

            CaO_bulk_wt           = out.SS_vec[out.SS_syms[:liq]].Comp_wt[findfirst(isequal("CaO"), out.oxides)]
            if CaO_fpat_wt*liq_wt  > CaO_bulk_wt
                @warn "Not enough CaO in the bulk composition to saturate in fapatite. Increasing the Sat_P2O5_liq to the available CaO content."
                factor            = CaO_bulk_wt / CaO_fpat_wt*liq_wt 
                Sat_P2O5_liq      = Sat_P2O5_liq + factor * (Cliq_P2O5 - Sat_P2O5_liq)
                fapt, CaO_fpat_wt = adjust_bulk_4_fapatite(Cliq_P2O5, Sat_P2O5_liq, liq_wt)
            end
        else
            fapt, CaO_fpat_wt = 0.0, 0.0
        end
    else
        C0_P2O5            =  C0[findfirst(KDs_database.element_name .== "P2O5")]
        fapt, CaO_fpat_wt  = adjust_bulk_4_fapatite(C0_P2O5, 0.0, 1.0)
    end

    bulk_cor_wt[findfirst(out.oxides .== "CaO")]  +=  CaO_fpat_wt

    return Sat_P2O5_liq, fapt, bulk_cor_wt
end


"""
    compute_CO2_sat_n_part(out, KDs_database, Cliq, bulk_cor_wt, C0, liq_wt; CO2Sat_model="SY26")

    Check CO₂ saturation and adjust the corrected bulk composition if the melt exceeds the CO₂ saturation limit.

    If CO₂ in the melt exceeds the saturation concentration, the excess degasses into a CO₂-bearing fluid and the
    corresponding weight is added back to the bulk CO₂ oxide entry in `bulk_cor_wt`.  Unlike mineral saturation
    phases, no stoichiometrically distinct oxide is returned — the excess CO₂ re-enters the CO₂ oxide budget directly.

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output.
    KDs_database : custom_KDs_database
        Trace element partitioning coefficient database (must include "CO2").
    Cliq : Vector{Float64}
        Current trace element concentrations in the melt [ppm].
    bulk_cor_wt : Vector{Float64}
        Corrected bulk oxide weight fractions (modified in place).
    C0 : Vector{Float64}
        Initial bulk trace element composition [ppm].
    liq_wt : Float64
        Melt weight fraction.
    CO2Sat_model : String, optional
        CO₂ saturation model (default: "SY26"). Passed to `co2_saturation`.

    Returns
    -------
    Sat_CO2_liq : Float64
        CO₂ saturation concentration in the melt [ppm].
    fl_CO2_wt : Float64
        Weight fraction of CO₂ fluid formed.
    bulk_cor_wt : Vector{Float64}
        Updated corrected bulk oxide weight fractions.
"""
function compute_CO2_sat_n_part(    out         :: MAGEMin_C.gmin_struct{Float64, Int64},
                                    KDs_database:: custom_KDs_database,
                                    Cliq, bulk_cor_wt, C0,
                                    liq_wt      :: Float64;
                                    CO2Sat_model :: String = "SY26")

    Sat_CO2_liq = NaN
    id_CO2      = findfirst(KDs_database.element_name .== "CO2")

    if liq_wt > 0.0
        Cliq_CO2  = Cliq[id_CO2]
        # Cliq_CO2    = C0[findfirst(KDs_database.element_name .== "CO2")] * (1.0/liq_wt)
        Sat_CO2_liq = co2_saturation(out; model = CO2Sat_model)
        if Cliq_CO2 > Sat_CO2_liq
            fl_CO2_wt, CO2_bulk_wt = adjust_bulk_4_fluid(Cliq_CO2, Sat_CO2_liq, liq_wt)
        else
            fl_CO2_wt, CO2_bulk_wt = 0.0, 0.0
        end
    else
        C0_CO2            =  C0[findfirst(KDs_database.element_name .== "CO2")]
        fl_CO2_wt, CO2_bulk_wt  = adjust_bulk_4_fluid(C0_CO2, 0.0, 1.0)
    end

    return Sat_CO2_liq, fl_CO2_wt, bulk_cor_wt
end


"""
    TE_prediction(out, C0, KDs_database, dtb; ZrSat_model="none", SSat_model="none", P2O5Sat_model="none", norm_TE=false)

    Perform trace element partitioning, optionally with zircon, sulfide, and/or apatite saturation corrections.

    Phases are classified via `mineral_classification`, elements are partitioned using the batch melting equation, and saturation corrections are applied when the corresponding model is not `"none"`. The corrected bulk composition (accounting for precipitated saturation phases) is also returned.

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output.
    C0 : Vector{Float64}
        Initial bulk trace element composition [ppm], in the order of `KDs_database.element_name`.
    KDs_database : custom_KDs_database
        Compiled trace element partitioning coefficient database.
    dtb : String
        Database identifier used for mineral classification (e.g., "ig", "mp").
    ZrSat_model : String, optional
        Zirconium saturation model — "none" disables Zr correction (default: "none"). Valid options: "CB", "W85", "BD92", "RZ93", "CZLD08".
    SSat_model : String, optional
        Sulfur saturation model — "none" disables S correction (default: "none"). Valid option: "1000ppm".
    P2O5Sat_model : String, optional
        Phosphate saturation model — "none" disables P₂O₅ correction (default: "none"). Valid options: "Klein26", "HWBea92", "Tollari06".
    norm_TE : Bool, optional
        Normalize phase fractions before computing KDs (default: false).

    Returns
    -------
    out_TE : out_tepm
        Structure containing melt/solid/mineral TE concentrations, saturation values, and corrected bulk compositions.
"""
function TE_prediction( out, C0, KDs_database, dtb;
                        ZrSat_model     :: String                       = "none",
                        SSat_model      :: String                       = "none",
                        P2O5Sat_model   :: String                       = "none",
                        CO2Sat_model    :: String                       = "none",
                        sat             :: Union{Nothing,SaturationConfig} = nothing,
                        norm_TE         :: Bool                         = false )

    # New-style SaturationConfig overrides individual keyword args when provided.
    # Old keyword args are kept for backward compatibility.
    if !isnothing(sat)
        ZrSat_model   = sat.Zr
        SSat_model    = sat.S
        P2O5Sat_model = sat.P2O5
        CO2Sat_model  = sat.CO2
        KDs_database  = _augment_KDs_for_saturation(KDs_database, sat)
    end

    # Initialize output variables
    Cliq, Csol, Cmin, ph_TE, ph_wt_norm, liq_wt_norm = NaN, NaN, NaN, nothing, NaN, NaN
    Sat_Zr_liq, zrc_wt, bulk_cor_wt       = NaN, NaN, NaN
    Sat_S_liq, sulf_wt                    = NaN, NaN, NaN
    Sat_P2O5_liq, fapt_wt,                = NaN, NaN, NaN
    Sat_CO2_liq, fl_CO2_wt                = NaN, NaN

    # input data
    liq_wt      = out.frac_M_wt
    sol_wt      = out.frac_S_wt

    status      = health_check_TE(C0, KDs_database)
    ph, ph_wt   =  mineral_classification(out, dtb)

    # first let's partition elements
    Cliq, Csol, Cmin, ph_TE, ph_wt_norm, liq_wt_norm, bulk_D = compute_TE_partitioning(     KDs_database,
                                                                                            out,
                                                                                            C0,
                                                                                            ph,
                                                                                            ph_wt, 
                                                                                            liq_wt,
                                                                                            sol_wt;
                                                                                            norm_TE = norm_TE)

    bulk_cor_wt = copy(out.bulk_wt); bulk_cor_wt .= 0.0;                                                                              

    if !isnothing(findfirst(KDs_database.element_name .== "Zr")) && ZrSat_model != "none"
        Sat_Zr_liq, zrc_wt, bulk_cor_wt = compute_Zr_sat_n_part(            out,
                                                                            KDs_database,
                                                                            Cliq, bulk_cor_wt, C0,
                                                                            liq_wt;
                                                                            ZrSat_model = ZrSat_model)
        push!(ph,"zrc")
        push!(ph_wt, zrc_wt)
    end

    if !isnothing(findfirst(KDs_database.element_name .== "S")) && SSat_model != "none"
        Sat_S_liq, sulf_wt, bulk_cor_wt = compute_S_sat_n_part(             out,
                                                                            KDs_database,
                                                                            Cliq, bulk_cor_wt, C0,
                                                                            liq_wt;
                                                                            SSat_model  = SSat_model)
        push!(ph,"sulf")
        push!(ph_wt, sulf_wt)
    end

    if !isnothing(findfirst(KDs_database.element_name .== "P2O5")) && P2O5Sat_model != "none"
        if P2O5Sat_model == "Tollari06" && out.oxides[findfirst(out.oxides .== "H2O")] != 0.0
            @warn "P2O5 saturation model 'Tollari06' is calibrated for dry systems. Use HWBea92 model for hydrous systems."
        end

        Sat_P2O5_liq, fapt_wt, bulk_cor_wt = compute_P2O5_sat_n_part(       out,
                                                                            KDs_database,
                                                                            Cliq, bulk_cor_wt, C0,
                                                                            liq_wt;
                                                                            P2O5Sat_model   = P2O5Sat_model)
        push!(ph,"fapt")
        push!(ph_wt, fapt_wt)
    end

    if !isnothing(findfirst(KDs_database.element_name .== "CO2")) && CO2Sat_model != "none"
        Sat_CO2_liq, fl_CO2_wt, bulk_cor_wt = compute_CO2_sat_n_part(       out,
                                                                            KDs_database,
                                                                            Cliq, bulk_cor_wt, C0,
                                                                            liq_wt;
                                                                            CO2Sat_model = CO2Sat_model)
        push!(ph,"flC")
        push!(ph_wt, fl_CO2_wt)
    end
    sum_wt       = sum(ph_wt)
    bulk_cor_wt .= bulk_cor_wt  ./ sum_wt
    ph_wt        = ph_wt        ./ sum_wt


    Cliq, Csol, Cmin, ph_TE, ph_wt_norm, liq_wt_norm, bulk_D = compute_TE_partitioning(     KDs_database,
                                                                                            out,
                                                                                            C0,
                                                                                            ph,
                                                                                            ph_wt, 
                                                                                            liq_wt,
                                                                                            sol_wt;
                                                                                            norm_TE = norm_TE)

    if liq_wt > 0.0
        if !isnothing(findfirst(KDs_database.element_name .== "Zr"))    && ZrSat_model != "none"
            id_Zr           = findfirst(KDs_database.element_name .== "Zr")
            Sat_Zr_liq      = zirconium_saturation(     out; 
                                                        model = ZrSat_model)   

            if Cliq[id_Zr] > Sat_Zr_liq
                Cliq[id_Zr] = Sat_Zr_liq
                Csol[id_Zr] = (C0[id_Zr] - Sat_Zr_liq*liq_wt_norm) / (1.0 - liq_wt_norm)
            else
                Csol[id_Zr] = 0.0
            end
        end
        if !isnothing(findfirst(KDs_database.element_name .== "S"))     && SSat_model != "none"
            id_S            = findfirst(KDs_database.element_name .== "S")
            Sat_S_liq       = sulfur_saturation(    out; 
                                                    model = SSat_model)  
            if Cliq[id_S] > Sat_S_liq
                Cliq[id_S] = Sat_S_liq
                Csol[id_S] = (C0[id_S] - Sat_S_liq*liq_wt_norm) / (1.0 - liq_wt_norm)
            else
                Csol[id_S] = 0.0
            end
        end
        if !isnothing(findfirst(KDs_database.element_name .== "P2O5"))  && P2O5Sat_model != "none"
            id_P2O5         = findfirst(KDs_database.element_name .== "P2O5")
            Sat_P2O5_liq    = phosphate_saturation(     out; 
                                                        model = P2O5Sat_model)  
            if Cliq[id_P2O5] > Sat_P2O5_liq
                Cliq[id_P2O5] = Sat_P2O5_liq
                Csol[id_P2O5] = (C0[id_P2O5] - Sat_P2O5_liq*liq_wt_norm) / (1.0 - liq_wt_norm)
            else
                Csol[id_P2O5] = 0.0
            end
        end
        if !isnothing(findfirst(KDs_database.element_name .== "CO2"))   && CO2Sat_model != "none"
            id_CO2          = findfirst(KDs_database.element_name .== "CO2")
            Sat_CO2_liq     = co2_saturation(   out;
                                                model = CO2Sat_model)
            if !isnan(Sat_CO2_liq) && Cliq[id_CO2] > Sat_CO2_liq
                Cliq[id_CO2] = Sat_CO2_liq
                Csol[id_CO2] = (C0[id_CO2] - Sat_CO2_liq*liq_wt_norm) / (1.0 - liq_wt_norm)
            else
                Csol[id_CO2] = 0.0
            end
        end
    end

    # compute corrected bulk molar composition
    n_sys           = [out.bulk_wt[i]       / get_molar_mass(out.oxides[i]) for i in eachindex(out.oxides) ]
    n_cor           = [bulk_cor_wt[i]       / get_molar_mass(out.oxides[i]) for i in eachindex(out.oxides) ]
    n_cor_fac       = n_cor./(n_sys .+ n_cor)
    bulk_cor_mol    = out.bulk .* n_cor_fac

    out_TE = out_tepm(  KDs_database.element_name, 
                        C0, Cliq, Csol, Cmin,
                        ph_TE, ph_wt_norm, liq_wt_norm, bulk_D, bulk_cor_wt, bulk_cor_mol,

                        Sat_Zr_liq,     zrc_wt/ sum_wt,
                        Sat_S_liq,      sulf_wt/ sum_wt,
                        Sat_P2O5_liq,   fapt_wt/ sum_wt,
                        Sat_CO2_liq,    fl_CO2_wt/ sum_wt)

    return out_TE
end


"""
    solve_with_saturation(P, T, data, X_mol, Xoxides, C0, KDs_dtb, dtb; sat, sys_in, tol, max_iter)

    Run the MAGEMin minimization + TE partitioning loop with saturation corrections
    until the corrected bulk composition converges.

    Parameters
    ----------
    P, T : Float64
        Pressure [kbar] and temperature [°C].
    data : MAGEMin data handle
        Initialised with `Initialize_MAGEMin`.
    X_mol : Vector{Float64}
        Normalised bulk molar composition (will not be mutated).
    Xoxides : Vector{String}
        Oxide names matching `X_mol`.
    C0 : Vector{Float64}
        Initial bulk trace-element composition [ppm].
    KDs_dtb : custom_KDs_database
        Trace-element partitioning database (real mineral phases only;
        saturation phases are appended automatically when `sat` is provided).
    dtb : String
        MAGEMin database identifier (e.g. "mp", "ig").
    sat : SaturationConfig, optional
        Saturation model selection. Defaults to all "none".
    sys_in : String, optional
        Composition input units for MAGEMin: "mol" (default) or "wt".
    tol : Float64, optional
        Convergence tolerance on the L2-norm of `bulk_cor_mol` (default 1e-6).
    max_iter : Int, optional
        Maximum number of iterations (default 32).

    Returns
    -------
    out : MAGEMin output at convergence.
    out_TE : out_tepm at convergence.
    converged : Bool — true if `tol` was reached.
    n_iter : Int — number of iterations performed.
"""
function solve_with_saturation( P           :: Float64,
                                T           :: Float64,
                                data,
                                X_mol       :: Vector{Float64},
                                Xoxides     :: Vector{String},
                                C0          :: Vector{Float64},
                                KDs_dtb     :: custom_KDs_database,
                                dtb         :: String;
                                sat         :: SaturationConfig = SaturationConfig(),
                                sys_in      :: String  = "mol",
                                tol         :: Float64 = 1e-6,
                                max_iter    :: Int     = 32      )

    X    = copy(X_mol)
    n0   = 0.0
    local out, out_TE
    for ite in 1:max_iter
        out    = single_point_minimization(P, T, data; X=X, Xoxides=Xoxides, sys_in=sys_in)
        out_TE = TE_prediction(out, C0, KDs_dtb, dtb; sat=sat)

        X   = X_mol .- out_TE.bulk_cor_mol
        res = abs(n0 - vec_norm(out_TE.bulk_cor_mol))
        n0  = vec_norm(out_TE.bulk_cor_mol)

        if res < tol
            return out, out_TE, true, ite
        end
    end

    @warn "solve_with_saturation: did not converge in $max_iter iterations"
    return out, out_TE, false, max_iter
end


# ---------------------------------------------------------------------------
# CO database — lattice strain models (Cornet 2017, TEPM v02.02)
# Requires lattice_strain.jl (D_cpx, D_gt, D_opx, D_pl, D_ol, D_amph, TE_names)
# to be included before this file is used with tedb="CO".
#
# make_comp_from_SS: extract wt% oxides from one SS_data entry into a NamedTuple
# that the D_* functions can consume directly.
# ---------------------------------------------------------------------------
function make_comp_from_SS(out::MAGEMin_C.gmin_struct{Float64,Int64}, ss_idx::Int)
    wt  = out.SS_vec[ss_idx].Comp_wt
    oxs = out.oxides
    g(name) = (i = findfirst(==(name), oxs); i === nothing ? 0.0 : wt[i])
    (SiO2=g("SiO2"), TiO2=g("TiO2"), Al2O3=g("Al2O3"), FeO=g("FeO"),
     MnO=g("MnO"),   MgO=g("MgO"),   CaO=g("CaO"),     Na2O=g("Na2O"),
     K2O=g("K2O"),   H2O=g("H2O"))
end

"""
    get_CO_KDs_database()

Build a `custom_KDs_database` using the lattice strain models of Cornet (2017)
with Thermocalc a-x solution model site fractions (cpx_G23, g_G23, opx_G23,
ol_H18, fsp_H21, hb_G16).

Requires `lattice_strain.jl` to be included before calling this function.

Element order (28 elements, fixed): Cs Rb K | Ba Sr | La…Lu Sc | Ti Hf Zr | U Th | Ta Nb

Phase names after `mineral_classification` that are handled:
  "cpx", "gt", "opx", "pl", "ol", "hb", "amp"

Returns a `custom_KDs_database` ready to pass directly to `TE_prediction`.
"""
function get_CO_KDs_database()
    n_el = length(TE_names)   # 28 — from lattice_strain.jl

    # (classified name, D_function, candidate original SS_syms keys)
    # "hb" and "amp" are two rows so either amphibole name returned by
    # mineral_classification will match via partition_TE's intersect.
    phase_info = [
        ("cpx", D_cpx,  [Symbol("cpx"), Symbol("pig"), Symbol("Na-cpx")]),
        ("gt",  D_gt,   [Symbol("gt")]),
        ("opx", D_opx,  [Symbol("opx")]),
        ("pl",  D_pl,   [Symbol("fsp")]),
        ("ol",  D_ol,   [Symbol("ol")]),
        ("hb",  D_amph, [Symbol("hb")]),
        ("amp", D_amph, [Symbol("amp"), Symbol("gl"), Symbol("act"),
                         Symbol("cumm"), Symbol("tr")]),
    ]

    n_ph     = length(phase_info)
    KDs_expr = Matrix{Function}(undef, n_ph, n_el)

    for (i, (_, D_func, orig_syms)) in enumerate(phase_info)
        # shared cache for the 28 closures in this row: (objectid(out), D_vec)
        cache = Ref{Tuple{UInt64,Vector{Float64}}}((UInt64(0), zeros(n_el)))

        for j in 1:n_el
            let D_f=D_func, syms=orig_syms, c=cache, j_cap=j, i_cap=i
                KDs_expr[i_cap, j_cap] = function(out)
                    liq_idx = get(out.SS_syms, :liq, 0)
                    liq_idx == 0 && return 0.0

                    min_idx = 0
                    for s in syms
                        min_idx = get(out.SS_syms, s, 0)
                        min_idx != 0 && break
                    end
                    min_idx == 0 && return 0.0

                    oid = objectid(out)
                    if c[][1] != oid
                        T_K   = out.T_C + 273.15
                        P_GPa = out.P_kbar / 10.0
                        melt  = make_comp_from_SS(out, liq_idx)
                        min_c = make_comp_from_SS(out, min_idx)
                        c[]   = (oid, D_f(T_K, P_GPa, melt, min_c))
                    end

                    return c[][2][j_cap]
                end
            end
        end
    end

    return custom_KDs_database(
        "CO — Lattice strain (Cornet 2017, TEPM v02.02) with TC a-x site fractions",
        copy(TE_names),
        [p[1] for p in phase_info],
        KDs_expr
    )
end

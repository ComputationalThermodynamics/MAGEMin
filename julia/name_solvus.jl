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
"""
    get_mineral_name(db, ss, SS_vec)
"""
# Function to return the mineral name from the solution phase name
function get_mineral_name(db, ss, SS_vec)

    mineral_name = ss
   
    if db == "ig" || db == "igad"
        x = SS_vec.compVariables
        if ss == "spl"
            if x[3] - 0.5 > 0.0;        mineral_name = "cm";
            elseif x[4] - 0.5 > 0.0;    mineral_name = "usp";
            elseif x[2] - 0.5 > 0.0;    mineral_name = "mgt";
            else                        mineral_name = "spl";    end
        elseif ss == "fsp"
            if x[2] - 0.5 > 0.0;       mineral_name = "afs";
            else                        mineral_name = "pl";    end
        elseif ss == "mu"
            if x[4] - 0.5 > 0.0;        mineral_name = "pat";
            else                        mineral_name = "mu";    end
        elseif ss == "amp"
            if x[3] - 0.5 > 0.0;        mineral_name = "gl";
            elseif -x[3] -x[4] + 0.2 > 0.0;   mineral_name = "act";
            else
                if x[6] < 0.1;          mineral_name = "cumm"; 
                elseif -1/2*x[4]+x[6]-x[7]-x[8]-x[2]+x[3]>0.5;      mineral_name = "tr";       
                else                    mineral_name = "amp";    end
            end  
        elseif ss == "ilm"
            if -x[1] + 0.5 > 0.0;       mineral_name = "hem";
            else                        mineral_name = "ilm";   end 
        elseif ss == "nph"
            if x[2] - 0.5 > 0.0;       mineral_name = "K-nph";
            else                        mineral_name = "nph";   end 
        elseif ss == "cpx"
            if x[3] - 0.6 > 0.0;        mineral_name = "pig";
            elseif x[4] - 0.5 > 0.0;    mineral_name = "Na-cpx";
            else                        mineral_name = "cpx";   end 
        end

    elseif db == "mp" || db == "mpe" || db == "mb" || db == "ume" || db == "mbe"
        x = SS_vec.compVariables
        if ss == "sp"
            if x[2] - 0.5 > 0.0;        mineral_name = "sp";
            else                        mineral_name = "mt";    end
        elseif ss == "spl"
            if x[3] - 0.5 > 0.0;        mineral_name = "cm";
            elseif x[2] - 0.5 > 0.0;    mineral_name = "mgt";
            else                        mineral_name = "sp";    end
        elseif ss == "fsp"
            if x[2] - 0.5 > 0.0;       mineral_name = "afs";
            else                        mineral_name = "pl";    end
        elseif ss == "mu"
            if x[4] - 0.5 > 0.0;        mineral_name = "pat";
            else                        mineral_name = "mu";    end
        elseif ss == "amp"
            if x[3] - 0.5 > 0.0;        mineral_name = "gl";
            elseif -x[3]-x[4]+0.2>0.0;  mineral_name = "act";
            else
                if x[6] < 0.1;          mineral_name = "cumm"; 
                elseif -1/2*x[4]+x[6]-x[7]-x[8]-x[2]+x[3]>0.5;      mineral_name = "tr";     
                else                    mineral_name = "amp";    end
            end  
        elseif ss == "ilmm"
            if x[1] - 0.5 > 0.0;        mineral_name = "ilmm";
            else                        mineral_name = "hemm";   end 
        elseif ss == "ilm"
            if 1.0 - x[1] > 0.5;        mineral_name = "hem";
            else                        mineral_name = "ilm";   end 
        elseif ss == "dio"
            if x[2] > 0.0 && x[2] <= 0.3;       mineral_name = "dio";
            elseif x[2] > 0.3 && x[2] <= 0.7;   mineral_name = "omph";
            else                                mineral_name = "jd";   end 
        elseif ss == "occm"
            if x[2] > 0.5;              mineral_name = "sid";
            elseif x[3] > 0.5;          mineral_name = "ank";  
            elseif x[1] > 0.25 && x[3] < 0.01;         mineral_name = "mag";  
            else                        mineral_name = "cc";   end
        elseif ss == "oamp"
            if x[2] < 0.3;              mineral_name = "anth";  #compositional variable y
            else                        mineral_name = "ged";   end
        end

    end

    return mineral_name
end

"""
    This function returns the solution phase name given the mineral name (handling solvus -> solution phase)
"""
function get_ss_from_mineral(db, mrl, mbCpx)

    ss = mrl
   
    if db =="ig" || db == "igad"

        if mrl == "cm" || mrl == "mgt" || mrl == "usp"
            ss = "spl"
        elseif mrl == "pat" || mrl == "mu"
            ss = "mu"
        elseif mrl == "afs" || mrl == "pl"
            ss = "fsp"
        elseif mrl == "gl" || mrl == "act" || mrl == "amp" || mrl == "cumm" || mrl == "tr"
            ss = "amp"
        elseif mrl == "hem" || mrl == "ilm"
            ss = "ilm"
        elseif mrl == "pig" || mrl == "Na-cpx"
            ss = "cpx"
        elseif mrl == "K-nph"
            ss = "nph"
        end

    elseif db == "mp" || db == "mpe" || db == "mb" || db == "ume" || db == "mbe"

        if mrl == "mt" || mrl == "sp"
            ss = "sp"
        elseif mrl == "cm" || mrl == "mgt" || mrl == "usp"
            ss = "spl"
        elseif mrl == "afs" || mrl == "pl"
            ss = "fsp"
        elseif mrl == "pat" || mrl == "mu"
            ss = "mu"
        elseif mrl == "gl" || mrl == "act" || mrl == "amp" || mrl == "cumm" || mrl == "tr"
            ss = "amp"
        elseif mrl == "hem" || mrl == "ilm"
            ss = "ilm"
        elseif mrl == "hemm" || mrl == "ilmm"
            ss = "ilmm"
        elseif mrl == "omph" || mrl == "dio" || mrl == "jd"
            if mbCpx == 0
                ss = "dio"
            else
                ss = "aug"
            end
        elseif mrl == "sid" || mrl == "mag" || mrl == "ank" || mrl == "cc"
            ss = "occm"
        elseif mrl == "anth" || mrl == "ged"
            ss = "oamp"
        end


    end

    return ss
end

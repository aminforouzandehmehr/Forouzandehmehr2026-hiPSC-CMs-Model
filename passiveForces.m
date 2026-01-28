
% Function to calculate the passive forces of the muscle
% The passive force is a function of the sarcomere length

function PF = passiveForces(SL,time,tDrugApplication,fmyo)

    SL_rest = 1.9;  % (um)
    PCon_titin = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(9)) * 1.0*0.002; %(normalised force) Def: 0.002
    PExp_titin = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(8)) * 1.0*10.0;  %(um-1)     Def: 10
    SL_collagen = 2.25; %(uM)
    PCon_collagen = 1*0.02; %(normalised force)
    PExp_collagen  = 1*70; %(um-1)
    a = 1.033947368421053;

    % Passive forces: Trabeculae: Titin and collagen, Single cells: Titin
    PF_titin = sign(SL-SL_rest)*PCon_titin*(exp(abs(SL-SL_rest)*PExp_titin)-1.0);
    
    if (SL>=SL_collagen)
        PF_collagen = PCon_collagen*(exp(PExp_collagen*(SL-SL_collagen))-1.0);
    else
        PF_collagen = 0;
    end
    
    
    PF = PF_titin + PF_collagen; %Trabeculae


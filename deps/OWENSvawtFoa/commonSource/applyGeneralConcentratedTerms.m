function [Kg,Mg,Cg] = applyGeneralConcentratedTerms(Kg,Mg,Cg,concStiff,concMass,concDamp)

if(~isempty(concStiff))
    for i=1:length(concStiff)
        nodeNum = concStiff(i).nodeNum;
        dof1 = concStiff(i).dof1;
        dof2 = concStiff(i).dof2;
        val = concStiff(i).val;
        gdof1 = (nodeNum-1)*6 + dof1;
        gdof2 = (nodeNum-1)*6 + dof2;
        
        Kg(gdof1,gdof2) = Kg(gdof1,gdof2) + val;
        
    end
end

if(~isempty(concMass))
    for i=1:length(concMass)
        nodeNum = concMass(i).nodeNum;
        dof1 = concMass(i).dof1;
        dof2 = concMass(i).dof2;
        val = concMass(i).val;
        gdof1 = (nodeNum-1)*6 + dof1;
        gdof2 = (nodeNum-1)*6 + dof2;
        
        Mg(gdof1,gdof2) = Mg(gdof1,gdof2) + val;
    end
end

if(~isempty(concDamp))
    for i=1:length(concDamp)
        nodeNum = concDamp(i).nodeNum;
        dof1 = concDamp(i).dof1;
        dof2 = concDamp(i).dof2;
        val = concDamp(i).val;
        gdof1 = (nodeNum-1)*6 + dof1;
        gdof2 = (nodeNum-1)*6 + dof2;
        
        Cg(gdof1,gdof2) = Cg(gdof1,gdof2) + val;
    end
end

end
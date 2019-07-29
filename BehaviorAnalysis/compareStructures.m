function [isDiff onlyIn1 onlyIn2] = compareStructures(struct1,struct2);
%compare structure fields
%isDiff- 1 if the two structs are different, 0 if same
%onlyIn1- fieldnames only in struct1
%onlyIn2- fieldnames only in struct2

    isDiff = 0;
    fields1 = fieldnames(struct1);
    fields2 = fieldnames(struct2);
    start1 = 1;
    onlyIn1 = [];
    for ifn = 1:length(fields1)
        if ~ismember(fields2,fields1(ifn))
            onlyIn1{start1} = fields1{ifn};
            start1 = start1+1;
        end
    end
    if start1>1
        isDiff = 1;
    end

    start2 = 1;
    onlyIn2 = [];
    for ifn = 1:length(fields2)
        if ~ismember(fields1,fields2(ifn))
            onlyIn2{start2} = fields2{ifn};
            start2 = start2+1;
        end
    end
    if start2>1
        isDiff = 1;
    end
end


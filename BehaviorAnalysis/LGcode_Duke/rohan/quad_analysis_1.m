cd('C:\Users\rohan\Documents\MATLAB\Repositories\BehaviorCode-Glickfeld-Hull-Master\BehaviorAnalysis')
%%
for i=1:length(input)
    vlength(i)=length(input{i}.tLeftTrial);
    tLeftTrial(i,1:vlength(i))=NaN;
   tLeftTrial(i,1:vlength(i))=double(celleqel2mat_padded(input{i}.tLeftTrial));
    
end
%%
% get quad times into array



for x=1:length(input.tLeftTrial)
     q_length(x)=length(input.quadratureTimesUs{x});
    end
tLeftTrial=input.tLeftTrial;    
%%
max_q=nanmax(q_length);
qtimes=nan(length(tLeftTrial), max_q);
%
%
for x=1:length(tLeftTrial)
    clear temp_qtime
    temp_qtime=double((input.quadratureTimesUs{x}));
    qtimes(x,1:length(temp_qtime)) = temp_qtime;
end
qtimesms=qtimes./1000
%%
% get quad values into array 
for x=1:length(tLeftTrial)
    q_length(x)=length(input.quadratureValues{x});
end
max_qval=nanmax(q_length);
qvals=nan(length(tLeftTrial), max_q);
for x=1:length(tLeftTrial)
    clear temp_qval
    temp_qval=double((input.quadratureValues{x}));
    qvals(x,1:length(temp_qval)) = temp_qval;
end
clear temp_qval; clear x 
%% plotting attempts
for x=1:length(tLeftTrial)
   
    p1(1,:,x)=1:length(qvals(x,:));
    p1(2,:,x)=(qvals);
    
end
hold off
%%
qfast=celleqel2mat_padded(input.qStartReact);
tRightTrial_idx=zeros(1,length(tLeftTrial));
tRightTrial_nums=find(~tLeftTrial);
tRightTrial(tRightTrial_nums)=1;
leftresponse=celleqel2mat_padded(input.tLeftResponse);
rightresponse=celleqel2mat_padded(input.tRightResponse);
ignores=intersect(find(~leftresponse), find(~rightresponse));
responded_idx=ones(1,length(tLeftTrial));
responded_idx(ignores)=0;
responded=find(responded_idx);
%%
responsevector=zeros(1,length(tLeftTrial)); %1 will mean left response, -1 equals right response)
for i=1:length(responded)
    c=responded(i)
    if qfast(c)>10
        responsevector(i)=1;
    elseif qfast(c)<-10
        responsevector(i)=-1;
    else
        responsevector(i)=0;
    end
end

early_rr=intersect(find(responsevector==-1), find(tRightTrial));
early_lr=intersect(find(responsevector==1), find(tLeftTrial));

early_rr_pct=length(early_rr)/sum(tRightTrial);
early_lr_pct=length(early_lr)/sum(tLeftTrial);
%%


   

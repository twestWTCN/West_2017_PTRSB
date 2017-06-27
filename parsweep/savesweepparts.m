function savesweepparts(appender,Nmn,MNxstore)
% This function saves each part of the parameter sweep to avoid memory
% problems

for n =1:Nmn
    for m = 1:Nmn
        tmpst = MNxstore{m,n};
            SAmax(n,m) = max(tmpst(5,:));
      [n m]
    save([pwd '\parsweep\saves\MNxstore_' num2str([n m]) appender],'tmpst')
%     save(['C:\Users\Tim\Documents\Work\OriginsofLife\cleanform\FinalGoAt\parsweep\saves\MNxstore_' num2str([n m]) appender],'tmpst')
    end
end

save('SAmax','SAmax')
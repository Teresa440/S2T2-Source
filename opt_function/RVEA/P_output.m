function P_output (Population,time,Algorithm,Problem,M)

FunctionValue = P_objective('value',Problem,M,Population);
if(strcmp(Algorithm, 'cRVEA'))
    FunctionValue = FunctionValue(:,1:end - 1);
end;
TrueValue = P_objective('true',Problem,M,1000);

NonDominated  = P_sort(FunctionValue,'first')==1;
Population    = Population(NonDominated,:);
FunctionValue = FunctionValue(NonDominated,:);

if(M == 2)
    Plot2D(TrueValue, FunctionValue, 'ro');
end;

if(M == 3)
    Plot3D(TrueValue, FunctionValue, 'ro');
end;

% eval(['save Data/',Algorithm,'/',Algorithm,'_',Problem,'_',num2str(M),'_',num2str(Run),' Population FunctionValue time'])
end



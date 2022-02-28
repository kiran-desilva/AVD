load test_weights_2

full_domain.FS = FS;
full_domain.TS = TS;
full_domain.HS = HS;
full_domain.N = N;
full_domain.weights = weights;

scatter_plot(full_domain);
multi_variable_plot(full_domain);



fail_domain = full_domain;

%filter for failure codes
nanidx = weights(:)>0;
fail_domain.FS(nanidx) = NaN;
fail_domain.TS(nanidx) = NaN;
fail_domain.HS(nanidx) = NaN;
fail_domain.N(nanidx) = NaN;
fail_domain.weights(nanidx) = NaN;

multi_variable_failure_plot(fail_domain);


function [fig] = multi_variable_plot(data)
    figure
    hold on

    X_new = [data.N(:),data.TS(:),data.HS(:),data.weights(:)];
    varNames = {'N stringers'; 'Stringer Thickness'; 'Stringer Height'; 'weights'};

    gplotmatrix(X_new,[],data.FS(:),[],[],15,[],[],varNames)
end

function [fig] = scatter_plot(data)
    figure
    hold on
    scatter3(data.TS(:),data.HS(:),data.N(:),100000*data.FS(:),data.weights(:));
    c = colorbar;
    zlabel('N Stringers')
    xlabel('Stringer Thickness [m]')
    ylabel('Stringer Height [m]')
    c.Label.String = ' Total Weight [Kg]'
end

function [fig] = multi_variable_failure_plot(data)
    figure
    hold on

    X_new = [data.N(:),data.TS(:),data.HS(:),data.FS(:)];
    varNames = {'N stringers'; 'Stringer Thickness'; 'Stringer Height'; 'Fuselage Thickness'};

    gplotmatrix(X_new,[],data.weights(:),[],[],15,[],[],varNames)
end

load test_weights

figure
hold on
% scatter3(FS(:),TS(:),HS(:),N(:)*5,weights(:))
scatter3(FS(:),TS(:),HS(:),exp(N(:)/3),weights(:))
c = colorbar;

xlabel('Fuselage Skin Thickness [m]')
ylabel('Stringer Thickness [m]')
zlabel('Stringer Height [m]')
c.Label.String = ' Total Weight [Kg]'

grid on

nanidx = isnan(weights(:));
FS(nanidx) = NaN;
TS(nanidx) = NaN;
HS(nanidx) = NaN;
N(nanidx) = NaN;

figure
hold on
bubblechart3(TS(:),HS(:),N(:),100000*FS(:),weights(:))
c = colorbar;
zlabel('N Stringers')
xlabel('Stringer Thickness [m]')
ylabel('Stringer Height [m]')
c.Label.String = ' Total Weight [Kg]'

figure
hold on

X_new = [N(:),TS(:),HS(:),weights(:)];
varNames = {'N stringers'; 'Stringer Thickness'; 'Stringer Height'; 'weights'};
% r_weights = round(weights,-2);
gplotmatrix(X_new,[],FS(:),[],[],15,[],[],varNames)

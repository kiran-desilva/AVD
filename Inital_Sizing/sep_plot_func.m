function [hObject,hCommon,hSpecific] = emplot(AC,h,n,v,AS_units,Ps_contours,r_contours,n_contours)
% EMPLOT generates an energy-maneuverability (E-M) plot.
%      [hObject,hCommon,hSpecific] = emplot(AC,h,n,v,Ps_contours,r_contours,n_contours)
%      plots the E-M characteristics of aircraft for simultaneous
%      comparison of capabilities during sustained maneuvers. 
%      
%      Inputs:
%           AC: 1xN (concatenated) structure array of aircraft properties
%               (N is number of aircraft concatenated)
%               LegendEntry: Aircraft name, string
%               DataColor: Color for limits and specific excess power contours
%               W: Gross weight, lb
%               BHP: Power*, BHP (put a 0 entry in either Power or Thrust)
%               T: Thrust*, lb (put a 0 entry in either Power or Thrust)
%               h_norm: Turbo-normalized density altitude, ft
%               eta: Engine efficiency (proportion of Power or Thrust for useful work)
%               b: Wingspan, ft
%               S: Planform area, ft^2
%               e: Oswald efficiency factor
%               n_max: Maximum positive load factor, g units
%               v_ne: Never exceed EAS, ft/s
%               v_s: Stall EAS, ft/s
%               CD0: Zero-lift drag coefficient
%           h: Density altitude, ft
%           n: Load factor vector, g units
%           v: EAS vector, ft/s
%           AS_units: Airspeed units for outputs (e.g. 'KEAS', 'KTAS', 'M')
%           Ps_contours: Specific excess power contour values vector, ft/min
%           r_contours: Turn radii contour values vector, ft
%           n_contours: Load factor contour values vector, g units
%
%      Outputs: 
%           hObject: Handle to figure
%           hCommon: Structure of common handles
%               n_contours
%               n_labels
%               r_contours
%               r_labels
%           hSpecific: 1xN (concatenated) structure array of aircraft-specific handles
%               Ps_contours
%               Ps_labels
%               limit_load
%               limit_ne
%               limit_stall
%               limit_power (hidden by default)
%
%      Furthermore, the function returns EAS for maximum turn rate (i.e.
%      maneuvering speed), EAS for maximum sustainable turn rate, and EAS
%      for maximum sustainable load factor for each aircraft evaluated.
%
%      Created by Robert Nordlund, USNA '12.

%% Unpack structure array of aircraft properties
N=length(AC); % Number of aircraft
for k=N:-1:1
    LegendEntry{k}=AC(k).LegendEntry; % Aircraft name, string
    W(k)=AC(k).W; % Gross weight, lb
    BHP(k)=AC(k).BHP; % Power, BHP (Put a 0 entry in either Power or Thrust)
    T(k)=AC(k).T; % Thrust, lb (Put a 0 entry in either Power or Thrust)
    h_norm(k)=AC(k).h_norm; % Turbo-normalized density altitude, ft
    eta(k)=AC(k).eta; % Engine efficiency (proportion of Power or Thrust for useful work)
    b(k)=AC(k).b; % Wingspan, ft
    S(k)=AC(k).S; % Planform area, ft^2
    e(k)=AC(k).e; % Oswald efficiency factor
    n_max(k)=AC(k).n_max; % Maximum positive load factor, g units
    v_ne(k)=AC(k).v_ne; % Never exceed EAS, ft/s
    v_s(k)=AC(k).v_s; % Stall EAS, ft/s
    CD0(k)=AC(k).CD0; % Zero-lift drag coefficient
end
%% Preliminaries (common)
g=32.174; % Gravitational acceleration, ft/s^2
OverlayColor=[0.7,0.7,0.7]; % Color for radius and load factor contours
[~,~,~,rho_ssl]=atmosisa(0); [~,a_h,~,rho_h]=atmosisa(h/3.2808); 
rho_ssl=rho_ssl/515.3788; rho_h=rho_h/515.3788;
a_h=a_h*3.28084; % Speed of sound at altitude, ft/s
sigma=rho_h/rho_ssl; % Density ratio
v=v/sqrt(sigma); % TAS vector, ft/s (from EAS)
v_ne=v_ne/sqrt(sigma); % Never exceed TAS, ft/s (from EAS)
v_s=v_s/sqrt(sigma); % Stall TAS, ft/s (from EAS)
%% Preliminaries (specific)
[~,~,~,rho_norm]=atmosisa(h_norm/3.2808); rho_norm=rho_norm/515.3788; % Density, slug/ft^3
%% Calculations (common)
for j=length(v):-1:1
    for i=length(n):-1:1
        omega(i,j)=g*(n(i)^2-1)^.5/v(j); % Turn rate array, rad/s
        r(i,j)=v(j)/omega(i,j); % Turn radius array, ft
    end
end
switch AS_units
    case 'KEAS'
        AS_factor=sqrt(sigma)*3600/6076;
    case 'KTAS'
        AS_factor = 3600/6076;
    case 'M'
        AS_factor = 1/a_h;
end
AS=repmat(reshape(v,1,[])*AS_factor,length(n),1); % AS array
n_array=repmat(reshape(n,[],1),1,length(v));
%% Calculations (specific)
R_a=b.^2./S; % Aspect ratio
CLmax=2*W./(rho_h*v_s.^2.*S); % Maximum lift coefficient
% Performance metrics K, A, & B reduced from Anderson's "Aircraft Performance and Design"
K=1./(pi*e.*R_a); % Combined proportionality constant for drag due to lift
A=(CD0*rho_h)./(2*W./S); % A
B=(2*K.*(W./S))/rho_h; % B
sigma_norm=rho_norm/rho_ssl;
BHP=reshape(BHP,[],1);
for k=N:-1:1
    n_stallu(k,:)=rho_h*v.^2*CLmax(k)*S(k)/(2*W(k)); % Stall load factor (upper)
    if T(k)~=0
        BHP(k,1:length(v))=T(k)*v/550; % Power, BHP (varies with TAS)
    else
        BHP(k,1:length(v))=BHP(k); % Power, BHP (invariant with TAS)
    end
    if h>h_norm(k)
        Pa(k,:)=BHP(k,:)*eta(k)*550/W(k)*sigma/sigma_norm(k); % Specific power available, ft/s
    else
        Pa(k,:)=BHP(k,:)*eta(k)*550/W(k); % Specific power available, ft/s
    end
    % Calculate specific excess power
    for j=length(v):-1:1
        for i=length(n):-1:1
            Pr=A(k)*(v(j))^3+B(k)*n(i)^2/v(j); % Specific power required, ft/s
            Ps(i,j,k)=Pa(k,j)-Pr; % Specific excess power, ft/s
        end
    end
    % Calculate limits
    n_s(k,:)=real(((Pa(k,:)-A(k)*v.^3).*v/B(k)).^.5); % Maximum sustainable load factor
    v_ns(k)=v(n_s(k,:)==max(n_s(k,:))); % TAS for maximum sustainable load factor
    omega_nm(k,:)=real(g*(n_max(k).^2-1).^.5./v); % Turn rate (maximum load factor), rad/s
        % For maximum sustainable load factor above, use max(n_s(k,:)) instead of n_max(k)
    omega_stall(k,:)=real(g*(n_stallu(k,:).^2-1).^.5./v); % Turn rate (stall), rad/s
    omega_s(k,:)=real(g*(n_s(k,:).^2-1).^.5./v); % Turn rate (sustainable), rad/s
    % Correct limits
    omega_s(k,omega_s(k,:)>=omega_stall(k,:))=NaN; omega_s(k,find(omega_s(k,:)==0,1)+1:end)=NaN;
    omega_nm(k,omega_nm(k,:)>omega_stall(k,:))=NaN; omega_nm(k,find(v>v_ne(k),1):end)=NaN;
    omega_stall(k,omega_stall(k,:)>omega_nm(k,:))=NaN; omega_stall(k,1:find(omega_stall(k,:),1)-2)=NaN; omega_stall(k,v>=v_ne(k))=NaN;
    %
    v_omegaA(k,:)=v(omega_stall(k,:)==max(omega_stall(k,:)));% TAS for maximum turn rate
    v_omegas(k,:)=v(omega_s(k,:)==max(omega_s(k,:))); % TAS for maximum sustainable turn rate
end
AS_units_BC='';
AS_units_AD=[AS_units, ' '];
if AS_units=='M'
    AS_units_BC=[AS_units, ' ']; AS_units_AD='';
end
fprintf(['V_A (',AS_units,' for maximum turn rate, i.e. maneuvering speed):\n'])
for k=1:N
    fprintf(['\tfor the ',char(LegendEntry(k)),': ',AS_units_BC,'%.2f ',AS_units_AD,'for %.1f deg/s and %.0f ft/min\n'],...
        [v_omegaA(k,:)*AS_factor,max(omega_nm(k,:))*180/pi,Ps(find(n>=n_max(k),1),v==v_omegaA(k,:),k)*60])
end
fprintf('\n')
fprintf(['V_omega(s) (',AS_units,' for maximum sustainable turn rate):\n'])
for k=1:N
    fprintf(['\tfor the ',char(LegendEntry(k)),': ',AS_units_BC,'%.2f ',AS_units_AD,'at %.1fg for %.1f deg/s\n'],...
        [v_omegas(k)*AS_factor,((v_omegas(k)*max(omega_s(k,:))/g)^2+1)^.5,max(omega_s(k,:))*180/pi])
end
fprintf('\n')
fprintf(['V_n(s) (',AS_units,' for maximum sustainable load factor):\n'])
for k=1:N
    fprintf(['\tfor the ',char(LegendEntry(k)),': ',AS_units_BC,'%.2f ',AS_units_AD,'at %.1fg for %.1f deg/s\n'],...
        [v_ns(k)*AS_factor,max(n_s(k,:)),(g*(max(n_s(k,:)).^2-1).^.5./v_ns(k))*180/pi])
end
fprintf('\n')
%% Determine axis limits
% Adapted from PlotEMDiagram.m (02 NOV 2011) by Assoc. Prof. D.S. Miklosovic
%%
% Determine x-axis based on 1-2-5 rule and set plot limits...
deltas=[500,200,100,50,20,10,5,2,1];
x_lims(1)=min(min(AS)); % Get the minimum from airpseed array
x_lims(2)=max(max(AS)); % Get the maximum from airspeed array
x_range=x_lims(2)-x_lims(1);
x_contours=x_range./deltas+1;
x_delta=deltas(find(x_contours<10,1,'last'));
x_lims(1)=x_delta*floor(x_lims(1)/x_delta); % set new airspeed lower limit for the x-axis
x_lims(2)=x_delta*ceil(x_lims(2)/x_delta); % set new airspeed upper limit for the x-axis
% Determine y-axis based on 1-2-5 rule and set plot limits...
y_lims(1)=0;
y_lims(2)=40; %max(max(omega_stall))*180/pi;
y_range=y_lims(2)-y_lims(1);
y_contours=y_range./deltas+1;
y_delta=deltas(find(y_contours<10, 1,'last'));
y_lims(1)=y_delta*floor(y_lims(1)/y_delta); % set new omega lower liit for the y-axis
y_lims(2)=y_delta*ceil(y_lims(2)/y_delta); % set new omega upper limit for the y-axis
hObject=figure;
axes('XLim',[x_lims(1),x_lims(2)],'XTick',x_lims(1):x_delta:x_lims(2),...
    'YLim',[y_lims(1),y_lims(2)],'YTick',y_lims(1):y_delta:y_lims(2));
set(gcf,'Name',sprintf('Energy-Maneuverability (DA=%.0f ft)',h),'Color','White',...
    'renderer','zbuffer','Colormap',lines(N)); % Alternate renderer enables dashed contours
cmap=get(gcf,'colormap'); 
title(get(gcf,'Name'));
hold on
xlabel(['Airspeed, V_\infty (',AS_units,')'])
ylabel('Turn Rate, \omega (deg/sec)')
%% Plots (common)
hold on
% Plot load factor contours
[hc_n,hh_n]=contour(AS,omega*180/pi,n_array,n_contours,':','LineColor',OverlayColor);
hl_n=clabel(hc_n,hh_n,'Color',OverlayColor,'BackgroundColor','White');
set(hl_n,{'String'},strcat(get(hl_n,'String'),' g'))
% Plot radius contours
[hc_r,hh_r]=contour(AS,omega*180/pi,r,r_contours,':','LineColor',OverlayColor);
hl_r=clabel(hc_r,hh_r,'Color',OverlayColor,'BackgroundColor','White');
set(hl_r,{'String'},strcat(get(hl_r,'String'),' ft'))
% Assemble common handles
hCommon.n_contours=hh_n;
hCommon.n_labels=hl_n;
hCommon.r_contours=hh_r;
hCommon.r_labels=hl_r;
%% Plots (specific)
% Specific excess power contours
for k=N:-1:1
    [hc_Ps,hh_Ps(k)]=contour(AS,omega*180/pi,Ps(:,:,k)*60,Ps_contours,'Color',cmap(k,:));
    hl_Ps=clabel(hc_Ps,hh_Ps(k),'Color',cmap(k,:),'BackgroundColor','White');
    set(hl_Ps,{'String'},strcat(get(hl_Ps,'String'),' ft/min'))
    % Begin code to change specific line
    ha_Ps=get(hh_Ps(k),'Children'); % get a handle to all the children of the contourgroup object
    for i=1:length(ha_Ps) % loop through all the children of the contourgroup object
        % Only consider the children that are patch objects
        % STRCMP compares two strings and return true if they are the same
        if strcmp(get(ha_Ps(i),'Type'),'patch') % 'UserData' indicates the elevation associated with a patch
            switch get(ha_Ps(i),'UserData')
                case 0 % Line at elevation 0
                    set(ha_Ps(i),'LineWidth',2.5);
            end
        end
    end
    hSpecific(k).Ps_contours=hh_Ps(k);
    hSpecific(k).Ps_labels=hl_Ps;
end
for k=N:-1:1
    % Plot limits
    hn(1)=line(v*AS_factor,omega_nm(k,:)*180/pi); % Load limit
    hn(2)=line(v*AS_factor,omega_stall(k,:)*180/pi); % Stall limit (upper)
    hn(3)=line([v_ne(k),v_ne(k)]*AS_factor,[0,omega_nm(k,find(v<=v_ne(k),1,'last'))]*180/pi); % Never exceed limit
    hn(4)=line(v*AS_factor,omega_s(k,:)*180/pi,'Visible','off'); % Power limit (Ps=0 contour)
    set(hn,'Color',cmap(k,:),'LineWidth',2.5)
    hSpecific(k).limit_load=hn(1);
    hSpecific(k).limit_stall=hn(2);
    hSpecific(k).limit_ne=hn(3);
    hSpecific(k).limit_power=hn(4);
end
legend(hh_Ps,char(LegendEntry))
box on

end
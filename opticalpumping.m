%% very general formulas used in multiple sections
pcterr= @(exp,theo) abs(exp-theo)/theo*100;
del_posquotient= @(a,del_a,b,del_b) (a./b).*((del_a./a)+(del_b./b));

%% 4B: calculate currents, fields, Lande' g-factors, and their ratios
R_bar=0.1639; %in meters
N=11;
f=8.991e-3;
BSwp= @(ISwp) f*ISwp*N/R_bar; %in gauss
Bres=BSwp(0.280);
B87=BSwp(0.585)-Bres;
B85=BSwp(0.799)-Bres;

mu0divbyh=1.3996; %MHz/gauss -- given in manual
gF= @(B,nu) (nu./B)/mu0divbyh; %(nu*h)/(mu0*B) is rearranged from the equation in the manual and uses the value of mu0/h for the units used here provided in the manual

RF=  [0.000,0.100,0.120,0.135,0.150,0.170,0.190,0.200]; %applied RF in MHz
ZFR= [0.282,0.283,0.282,0.282,0.282,0.281,0.282,0.281]; %voltage (in volts) at which the ZFR occurs for each RF
Rb87=[0.282,0.516,0.558,0.591,0.628,0.679,0.729,0.751]; %voltage (in volts) at which the 87Rb resonance occurs for each RF
Rb85=[0.282,0.632,0.696,0.746,0.801,0.875,0.953,NaN  ]; %voltage (in volts) at which the 85Rb resonance occurs for each RF
gF87=gF(BSwp(Rb87)-Bres,RF); %!!discard 1st element b/c it doesn't correspond to anything physical -- the number there is just a placeholder
gF85=gF(BSwp(Rb85)-Bres,RF); %!!discard 1st and last elements b/c they don't correspond to anything physical -- these numbers are just placeholders
pcterrgF87=pcterr(gF87,1/2);
pcterrgF85=pcterr(gF85,1/3);
landeratio=gF87./gF85;

%% 4B: linear regression on current - RF data and ratio of g-factors calculated this way
mdl87=fitlm(Rb87,RF);
mdl85=fitlm(Rb85,RF);

Rb87int=mdl87.Coefficients{1,1};
Rb85int=mdl85.Coefficients{1,1};
Rb87slope=mdl87.Coefficients{2,1};
Rb85slope=mdl85.Coefficients{2,1};
RB87slUnc=mdl87.Coefficients{2,2};
RB85slUnc=mdl85.Coefficients{2,2};
rb87rsq=mdl87.Rsquared.Ordinary;
rb85rsq=mdl85.Rsquared.Ordinary;
slopesratio=Rb87slope/Rb85slope;

Rb87eqn=['87Rb: y=',num2str(Rb87slope,3),'x',num2str(Rb87int,3),'; R^2=',num2str(rb87rsq,3),'.00']; %hard-code the significant zeroes that MATLAB leaves implicit
Rb85eqn=['85Rb: y=',num2str(Rb85slope,3),'x',num2str(Rb85int,3),'; R^2=',num2str(rb85rsq,3)];

%% 4B: propagated uncertainties for all derived quantities in this section
%%DEFINING ANONYMOUS FUNCTIONS FOR ALL THE ERRORS THAT NEED TO BE PROPAGATED IN THIS SECTION
del_RF=0.005; %MHz
del_Rmon=0.01;
del_V=0.009;
del_I= @(V) del_V+V*del_Rmon;
del_Bnet= @(V_ZF,V_iso) ((f*N)/(R_bar))*(2*del_V+(V_iso-V_ZF)*del_Rmon); %should spit out a 1x8 where none of the values are NaN ... but, discard the first entry, because the original expression evaluates to zero at zero RF b/c there are no isotopes, so the net is just ZF minus ZF
del_gF= @(V,nu) ((R_bar)./(f*N*V*mu0divbyh)).*(del_RF+nu*del_Rmon+nu.*del_V./V); %should spit out a 1x8 where only the first (and if 85, also last) value is NaN

del_Bres= @(V_ZF) ((f*N)/(R_bar))*(del_V+V_ZF*del_Rmon); %really only need to call this once -- same result for both isotopes
del_Bres_both=del_Bres(ZFR(1)); %just need this to populate the RF-0 row for both isotopes

%%ERROR PROPAGATION CALCULATIONS FOR 87RB -- for all three resulting uncertainty vectors, discard 1st element b/c it doesn't correspond to anything physical -- the number there is just a placeholder
del_I_87=del_I(Rb87); % currents - 87 
del_Bnet_87=del_Bnet(ZFR,Rb87); % net sweep field - 87
del_gF_87=del_gF(Rb87,RF); % Lande' g-factor - 87

%%ERROR PROPAGATION CALCULATIONS FOR 85RB -- for all three resulting uncertainty vectors, ignore the first and last entries because they are just placeholders and don't represent anything physical
del_I_85=del_I(Rb85); % currents - 85 
del_Bnet_85=del_Bnet(ZFR,Rb85); % net sweep field - 85
del_gF_85=del_gF(Rb85,RF); % Lande' g-factor - 85

%Lande' g-factor ratio
del_Landeratio=del_posquotient(gF87,del_gF_87,gF85,del_gF_85);

%uncertainty in the ratio of slopes of the current-frequency data for this section 
del_slopesratio=del_posquotient(Rb87slope,RB87slUnc,Rb85slope,RB85slUnc);


%% 4B: low-field Zeeman effect plot using the linear regression calculated above
%piecing together vectors for error bars
errI87=[0,del_I_87(2:8)]; %discard 1st element b/c it is not physically meaningful
errI85=[0,del_I_85(2:7),0]; %discard 1st and last elements b/c they are meaningless placeholders
errRF=del_RF*ones(8,1);

%actually making the figure
figure 
hold on
plot(mdl87,'MarkerSize',3,'Marker','o','MarkerFaceColor','red')
han85=plot(mdl85,'MarkerSize',3,'Marker','o','MarkerFaceColor','blue');
fitHandle = findobj(han85,'DisplayName','Fit');
cbHandles = findobj(han85,'DisplayName','Confidence bounds');
cbHandles = findobj(han85,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
fitHandle.Color = 'blue';
set(cbHandles, 'Color', 'b', 'LineWidth', 1)
xlabel('sweep coil current (Amperes)')
ylabel('transition frequency (MHz)')
title('low-field Zeeman effect for 87Rb and 85Rb')
axis([0.2,1.0,0,0.250])
text(Rb87(3),RF(4),Rb87eqn,'Position',[0.3,0.14,0])
text(Rb85(3),RF(4),Rb85eqn,'Position',[0.7,0.11,0])
errorbar(Rb87,RF,-errRF,errRF,-errI87,errI87,'.','Color','red')
errorbar(Rb85,RF,-errRF,errRF,-errI85,errI85,'.','Color','blue')
legend('measured data','linear regression','confidence bounds','Location','southeast')
hold off

%% 4D: plot period as a function of RF amplitude (1/x relationship)
RF=  [1.0, 1.5, 2.0, 2.5, 3.0, 3.5];
Rb87=[5.0, 3.0, 2.5, 2,   1.8, 1.6]; %from half-increment convergent estimate
Rb85=[8,   5.4, 4.2, 2.4, 2.8, 2.4]; %slopes ratio = 1.3697

perioderr85=2*[0.5, 0.2, 0.1, 0.2,0.1,0.1]; %propagate the half-period errors into full-period errors
perioderr87=2*[0.25,0.25,0.25,0.5,0.1,0.1];
invpererr85=perioderr85./(Rb85.^2);
invpererr87=perioderr87./(Rb87.^2);
rabierrRF=0.0001*ones(6,1);

figure
plot(RF,Rb85,'o','MarkerFaceColor','blue')
hold on
plot(RF,Rb87,'o','MarkerFaceColor','red')
title('Rabi periods as a function of RF amplitude')
errorbar(RF,Rb85,-perioderr85,perioderr85,-rabierrRF,rabierrRF,'.','Color','blue')
errorbar(RF,Rb87,-perioderr87,perioderr87,-rabierrRF,rabierrRF,'.','Color','red')
legend('Rb85','Rb87')
xlabel('RF amplitude (V)')
ylabel('period of Rabi oscillations (ms)')
hold off

%% 4D: linear regression on and plot of inverse Rabi period as a function of RF amplitude
mdl85=fitlm(RF,Rb85.^-1);
mdl87=fitlm(RF,Rb87.^-1);

rabislope85=mdl85.Coefficients{2,1};
rabislope87=mdl87.Coefficients{2,1};
rabislunc85=mdl85.Coefficients{2,2};
rabislunc87=mdl87.Coefficients{2,2};
rabiint85=mdl85.Coefficients{1,1};
rabiint87=mdl87.Coefficients{1,1};
rb85rsq=mdl85.Rsquared.Ordinary;
rb87rsq=mdl87.Rsquared.Ordinary;

Rb85eqn=['87Rb: y=',num2str(rabislope87,3),'x+',num2str(rabiint87,3),'; R^2=',num2str(rb85rsq,3)];
Rb87eqn=['85Rb: y=',num2str(rabislope85,3),'x+',num2str(rabiint85,3),'; R^2=',num2str(rb87rsq,3)];

figure
han85=plot(mdl85,'MarkerSize',3,'Marker','o','MarkerFaceColor','blue'); %in this block, I referenced https://www.mathworks.com/matlabcentral/answers/538543-how-to-change-the-color-of-different-parts-of-regression-plots-using-fitlm-function
dataHandle = findobj(han85,'DisplayName','Data');
fitHandle = findobj(han85,'DisplayName','Fit');
cbHandles = findobj(han85,'DisplayName','Confidence bounds');
cbHandles = findobj(han85,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
fitHandle.Color = 'blue';
set(cbHandles, 'Color', 'b', 'LineWidth', 1)
hold on
han87=plot(mdl87,'MarkerSize',3,'Marker','o','MarkerFaceColor','red');
title('linearized 87Rb and 85Rb Rabi period data as a function of RF amplitude')
text(Rb87(5),RF(5),Rb87eqn,'Position',[0.8,0.45,0])
text(Rb85(5),RF(5),Rb85eqn,'Position',[2.4,0.23,0])
xlabel('RF amplitude (V)')
ylabel('inverse Rabi period (ms)^-1')
errorbar(RF,Rb85.^-1,-invpererr85,invpererr85,-rabierrRF,rabierrRF,'.','Color','blue')
errorbar(RF,Rb87.^-1,-invpererr87,invpererr87,-rabierrRF,rabierrRF,'.','Color','red')
legend('measured data','linear regression','confidence bounds','Location','southeast')
hold off

%% 4D: ratio of slopes for linear regression and comparison to ratio of g_F
rabislopesratio=rabislope87/rabislope85;
rabislopesratiopcterr=pcterr(rabislopesratio,1.5);

%% 4D: propagated uncertainty in the ratio of slopes [use standard error from lmfit()]
del_posquotient= @(a,del_a,b,del_b) (a./b).*((del_a./a)+(del_b./b));
rabislopesratioerr=del_posquotient(rabislope87,rabislunc87,rabislope85,rabislunc85);
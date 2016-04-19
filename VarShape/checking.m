clear all;
close all;

alpha=2e-4
W=4.5;%define the work function
d=2e4;

debug=0;
maxiterR=10;
maxiterF=15;
Rtol=0.01;
Ftol=0.001;

Fmaximum=W^2/1.44;
Fminimum=Fmaximum/20;
Rminimum=1;
Rmaximum=30;

%This script fits the generalized FN equation to experimental data in order
%to extract the following parameters: enhancementa factor beta, radius of
%curvature R and prefactor sigma*Aefficient

%% give the experimental data straight from mat file
%There must be two vectors in the mat file: xexp and yexp. They contain the
%I-V data in the corresponding FN coordinates. xexp=1/V, yexp=log(I/V^2).
%Experimental data must be very smooth and accurate. Otherwise the method faces
%instabilities and it would be much better to fit a polynomial to the
%experimental data and give the data in polynomial form. In this case
%activate lines 26-30 in the next section.

% load ExpData1.mat;
% xexp=xexp1;
% yexp=yexp1;
fromCSV=csvread('./Exp_data/2003_CERNpaper.csv');
xexp=fromCSV(:,1);
yexp=fromCSV(:,2);

%% give experimental data in polynomial form
% fit a polynomial to the data of the form yexp=p1*xexp^2+p2*xexp+p3
% activate this section if the data is noisy and you choose to fit a
% polynomial first
[f1,g1]=fit(xexp,yexp,'poly2');%fit polynomial
xstart=min(xexp); %give starting xexp
xend=max(xexp); %give last xexp  
xe=(xstart:(xend-xstart)/50:xend)';
ye=f1(xe);%calculate values

%% main section

hoRofF= @(Floc,Fmac,alpha) (Floc-2*Fmac)./(Fmac+Fmac.*alpha.*Floc.^2);
Ie=1e6*exp(ye)./xe.^2;
ye=log(Ie);
Fmac=1./xe;
Ra=Rminimum;%extreme values for R
Rb=Rmaximum;
for i=1:maxiterR
    R=0.5*(Ra+Rb);%bisection for R
    Fa=Fminimum;%extreme maximum fields
    Fb=Fmaximum;
    for j=1:maxiterF
        Fmax=0.5*(Fa+Fb);%bisection for maximum field
        hoR=hoRofF(Fmax,max(Fmac),alpha);
        F=arrayfun(@(Fm) FlocofFmac(Fm,hoR,alpha),Fmac);
        JofF=@(f) J_sph_approx(f,W,R);
        I=arrayfun(JofF,F);%current density
        yanal=log(I);%FNplot vertical variable
        sigma=max(yanal)-max(ye);%FN plot vertical offset to match
        yanal=yanal-sigma;%match plots at the left by appropriate prefactor
        if min(yanal)-min(ye)>Ftol%
            Fb=Fmax;%bisection to the correct side, match plots at the right
        elseif min(ye)-min(yanal)>Ftol%this way beta or Fmax is defined
            Fa=Fmax;
        else
            break;%if they match with 1% accuracy break the loop
        end
        if debug
            plot(xe,log(Ie),'*')
            hold on
            plot(xe,yanal,'r')
            pause(0.01)
        end
    end
    if mean(yanal)-mean(ye)>Rtol
        Ra=R;%bisection at the correct side
    elseif mean(ye)-mean(yanal)>Rtol
        Rb=R;%once left-right endpoints are matched, R is tuned to match the middle
    else
        break;%if the error at the mean is smaller than 0.001% break
    end
end
sigmaAeff=exp(-sigma);%correcting prefactor
V=Fmac*d;

zeta=(I./F.^2).*sqrt(V)*1.904e5;
yy=1.44*F/W^2;

%%   plotting the results

figure;
plot(xe,ye,'*');%plot data with V in the free variable
hold on
plot(xe,yanal,'r');
xlabel('1/F_{mac} (nm/V)');
ylabel('log(I) A');
disp(['The radius of curvature is R=' num2str(R) ' nm']);%show results of R, beta and sigmaAeff
disp(['The hoR is =' num2str(hoR)]);
disp(['The pre-factor is sigmaAeff=' num2str(sigmaAeff) ' nm^2']);

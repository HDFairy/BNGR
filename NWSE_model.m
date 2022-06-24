tic ;

%% Load SAR-SFMR dataset 
Data = load('SAR_SFMRWinds_Dataset.mat') ;
WSpd = Data.BNGR_SFMR_WSpd;WSpd=WSpd'; % SFMR wind speed
NRCS_VH = Data.BNGR_NRCS_VH ;NRCS_VH=NRCS_VH' ;   % SAR VH-NRCS
NRCS_VV = Data.BNGR_NRCS_VV ;NRCS_VV=NRCS_VV' ;    % SAR VV-NRCS
Angle = Data.BNGR_Angle ;Angle=Angle' ;
Angle = Angle/180 * pi;        % SAR incidence angle, need to be converted to rad unit

N = length(WSpd) ;  % The size of dataset

%% Dataset and parameters for selected BNGR model
X_SAR_SFMR = [NRCS_VH; Angle];
Y_SAR_SFMR = WSpd;

%% Parameters;
s1 = 1.1615e-04 ;
s2 = 1.5075 ;

%% SAR Prediction: e.g., Irma TC on 10:30 UTC, September 7, 2017 
SAR = load('IRMA_20170907_S1A.mat') ;
SAR_NRCS_VH = SAR.NRCS_VH_3KM;
SAR_NRCS_VV = SAR.NRCS_VV_3KM;
SAR_ANGLE = SAR.Angle_3KM ;
SAR_ANGLE = SAR_ANGLE/180 * pi;
SAR_Lon = SAR.Lon_3KM ;
SAR_Lat = SAR.Lat_3KM ;


for i = 1 : length(SAR_ANGLE(:,1))
    for j = 1 : length(SAR_ANGLE(1,:))
        if isnan(SAR_ANGLE(i,j))==1
            SAR_ANGLE(i,j)=nanmean(SAR_ANGLE(:,j));
        end
    end
end

BNGR_WSpd = SAR_NRCS_VH;
for m = 1 : size(SAR_NRCS_VH, 1)
    for n = 1 : size(SAR_NRCS_VH, 2)
        if isnan(SAR_NRCS_VH(m, n))==0
            sum1 = 0 ;
            sum2 = 0 ;
            X_SAR = [SAR_NRCS_VH(m, n); SAR_ANGLE(m, n)];

            for j = 1 : N   
                temp1 = (X_SAR - X_SAR_SFMR(:, j))' * (X_SAR - X_SAR_SFMR(:, j)) ;
                temp2 = exp(-2 * ((X_SAR - X_SAR_SFMR(:, j))' * (X_SAR - X_SAR_SFMR(:, j)))) ;
                sum1 = sum1 + temp1 ;
                sum2 = sum2 + temp2 ;
            end
            SAR_delta1 = (s1 / N) * sum1 ;
            SAR_delta2= s2 / sum2 ;
            
            sum3 = 0 ;
            sum4 = 0 ;
            for j = 1 : N
                temp3 = Y_SAR_SFMR(j) * exp(-((X_SAR - X_SAR_SFMR(:, j))' * (X_SAR - X_SAR_SFMR(:, j))) / (2 * SAR_delta1)) ;
                temp4 = exp(-((X_SAR - X_SAR_SFMR(:, j))' * (X_SAR - X_SAR_SFMR(:, j))) / (2 * SAR_delta1)) ;
                sum3 = sum3 + temp3 ;
                sum4 = sum4 + temp4 ;
            end
            BNGR_WSpd(m, n) = sum3 / sum4 ;
        end
    end
end

%%
mycolor = load('mycolor.mat') ;
Lms = [min(min(SAR_Lon)) max(max(SAR_Lon)) min(min(SAR_Lat)) max(max(SAR_Lat))];
figure('Visible','on')
set(gcf, 'Units', 'centimeter', 'Position', [5 5 40 28]) ;
m_proj('miller','lon',[Lms(1) Lms(2)],'lat',[Lms(3) Lms(4)])
m_pcolor(SAR_Lon,SAR_Lat,BNGR_WSpd)
m_grid('tickdir', 'in', 'linewi', 0.5, 'linestyl', '-.', 'fontsize', 16) ;
m_coast('patch', [.7 .7 .7], 'edgecolor', 'none') ;
colormap(mymap('jet')) ;
hold on
shading interp
caxis([0 75])
title('NWSE', 'fontsize', 20) ;
h=colorbar('eastoutside') ;
set(get(h, 'ylabel'), 'String', 'Wind Speed (m/s)', 'fontsize', 20) ;
set(gcf, 'color', 'w') ;
set(gca, 'FontSize', 20) ;


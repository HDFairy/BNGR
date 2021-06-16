tic ;

%% Load SAR-SFMR dataset N = 519
Data = load('SAR-SFMR Dataset.mat') ;
WSpd = Data.Wind_Speed';               % Measured wind speed
NRCS_VH = Data.VH_NRCS' ;   % SAR VH-NRCS
NRCS_VV = Data.VV_NRCS' ;    % SAR VV-NRCS
Angle = Data.Incidence_Angle' ;
Angle = Angle/180 * pi;        % SAR incidence angle, need to be converted to rad unit

N = length(WSpd) ;  % The size of dataset

%% Dataset and parameters for selected BNGR model
X_SAR_Winds = [NRCS_VH; Angle];
Y_SAR_Winds = WSpd;

%% Parameters
s1 = 6.7103e-04;
s2 = 1.8247;

%% SAR Prediction: e.g., Irma TC on 10:30 UTC, September 7, 2017 
Irma_data = load('IRMA_20170907_S1A.mat') ;
Irma_VH_NRCS = Irma_data.Irma_VH_NRCS;
Irma_Inc_angle = Irma_data.Irma_Inc_angle ;
Irma_Inc_angle = Irma_Inc_angle/180 * pi;
Irma_Lon = Irma_data.Irma_Lon ;
Irma_Lat = Irma_data.Irma_Lat ;

%% Estimation SAR wind speed with Equation (6)
BNGR_WSpd = zeros(size(Irma_VH_NRCS,1),size(Irma_VH_NRCS,2));
for m = 1 : size(Irma_VH_NRCS, 1)
    for n = 1 : size(Irma_VH_NRCS, 2)
        if isnan(Irma_VH_NRCS(m, n))==0
            sum1 = 0 ;
            sum2 = 0 ;
            X_SAR = [Irma_VH_NRCS(m, n); Irma_Inc_angle(m, n)];
            %% Computer the smoothing and prediction-error parameters with Equation (21) (22)
            for j = 1 : N   
                temp1 = (X_SAR - X_SAR_Winds(:, j))' * (X_SAR - X_SAR_Winds(:, j)) ;
                temp2 = exp(-2 * ((X_SAR - X_SAR_Winds(:, j))' * (X_SAR - X_SAR_Winds(:, j)))) ;
                sum1 = sum1 + temp1 ;
                sum2 = sum2 + temp2 ;
            end
            Delta = (s1 / N) * sum1 ;
            Gamma= s2 / sum2 ;
            
            sum3 = 0 ;
            sum4 = 0 ;
            for j = 1 : N
                temp3 = Y_SAR_Winds(j) * exp(-((X_SAR - X_SAR_Winds(:, j))' * (X_SAR - X_SAR_Winds(:, j))) / (2 * Delta)) ;
                temp4 = exp(-((X_SAR - X_SAR_Winds(:, j))' * (X_SAR - X_SAR_Winds(:, j))) / (2 * Delta)) ;
                sum3 = sum3 + temp3 ;
                sum4 = sum4 + temp4 ;
            end
            BNGR_WSpd(m, n) = sum3 / sum4 ;
        end
    end
end

%% Drawing
figure('Visible','on')
pcolor(Irma_Lon, Irma_Lat, BNGR_WSpd)
shading interp
colormap(jet)
xlabel('Longitude');
ylabel('Latitude');
ylabel(colorbar, 'Wind Speed (m/s)');
title('Irma (2017)')
set(gca,'FontSize',14,'FontWeight','bold','Fontname', 'Times New Roman');

toc;

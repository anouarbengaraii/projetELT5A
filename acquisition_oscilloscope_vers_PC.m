% Programme de réception
%  Oscilloscope AGILENT DSO/MSO
%

%Programme:
%%
global plot1  plot2 plot3 plot4 Y2 f2 dist Fb Distance psor fbeat line dist2 Yref;  
global acquisitionState;

%Paramètres:
global permittiviteR
permittiviteR = 1; % Permitivité diélectrique relative du Tflon (PTFE). Subtrat diélectrique utilisé pour les cables coaxiaux SMA de l'expérience. 
global c
c = 3e8/sqrt(permittiviteR);
global Tm 
Tm = 200e-3;
global B
B = 2e8;
global k 
k = B/Tm;

% Initialisation
[visaObj, NbPts, NbAcq] = initializeScope();

digitize(visaObj);
[tps, ch1] = captureChannel(visaObj, 'CHAN1', NbPts);
[tps, ch2] = captureChannel(visaObj, 'CHAN2', NbPts);

% Afficher le graphique: 
[fig, Y1, f1] = AfficherFigure(tps, ch1, ch2)
afficherControls(fig, Y1, f1);% Afficher le graphique: 

%Attendre 2sec
pause(2);

%Indiquer à l'utilisateur que la mesure démarre dans 5sec:
updateAcquisitionState('Acquisition du signal de rérérence...');

%Mettre à jour le signal de référence
updateYref(visaObj, NbPts);

%Attendre 2sec
pause(1);

%Indiquer à l'utilisateur que la mesure démarre dans 5sec:
updateAcquisitionState('Acquisition des signaux...');

%Attendre 2sec
pause(1);

% Capture des signaux en continue
while true

    %Capturer les signaux de l'oscilloscope
    digitize(visaObj);
    [tps, ch1] = captureChannel(visaObj, 'CHAN1', NbPts);
    [tps, ch2] = captureChannel(visaObj, 'CHAN2', NbPts);

    %Filtrer le signal parasite
    ch1 = FiltreParasite(ch1);
    
    %Moyennage des chirps:
    [ch1_temp, ch2_temp, tps_temp, Y1, Y2, f1, f2, psorM, fbeatM, DistanceM] = MoyennageChirps(ch1, ch2, tps);
       
    %Mettre à jour l'UI
    updateWaveform(ch1_temp, ch2_temp, tps_temp, Y1, Y2, f1, f2);
    updateControls(DistanceM, fbeatM, psorM);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fig, Y1, f1] = AfficherFigure(tps, ch1, ch2)

global plot1 plot2 plot3 plot4;

% Ouvrir une fentre pour afficher les graphiques:
fig = figure(name="Acquisition des signaux de EV-RADAR MMIC2");
fig.WindowState = 'maximized'; % Met la fenêtre en mode plein écran

%% Temporel:
% CHANNEL1:
subplot1 = subplot(2, 2, 1);
plot1 = plot(tps, ch1, "r");
title("Signal de Battement");
xlabel('Temps (s)');
ylabel('Amplitude (mV)');
legend("CHANNEL1");
subplot1.Position = subplot1.Position + [0.05 0 0 0]

% CHANNEL2:
subplot2 = subplot(2, 2, 2);
plot2 = plot(tps, ch2, "g");
title("Tension de Commande du VCO");
xlabel('Temps (s)');
ylabel('Amplitude');
legend("CHANNEL2");
subplot2.Position = subplot2.Position + [0.05 0 0 0]

% Calcul des FFT:
[Y1, f1] = calculerFFT(ch1, tps);
[Y2, f2] = calculerFFT(ch2, tps);

%% Fréquentiel:
% CHANNEL1:
subplot3 = subplot(2, 2, 3);
plot3 = plot(f1, Y1, "r", 'LineWidth', 2);
title("Fréquence de Battement");
xlabel('Fréquence (Hz)');
ylabel('Amplitude');
xlim([0 100])
AfficherDistance(Y1, f1);
legend("CHANNEL1", "");
subplot3.Position = subplot3.Position + [0.05 0 0 0]

% CHANNEL2:
subplot4 = subplot(2, 2, 4);
plot4 = plot(f2, Y2, "g", 'LineWidth', 2);
title("FFT de Tension de Commande du VCO");
xlabel('Fréquence (Hz)');
ylabel('Amplitude');
xlim([0 1000]);
legend("CHANNEL2");
subplot4.Position = subplot4.Position + [0.05 0 0 0]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AfficherDistance(Y, f)
    global c k;
    global line dist2 Fb2;

    [psor, fbeat, Distance] = Mesure(Y, f, c, k)

    hold on;

    dist2 = text(fbeat + 30,psor-0.01, "Distance = " + Distance + "m", 'Color', 'b');

    Fb2 = text(fbeat + 10 ,0.002, "Fbeat = " + floor(fbeat) + "Hz", 'Color', 'b');

    %Afficher un trait pour localiser la fréquence de battement
    line = plot([fbeat, fbeat], [0, psor], 'b--', 'LineWidth', 0.1); % Fréquence de battement en rouge

    hold off;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psor, fbeat, Distance] = Mesure(Y, f, c, k)

    % Récupérer le maximum du pic (X, Y):
    [psor,fbeat] = findpeaks(Y,f,'SortStr','descend', 'NPeaks', 1);

    % Afficher la distance à partir de la formule du chirp(Formule: fbeat = 2 x k x distance/c)
    Distance = (fbeat * c) / (2 * k);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y, f] = calculerFFT(ch1, tps)
    % Calcul de la FFT du signal
    N = length(ch1); % Longueur de la séquence du signal
    Y = fft(ch1); % Calcul de la FFT
    
    % Calcul du vecteur de fréquences
    T = tps(2) - tps(1); % Calcul de l'intervalle de temps (supposé constant)
    Fs = 1/T; % Fréquence d'échantillonnage
    f = Fs*(0:(N/2))/N; % Vecteur de fréquences
    
    % Normalisation de la FFT et extraction de la moitié (symétrie)
    Y = abs(Y/N); % Normalisation
    Y = Y(1:N/2+1); % Extraction de la moitié de la FFT pour symétrie
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function afficherControls(fig, Y1, f1)

    global c B Tm k permittiviteR; % Constantes
    global dist Fb fbeat acquisitionState; % Composant UI

% Ajouter un uipanel à gauche pour les paramètres
panelParam1 = uipanel('Parent', fig, 'Title', 'Paramètres', 'Position', [0.05, 0.6, 0.05, 0.3], 'BackgroundColor',[1 0.5 0.1]); % Position [x y largeur hauteur], normalisée à la taille de la figure

% Ajouter des éléments uicontrol dans le panel pour afficher les paramètres

% Excursion coefficient de pente:
uicontrol('Parent', panelParam1, 'Style', 'text', 'String', "εr = " + permittiviteR, 'Position', [10 180 60 20], 'BackgroundColor',[1 0.5 0.1], 'HorizontalAlignment','left');

% Vitesse de l'onde dans le milieu:
uicontrol('Parent', panelParam1, 'Style', 'text', 'String', "C = " + round(c/1e6)/1e2 + "e8", 'Position', [10 140 60 20], 'BackgroundColor',[1 0.5 0.1], 'HorizontalAlignment','left');

% Excursion en fréquence:
uicontrol('Parent', panelParam1, 'Style', 'text', 'String', "B = " + B/1e6 + "Mhz", 'Position', [10 100 60 20], 'BackgroundColor',[1 0.5 0.1], 'HorizontalAlignment','left');

% Temps de montée:
uicontrol('Parent', panelParam1, 'Style', 'text', 'String', "Tm = " + Tm * 1e3 + "ms", 'Position', [10 60 60 20], 'BackgroundColor',[1 0.5 0.1], 'HorizontalAlignment','left');

% coefficient de pente:
uicontrol('Parent', panelParam1, 'Style', 'text', 'String', "K = " + k/1e9 + "GHz²", 'Position', [10 20 60 20], 'BackgroundColor',[1 0.5 0.1], 'HorizontalAlignment','left');

% Ajouter un uipanel à gauche pour la formule
panelParam2 = uipanel('Parent', fig, 'Title', 'Formule', 'Position', [0.035, 0.4, 0.1, 0.1], 'BackgroundColor',[0.1 0.5 1]); % Position [x y largeur hauteur], normalisée à la taille de la figure

uicontrol('Parent', panelParam2, 'Style', 'text', 'String', 'Distance = (Fbeat x C) / 2K', 'Position', [10 20 150 20], 'BackgroundColor', [0.1 0.5 1], 'HorizontalAlignment', 'left');

% Ajouter un uipanel à gauche pour les mesures
panelParam3 = uipanel('Parent', fig, 'Title', 'Mesure', 'Position', [0.035, 0.2, 0.1, 0.1], 'BackgroundColor',[0.1 1 0.5]); % Position [x y largeur hauteur], normalisée à la taille de la figure

[psor, fbeat, Distance] = Mesure(Y1, f1, c, k)

dist = uicontrol('Parent', panelParam3, 'Style', 'text', 'String', 'Distance = ' + round(Distance), 'Position', [10 30 150 20], 'BackgroundColor', [0.1 1 0.5], 'HorizontalAlignment', 'left');

Fb = uicontrol('Parent', panelParam3, 'Style', 'text', 'String', 'Fbeat = ' + round(fbeat), 'Position', [10 10 150 20], 'BackgroundColor', [0.1 1 0.5], 'HorizontalAlignment', 'left');

%Afficher l'état de l'expérience:
acquisitionState = annotation('textbox', [0.5, 0.9, 0.1, 0.1], 'String', 'Etalonnage...');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [visaObj, NbPts, NbAcq] = initializeScope()
    persistent vObj  % Utilisez un nouveau nom pour éviter les confusions
    if isempty(vObj)
        % Initialiser l'objet VISA seulement si ce n'est pas déjà fait
        vObj = visa('keysight', 'USB0::10893::5990::MY58493285::INSTR');
        vObj.InputBufferSize = 1000000;
        vObj.Timeout = 100;
        vObj.ByteOrder = 'littleEndian';
        fopen(vObj);
        
        % Configuration de l'oscilloscope seulement lors de la première initialisation
        fprintf(vObj, ':TIMEBASE:MODE MAIN');
        fprintf(vObj, ':ACQUIRE:TYPE AVERAGE');
        fprintf(vObj, ':ACQUIRE:AVERAGES 2');
        %fprintf(vObj, ':ACQUIRE:COUNT 4');
        fprintf(vObj, ':WAV:POINTS:MODE MAX');
        fprintf(vObj, ':WAV:POINTS 500000');
        fprintf(vObj, ':STOP');
    end
    visaObj = vObj;  % Retourne l'objet VISA persistant
    
    % Paramètres de capture
    NbPts = 500000; % Nombre de points de capture
    NbAcq = 2;    % Nombre d'acquisitions
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fonction qui envoie l'ordre à l'oscilloscope de convertir en numérique les
%signaux analogiques de tous les channels

function digitize (visaObj)

    fprintf(visaObj,':DIGITIZE');
    operationComplete = str2double(query(visaObj,'*OPC?'));
    while ~operationComplete
        operationComplete = str2double(query(visaObj,'*OPC?'));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tps, chData] = captureChannel(visaObj, channel, NbPts)
    
    fprintf(visaObj, sprintf(':WAVEFORM:SOURCE %s', channel));
    fprintf(visaObj, ':WAVEFORM:FORMAT WORD');
    fprintf(visaObj, ':WAVEFORM:BYTEORDER LSBFirst');
    
    preambleBlock = query(visaObj, ':WAVEFORM:PREAMBLE?');
    fprintf(visaObj, ':WAV:DATA?');
    rawData = binblockread(visaObj, 'uint16'); fread(visaObj, 1);
    
    waveform = processWaveform(preambleBlock, rawData);
    tps = waveform.XData;
    chData = waveform.YData;
    
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function waveform = processWaveform(preambleBlock, rawData)
    preambleBlock = regexp(preambleBlock, ',', 'split');
    maxVal = 2^16;
    
    % Extraction des paramètres de la forme d'onde
    waveform.Points = str2double(preambleBlock{3});
    waveform.XIncrement = str2double(preambleBlock{5}); 
    waveform.XOrigin = str2double(preambleBlock{6});   
    waveform.XReference = str2double(preambleBlock{7});
    waveform.YIncrement = str2double(preambleBlock{8});
    waveform.YOrigin = str2double(preambleBlock{9});
    waveform.YReference = str2double(preambleBlock{10});
    
    waveform.XData = (0:(waveform.Points-1)) * waveform.XIncrement + waveform.XOrigin;
    waveform.YData = (double(rawData) - waveform.YReference) * waveform.YIncrement + waveform.YOrigin;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateWaveform(ch1, ch2, tps, fft1, fft2, f1, f2)

% Mettre à jour les données en X et Y du tracé
    global plot1 plot2 plot3 plot4;

    set(plot1, 'YData', ch1);
    set(plot1, 'XData', tps); 
    drawnow;
    set(plot2, 'YData', ch2); 
    set(plot2, 'XData', tps); 
    drawnow;
    set(plot3, 'YData', fft1);
    set(plot3, 'XData', f1); 
    drawnow;
    set(plot4, 'YData', fft2);
    set(plot4, 'XData', f2);
    drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateAcquisitionState(state)

% Mettre à jour les données en X et Y du tracé
    global acquisitionState;

    set(acquisitionState, 'String', state);

    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateControls(Distance, fbeat, psor)
    global dist Fb line dist2 Fb2;

    % Convertit 'Distance' en chaîne et concatène
    distanceStr = ['Distance = ', num2str(Distance), ' m'];

    % Met à jour le premier contrôle UI pour afficher la distance
    set(dist, 'String', distanceStr);
    drawnow;

    % Met à jour le deuxième contrôle UI pour afficher la distance
    set(dist2, 'String', distanceStr);
    drawnow;

    % Convertit 'fbeat' en chaîne et concatène
    FbStr = ['Fbeat = ', num2str(round(fbeat)), ' Hz'];

    % Met à jour le contrôle UI pour la valeur de Fbeat
    set(Fb, 'String', FbStr);
    drawnow;

    % Met à jour le contrôle UI pour la valeur de Fbeat
    set(Fb2, 'String', FbStr);
    drawnow;

    % Met à jour la ligne verticale sur le graphe
    set(line, 'XData', [fbeat, fbeat], 'YData', [0, psor]);
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function signalFLT = FiltreParasite(Y1)
    global Yref;
    signalFLT = Y1 - Yref;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateYref(visaObj, NbPts)
    global Yref;
    
    digitize(visaObj);
    [tps, Yref] = captureChannel(visaObj, 'CHAN1', NbPts);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction qui fournit les valeurs minimum d'un signal
%Rq: On utilise Findpeaks de matlab au lieu de la recherche par la méthode de
%la dérivée car le signal est bruité.

function mins = findMins(tps, Y2, hauteurSeuil, distanceSeuil)

%Inverser le signal
invertedSignal = -Y2;

% Utilisation de findpeaks avec un seuil de hauteur en %
[peaks, mins] = findpeaks(invertedSignal, tps, 'MinPeakProminence', hauteurSeuil, 'MinPeakDistance', distanceSeuil);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fonction qui fournit les valeurs maximum d'un signal
%Rq: On utilise Findpeaks de matlab au lieu de la recherche par la méthode de
%la dérivée car le signal est bruité.

function maxs = findMaxs(tps, Y2, hauteurSeuil, distanceSeuil)

% Utilisation de findpeaks avec un seuil de hauteur et de distance entre
% les pics
[peaks, maxs] = findpeaks(Y2, tps, 'MinPeakProminence', hauteurSeuil, 'MinPeakDistance', distanceSeuil);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fonction qui renvoie une matrice contenant des couples de minimum et
% maximum pour un signal
%Rq: Les parametres optimaux hauteurSeuil et distanceSeuil de findMin et findMax ont été obtenu en
%expérimentant
function [MinMax, MinMaxSize] = findCorrespondingMinMax(Y2, tps)

% Récupérer les minimums et maximums
    maxs = findMaxs(tps, Y2, 0.5, 0.18);
    mins = findMins(tps, Y2, 0.5, 0.18);

% Obtenir le couple de valeurs
MinMax = findClosestPairs(mins, maxs);
MinMaxSize = size(MinMax, 1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fonction qui trouve les paires de valeurs les plus proches entre deux vecteurs.
%
% Arguments:
%   vec1 - Premier vecteur d'entrée.
%   vec2 - Deuxième vecteur d'entrée.
%
% Retourne:
%   closestPairs - Une matrice N x 2 contenant les paires de valeurs les plus proches,
%                  avec N étant le nombre de couples. Chaque ligne [a, b]
%                  représente la valeur la plus proche trouvée : a de vec1
%                  et b de vec2, ou vice versa. b>a!

function pairs = findClosestPairs(vec1, vec2)
    % Combine et trie les vecteurs
    vec3 = sort([vec1, vec2]);
    
    % Initialiser un tableau pour stocker les paires temporaires
    tempPairs = [];
    
    % Parcourir chaque élément de vec3 et trouver les paires
    for i = 1:length(vec3)
        value = vec3(i);
        if ismember(value, vec1)
            % Si la valeur est dans vec1, trouver la valeur supérieure la plus proche dans vec2
            possibleVals = vec2(vec2 > value);
            if ~isempty(possibleVals)
                closestVal = min(possibleVals);
                tempPairs = [tempPairs; value, closestVal];
            end
        else
            % Si la valeur est dans vec2, trouver la valeur supérieure la plus proche dans vec1
            possibleVals = vec1(vec1 > value);
            if ~isempty(possibleVals)
                closestVal = min(possibleVals);
                tempPairs = [tempPairs; value, closestVal];
            end
        end
    end
    
    % Retirer les paires répétées et trier
    pairs = unique(tempPairs, 'rows');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ch1_temp, ch2_temp, tps_temp, Y1Moyen, Y2, f1_temp, f2, psor, fbeat, Distance] = MoyennageChirps(ch1, ch2, tps)
    % Fonction qui redimensionne les signaux tps, CHANNEL1 et CHANNEL2, calcule
    % pour chaque plage la fréquence de battement et retourne la moyenne de la
    % distance estimée
    % La parité de l'indice permet de sélectionner la pente du chirp (montant ou descendant) que l'on veut capturer

    global c k;

    % Récupération du couple (a, b) pour la plage d'acquisition, b>a
    [X1X2, MinMaxSize] = findCorrespondingMinMax(ch2, tps);
    Y1Cumul = [];

    tailleCumul = 0;

    for indice = 1:MinMaxSize
        X1 = X1X2(indice, 1);
        X2 = X1X2(indice, 2);

        % Indexation logique
        validIndices = tps > X1 & tps < X2; % validIndices est un vecteur logique qui contient true pour les éléments de tps qui se trouvent dans l'intervalle (Min, Max)

        % Filtrer les vecteurs
        tps_temp = tps(validIndices);
        ch1_temp = ch1(validIndices);
        ch2_temp = ch2(validIndices);

        % Calculer les FFT
        [Y1, f1] = calculerFFT(ch1_temp, tps_temp); 
        [Y2, f2] = calculerFFT(ch2_temp, tps_temp);

       % Initialiser Y1Cumul avec la première FFT si elle est vide
        if isempty(Y1Cumul)
            Y1Cumul = Y1;
        else
            % Prendre le minimum des tailles pour l'addition
            minLength = min(length(Y1Cumul), length(Y1));
            Y1Cumul(1:minLength) = Y1Cumul(1:minLength) + Y1(1:minLength);
        end
        tailleCumul = tailleCumul + 1;
    end

    % Calculer la moyenne des FFT sur tous les segments
    Y1Moyen = Y1Cumul / tailleCumul;

    % Calcul de la fréquence d'échantillonnage
    Fs = 1 / mean(diff(tps_temp)); % Supposant tps_temp régulièrement espacé

    % Recalculer le vecteur de fréquences ajusté
    f1_temp = ajusterVecteurFrequences(Fs, length(Y1Moyen));

    % Appliquer la fonction Mesure sur le signal moyen Y1Moyen
    [psor, fbeat, Distance] = Mesure(Y1Moyen, f1_temp, c, k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f1_temp = ajusterVecteurFrequences(Fs, N)
    % Génère un vecteur de fréquences pour la moitié positive du spectre, incluant la fréquence Nyquist
    f1_temp = linspace(0, Fs/2, N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

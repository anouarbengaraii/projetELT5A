%% Simulation multiplication de deux chirps de type montant

%% Paramètres de simulation

%% Directes:
Tm = 5e-3;   % Temps de montée du chirp
B = 200e6;   % Excursion en fréquence du chirp
f0 = 0;   % Fréquence initiale du chirp1 (en général 30Mhz)
distance=200  ;% Distance aller radar-cible
c = 3e8;      % Célérité de la lumière
N = 1e5;     % Nombre de points pour la simulation (lié à la fréquence d'échantillonnage)
n = 4;  % Ordre du filtre Par exemple, un filtre d'ordre 4


%% Indirectes
retard = 2*distance/c;  % Retard entre le chirp émis et reçu
k = B / Tm;             % Coefficient de pente
fbeat=2*k*distance/c; % calcul de la frequence de battement fbeat à partir des paramètres

%% 1) Synthèse du premier chirp
t = linspace(0, Tm, N); % Vecteur de temps t
chirp1 = cos(2 * pi * (f0 * t + 0.5 * k * t .^ 2));  % formule du chirp1

%% 2) Synthèse du deuxième chirp (effet doppler négligé)
t2 = t - retard; % Correction du temps pour le deuxième chirp
chirp2 = cos(2 * pi * (f0 * t2 + 0.5 * k * t2 .^ 2));  % formule du chirp2

%% 3) Multiplication des deux chirps
Smult = chirp1 .* chirp2; % Multiplication du chirp1 et du chirp2
Fs = 1 / (t(2) - t(1)); % Fréquence d'échantillonnage calculée à partir du vecteur de temps

%conception du filtre passe-bande
fc_low=0.95*fbeat;
fc_high=1.05*fbeat;

Smult_filtered=filtrePasseBandeIdeal(Smult,Fs,fc_low,fc_high);

%% 4) FFT de Smult et préparation pour le tracé
fft_result = fft(Smult, N);
fft_result_shifted = fftshift(fft_result); % Décalage du spectre FFT pour centrage

fft_result2 = fft(Smult_filtered, N);
fft_result_shifted2 = fftshift(fft_result2);

f = linspace(-Fs/2, Fs/2, N); % Vecteur de fréquences centré autour de 0

%% 5) Synthèse en fréquence des deux chirps:
F1 = f0 + k * t;
F2 = f0 + k * t2;

%% Affichage des graphiques
fig=figure(Name="Simulation de la multiplication de deux chirps");
fig.WindowState="maximized";

%% Tracé des deux Chirp fréquence = f(temps)
subplot(2, 2, 1);
hold on;
plot(t, F1, 'b-');
plot(t, F2, 'r-');
xlabel('Temps (s)');
ylabel('Fréquence (Hz)');
title('Signal émis et reçu');
ylim([0 2e8])
legend("Chirp émis", "Chirp reçu")

%% Tracé dans le domaine temporel de Smult
subplot(2, 2, 2);
plot(t, Smult);
title('Signal de battement avant filtrage');
xlim([2e-4, 2.3e-4]); %à adapter en fonction de la plage voulue

subplot(2, 2, 3);
plot(t, Smult_filtered);
xlabel('Temps (s)');
ylabel('Amplitude');
title('Signal de battement après filtrage');
xlim([1e-3, 2e-3]); %à adapter en fonction de la plage voulue

%% Tracé dans le domaine fréquentiel de Smult

% Sélection du quatrième sous-graphique pour le dessin
subplot(2, 2, 4);

% Tracé du premier résultat de FFT
hold on;
h1=plot(f, abs(fft_result_shifted2),'DisplayName','Après filtrage');
 % Permet de superposer le prochain graphique sans effacer le précédent

% Tracé du deuxième résultat de FFT
h2=plot(f, abs(fft_result_shifted),'DisplayName','Avant filtrage');
legend([h1,h2]);

% Configuration des étiquettes et du titre
xlabel('Fréquence (Hz)');
ylabel('Amplitude');
title('Fréquence de battement');
% Ajustement de la limite de l'axe des x pour une meilleure visualisation
xlim([0.7*fbeat, 1.3*fbeat]);

% Activation de la grille pour une meilleure lisibilité
grid on;

hold on;
%tracé de la fréquence de battement et la distance reels
[psor, pics] = findpeaks(abs(fft_result_shifted), f, 'SortStr', 'descend', 'Threshold', 1e-3);%recherche des pics de la FFT
fbeat_real = pics(2);  %la frequence de battement correspond au deuxieme pic
distance_reel=fbeat_real*c/(2*k); %fbeat=2*k*distance/c;

plot([fbeat_real, fbeat_real], [0, max(abs(fft_result_shifted))], 'y--', 'LineWidth', 2);
point_y = max(abs(fft_result_shifted)) / 3;
text(Fs/8, point_y, sprintf('Distance = %.2f m', distance_reel), 'HorizontalAlignment', 'center', 'Color', 'green');



hold off; % Désactive le mode superposition pour les futurs graphiques

afficherControls(fig, distance, distance_reel,fbeat,fbeat_real)

%fonction qui realise un filtre passe bas ideal
function y = filtrePasseBandeIdeal(x, Fs, Fpass1, Fpass2)
    % x : signal d'entrée
    % Fs : fréquence d'échantillonnage du signal
    % Fpass1 : fréquence de coupure basse du filtre passe-bande
    % Fpass2 : fréquence de coupure haute du filtre passe-bande
    
    N = length(x); % Nombre de points dans le signal
    f = (-N/2:N/2-1)*(Fs/N); % Vecteur de fréquences
    
    % Création du filtre passe-bande idéal
    H = double((f >= Fpass1) & (f <= Fpass2)); % Fonction rectangulaire
    H = fftshift(H); % Décalage du filtre pour l'aligner avec la FFT du signal
    
    % Application du filtre
    X = fft(x); % Transformée de Fourier du signal
    Y = H .* X; % Filtrage dans le domaine fréquentiel
    y = ifft(Y, 'symmetric'); % Transformée de Fourier inverse pour revenir au domaine temporel
    
    % y est le signal filtré
end

%fonction qui affiche les paramètres
% function that displays the parameters
function afficherControls(fig, distance, distance_reel, fbeat, fbeat_real)
    % Ajouter un uipanel à gauche pour les paramètres théoriques
    panelParam1 = uipanel('Parent', fig, 'Title', 'Paramètres', 'Position', [0, 0.55, 0.1, 0.26], 'BackgroundColor',[0.9 0.9 0.7]);

    % Ajouter des éléments uicontrol dans le panel pour afficher les paramètres théoriques

    % Use handles to store references to the uicontrols
    uicontrol('Parent', panelParam1, 'Style', 'text', 'String', ['Distance théorique =  ' num2str(distance)], 'Position', [0 50 300 20], 'BackgroundColor',[0.9 0.9 0.9], 'HorizontalAlignment','left', 'FontSize', 7,'FontWeight', 'bold');
    uicontrol('Parent', panelParam1, 'Style', 'text', 'String', ['fbeat théorique = ' num2str(fbeat)], 'Position', [0 10 300 20], 'BackgroundColor',[0.9 0.9 0.9], 'HorizontalAlignment','left', 'FontSize', 7,'FontWeight', 'bold');
    
    % Ajouter un uipanel à gauche pour les paramètres obtenus
    panelParam2 = uipanel('Parent', fig, 'Title', 'Résultats', 'Position', [0, 0.2, 0.1, 0.3], 'BackgroundColor',[0.9 0.9 0.5]);

    % Ajouter des éléments uicontrol dans le panel pour afficher les paramètres obtenus

    uicontrol('Parent', panelParam2, 'Style', 'text', 'String', ['Distance obtenue = ' num2str(distance_reel)], 'Position', [0 120 300 20], 'BackgroundColor',[0.9 0.9 0.9], 'HorizontalAlignment','left', 'FontSize', 7,'FontWeight', 'bold');
    uicontrol('Parent', panelParam2, 'Style', 'text', 'String', ['fbeat obtenue =  ' num2str(fbeat_real)], 'Position', [0 80 300 20], 'BackgroundColor',[0.9 0.9 0.9], 'HorizontalAlignment','left', 'FontSize', 7,'FontWeight', 'bold');

end

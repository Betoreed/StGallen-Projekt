clear all; clc; close all;

%% 1. Vorgaben und Initialisierung laut Projektauftrag
% Punkte im Raum [x, y, z] in Metern
P1 = [0.0, 0.3, 0.2];
P2 = [0.3, 0.0, 0.1];
P3 = [0.3, 0.3, 0.25];
P4 = [0.4, 0.1, 0.15];

% Definiere die Gesamtbahn: P1 -> P2 -> P3 -> P4 -> P1
P_seq = [P1; P2; P3; P4; P1];
num_segments = size(P_seq, 1) - 1;

v_soll = 1;      % Bahngeschwindigkeit [m/s]
a_soll = 10;     % Bahnbeschleunigung [m/s^2]
ts = 1e-3;       % Abtastzeit 1ms

% Initialisierung der Sammel-Arrays für die gesamte Trajektorie
BahnkoordinatenGesamt = [];
t_gesamt = [];
v_gesamt = [];
a_gesamt = [];

t_aktuell = 0; % Mitlaufende Zeit

%% 2. Berechnung der Bahnpunkte (Trapezförmiges Profil)
for k = 1:num_segments
    StartP = P_seq(k, :);
    EndP   = P_seq(k+1, :);
    
    Bahn = EndP - StartP;
    Strecke = norm(Bahn);
    Richtung = Bahn ./ Strecke;
    
    % Berechnung der Zeiten des trapezförmigen Geschwindigkeitsverlaufs
    ta = v_soll / a_soll;                   % Theoretische Beschleunigungszeit
    tv = (Strecke - v_soll * ta) / v_soll;  % Theoretische Zeit für Konstantfahrt
    
    % Gesamtzeit auf ein ganzzahliges Vielfaches von ts aufrunden
    tges = ceil((2 * ta + tv) / ts) * ts;
    
    % Leicht angepasste Beschleunigung, um genau im ts-Raster zu landen
    % (Verhindert Positionsfehler am Ende des Segments)
    a_neu = (v_soll^2) / (v_soll * tges - Strecke);
    ta_neu = v_soll / a_neu;
    
    % Zeitvektor für dieses Segment generieren
    t_seg = (ts : ts : tges)';
    
    % Profil-Arrays für das aktuelle Segment
    s_seg = zeros(length(t_seg), 1);
    v_seg = zeros(length(t_seg), 1);
    a_seg = zeros(length(t_seg), 1);
    
    % Weg-, Geschwindigkeits- und Beschleunigungsgesetz (Trapezprofil)
    for i = 1:length(t_seg)
        t = t_seg(i);
        if t <= ta_neu % Beschleunigungsphase
            a_seg(i) = a_neu;
            v_seg(i) = a_neu * t;
            s_seg(i) = 0.5 * a_neu * t^2;
        elseif t <= (tges - ta_neu) % Konstantfahrtphase
            a_seg(i) = 0;
            v_seg(i) = v_soll;
            s_seg(i) = 0.5 * a_neu * ta_neu^2 + v_soll * (t - ta_neu);
        else % Verzögerungsphase
            dt = t - (tges - ta_neu);
            a_seg(i) = -a_neu;
            v_seg(i) = v_soll - a_neu * dt;
            s_seg(i) = (0.5 * a_neu * ta_neu^2) + (v_soll * (tges - 2*ta_neu)) + (v_soll * dt - 0.5 * a_neu * dt^2);
        end
    end
    
    % Koordinaten im Raum berechnen (Lineare Interpolation)
    Bahnkoordinaten = StartP + s_seg * Richtung;
    
    % Anfügen an die Gesamttrajektorie
    BahnkoordinatenGesamt = [BahnkoordinatenGesamt; Bahnkoordinaten];
    t_gesamt = [t_gesamt; t_seg + t_aktuell];
    v_gesamt = [v_gesamt; v_seg];
    a_gesamt = [a_gesamt; a_seg];
    
    t_aktuell = t_aktuell + tges;
end

%% 3. Analytische Inverse Kinematik (Optimiert für MATLAB 2025b, ohne 'syms')
l1 = 0.200; % Armlänge 1 in m
l2 = 0.250; % Armlänge 2 in m

x = BahnkoordinatenGesamt(:, 1);
y = BahnkoordinatenGesamt(:, 2);

% Kosinussatz für phi2
cos_phi2 = (x.^2 + y.^2 - l1^2 - l2^2) ./ (2 * l1 * l2);
% Begrenzung auf [-1, 1] zur Vermeidung von numerischen Fehlern
cos_phi2 = max(min(cos_phi2, 1), -1);

% Winkel berechnen
vecPhi2 = acos(cos_phi2); 
vecPhi1 = atan2(y, x) - atan2(l2 * sin(vecPhi2), l1 + l2 * cos_phi2);

% Verschiebung des Koordinatensystems (wie in deinem Ursprungsskript)
vecPhi1 = vecPhi1 - pi/4;

% Export für SimulationX
tempVec1 = [t_gesamt, vecPhi1]; 
tempVec2 = [t_gesamt, vecPhi2]; 
dlmwrite('Winkel_Phi1_lang.txt', tempVec1, 'precision', 10);
dlmwrite('Winkel_Phi2_lang.txt', tempVec2, 'precision', 10);

%% 4. Numerische Ableitung (Geschwindigkeit)
% Da wir ein festes Raster ts haben, können wir diff nutzen
vecPhi1_punkt = [0; diff(vecPhi1) / ts];
vecPhi2_punkt = [0; diff(vecPhi2) / ts];

%% 5. Plots
figure('Name', 'Bahn im Raum');
plot3(BahnkoordinatenGesamt(:,1), BahnkoordinatenGesamt(:,2), BahnkoordinatenGesamt(:,3), 'b-', 'LineWidth', 1.5)
hold on
scatter3(P_seq(:,1), P_seq(:,2), P_seq(:,3), 50, 'k', 'filled')
text(P1(1)+0.01, P1(2), P1(3), ' P1 (Start/Ende)', 'FontWeight', 'bold')
text(P2(1)+0.01, P2(2), P2(3), ' P2')
text(P3(1)+0.01, P3(2), P3(3), ' P3')
text(P4(1)+0.01, P4(2), P4(3), ' P4')
grid on
title('Sollbahn (P1 -> P2 -> P3 -> P4 -> P1)')
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]')

figure('Name', 'Dynamik der Bahn');
subplot(2,1,1)
plot(t_gesamt, v_gesamt, 'LineWidth', 1.5); grid on;
title('Bahngeschwindigkeit')
ylabel('Geschwindigkeit [m/s]')

subplot(2,1,2)
plot(t_gesamt, a_gesamt, 'LineWidth', 1.5); grid on;
title('Bahnbeschleunigung')
xlabel('Zeit [s]'); ylabel('Beschleunigung [m/s^2]')

figure('Name', 'Achswinkel');
subplot(2,1,1)
plot(t_gesamt, vecPhi1*180/pi, t_gesamt, vecPhi2*180/pi, 'LineWidth', 1.5); grid on;
title('Winkel \phi_1 und \phi_2')
ylabel('Winkel [°]')
legend({'\Phi_1', '\Phi_2'}, 'Location', 'best');

subplot(2,1,2)
plot(t_gesamt, vecPhi1_punkt, t_gesamt, vecPhi2_punkt, 'LineWidth', 1.5); grid on;
title('Winkelgeschwindigkeit \omega_1 und \omega_2')
xlabel('Zeit [s]'); ylabel('Geschwindigkeit [rad/s]')
legend({'\omega_1', '\omega_2'}, 'Location', 'best');

%% Überprüfung per Handrechnung

% Test und Daten für die Handrechnung
% Winkelgeschwindigkeiten mittels zentraler Differenz
vecPhi1_punkt = gradient(vecPhi1) / ts; 
vecPhi2_punkt = gradient(vecPhi2) / ts; 

testpos = 900; % Entspricht t = 0.9s bei ts = 1ms

% Extraktion der exakten Werte an der Testposition für den Bericht
phi1_test = vecPhi1(testpos);
phi2_test = vecPhi2(testpos);
omega1_test = vecPhi1_punkt(testpos);
omega2_test = vecPhi2_punkt(testpos);

% Deine externe Funktion zur Kontrolle (alternativ hier direkt die Formeln eintragen)
[v_x, v_y, x_tcp, y_tcp] = IdealWert(vecPhi1, vecPhi2, vecPhi1_punkt, vecPhi2_punkt, testpos);

% Visualisierung
figure('Name', 'Winkelgeschwindigkeiten')
plot(t, vecPhi1_punkt, 'LineWidth', 1.5);
hold on; grid on;
plot(t, vecPhi2_punkt, 'LineWidth', 1.5);
title('Winkelgeschwindigkeit über der Zeit \omega_1 und \omega_2')
xlabel('Zeit [s]')
ylabel('Winkelgeschwindigkeit [rad/s]')
legend({'\omega_1','\omega_2'}, 'Location', 'best');

%% LOKALE FUNKTIONEN

function [v_x, v_y, x_tcp, y_tcp] = IdealWert(vecPhi1, vecPhi2, vecPhi1_punkt, vecPhi2_punkt, testpos)
    % 1. Konstanten laut Datenblatt (Armlängen in m)
    l1 = 0.200;
    l2 = 0.250;

    % 2. Werte für den gewählten Zeitpunkt (testpos) extrahieren
    % WICHTIG: Die -pi/4 Verschiebung aus der inversen Kinematik muss 
    % hier für die globale Weltkoordinate wieder rückgängig gemacht werden (+pi/4)!
    phi1 = vecPhi1(testpos) + pi/4; 
    phi2 = vecPhi2(testpos);
    
    omega1 = vecPhi1_punkt(testpos);
    omega2 = vecPhi2_punkt(testpos);

    % 3. Vorwärtskinematik (Position berechnen)
    x_tcp = l1 * cos(phi1) + l2 * cos(phi1 + phi2);
    y_tcp = l1 * sin(phi1) + l2 * sin(phi1 + phi2);

    % 4. Differenzielle Kinematik (Geschwindigkeiten über Jacobi-Matrix berechnen)
    % J11 = -l1*sin(phi1) - l2*sin(phi1+phi2)
    % J12 = -l2*sin(phi1+phi2)
    % J21 = l1*cos(phi1) + l2*cos(phi1+phi2)
    % J22 = l2*cos(phi1+phi2)
    
    v_x = (-l1 * sin(phi1) - l2 * sin(phi1 + phi2)) * omega1 + (-l2 * sin(phi1 + phi2)) * omega2;
    v_y = (l1 * cos(phi1) + l2 * cos(phi1 + phi2)) * omega1 + (l2 * cos(phi1 + phi2)) * omega2;
    
    % Konsolenausgabe für die Dokumentation / Handrechnung
    fprintf('\n--- ERGEBNISSE DER HANDRECHNUNG BEI INDEX %d ---\n', testpos);
    fprintf('Eingangswerte:\n');
    fprintf('  Phi_1:    %.4f rad\n', phi1);
    fprintf('  Phi_2:    %.4f rad\n', phi2);
    fprintf('  Omega_1:  %.4f rad/s\n', omega1);
    fprintf('  Omega_2:  %.4f rad/s\n', omega2);
    fprintf('Berechnete kartesische Werte:\n');
    fprintf('  X_tcp:    %.4f m\n', x_tcp);
    fprintf('  Y_tcp:    %.4f m\n', y_tcp);
    fprintf('  v_x:      %.4f m/s\n', v_x);
    fprintf('  v_y:      %.4f m/s\n\n', v_y);
end
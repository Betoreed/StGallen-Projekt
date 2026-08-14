clear all; clc; close all;

%% 1. Vorgaben und Initialisierung laut Projektauftrag
P1 = [0.0, 0.3, 0.2];
P2 = [0.3, 0.0, 0.1];
P3 = [-0.3, 0.3, 0.25];
P4 = [-0.4, -0.1, 0.15];

% Definiere die Gesamtbahn: P1 -> P2 -> P3 -> P4 -> P1
P_seq = [P1; P2; P3; P4; P1];
num_segments = size(P_seq, 1) - 1;

v_soll = 1;      % Bahngeschwindigkeit [m/s]
a_soll = 10;     % Bahnbeschleunigung [m/s^2]
ts = 1e-3;       % Abtastzeit 1ms

BahnkoordinatenGesamt = [];
t_gesamt = [];
v_gesamt = [];
a_gesamt = [];
t_aktuell = 0;

%% 2. Berechnung der Bahnpunkte (Trapezförmiges Profil)
for k = 1:num_segments
    StartP = P_seq(k, :);
    EndP   = P_seq(k+1, :);
    
    Bahn = EndP - StartP;
    Strecke = norm(Bahn);
    Richtung = Bahn ./ Strecke;
    
    ta = v_soll / a_soll;                   
    tv = (Strecke - v_soll * ta) / v_soll;  
    
    tges = ceil((2 * ta + tv) / ts) * ts;
    
    a_neu = (v_soll^2) / (v_soll * tges - Strecke);
    ta_neu = v_soll / a_neu;
    
    t_seg = (ts : ts : tges)';
    
    s_seg = zeros(length(t_seg), 1);
    v_seg = zeros(length(t_seg), 1);
    a_seg = zeros(length(t_seg), 1);
    
    for i = 1:length(t_seg)
        t = t_seg(i);
        if t <= ta_neu
            a_seg(i) = a_neu;
            v_seg(i) = a_neu * t;
            s_seg(i) = 0.5 * a_neu * t^2;
        elseif t <= (tges - ta_neu)
            a_seg(i) = 0;
            v_seg(i) = v_soll;
            s_seg(i) = 0.5 * a_neu * ta_neu^2 + v_soll * (t - ta_neu);
        else
            dt = t - (tges - ta_neu);
            a_seg(i) = -a_neu;
            v_seg(i) = v_soll - a_neu * dt;
            s_seg(i) = (0.5 * a_neu * ta_neu^2) + (v_soll * (tges - 2*ta_neu)) + (v_soll * dt - 0.5 * a_neu * dt^2);
        end
    end
    
    Bahnkoordinaten = StartP + s_seg * Richtung;
    
    BahnkoordinatenGesamt = [BahnkoordinatenGesamt; Bahnkoordinaten];
    t_gesamt = [t_gesamt; t_seg + t_aktuell];
    v_gesamt = [v_gesamt; v_seg];
    a_gesamt = [a_gesamt; a_seg];
    
    t_aktuell = t_aktuell + tges;
end

%% 3. Analytische Inverse Kinematik
l1 = 0.200; % Armlänge 1 in m
l2 = 0.250; % Armlänge 2 in m
x = BahnkoordinatenGesamt(:, 1);
y = BahnkoordinatenGesamt(:, 2);

cos_phi2 = (x.^2 + y.^2 - l1^2 - l2^2) ./ (2 * l1 * l2);
cos_phi2 = max(min(cos_phi2, 1), -1);

vecPhi2 = acos(cos_phi2); 
vecPhi1 = atan2(y, x) - atan2(l2 * sin(vecPhi2), l1 + l2 * cos_phi2);

% Winkelsprünge entfernen
vecPhi1 = unwrap(vecPhi1);
vecPhi2 = unwrap(vecPhi2);

% Einheitlicher Offset (pi/4)
phi_offset = pi/4; 
vecPhi1 = vecPhi1 - phi_offset;

% Export für SimulationX
tempVec1 = [t_gesamt, vecPhi1]; 
tempVec2 = [t_gesamt, vecPhi2]; 
dlmwrite('Phi1.txt', tempVec1, 'precision', 10);
dlmwrite('Phi2.txt', tempVec2, 'precision', 10);

%% 4. Numerische Ableitung (Geschwindigkeit)
vecPhi1_punkt = gradient(vecPhi1) / ts; 
vecPhi2_punkt = gradient(vecPhi2) / ts; 

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

%% 6. Überprüfung per Handrechnung
testpos = 1000; % anpassen für verschiedene Testpositionen

[v_x, v_y, x_tcp, y_tcp] = IdealWert(vecPhi1, vecPhi2, vecPhi1_punkt, vecPhi2_punkt, testpos, phi_offset);

%% LOKALE FUNKTIONEN
function [v_x, v_y, x_tcp, y_tcp] = IdealWert(vecPhi1, vecPhi2, vecPhi1_punkt, vecPhi2_punkt, testpos, phi_offset)
    l1 = 0.200;
    l2 = 0.250;
    
    % Winkel im Weltkoordinatensystem wiederherstellen (+ phi_offset)
    phi1 = vecPhi1(testpos) + phi_offset; 
    phi2 = vecPhi2(testpos);
    
    omega1 = vecPhi1_punkt(testpos);
    omega2 = vecPhi2_punkt(testpos);
    
    % Vorwärtskinematik
    x_tcp = l1 * cos(phi1) + l2 * cos(phi1 + phi2);
    y_tcp = l1 * sin(phi1) + l2 * sin(phi1 + phi2);
    
    % Jacobi-Matrix Geschwindigkeiten
    v_x = (-l1 * sin(phi1) - l2 * sin(phi1 + phi2)) * omega1 + (-l2 * sin(phi1 + phi2)) * omega2;
    v_y = (l1 * cos(phi1) + l2 * cos(phi1 + phi2)) * omega1 + (l2 * cos(phi1 + phi2)) * omega2;
    
    % Konsolenausgabe
    fprintf('\n--- ERGEBNISSE DER HANDRECHNUNG BEI t = %.3f s (INDEX %d) ---\n', testpos*1e-3, testpos);
    fprintf('Eingangswerte:\n');
    fprintf('  Phi_1:    %.5f rad (%.2f°)\n', phi1, phi1*180/pi);
    fprintf('  Phi_2:    %.5f rad (%.2f°)\n', phi2, phi2*180/pi);
    fprintf('  Omega_1:  %.5f rad/s\n', omega1);
    fprintf('  Omega_2:  %.5f rad/s\n', omega2);
    fprintf('Berechnete kartesische Werte:\n');
    fprintf('  X_tcp:    %.5f m  \n', x_tcp);
    fprintf('  Y_tcp:    %.5f m  \n', y_tcp);
    fprintf('  v_x:      %.5f m/s \n', v_x);
    fprintf('  v_y:      %.5f m/s \n\n', v_y);
end
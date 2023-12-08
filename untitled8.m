clc;
clear all;

K_values = [0, 10, 100]; % Valores de K a evaluar
thresholds = [3, 5, 10, 20]; % Umbrales de potencia en dB

LCR_values = zeros(length(thresholds), length(K_values)); % Inicializar matriz para almacenar LCR
AFD_values = zeros(length(thresholds), length(K_values)); % Inicializar matriz para almacenar AFD

figure; % Crear una sola figura para los gráficos de Pnorm_db y su CDF

% Gráfico de Pnorm_db y CDF para diferentes valores de K
subplot(2, 1, 1);
hold on; % Mantener el gráfico actual para agregar más curvas

for i = 1:length(K_values)
    K = K_values(i);
    
    A = sqrt(2 * K); % Calcular A basado en el valor de K
    
    c = A + randn(1000, 1) + 1i * randn(1000, 1); % Generar señal con componente A
    r = abs(c);
    
    Pveces = r.^2;
    Pnorm_veces = Pveces / mean(Pveces);
    Pnorm_db = 10 * log10(Pnorm_veces);
    
    plot(Pnorm_db, 'LineWidth', 1.5, 'DisplayName', ['K = ', num2str(K)]); % Graficar Pnorm_db con distintos colores para cada K
end

hold off; % Terminar de agregar curvas al gráfico
title('Gráfico de dB vs muestras');
legend('Location', 'best');
xlabel('Muestras');
ylabel('Valor en dB');

% Gráfico de CDF para diferentes valores de K
subplot(2, 1, 2);
hold on; % Mantener el gráfico actual para agregar más curvas

for i = 1:length(K_values)
    K = K_values(i);
    
    A = sqrt(2 * K); % Calcular A basado en el valor de K
    
    c = A + randn(1000, 1) + 1i * randn(1000, 1); % Generar señal con componente A
    r = abs(c);
    
    Pveces = r.^2;
    Pnorm_veces = Pveces / mean(Pveces);
    Pnorm_db = 10 * log10(Pnorm_veces);
    
    % Calcular la CDF de Pnorm_db y graficarla
    [f, x] = ecdf(Pnorm_db);
    plot(x, f, 'LineWidth', 1.5, 'DisplayName', ['K = ', num2str(K)]); % Graficar la CDF
end

hold off; % Terminar de agregar curvas al gráfico
title('CDF Empírico');
legend('Location', 'best');
xlabel('Valor en dB');
ylabel('Probabilidad acumulada');

% Gráfico de LCR y AFD vs. factor K
figure; % Crear una nueva figura para el gráfico de LCR y AFD vs. factor K

for i = 1:length(K_values)
    K = K_values(i);
    
    A = sqrt(2 * K); % Calcular A basado en el valor de K
    
    c = A + randn(1000, 1) + 1i * randn(1000, 1); % Generar señal con componente A
    r = abs(c);
    
    Pveces = r.^2;
    Pnorm_veces = Pveces / mean(Pveces);
    Pnorm_db = 10 * log10(Pnorm_veces);
    
    for j = 1:length(thresholds)
        threshold = thresholds(j);
        
        % Encontrar cruces por cero para el umbral de potencia actual
        crossings = nnz(diff(Pnorm_db > threshold)); % Contar cruces por cero
        
        % Calcular el LCR para el umbral de potencia actual y el factor K actual
        LCR_values(j, i) = crossings / (length(Pnorm_db) - 1); % Calcular LCR
        
        % Calcular el AFD para el umbral de potencia actual y el factor K actual
        ms_0 = 0;
        for x = 1:length(Pnorm_db)
            if Pnorm_db(x) < threshold
                ms_0 = ms_0 + 1;
            end
        end
        AFD_values(j, i) = ms_0 / length(Pnorm_db) * 100; % Calcular AFD en porcentaje
    end
end

% Graficar LCR vs. factor K
subplot(2, 1, 1);
plot(K_values, LCR_values, 'o-', 'LineWidth', 1.5);
title('LCR vs. Factor K');
legend('3 dB', '5 dB', '10 dB', '20 dB', 'Location', 'best');
xlabel('Factor K');
ylabel('LCR');

% Graficar AFD vs. factor K
subplot(2, 1, 2);
plot(K_values, AFD_values, 'o-', 'LineWidth', 1.5);
title('AFD vs. Factor K');
legend('3 dB', '5 dB', '10 dB', '20 dB', 'Location', 'best');
xlabel('Factor K');
ylabel('AFD (%)');

function [P, M, N, navSolutions] = initNavSolutions(channelNr, settings, measNrSum)
    P = channelNr;              % Número de PRNs usados
    M = settings.numberOfChannels; % Número de canais GNSS
    N = measNrSum;              % Número de épocas
    
    % Estrutura para agrupar as matrizes das soluções de navegação
    navSolutions.PRNs          = zeros(P, N);
    navSolutions.correctedP   = zeros(P, N);
    navSolutions.DOP          = zeros(P + 1, N);
    navSolutions.el           = zeros(M, N);
    navSolutions.az           = zeros(M, N);
    navSolutions.transmitTime = zeros(M, N);
    navSolutions.satClkCorr   = zeros(M, N);
    navSolutions.rawP         = zeros(M, N);
    navSolutions.X            = zeros(1, N);
    navSolutions.Y            = zeros(1, N);
    navSolutions.Z            = zeros(1, N);
    navSolutions.dt           = zeros(1, N);
    navSolutions.currMeasSample = zeros(1, N);
    navSolutions.latitude     = zeros(1, N);
    navSolutions.longitude    = zeros(1, N);
    navSolutions.height       = zeros(1, N);
    navSolutions.utmZone      = 0;
    navSolutions.E            = zeros(1, N);
    navSolutions.N            = zeros(1, N);
    navSolutions.U            = zeros(1, N);
end
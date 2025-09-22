function y = wrap_to_pm(x, halfspan)
% mapeia x para o intervalo [-halfspan, +halfspan]
y = mod(x + halfspan, 2*halfspan) - halfspan;
end
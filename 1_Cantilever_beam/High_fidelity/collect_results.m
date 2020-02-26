close all
clear all


nsim = 3000;

for isim = 1:nsim
    disp(isim)    
    % Generate the new directory
    dir_name = ['out' num2str(isim)];
    Uf(:,isim)  = load(['./' dir_name '/solution_top.txt']);   
   % Uf(:,isim) = U;
end

plot(Uf)

save Uf.mat Uf

function T=tempo(T_p)
%%dato il periodo dell'orbita in secondi, lo trasforma in giorni,ore,min e
%%secondi
secondi=T_p;
giorni=floor(secondi/(24*3600));
decimale_g=(secondi/(24*3600))-giorni;
ore=floor(decimale_g*24);
decimale_o=decimale_g*24-ore;
minuti=floor(decimale_o*60);
decimale_m=decimale_o*60-minuti;
secondi=floor(decimale_m*60);

T=[giorni ore minuti secondi];
%fprintf('\nil tempo di %f s equivale a:\ngiorni=%f\nore=%f\nminuti=%f\nsecondi=%f\n',T_p,giorni,ore,minuti,secondi);
return

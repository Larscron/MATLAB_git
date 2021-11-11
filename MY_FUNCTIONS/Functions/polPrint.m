function [] = polPrint(recNum)
% polPrint Pol prints takes in a complex number in rectangular form 
% and prints it in polar form 

mag=num2str(abs(recNum)); %get magnitude and convert to string
ang=num2str(angle(recNum)*180/pi); %get angle and convert to string
str=append(mag,'<',ang,' deg\n'); %combine strings
fprintf(str); %print strings
end


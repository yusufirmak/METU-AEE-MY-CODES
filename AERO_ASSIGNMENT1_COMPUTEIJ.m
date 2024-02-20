function [I,J] = COMPUTE_IJ(xc,yc,xbp,ybp,gama,S)

num_panel= length(xc)

I= zeros(num_panel,num_panel);
J= zeros(num_panel,num_panel);

for i=1:1:num_panel
    for j=1:1:num_panel
        if(j ~= i)
            A= -(xc(i)-xbp(j))*cosd(gama(j))-(yc(i)-ybp(j))*sind(gama(j));
            B= (xc(i)-xbp(j))^2 + (yc(i)-ybp(j))^2;
            Cn= sind(gama(i)-gama(j)); 
            Dn = -(xc(i)-xbp(j))*sind(gama(i))+(yc(i)-ybp(j))*cosd(gama(i));
            Ct = -cosd(gama(i)-gama(j));
            Dt = (xc(i)-xbp(j))*cosd(gama(i))+(yc(i)-ybp(j))*sind(gama(i)); 
            E  = sqrt(B-A^2); 
            if (~isreal(E))
                E = 0
            end

            %computation of I
            term1  = 0.5*Cn*log((S(j)^2+2*A*S(j)+B)/B);
            term2  = ((Dn-A*Cn)/E)*(atan2((S(j)+A),E) - atan2(A,E));  
            I(i,j) = term1 + term2;
            %computation of J
            term1  = 0.5*Ct*log((S(j)^2+2*A*S(j)+B)/B);
            term2  = ((Dt-A*Ct)/E)*(atan2((S(j)+A),E) - atan2(A,E)); 
            J(i,j) = term1 + term2;
        end

    end
end        

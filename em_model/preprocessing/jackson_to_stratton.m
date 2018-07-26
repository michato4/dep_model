function [s] = jackson_to_stratton(j)
%   Convert Jackson's type of cartesian multipoles to Stratton's type
	if numel(j) <= 3^1
        s = j;
    elseif numel(j) == 3^2
        % qudrupole
        p2 = j;
%         aux = 1/3*p2;
%         aux(3,3) = 0; % libovolné èíslo (dost malé kvùli numerice)
%         aux(2,2) = 1/3*(2*p2(2,2)+p2(1,1)+3*aux(3,3));
%         aux(1,1) = 1/2*(p2(1,1)+aux(2,2)+aux(3,3));
%         p2 = aux;
%         s = p2;
        s = p2/3;
    elseif numel(j) == 3^3
        p3 = j;
%         % octupole
%         O111S = 0;
%         O112S = (4*p3(1,1,2))/45 + p3(3,3,2)/45;
%         O113S = (4*p3(1,1,3))/45 + p3(3,2,2)/45;
%         O222S = 0;
%         O123S = p3(1,2,3)/15;
%         O333S = 0;
%         O122S = (4*p3(1,2,2))/45 + p3(1,3,3)/45;
%         O133S = p3(1,2,2)/45 + (4*p3(1,3,3))/45;
%         O322S = p3(1,1,3)/45 + (4*p3(3,2,2))/45;
%         O332S = p3(1,1,2)/45 + (4*p3(3,3,2))/45;
%         p3(:,:,1) = [O111S O112S O113S; O112S O122S O123S; O113S O123S O133S];
%         p3(:,:,2) = [O112S O122S O123S; O122S O222S O322S; O123S O322S O332S];
%         p3(:,:,3) = [O113S O123S O133S; O123S O322S O332S; O133S O332S O333S];
%         s = p3;
         s = p3/15;
    elseif numel(j) == 3^4
         p4 = j;
         s = p4/105;
    elseif numel(j) == 3^5
         p5 = j;
         s = p5/(945);
    end
end


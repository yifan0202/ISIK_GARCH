function [u1, u2, u3, u4] = Transf(s,k)
% Input: s, k
% Output: u1, u2, u3, u4
    Gammat = 1+s/6+(k-3).^2/24;
    m1 = s.*(k-3)./3./Gammat;
    m2 = (1+7/6*s.^2+3/8*(k-3).^2)./Gammat;
    m3 = (2*s+4*s.*(k-3))./Gammat;
    m4 = (3+2*(k-3)+25/2*s.^2+41/8*(k-3))./Gammat;
    
    u1 = m1;
    u2 = m2-m1.^2;
    u3 = m3-3*m2.*m1+2*m1.^3;
    u4 = m4-4*m3.*m1+6*m2.*m1.^2-3*m1.^4;
end
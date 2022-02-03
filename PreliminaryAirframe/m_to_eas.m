function [EAS] = m_to_eas(m,alt)
    [~,a0,p0,~] = atmosisa(0);
    [~,~,p1,~] = atmosisa(alt);
    EAS = a0*m*sqrt(p1/p0);
end

    
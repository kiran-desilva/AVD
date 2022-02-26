function K=d_cell_k_s(a,b,R1,t1)
persistent cd
persistent cdi

if a>=b
    x=b/sqrt(R1*t1);
    z=a/b;
    load a_b_int.mat
    
    
    
    if isempty(cd)
        load a_b_int.mat;
        cd = a_b_fit_arr;
    end
    if isempty(cdi)
        load a_b_int.mat a_b_Int3D;
        cdi = a_b_Int3D;
    end
    
    
    %select the right curve
    %     if use3D
%     k_s = cdi(x,z);
%         else
            [~, idx] = min(abs(a_b_range - z));
    
            k_s = a_b_fit_arr{idx}(x); % Read off of catchpole boi
            % 	end
    
    
%         end
    
    if b>a
        x=a/sqrt(R1*t1);
        z=b/a;
        
        load b_a_int.mat
        if isempty(cd)
            load b_a_int.mat;
            cd = b_a_fit_arr;
        end
        if isempty(cdi)
            load b_a_int.mat b_a_Int3D;
            cdi = b_a_Int3D;
        end
        
        
%         select the right curve
%         if use3D
%             k_s = cdi(x,z);
%         else
            [~, idx] = min(abs(cd.b_a_range - z));
            
            k_s = cd.b_a_fit_arr{idx}(x); % Read off of catchpole boi
%         end
        
        
    end
    
    K=k_s(end);
end

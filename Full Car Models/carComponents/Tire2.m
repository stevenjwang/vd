classdef Tire2  
    properties 
        gamma 
        p_i 
        Fx_parameters
        Fy_parameters       
        friction_scaling_factor
        
        % Fast compiled interpolants
        camber4_interp
        camber2_interp
        camberM2_interp
        camberM4_interp
    end
           
    methods
        function obj = Tire2(p_i,Fx_parameters,Fy_parameters,friction_scaling_factor)
            obj.p_i = p_i;
            obj.Fx_parameters = Fx_parameters;
            obj.Fy_parameters = Fy_parameters;      
            obj.friction_scaling_factor = friction_scaling_factor;
            
            % Load camber data
            load("camberratiossmoothed.mat", "camber2indices", "camber4indices", "camber2ratio", "camber4ratio");
            
            % --- FIXED: Ensuring ascending order for negative indices ---
            obj.camber4_interp  = griddedInterpolant(camber4indices, camber4ratio, 'linear', 'nearest');
            obj.camber2_interp  = griddedInterpolant(camber2indices, camber2ratio, 'linear', 'nearest');
            
            % Flip both the negated indices and the ratios so they go from (e.g.) -10 to 0
            obj.camberM2_interp = griddedInterpolant(flip(-camber2indices), flip(camber2ratio), 'linear', 'nearest');
            obj.camberM4_interp = griddedInterpolant(flip(-camber4indices), flip(camber4ratio), 'linear', 'nearest');
        end
        
        function out = F_y(obj,alpha,kappa,F_z,gamma)
            
            % Calculate Camber Multipliers using fast interpolants
            cambermultiplier4 = obj.camber4_interp(alpha);
            cambermultiplier2 = obj.camber2_interp(alpha);
            cambermultiplier0 = ones(size(alpha)); % Fast array of 1s
            cambermultiplierminus2 = obj.camberM2_interp(alpha);
            cambermultiplierminus4 = obj.camberM4_interp(alpha);
            
            % 2D interpolation for the specific gamma value
            % If gamma is an array (4 tires), we loop it fast.
            cambermultiplier = zeros(size(gamma));
            for i = 1:numel(gamma)
                y_vals = [cambermultiplier4(i), cambermultiplier2(i), cambermultiplier0(i), ...
                          cambermultiplierminus2(i), cambermultiplierminus4(i)];
                % A single fast interp1 per tire is acceptable here, 
                % but if gamma is locked static in Car.m, this evaluates instantly.
                cambermultiplier(i) = interp1([4, 2, 0, -2, -4], y_vals, gamma(i), 'linear', 'extrap');
            end
            
            cambershiftMod = 16.125 .* gamma .* (F_z ./ 250);
            
            % Convert units
            gamma = 0; % Enforced from original script logic
            alpha = alpha .* 0.0174533; % deg to rad
            F_z = abs(F_z .* 0.224809); % N to lbf
            
            % --- OPTIMIZATION 2: SCALAR BROADCASTING & MATH STABILITY ---
            F_z0 = 200; % nominal load
            p_i0 = 13;  % nominal pressure
            
            % Standard Pacejka uses F_z0 in denominator to prevent NaN crashes on wheel lift
            df_z = (F_z - F_z0) ./ F_z0; 
            dp_i = (obj.p_i - p_i0) ./ obj.p_i;
            
            % Parameters mapping (Kept original array mapping as it is fast enough)
            p_cy1 = obj.Fy_parameters(1); p_dy1 = obj.Fy_parameters(2); p_dy2 = obj.Fy_parameters(3);
            p_dy3 = obj.Fy_parameters(4); p_ey1 = obj.Fy_parameters(5); p_ey2 = obj.Fy_parameters(6);
            p_ey3 = obj.Fy_parameters(7); p_ey4 = obj.Fy_parameters(8); p_ey5 = obj.Fy_parameters(9);
            p_hy1 = obj.Fy_parameters(10); p_hy2 = obj.Fy_parameters(11); p_ky1 = obj.Fy_parameters(12);
            p_ky2 = obj.Fy_parameters(13); p_ky3 = obj.Fy_parameters(14); p_ky4 = obj.Fy_parameters(15);
            p_ky5 = obj.Fy_parameters(16); p_ky6 = obj.Fy_parameters(17); p_ky7 = obj.Fy_parameters(18);
            p_py1 = obj.Fy_parameters(19); p_py2 = obj.Fy_parameters(20); p_py3 = obj.Fy_parameters(21);
            p_py4 = obj.Fy_parameters(22); p_py5 = obj.Fy_parameters(23); p_vy1 = obj.Fy_parameters(24);
            p_vy2 = obj.Fy_parameters(25); p_vy3 = obj.Fy_parameters(26); p_vy4 = obj.Fy_parameters(27);
            r_by1 = obj.Fy_parameters(28); r_by2 = obj.Fy_parameters(29); r_by3 = obj.Fy_parameters(30);
            r_by4 = obj.Fy_parameters(31); r_cy1 = obj.Fy_parameters(32); r_ey1 = obj.Fy_parameters(33);
            r_ey2 = obj.Fy_parameters(34); r_hy1 = obj.Fy_parameters(35); r_hy2 = obj.Fy_parameters(36);
            r_vy1 = obj.Fy_parameters(37); r_vy2 = obj.Fy_parameters(38); r_vy3 = obj.Fy_parameters(39);
            r_vy4 = obj.Fy_parameters(40); r_vy5 = obj.Fy_parameters(41); r_vy6 = obj.Fy_parameters(42);
     
            lambda_cy = 1; lambda_ey = 1; lambda_hy = 0; lambda_kyalpha = 1;  
            lambda_kygamma = 1; lambda_muy = 1; lambda_vy = 0; lambda_ykappa = 1; lambda_vykappa = 0;  
            
            % Magic Formula Equations
            K_yalpha = p_ky1.*F_z0.*(1+p_py1.*dp_i).*sin(p_ky4.*atan(F_z./...
                ((p_ky2 + p_ky5.*gamma.^2).*(1+p_py2.*dp_i).*F_z0))).*...
                (1-p_ky3.*abs(gamma)).*lambda_kyalpha;
            K_ygamma = (p_ky6+p_ky7.*df_z).*(1+p_py5.*dp_i).*F_z.*lambda_kygamma;
            S_vy0 = F_z.*(p_vy1+p_vy2.*df_z).*lambda_vy.*lambda_muy;
            S_vygamma = F_z.*(p_vy3+p_vy4.*df_z).*gamma.*lambda_kygamma.*lambda_muy;
            S_vy = S_vy0+S_vygamma;
            S_hy0 = (p_hy1+p_hy2.*df_z).*lambda_hy;
            S_hygamma = (K_ygamma.*gamma-S_vygamma)./K_yalpha;
            S_hy = S_hy0+S_hygamma;
            alpha_y = alpha+S_hy;
            mu_y = (p_dy1+p_dy2.*df_z).*(1-p_dy3.*gamma.^2).*...
                (1+p_py3.*dp_i+p_py4.*dp_i.^2).*lambda_muy;
            
            C_y = p_cy1.*lambda_cy;
            D_y = mu_y.*F_z;
            E_y = (p_ey1+p_ey2.*df_z).*(1+p_ey5.*gamma.^2-(p_ey3+p_ey4.*gamma).*sign(alpha_y)).*lambda_ey;
            B_y = K_yalpha./C_y./D_y;
            
            % Combined Slip
            B_ykappa = (r_by1 + r_by4.*gamma.^2).*cos(atan(r_by2.*(alpha - r_by3))).*lambda_ykappa;
            C_ykappa = r_cy1;
            E_ykappa = r_ey1 + r_ey2.*df_z;
            S_hykappa = r_hy1 + r_hy2.*df_z; 
            kappa_s = kappa + S_hykappa;
            D_vykappa = mu_y.*F_z.*(r_vy1 + r_vy2.*df_z + r_vy3.*gamma).*cos(atan(r_vy4.*alpha));
            S_vykappa = D_vykappa.*sin(r_vy5.*atan(r_vy6.*kappa)).*lambda_vykappa;
            G_ykappa = cos(C_ykappa.*atan(B_ykappa.*kappa_s - E_ykappa.*(B_ykappa.*kappa_s - ...
                atan(B_ykappa.*kappa_s))))./cos(C_ykappa.*atan(B_ykappa.*S_hykappa - E_ykappa...
                .*(B_ykappa.*S_hykappa - atan(B_ykappa.*S_hykappa))));
            
            % Lateral Force
            F_y = D_y.*sin(C_y.*atan(B_y.*alpha_y-E_y.*(B_y.*alpha_y-atan(B_y.*alpha_y))))+S_vy;
            
            % Transpose is slow, keeping matrices aligned is better.
            % G_ykappa is 1x4, F_y is 1x4, S_vykappa is 1x4
            F_y = (G_ykappa .* F_y) + S_vykappa; 
            
            F_y(F_z == 0) = 0; % zero load
            
            % Scale back to N
            F_y2 = F_y .* 4.44822 .* obj.friction_scaling_factor; 
            out = (F_y2 .* cambermultiplier) + cambershiftMod;
            
            % Transpose final output to match expected Car.m orientation
            out = out.'; 
        end        
        
        function out = F_x(obj,alpha,kappa,F_z, gamma)
            
            gamma = gamma .* 0.0174533; % deg to rad
            alpha_f = alpha .* 0.0174533; % deg to rad
            F_z = abs(F_z .* 0.224809); % N to lbf
            
            % --- SCALAR BROADCASTING OPTIMIZATION ---
            F_z0 = 200; 
            p_i0 = 13; 
            
            df_z = (F_z - F_z0) ./ F_z0; % Fixed Pacejka denominator
            dp_i = (obj.p_i - p_i0) ./ obj.p_i;
            
            % Parameters
            p_cx1 = obj.Fx_parameters(1); p_dx1 = obj.Fx_parameters(2); p_dx2 = obj.Fx_parameters(3);
            p_dx3 = obj.Fx_parameters(4); p_ex1 = obj.Fx_parameters(5); p_ex2 = obj.Fx_parameters(6);
            p_ex3 = obj.Fx_parameters(7); p_ex4 = obj.Fx_parameters(8); p_hx1 = obj.Fx_parameters(9);
            p_hx2 = obj.Fx_parameters(10); p_kx1 = obj.Fx_parameters(11); p_kx2 = obj.Fx_parameters(12);
            p_kx3 = obj.Fx_parameters(13); p_px1 = obj.Fx_parameters(14); p_px2 = obj.Fx_parameters(15);
            p_px3 = obj.Fx_parameters(16); p_px4 = obj.Fx_parameters(17); p_vx1 = obj.Fx_parameters(18);
            p_vx2 = obj.Fx_parameters(19); r_bx1 = obj.Fx_parameters(20); r_bx2 = obj.Fx_parameters(21);
            r_bx4 = obj.Fx_parameters(22); r_cx1 = obj.Fx_parameters(23); r_ex1 = obj.Fx_parameters(24);
            r_ex2 = obj.Fx_parameters(25); r_hx1 = obj.Fx_parameters(26);
            
            lambda_cx = 1; lambda_ex = 1; lambda_hx = 0; lambda_kxkappa = 1; 
            lambda_mux = 1; lambda_vx = 0; lambda_xalpha = 1; 
            
            % Magic Formula Equations
            mu_x = (p_dx1+p_dx2.*df_z).*(1-p_dx3.*gamma.^2).*(1+p_px3.*dp_i+p_px4.*dp_i.^2).*lambda_mux;
            
            K_xkappa = (p_kx1+p_kx2.*df_z).*exp(p_kx3.*df_z).*(1+p_px1.*dp_i+p_px2.*dp_i.^2).*F_z.*lambda_kxkappa;
            K_xkappa(~isfinite(K_xkappa)) = 1e200; % Vectorized finiteness check
            
            S_hx = (p_hx1+p_hx2.*df_z).*lambda_hx;
            S_vx = (p_vx1+p_vx2.*df_z).*F_z.*lambda_vx.*lambda_mux;
            kappa_x = kappa + S_hx;
            C_x = p_cx1.*lambda_cx;
            D_x = mu_x.*F_z;
            E_x = (p_ex1+p_ex2.*df_z+p_ex3.*df_z.^2).*(1-p_ex4.*sign(kappa_x)).*lambda_ex;
            B_x = K_xkappa./(C_x.*D_x);
            
            % Combined slip
            B_xalpha = (r_bx1+r_bx4.*gamma.^2).*cos(atan(r_bx2.*kappa)).*lambda_xalpha;
            C_xalpha = r_cx1;
            E_xalpha = r_ex1 + r_ex2.*df_z;
            S_hxalpha = r_hx1;
            alpha_s = alpha_f + S_hxalpha;
            G_xalpha = cos(C_xalpha.*atan(B_xalpha.*alpha_s - E_xalpha.*(B_xalpha.*alpha_s - ...
                atan(B_xalpha.*alpha_s))))./cos(C_xalpha.*atan(B_xalpha.*S_hxalpha - E_xalpha...
                .*(B_xalpha.*S_hxalpha - atan(B_xalpha.*S_hxalpha))));
            
            % Longitudinal Force
            F_x = (D_x.*sin(C_x.*atan(B_x.*kappa_x-E_x.*(B_x.*kappa_x-atan(B_x.*kappa_x))))+S_vx).*G_xalpha;
            F_x(F_z == 0) = 0; % zero load
            
            out = (F_x .* 4.44822 .* obj.friction_scaling_factor).'; % Scale and transpose to match Car.m
        end
        
    end
end
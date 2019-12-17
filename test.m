%# common variables
wb = 2*pi*50;
s = sym('s');
I = eye(2);
w = logspace(-2,5,500)*2*pi;
%w = [-flip(w),w];

%%
%# default parameters
para1.J = 3.5*2/(2*pi*50)^2;
para1.D = 1/(2*pi*50)^2;
para1.L = 0.05/(2*pi*50);
para1.R = 0.01;

para2.V_dc = 2.5;
para2.C_dc = 2*0.1*para2.V_dc^2;
para2.kp_v_dc = para2.V_dc*para2.C_dc*(10*2*pi);
para2.ki_v_dc = para2.kp_v_dc*(10*2*pi)/4;
para2.kp_pll = 2*2*pi;
para2.ki_pll = para2.kp_pll * (2*2*pi)/4; 
para2.tau_pll = 1/(2*pi*200);
para2.k_pf = 0;
para2.L = 0.05/(2*pi*50);
para2.R = 0.01;
para2.kp_i_dq = para2.L * (500*2*pi);
para2.ki_i_dq = para2.kp_i_dq *(500*2*pi)/4;

%%
%# Test B: low inertia interact with pll
if 1
    
    layout = 3;
    
    if layout == 1      
        %2 buses : test the match between symbolic and tf models     
        %-----------------------------------------------------------
        %         |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi | Qmin | Qmax |
        Bus     = [ 1     1     1.0     0     1.0    0    0.0    0     -1     1;
                    2     3     1.0     0     1.0  -0.2   0.0    0     -1     1];
        %-----------------------------------------------------------
        %         |  From |  To   |   R   |   L   |   C   |   G   |
        %         |  Bus  |  Bus  |       |       |       |       |
        Line    = [  1       2      0.00     0.3    0.00      inf;
                     1       1        0       0     1e-1      2.0;
                     2       2        0       0     1e-1      0.0];     
    elseif layout == 2      
        %3 buses : test generator interaction        
        %-----------------------------------------------------------
        %         |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi | Qmin | Qmax |
        Bus     = [ 1     1     1.0     0     0.5    0    0.0    0     -1     1;
                    2     2     1.0     0     0.5    0    0.0    0     -1     1;
                    3     2     1.0     0     1.0    0    0.0    0     -1     1];            
        %-----------------------------------------------------------
        %         |  From |  To   |   R   |   L   |   C   |   G   |
        %         |  Bus  |  Bus  |       |       |       |       |
        Line    = [  1       2      0.01     0.3      0      inf;
                     2       3      0.01     0.3      0      inf;
                     3       1      0.01     0.3      0      inf;
                     1       1        0       0     2e-2     1.0;
                     2       2        0       0     2e-2     1.0;
                     3       3        0       0     2e-2     0.0]; 
    elseif layout == 3      
        %4 buses : test generator converter interaction  
        %-----------------------------------------------------------
        %         |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi | Qmin | Qmax |
        Bus     = [ 1     1     1.0     0     0.5    0    0.0    0     -1     1;
                    2     2     1.0     0     0.5    0    0.0    0     -1     1;
                    3     2     1.0     0     0.5    0    0.0    0     -1     1;
                    4     3     1.0     0     0.5  -0.2   0.0    0     -1     1];
        %-----------------------------------------------------------
        %         |  From |  To   |   R   |   L   |   C   |   G   |
        %         |  Bus  |  Bus  |       |       |       |       |
        Line    = [  1       2      0.01     0.5      0      inf;
                     2       3      0.01     0.5      0      inf;
                     3       1      0.01     0.5      0      inf;
                     3       4      0.01     0.5      0      inf;
                     1       1        0       0     2e-2     1.0;
                     2       2        0       0     2e-2     1.0;
                     3       3        0       0     2e-2     0.1;
                     4       4        0       0     2e-2     0.1];           
    end

    [~,~,Ang0,P0,Q0,V0]=PowerFlow(Bus,Line);
    
    nbus = max(Bus(:,1));
    
    [Ytf1,Yb1] = YbusCalcTF(Line(1:(end-nbus),:),wb);
    [Ytf2,Yb2] = YbusCalcTF(Line((end-nbus+1):end,:),wb);
    Yb1 = minreal(Yb1);
    Zb2 = inv(Yb2);
    Zb2 = minreal(Zb2);
    Zb = feedback(Zb2,Yb1);
            
    for gain = logspace(log10(1),log10(20),10)
        
        para1_ = para1;
        para1_.J = para1_.J/10;
        para2_ = para2;
        para2_.kp_pll = 2*2*pi/2 *gain/2;
        para2_.ki_pll = para2_.kp_pll * (2*2*pi)/4/2 *gain/2;
        para2_.kp_v_dc = para2_.V_dc*para2_.C_dc*(10*2*pi) /10 *gain/2;
        para2_.ki_v_dc = para2_.kp_v_dc*(10*2*pi)/4 /10 *gain/2;
          
        if layout == 1
            type = {0,10};
            para = {para1_,para2_};
        elseif layout == 2
            % three generators
            type = {0,0,0};
            para = {para1,para1,para1};
        elseif layout == 3
            % four generators
            %type = {0,0,0,0};
            %para = {para1,para1,para1_,para1};
            % three generators and one wind farm
            type = {0,0,0,10};
            para = {para1,para1,para1_,para2_};
        end
        
        Gm = cell(1,nbus);
        Gc = cell(1,nbus);
        for n = 1:nbus
            [~,Gm{n},Gc{n},~] = MdlCreate('type', type{n} ,'flow',[-P0(n) -Q0(n) V0(n) Ang0(n) wb],'para',para{n});
        end
        
        if layout == 1
            Ys1 = Gc{1}(2:3,2:3);
            Ys2 = Gc{2}(2:3,2:3);
            Ys1 = Ys1 + eye(2)*Line(2,6);
            Ys1 = Ys1 + [(1j + 1/wb*s) 0;0 (-1j + 1/wb*s)]*Line(2,5);
            Ys2 = Ys2 + [(1j + 1/wb*s) 0;0 (-1j + 1/wb*s)]*Line(3,5);                       
            Zs1 = Ys1^(-1);
            Zs1 = Zs1 + [(1j + 1/wb*s) 0;0 (-1j + 1/wb*s)]*Line(1,4);
            Ys1 = Zs1^(-1);
            Y = Ys1 + Ys2;
            Z = Y^(-1);
            pc = vpasolve(1/Z(1,1))/2/pi;
        end
            
        Gm = MdlLink(Gm);
        Gsys = feedback(Gm,Zb,(nbus+1):(3*nbus),(nbus+1):(3*nbus));
        
        psys = pole(Gsys)/2/pi;
        figure(layout);
        scatter(real(psys),imag(psys),'x','LineWidth',1.5);
        hold on; grid on;
        axis([-6,2,-10,10]);
        
        if layout == 1
            scatter(real(pc),imag(pc),'o','LineWidth',1.5);            
        end

    end
end

print(gcf,'fig.png','-dpng','-r600');
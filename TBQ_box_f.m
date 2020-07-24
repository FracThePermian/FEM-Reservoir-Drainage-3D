function [ T,Q,B ] = TBQ_box_f(BC,P_B,para_wells,q_wells)
global phi k
[para] = reservoir; N = para.N; NX = para.NX; NY = para.NY; mu = para.mu; Bw = para.Bw; dx = para.dx; dy = para.dy; Ax = para.Ax; Ay = para.Ay;
h = para.h; ct = para.ct; Bw = para.Bw
[res,block,nums,fluid]=res_parameters;
T = sparse(N,N); B = sparse(N,N); Q = sparse(N,1); Q(para_wells) = Q(para_wells) + q_wells;
 
for i = 1:N
    if i+NX <= N 
        T_diag = T_frac(i,i+NX,Ay,mu,Bw,dy);
        T(i,i+NX) = T(i,i+NX)-T_diag;
        T(i+NX,i) = T(i+NX,i)-T_diag;
    end
    
    % This is for diagonals toward the inside of sparsity
    if (mod(i,NX) ~= 0) && (i+1 <= N)  
        T_diag = T_frac(i,i+1,Ax,mu,Bw,dx);
        T(i,i+1) = T(i,i+1)-T_diag;
        T(i+1,i) = T(i+1,i)-T_diag;
    end
    T(i,i) = abs(sum(T(i,:)));%Sum up every row on every iteration
    if BC(i) ~= 0 %if and only if the boundary condition exists
        T(i,i) = T(i,i)+BC(i)*2*T_frac(i,i,Ax,mu,Bw,dx);
    end
B(i,i) = dx*dy*h*phi(i)*ct/Bw; %phi as such for next project '2'
Q(i) = Q(i) + 2*BC(i)*P_B(i)*T_frac(i,i,Ax,mu,Bw,dx);
end
end





s=second_pass;

K=length(tt);
T=round(min(tt)):round(max(tt));
N=length(T);

dt=1; % 1 hour on each side
gamma=2;
sigma_h=1;
std_thresh=0.1;

rhs=zeros(N,1);
v=zeros(N,1);
sw=zeros(N,1);
for k=1:N
    
    ind_s=find(tt>=T(k)-dt & tt<=T(k)+dt & isfinite(sst)==1);
    N_s=length(ind_s);
    
    ind_a=find(tt>=T(k)-dt & tt<=T(k)+dt & isfinite(a)==1);
    N_a=length(ind_a);
    
    if (N_a==0)
        v(k)=0;
        rhs(k)=0;
        sw(k)=gamma;
    else if (N_a==1) 
        v(k)=1;
        rhs(k)=a(ind_a);
        sw(k)=gamma;
    else if (N_a==2)
            sw(k)=gamma; 
            w=exp(-(tt(ind_a)-T(k)).^2/(2*(sigma_h)^2));
            v(k)=sum(w);
            rhs(k)=sum(w.*a(ind_a));
%             v(k)=1;
%             rhs(k)=max(a(ind_a));
        else 
            x0=ind_a(1);
            x=(ind_a-x0)';
            y0=a(ind_a(1));
            y=a(ind_a)-y0;
            A=[ones(N_a,1) x x.^2];
            p=A\y';
%             g=A*p;
%             sw(k)=std(g-y');
            x=(ind_s-x0)';
            A=[ones(N_s,1) x x.^2];
            g=A*p;
            sw(k)=std(g-(sst(ind_s)-y0)');
            w=exp(-(tt(ind_a)-T(k)).^2/(2*(sigma_h)^2));
            v(k)=sum(w);
            rhs(k)=sum(w.*a(ind_a));
            if sw(k)>std_thresh
                sw(k)=gamma+sw(k);
            end
        end
        end
    end
end

W=zeros(N,N); 
W(1,1)=v(1)+sw(1); W(1,2)=-sw(1);
W(N,N)=v(N)+sw(N); W(N,N-1)=-sw(N);
for k=2:N-1
    W(k,k-1)=-sw(k-1);
    W(k,k)=v(k)+sw(k-1)+sw(k);
    W(k,k+1)=-sw(k);
end
c=W\rhs;
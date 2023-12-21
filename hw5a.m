function monte_carlo_a(n_methane)
%%Open MFI Data
fileId=fopen('MFI.txt','r');
MFI=fscanf(fileId,'%g %g %g',[3 inf]);
MFI=MFI';
fclose(fileId);

% Input Values
n_zeolite=192;
N=200000;
T=300;
kb=1;
Lx=20.09;
Ly=19.738;
Lz=13.142;
ep1=148.0;
sig1=3.73;
ep2=115.0;
sig2=3.47;
R_cut=6.5;

%Initialization of arrays
total_potential_energy=zeros(N,1);
xz=zeros(n_zeolite,1);
yz=zeros(n_zeolite,1);
zz=zeros(n_zeolite,1);
xm=zeros(n_methane,1);
ym=zeros(n_methane,1);
zm=zeros(n_methane,1);

%% Molecule allocation
for k=1:n_zeolite
    xz(k)=MFI(k,1);
    yz(k)=MFI(k,2);
    zz(k)=MFI(k,3);
end

for k=1:n_methane
    xm(k)=Lx*rand;
    ym(k)=Ly*rand;
    zm(k)=Lz*rand;
end

%% Initial potential energy
% Methane-Methane
U_mm=0;
for k=1:n_methane
    for kk=k+1:n_methane
        dx=xm(k)-xm(kk);
        dy=ym(k)-ym(kk);
        dz=zm(k)-zm(kk);
        r=sqrt(dx^2+dy^2+dz^2);
        if r>R_cut
            U_pair=0;
        else
            U_pair=4*ep1*(sig1^12/r^12-sig1^6/r^6);
        end
        U_mm=U_mm+U_pair;
    end
end

% Methane-Zeolite
U_mz=0;
for k=1:n_zeolite
    for kk=1:n_methane
        dx=xz(k)-xm(kk);
        dy=yz(k)-ym(kk);
        dz=zz(k)-zm(kk);
        r=sqrt(dx^2+dy^2+dz^2);
        if r>R_cut
            U_pair=0;
        else
            U_pair=4*ep2*(sig2^12/r^12-sig2^6/r^6);
        end
        U_mz=U_mz+U_pair;
    end
end

U_total=U_mm+U_mz;

%% Monte Carlo Loop
for j=1:N
    rand_idx=randi([1,n_methane]);
    % Compute old potential energy of Methane-Methane
    U_mm_old=0;
    for kk=1:n_methane
        if kk~=rand_idx
            dx=xm(rand_idx)-xm(kk);
            dy=ym(rand_idx)-ym(kk);
            dz=zm(rand_idx)-zm(kk);
            r=sqrt(dx^2+dy^2+dz^2);
            if r>R_cut
                U_pair=0;
            else
                U_pair=4*ep1*(sig1^12/r^12-sig1^6/r^6);
            end
            U_mm_old=U_mm_old+U_pair;
        end
    end

    % Compute dld total potential energy of Methane-Zeolite
    U_mz_old=0;
    for kk=1:n_zeolite
        dx=xm(rand_idx)-xz(kk);
        dy=ym(rand_idx)-yz(kk);
        dz=zm(rand_idx)-zz(kk);
        r=sqrt(dx^2+dy^2+dz^2);
        if r>R_cut
            U_pair=0;
        else
            U_pair=4*ep2*(sig2^12/r^12-sig2^6/r^6);
        end
        U_mz_old=U_mz_old+U_pair;
    end
    U_old=U_mm_old+U_mz_old;
    
    %Move a molecule at random
    xnew=xm(rand_idx)+(4*rand-2);
    ynew=ym(rand_idx)+(4*rand-2);
    znew=zm(rand_idx)+(4*rand-2);

    %Compute new potential energy of Methane-Methane
    U_mm_new=0;
    for kk=1:n_methane
        if kk~=rand_idx
            dx=xnew-xm(kk);
            dy=ynew-ym(kk);
            dz=znew-zm(kk);
            r=sqrt(dx^2+dy^2+dz^2);
            if r>R_cut
                U_pair=0;
            else
                U_pair=4*ep1*(sig1^12/r^12-sig1^6/r^6);
            end
            U_mm_new=U_mm_new+U_pair;
        end
    end
        
    %Compute new potential energy of Methane-Zeolite
    U_mz_new=0;
    for kk=1:n_zeolite
        dx=xnew-xz(kk);
        dy=ynew-yz(kk);
        dz=znew-zz(kk);
        r=sqrt(dx^2+dy^2+dz^2);
        if r>R_cut
            U_pair=0;
        else
            U_pair=4*ep2*(sig2^12/r^12-sig2^6/r^6);
        end
        U_mz_new=U_mz_new+U_pair;
    end
    U_new=U_mm_new+U_mz_new;
    
    %%Metropolis-Hastings algorithm
    metro=min(1,exp(-(U_new-U_old)/(kb*T)));
    if ~((xnew<0)||(xnew>Lx)||(ynew<0)||(ynew>Ly)||(znew<0)||(znew>Lz))
        if metro>=rand
            xm(rand_idx)=xnew;
            ym(rand_idx)=ynew;
            zm(rand_idx)=znew;
            U_total=U_total+(U_new-U_old);
        end
    end
    total_potential_energy(j)=U_total;
end

average_potential_energy=sum(total_potential_energy(0.1*N+1:N))/(0.9*N);
fprintf('average_potential_energy=%.4e Kelvin \n',average_potential_energy);
end

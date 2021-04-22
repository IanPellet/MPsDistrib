Viteglobal dt dx
Sauvegarde=true;
FichConcentration = 'Test1';

tf= 100*86400; %duree 
dtmax=100; %pas de temps max (/!\ la particule ne doit pas parcourir plus d'une maille par dt)
Tdes=60*60; %pas dessin concentration
Nsauv=floor(tf/Tdes)+1;


%CalculVitesses
for indNom=1:size(Nom,1)
    Ws(:,:,indNom)=eval(['Vitesse' cell2mat(Nom(indNom)) '(d,S,D_);']);
end


% schema='Centre1';
% schema='CentreVF';
schema='UpWind';
% schema='Ordre4';



% %_______________ Une Vitesse______________________________________________%

Concentration(1,:)=C;

    if rop<row
        u0=-(Ws(:,:,3));
    else
        u0=(Ws(:,:,3));
    end

u0_=max(u0);Nu0_=max(Nu);
       if u0_~=0 & Nu0_~=0; 
           dt=min(dx/abs(u0_)*0.5,dx*dx/(2*Nu0_)*0.5); 
       elseif u0_==0 & Nu0_~=0; 
           dt=dx*dx/(2*Nu0_)*0.5;
       elseif u0_~=0 & Nu0_==0; 
           dt=min(dx/abs(u0_)*0.5,dx*dx/(2*Nu0_)*0.5); 
       else, 
           dt=dtmax;
       end

    u= u0_*ones(1,size(C,2)+1);u(1)=0;u(end)=0;
    C=Concentration(1,:);

    for t= dt:dt:tf
        C=StepTransport (u,Nu,C,schema); 
        %disp(mod(t,Tdes))

        if (mod(t,Tdes)<=dt/2 | Tdes-mod(t,Tdes)<=dt/2 )
           index=round(t/Tdes)+1;
           disp([' Temps : ' num2str(t) 's -' ...
                 ' Numero de sauvegarde : ' num2str(index)  ...
                 ' - Concentration Totale : ' num2str(sum(C)) ])

           Concentration(index,:)=C;

           AffichageConcentration

        end
        
           

        
        
    end
    
Nsauv_=min(size(Concentration,1),Nsauv);
Csaved=Concentration(1:Nsauv_,:);
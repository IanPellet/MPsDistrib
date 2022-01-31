% Fit a power law to data to compare the probability p_c to this found for
% the two fragility conditions : pillar of strength and plane of weakness.
%
% Theorical data :
% Turcotte, D. L. (1986). Fractals and fragmentation. 
% Journal of Geophysical Research, 91(B2), 1921–1926. https://doi.org/10.1029/JB091iB02p01921
%
% Réunion biblio 210621

% global sizes nSizes
data = load_MPs_data('../Data/data_mps.txt');
data = data(data.date == datetime(2021,03,18),:);
sizes = unique(data.med_size);
nSizes = NaN(size(sizes));
for i = 1:length(nSizes)
    nSizes(i) = sum(data(data.med_size==sizes(i),:).n);
end, clear i,
bins = [50 200 250 500 1000:500:5000 30000];

nSizes = nSizes./diff(bins');

weak = {'Pillar', 'Plane'};

nmin = 3;
delta = NaN(length(sizes)-nmin,length(weak));
rmse = NaN(length(sizes)-nmin,length(weak));

for j=1:length(weak)
    for i = 1:length(sizes)-nmin
        exclude = i;
        weakness = weak{j};
        
        curvefit = fit(sizes(exclude:end),nSizes(exclude:end),'Power1');
        D = -curvefit.b;

        [pcData, ~, pc] = fragility(D, weakness);

        delta(i,j) = abs(pcData-pc);
        rmse(i,j) = sqrt(mean((curvefit(sizes)-nSizes).^2));
        
    end, clear i,
end, clear j, 

[I,J] = find(delta==min(delta));

for n=1:length(weak)
    exclude = I(n);
    weakness = weak{J(n)};
    
    curvefit = fit(sizes(exclude:end),nSizes(exclude:end),'Power1');
    D = -curvefit.b;

    [pcData, ~, pc] = fragility(D, weakness);
    
    figure(n), clf
    loglog(sizes,nSizes,'+', 'DisplayName', 'Data')
    hold on
    loglog(sizes, curvefit(sizes), 'DisplayName', ['ax^b, a = ' num2str(curvefit.a) ', b = ' num2str(curvefit.b)])
    hold off
    legend('Location', 'best')
    title([weakness ' - p_c = ' num2str(pcData)] , [ 'Delta p_c = ' num2str(delta(I(n),J(n))) ' - Excluded the first ' num2str(I(n)-1) ' data points'])
    xlabel('Size (µm)')
    ylabel('Normalised number of particles')
end, clear n,

% x = fminbnd(@powerLawDelta,-1e10,1e10);
% figure, clf,
% loglog(sizes,nSizes,'+')
% hold on
% loglog(sizes, x.*sizes.^(-log(8*0.4901)/log(2)))
% hold off


function [pcData, pn, pc] = fragility(D, weakness)
pn1 = 0:1e-6:1;
if strcmpi(weakness, "pillar")
    pn = pn1.^4.*(3.*pn1.^4-8.*pn1.^3+4.*pn1.^2+2);
%     pcData = 2^D/8;
elseif strcmpi(weakness, "plane")
%     pn = 3.*pn1.^8 - 32.*pn1.^7 + 88.*pn1.^6 - 69.*pn1.^5 + 38.*pn1.^4;
    pn = pn1.^8 + 8.*pn1.^7.*(1-pn1) + 28.*pn1.^6.*(1-pn1).^2 ...
         + 56.*pn1.^5.*(1-pn1).^3 + 38.*pn1.^4.*(1-pn1).^4;
%     pcData = sqrt(4^D/64);
end
pcData = 2^D/8;
pc = pn(abs(pn-pn1) == min(abs(pn(2:end-1)-pn1(2:end-1))));
end

% function [res] = powerLawDelta(a)
% global sizes nSizes
% pc = 0.4901;
% x = sizes(2:end);
% y = nSizes(2:end);
% b = log(8*pc)/log(2);
% powLaw = a.*x.^(-b);
% res = sqrt(mean((powLaw - y).^2));
% end
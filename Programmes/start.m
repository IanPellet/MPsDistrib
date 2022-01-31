% add curent folder and its subdirectories 
currDir = genpath(pwd);
addpath(currDir) 

% add path to data
addpath('../Data/')

% save added paths
savepath pathdef.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Add the paths of the different toolboxes
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Patrick Marchesiello and Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated    10-Sep-2006 by Pierrick Penven
%  Updated    22-Sep-2006 by Pierrick Penven (64 bits test)
%  Updated    24-Oct-2006 by Pierrick Penven (mask added)
%  Updated    16-jan-2007 by Pierrick Penven (quikscat added)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath('Programmes')
%addpath('Dessins')
% addpath('/data_model/OUTPUTS/RHOMA/RHOMA400m/RHOMA003')
% addpath('/data_model/OUTPUTS/glxl/2018/20180930/')

% addpath /home/mlaval/ProgrammesTraitementDonnees/ProgCommun/smartquad
% addpath /home/mlaval/ProgrammesTraitementDonnees/ProgCommun/DessinRobert
% addpath /home/mlaval/ProgrammesTraitementDonnees/ProgCommun/TraitementTransect
% addpath /home/mlaval/ProgrammesTraitementDonnees/ProgCommun
% addpath /home/mlaval/ProgrammesTraitementDonnees/Foreman
% 
% addpath('/home/mlaval/ProgrammesTraitementDonnees/DiversPetitProg')

% addpath /media/ian/Elements/Ian_Plastique/LibrairiesMatlab/ProgrammesTraitementDonnees/ProgCommun/smartquad
% addpath /media/ian/Elements/Ian_Plastique/LibrairiesMatlab/ProgrammesTraitementDonnees/ProgCommun/DessinRobert
% addpath /media/ian/Elements/Ian_Plastique/LibrairiesMatlab/ProgrammesTraitementDonnees/ProgCommun/TraitementTransect
% addpath /media/ian/Elements/Ian_Plastique/LibrairiesMatlab/ProgrammesTraitementDonnees/ProgCommun
% addpath /media/ian/Elements/Ian_Plastique/LibrairiesMatlab/ProgrammesTraitementDonnees/Foreman
% addpath /media/ian/Elements/Ian_Plastique/LibrairiesMatlab/loadods

% addpath('/media/ian/Elements/Ian_Plastique/LibrairiesMatlab/ProgrammesTraitementDonnees/DiversPetitProg')



disp(['Add the paths of the different toolboxes'])
NetCdfPath=...
'./LibrairiesMatlab/LibrairiesNetcdf/';
%-------------------------------------------------------
%
% Get the path to the mexcdf (it depends on the architecture)
% Comment  all these lines if you don't want to pass in these tests
!uname -m > .mysystem
fid=fopen('.mysystem');
mysystem=fscanf(fid,'%s');

% if ( strcmp(mysystem(end-1:end),'86') )
 mysystem2='32';
% elseif ( strcmp(mysystem(end-1:end),'64') )
%  mysystem2='64';
% end

fclose(fid); 
matversion=version('-release');
myversion=str2num(matversion(1:2));
!rm -f .mysystem
disp(['Arch : ',mysystem,' - Matlab version : ',matversion])


if ((myversion > 13)    )
  disp(['Use of mexnc and loaddap in ',mysystem2,' bits.'])
  addpath([NetCdfPath,'mexcdf/mexnc'])   % 32 and 64 bits version of mexnc 
%
% - If these directories are already in your matlab native path, 
% you can comment these lines
addpath([NetCdfPath,'mexcdf/netcdf_toolbox/netcdf'])
addpath([NetCdfPath,'mexcdf/netcdf_toolbox/netcdf/ncsource'])
addpath([NetCdfPath,'mexcdf/netcdf_toolbox/netcdf/nctype'])
addpath([NetCdfPath,'mexcdf/netcdf_toolbox/netcdf/ncutility'])
%
%-------------------------------------------------------
elseif (myversion <= 13)
  disp('Use of mex60 and loaddap in 32 bits.')
  addpath([NetCdfPath,'mex60'])         % Older/32 bits version of mexcdf

% - If these directories are already in your matlab native path, 
% you can comment these lines
% - In this case, if problems with subsrefs.m ans subsasign.m,
% it is because there is a conflict with another native subs.m routines in the
% symbolic native toolbox
addpath([NetCdfPath,'netcdf_matlab_60'])
addpath([NetCdfPath,'netcdf_matlab_60/ncfiles'])
addpath([NetCdfPath,'netcdf_matlab_60/nctype'])
addpath([NetCdfPath,'netcdf_matlab_60/ncutility'])

else
  disp(['Arch : ',mysystem,...
       ' you should provide the paths of your own loaddap and mexcdf directories'])
end

clear
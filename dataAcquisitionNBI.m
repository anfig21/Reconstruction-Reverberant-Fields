function Data = dataAcquisitionNBI(Data)
%Data = dataAcquisitionNBI(Data) Read eigenmike dataset.
%   Input:
%       - Data      : data structure.
%
% Author: Antonio Figueroa Durán
% Date: August 2023

%% DATASETS
em32 = 'rir_NBI_em32.h5';
line = 'rir_NBI_line.h5';

%% SOURCE
Data.Source.pos = h5read(em32,'/source/position');

%% EIGENMIKES
Toff = 2;
Data.Mic.pos = nan(Data.Mic.nM,3,Data.Mic.nArrays);
Data.Mic.h = nan(Toff*Data.Fs,Data.Mic.nM,Data.Mic.nArrays);
Data.Mic.R0 = nan(Data.Mic.nArrays,3);
for ii = 1:Data.Mic.nArrays
    Data.Mic.pos(:,:,ii) = h5read(em32,['/eigenmike' num2str(Data.Mic.IdxArrays(ii),'%02d') '/position_microphones']);
    Data.Mic.h(:,:,ii) = h5read(em32,['/eigenmike' num2str(Data.Mic.IdxArrays(ii),'%02d') '/total_impulse_response']);
    Data.Mic.R0(ii,:) = h5read(em32,['/eigenmike' num2str(Data.Mic.IdxArrays(ii),'%02d') '/position_centre']);
end

% Mic Orientations
Data.Mic.Omega = nan(Data.Mic.nM,2,Data.Mic.nArrays);
for ii = 1:Data.Mic.nArrays
    posLocal = Data.Mic.pos(:,:,ii)-Data.Mic.R0(ii,:);
    [~,Data.Mic.Omega(:,1,ii),Data.Mic.Omega(:,2,ii)] = cart2sph2(posLocal(:,1),posLocal(:,2),posLocal(:,3));
end

%% REFERENCE LINE
% Line 1
Data.Line1.pos = h5read(line,'/line01/posRIR');
Data.Line1.h = h5read(line,'/line01/impulse_response');

% Line 2
Data.Line2.pos = h5read(line,'/line02/posRIR');
Data.Line2.h = h5read(line,'/line02/impulse_response');

% Line 3
Data.Line3.pos = h5read(line,'/line03/posRIR');
Data.Line3.h = h5read(line,'/line03/impulse_response');

disp('Reading data... OK')

end

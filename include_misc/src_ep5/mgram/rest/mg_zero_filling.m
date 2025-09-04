function images = mg_zero_filling(images, zero_params)

% Version: Maximilian Gram, 21.03.2024

%% single image case
if ndims(images)==2
    temp = images; clear images;
    images(1,:,:) = temp(:,:);
    clear temp;
end

%% check matrix size
[Nimages, Ny, Nx] = size(images);

%% create tukey window
N     = max([Ny,Nx]);
tukey = tukeywin(N, zero_params.radius);
tukey(1:round(numel(tukey)/2)) = [];
rq    = 0 : 0.1 : numel(tukey)-1;
tukey = interp1( (0:numel(tukey)-1), tukey', rq);
[x,y] = meshgrid(1:N, 1:N);
r     = round(sqrt((x-floor(N/2+1)).^2+(y-floor(N/2+1)).^2),1);

tukey2D = zeros(N, N);
for y=1:N
for x=1:N
    if r(y,x) <= max(rq)
        ind          = find( abs(rq-r(y,x)) == min(abs(rq-r(y,x))) );
        tukey2D(y,x) = tukey(ind(1));
    end
end
end
clear tukey x y r rq ind;

if Ny ~= Nx
    [X,  Y]  = meshgrid(1:N,  1:N);
    [Xq, Yq] = meshgrid(1:Nx, 1:Ny);
    X        = X  / max(max(X));
    Y        = Y  / max(max(Y));
    Xq       = Xq / max(max(Xq));
    Yq       = Yq / max(max(Yq));
    tukey2D  = interp2(X, Y, tukey2D, Xq, Yq);
    tukey2D(isnan(tukey2D)) = 0;
end
clear N X Y Xq Yq;

tukey3D(1,:,:) = tukey2D(:,:);
tukey3D = repmat(tukey3D, [Nimages, 1, 1]);

%% ksapce filter: tukey multiplication
kspace = image2kspace(images) .* tukey3D;

%% zero interpolation filling
Nx0     = floor(zero_params.factor*Nx);
Ny0     = floor(zero_params.factor*Ny);
d1x     = floor((Nx0-Nx)/2+1);
d2x     = d1x+Nx-1;
d1y     = floor((Ny0-Ny)/2+1);
d2y     = d1y+Ny-1;
kspace0 = zeros(Nimages, Ny0, Nx0);
kspace0(:,d1y:d2y, d1x:d2x) = kspace(:,:,:);

%% inverse fourier transform
images = kspace2image(kspace0);

%% single image case
images = squeeze(images);

end


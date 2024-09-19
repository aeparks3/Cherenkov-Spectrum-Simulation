function simulationCode(wavelength) 
    MCmatlab.closeMCmatlabFigures();
    model = MCmatlab.model;

    model.G.nx                = 461; % Number of bins in the x direction
    model.G.ny                = 461; % Number of bins in the y direction
    model.G.nz                = 301; % Number of bins in the z direction
    model.G.Lx                = 20; % [cm] x size of simulation cuboid
    model.G.Ly                = 20; % [cm] y size of simulation cuboid
    model.G.Lz                = 30.0; % [cm] z size of simulation cuboid
    
    model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
    model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
    model = plot(model,'G');

    %%%%%%%% 

    model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
    model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
    model.MC.boundaryType             = 2; % 2: Top cuboid boundary only is escaping
    model.MC.nPhotonsRequested        = round((1/wavelength)^2 * (1e10)); % # of photons to launch, variable for beam energy and field size
    model.MC.wavelength = wavelength;     
    [X,Y,Z,lambda] = ndgrid(model.G.x, model.G.y, model.G.z, model.MC.wavelength);
    model.MC.lightSource.sourceType   = 0; 

    myArray = h5read("10x10.h5", '/array_data'); %Replace with 3D matrix from TPS 
    myArray2 = replace_nan_with_zero(myArray);
    nx = 461; 
    ny = 461;  
    nz = 301; 
    xval = double(1):double(nx); 
    yval = double(1):double(ny); 
    zval = double(1):double(nz); 
    
    [tempx, tempy, tempz] = ndgrid(model.G.x, model.G.y, model.G.z);

    F = griddedInterpolant(tempx, tempy, tempz, myArray2);
    
    % Create query points using ndgrid based on the model grid
    [queryX, queryY, queryZ] = ndgrid(model.G.x, model.G.y, model.G.z);
    
    % Evaluate the interpolant at the query points
    model.MC.sourceDistribution = F(queryX, queryY, queryZ).*(Z>3);


    model = runMonteCarlo(model);
    model = plot(model,'MC');
end

function array_without_nan = replace_nan_with_zero(input_array)
    % Replace NaN values with 0 in the input array
    array_without_nan = input_array;
    nan_indices = isnan(array_without_nan);
    array_without_nan(nan_indices) = 0;
end

%% Geometry function(s) 
function M = geometryDefinition(X,Y,Z,parameters)
  zSurface = 3; 
  M = ones(size(X)); % Air
  M(Z > zSurface) = 2; % "Standard" tissue
end

%% Media Properties function 
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  j=1;
  mediaProperties(j).name  = 'air';
  mediaProperties(j).mua   = 1e-25; % Absorption coefficient [cm^-1]
  mediaProperties(j).mus   = 1e-8; % Scattering coefficient [cm^-1]
  mediaProperties(j).g     = 1; % Henyey-Greenstein scattering anisotropy

  j=2;
  mediaProperties(j).name  = 'breast tissue';
  mediaProperties(j).mua = @func_mua; %absorption coefficient 
  function mua = func_mua(wavelength)
    B = 0.00976; % Blood content
    S = 0.6802; % Blood oxygen saturation
    W = 0.16025; % Water content
    M = 0; % Melanin content
    F = 0.62325; % Fat content 
    mua = calc_mua(wavelength,S,B,W,F,M); % Jacques "Optical properties of biological tissues: a review" eq. 12
  end

  mediaProperties(j).mus = @func_mus;%scattering coefficient 
  function mus = func_mus(wavelength)
    aPrime = 18.7; 
    fRay = 0.288; % Fraction of scattering due to Rayleigh scattering
    bMie = 0.685; % Scattering power for Mie scattering
    g = 0.9; % Scattering anisotropy
    mus = calc_mus(wavelength,aPrime,fRay,bMie,g); % Jacques "Optical properties of biological tissues: a review" eq. 2
  end
  mediaProperties(j).g     = 0.9; % Henyey-Greenstein scattering anisotropy
end

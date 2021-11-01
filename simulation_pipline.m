
%% Recurrent model


%parameter
tauS = 0.8;
nreps=15;
tauI = 0.05; tauE = 0.05;
tau=0.1;
recurrent_connectivity=0.1;
claustrum2otherarea_connectivity=0.0;
other2claustrum_connectivity=0.0; 
otherarea_recurrent_connectivity=0.0;

gInh=3.1;
g=1;

activation_power=g/2;
halo_power=0; 
rstore=Recurent_simulation_function(tau,tauE,tauS,tauI,halo_power,g,gInh,nreps,...
    activation_power,recurrent_connectivity,otherarea_recurrent_connectivity,...
    claustrum2otherarea_connectivity,other2claustrum_connectivity);

%inhibition
halo_power=-5*g;
rstore_inhibition=Recurent_simulation_function(tau,tauE,tauS,tauI,halo_power,g,gInh,nreps,...
    activation_power,recurrent_connectivity,otherarea_recurrent_connectivity,...
    claustrum2otherarea_connectivity,other2claustrum_connectivity);


%% Ext input model

tauS = 0.8;
nreps=15;
tauI = 0.05; tauE = 0.05;
tau=0.1;
recurrent_connectivity=0.0001;
claustrum2otherarea_connectivity=0.01;
other2claustrum_connectivity=0.05; 
otherarea_recurrent_connectivity=0.1;

gInh=4.5;
g=1;

activation_power=g/2;

halo_power=0; 
rstore=Recurent_simulation_function(tau,tauE,tauS,tauI,halo_power,g,gInh,nreps,...
    activation_power,recurrent_connectivity,otherarea_recurrent_connectivity,...
    claustrum2otherarea_connectivity,other2claustrum_connectivity);

%inhibition
halo_power=-5*g;
rstore_inhibition=Recurent_simulation_function(tau,tauE,tauS,tauI,halo_power,g,gInh,nreps,...
    activation_power,recurrent_connectivity,otherarea_recurrent_connectivity,...
    claustrum2otherarea_connectivity,other2claustrum_connectivity);
## Simulation test

Application of local KDE classification of referfence and admixed simulated genotypes.

### Evolutionary tree:

The R package `ms` was used to simulate three populations coalescing along the follwing scheme:

`
         |                   
         |                   
         |                   
        / \        <--- 0.2  
       /   \                 
      /     \                
     /       \               
    /\ <- - - \ - - - - 0.15  
   /  \        \             
  /    \        \            
 P1    P2       P3           
                =X           

`

The number of initial haplotypes was set to 1e4, the effective size of size of populations P1, P2 and P3 was set to 
600, 600 and 300. Mutation rate was kept constant at 1e9 per generation.

### Admixture simulations.

An in-house script was used to generate three admixed populations as defined 
by their admixture proportions:

`
    P1 P2 P3
A1: 300 0   300
A2: 400 200 0
A3: 200 200 200
`

Admixed populations were allowed 17 generations of random crossing followed
by 8 of selfing. 

*Note* the software used to simulate an admixture framework with optional
reproduction system and generation number is under development. Contact us 
for details and access to the latest version.

## Running locally.

To run the kde classifier on this data, the following modules need to be installed:

`
R: version 3.4.3
python: version 3.4.3
plink: version 1.90b3v
sNMF_CL: version 1.2
`

Run the KDE and Ideogram scripts from within the directory, using the 
[example commands] (https://github.com/SantosJGND/Galaxy_KDE_classifier/blob/master/Simulation_related/Commands_sim.txt).




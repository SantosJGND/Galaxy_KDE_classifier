import numpy as np
from scipy.stats import beta
from sklearn.decomposition import PCA

def generate_vectors_Beta(L,n,rangeA= [1,2],rangeB = [.1,4],steps= 20,n_comp = 100):
    '''
    Starshape structure.
    '''
    # Simulate frequency vectors. 
    # We must first define the number of populations, the length of the haplotypes desired, and their respective population sizes

    import itertools as it

    # Vary a (beta distribution parameter).
    a_range= np.linspace(rangeA[0],rangeA[1],steps)
    a_set= [i for i in a_range for _ in range(n)]

    # vary b.
    b_range= np.linspace(rangeB[0],rangeB[1],steps)
    b_set= [i for i in b_range for _ in range(n)]

    ## length of haplotypes to extract.
    L_set= [L] * n * steps


    background_1= np.array([a_set,b_set,L_set]).T

    vector_lib= []
    for k in range(background_1.shape[0]):

        probs= beta.rvs(background_1[k,0], background_1[k,1], size=int(background_1[k,2]))
        probs[(probs > 1)]= 1


        vector_lib.append(probs)

    vector_lib= np.array(vector_lib)

    n_comp = 100

    
    return vector_lib

    
def generate_Branches_Beta(Nbranches,density,L,n,rangeA= [1,2],rangeB = [.1,4],steps= 20,n_comp = 100):
    
    vector_lib= generate_vectors_Beta(L,n,rangeA,rangeB,steps,n_comp)
    
    pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized').fit(vector_lib)
    features = pca.transform(vector_lib)# * pca.explained_variance_ratio_
    
    Iter= density
    target= [0,1]
    stairs= Nbranches

    MRCA= np.random.choice(range(vector_lib.shape[0]),1)
    calypso= []
    feat= []

    for inter in range(stairs):
        Pair= np.random.choice(range(vector_lib.shape[0]),2,replace= False)
        Pair[1]= MRCA
        print(Pair)

        coords= features[Pair,:]

        vector2= coords[target[1]] - coords[target[0]]
        for angle in np.linspace(-20,20,Iter):
            new_guy = coords[target[0]] + [angle / 10 * x for x in vector2]

            feat.append(new_guy)

            new_guy= pca.inverse_transform(new_guy)
            new_guy[new_guy < 0]= 1e-4
            new_guy[new_guy > 1]= 1

            calypso.append(new_guy)

    features= np.array(feat)
    vector_lib= np.array(calypso)

    return features, vector_lib


from utils import *
import pandas as pd
import sys
import os
import lzma

#demog_file = sys.argv[1]#'../demographies/generic_models.csv'
out_path = sys.argv[1]#'../results/generic/'
replicate = str(sys.argv[2]) #replicate name


#demographies = pd.read_csv(demog_file)
#models = np.array(demographies.columns[:-1])
#print(models)
Nstart = 10000

recregion =[fp11.Region(i,i+1,1, coupled=True) for i in range(21)]

################## simulate neutral ############################

sregion,nregion,rates = regions_dfe(species='test',neutral=True)

# constant size for 10 N generations
burnin=np.array([Nstart]*int(10*Nstart),dtype=np.uint32)

mypop =  fp11.SlocusPop(Nstart)

#prepare random number gernerator
rng2 = fp11.GSLrng(np.random.randint(420000))


[print(i) for i in recregion]
p = {'nregions':nregion,
'sregions': sregion,
'recregions':recregion,
'rates':rates,
'demography':burnin,
}
params = fp11.model_params.SlocusParams(**p)

burn_rec = track_burnin()
# simulate until equilibrium
wf.evolve(rng2, mypop,params,burn_rec)




#ppop = pickle.dumps(mypop,-1)
# pickle equilibirum population
burnin_name = out_path + "burnins/burnin_neut_%s.lzma9" % replicate
with lzma.open(burnin_name, "wb", preset=9) as f:
    pickle.dump(mypop, f, -1)

print('burnin done')
print('Generation',mypop.generation)

model_path = out_path + 'maize/'
if not os.path.exists(model_path):
    os.makedirs(model_path)
demog = np.array([Nstart]*int(0.1*Nstart),dtype=np.uint32)

#Unpickle to create a new pop:
#pop2 = pickle.loads(ppop)

with lzma.open(burnin_name, 'rb') as f:
    pop2 = pickle.load(f)
print(mypop==pop2)
print(pop2.generation)

p['demography'] = demog
params = fp11.model_params.SlocusParams(**p)

# add recorder that records pi, singletons and tajimas D
set_gen = (10 * Nstart) + 200  # adjust generation labels without burnin and start
rec1 = neutral_div(set_gen, final=pop2.generation + len(demog) + 200, Nstart=Nstart)

wf.evolve(rng2, pop2, params, rec1)
print('Generation', pop2.generation)

# write output
write_output(rec1, model_path, 'neutral', replicate)



################## simulate BGS ############################

sregion,nregion,rates = regions_dfe(species='test',neutral=False)

# constant size for 10 N generations
burnin=np.array([Nstart]*int(10*Nstart),dtype=np.uint32)

mypop =  fp11.SlocusPop(Nstart)

#prepare random number gernerator
rng2 = fp11.GSLrng(np.random.randint(420000))


p = {'nregions':nregion,
'sregions': sregion,
'recregions':recregion,
'rates':rates,
'demography':burnin,
}

params = fp11.model_params.SlocusParams(**p)

burn_rec = track_burnin()
print('sregions')
[print(i) for i in p['sregions']]
print('nregions')
[print(i) for i in p['nregions']]
print('Recomb region')
[print(i) for i in p['recregions']]
print(p['rates'])

# simulate until equilibrium


wf.evolve(rng2, mypop,params,burn_rec)
print('burnin done')
print('Generation',mypop.generation)


#ppop = pickle.dumps(mypop,-1)
# pickle equilibirum population
burnin_name = out_path + "burnins/burnin_bgs_%s.lzma9" % replicate
with lzma.open(burnin_name, "wb", preset=9) as f:
    pickle.dump(mypop, f, -1)

model_path = out_path + 'maize/'
if not os.path.exists(model_path):
    os.makedirs(model_path)
demog = np.array([Nstart]*int(0.1*Nstart),dtype=np.uint32)

with lzma.open(burnin_name, 'rb') as f:
    pop2 = pickle.load(f)
print(pop2.generation)
#Unpickle to create a new pop:
#    pop2 = pickle.loads(ppop)
#    print(mypop==pop2)

p['demography'] = demog
params = fp11.model_params.SlocusParams(**p)

# add recorder that records pi, singletons and tajimas D
set_gen = (10 * Nstart) + 200  # adjust generation labels without burnin and start
rec1 = neutral_div(set_gen, final=pop2.generation + len(demog) + 200, Nstart=Nstart)

wf.evolve(rng2, pop2, params, rec1)
print('Generation', pop2.generation)

# write output
write_output(rec1, model_path, 'bgs', replicate)
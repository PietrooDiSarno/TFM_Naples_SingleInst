import PMOT as pm
import sys
import spiceypy as spice

if len(sys.argv) > 1:
    print(sys.argv[1])
    pathToData = sys.argv[1]
    caseName = sys.argv[2]

    GApop = input('Population: ')
    GAgen = input('Generations: ')
    nm = input('Mutants: ')
    nd = input('Descendants: ')
    ne = input('Elites: ')
    ncm = input('Can mutate: ')
    savingIte = input('Start saving at generation: ')
    loadingIte = input('Start loading at generation: ')
    info = input('Info: ')
    nAgents = input('Number of agents: ')
else:
    GApop = input('Population: ')
    GAgen = input('Generations: ')
    nm = input('Mutants: ')
    nd = input('Descendants: ')
    ne = input('Elites: ')
    ncm = input('Can mutate: ')
    savingIte = input('Start saving at generation: ')
    loadingIte = input('Start loading at generation: ')
    info = input('Info: ')
    nAgents = input('Number of agents: ')
    pathToData = 'JUICE_multiAgent'
    caseName = 'JUICE'

ads = pm.agentDataSharing(True, pathToData, caseName)
lastAgentFound = ads.getLastKnownAgent()
print('lastAgentFound=', lastAgentFound)
arguments = [GApop, GAgen, nm, nd, ne, ncm, savingIte, loadingIte, info]
firstAgent = lastAgentFound + 1
if True:
    sp = pm.spawnAgents('JUICE_schedule_spawn', range(firstAgent, firstAgent + int(nAgents)), None, pathToData, caseName,
                        arguments=arguments)
    sp.spawn()

agent, generation, fit, ind = ads.findBest()
print('Best result obtained: agent =', agent, 'generation =', generation, 'fit =', fit, 'ind =', ind)
ads.plotAllHistory()

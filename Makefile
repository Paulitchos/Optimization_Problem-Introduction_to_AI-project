evolutionary:
	gcc ./EvolutionaryAlgorithm/*.c -o ./dist/evolutionary -lm && ./dist/evolutionary

evolutionaryrefe:
	gcc ./EvolutionaryAlgorithm/refe/*.c -o ./dist/evolutionaryrefe -lm && ./dist/evolutionaryrefe

localsearch:
	gcc ./LocalSearchAlgorithm/*.c -o ./dist/localsearch -lm && ./dist/localsearch

hybrid:
	gcc ./HybridMethod/*.c -o ./dist/hybrid -lm && ./dist/hybrid
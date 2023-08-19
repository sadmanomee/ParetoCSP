import math

from pymoo.core.infill import InfillCriterion


class Mating(InfillCriterion):

    def __init__(self,
                 selection,
                 crossover,
                 mutation,
                 **kwargs):

        super().__init__(**kwargs)
        self.selection = selection
        self.crossover = crossover
        self.mutation = mutation

    def _do(self, problem, pop, n_offsprings, parents=None, **kwargs):
        #print('entering do method of core.Mating ...........................')

        # how many parents need to be select for the mating - depending on number of offsprings remaining
        n_matings = math.ceil(n_offsprings / self.crossover.n_offsprings)
        #print(f'n_offsprings: {n_offsprings} ***********************************')
        #print(f'Number of matings: {n_matings} [each mating includes 2 parents by default] *******************************')

        # if the parents for the mating are not provided directly - usually selection will be used
        if parents is None:

            # select the parents for the mating - just an index array
            parents = self.selection(problem, pop, n_matings, n_parents=self.crossover.n_parents, **kwargs)
            #print(f'shape of parents: {parents.shape} ***********************')

        for i, _parents in enumerate(parents):
            parents[i][0].age += 1
            parents[i][1].age += 1
        # do the crossover using the parents index and the population - additional data provided if necessary
        off = self.crossover(problem, parents, **kwargs)
        #print(f'offspring shape after crossovers: {off.shape} **********************')
        # print('########## 13')
        # for i, _off in enumerate(off):
        #     print(off[i].age)
        # print('########## 14')

        # do the mutation on the offsprings created through crossover
        off = self.mutation(problem, off, **kwargs)
        #print(f'offspring shape after mutation: {off.shape} ***********************')

        #print('exiting do method of core.Mating ...........................')
        return off




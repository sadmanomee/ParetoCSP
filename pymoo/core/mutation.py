from copy import deepcopy

import numpy as np

from pymoo.core.operator import Operator
from pymoo.core.variable import Real, get


class Mutation(Operator):

    def __init__(self, prob=1.0, prob_var=None, **kwargs) -> None:
        super().__init__(**kwargs)
        self.prob = Real(prob, bounds=(0.7, 1.0), strict=(0.0, 1.0))
        self.prob_var = Real(prob_var, bounds=(0.0, 0.25), strict=(0.0, 1.0)) if prob_var is not None else None

    def do(self, problem, pop, inplace=True, **kwargs):

        # if not inplace copy the population first
        if not inplace:
            pop = deepcopy(pop)

        n_mut = len(pop)

        # get the variables to be mutated
        X = pop.get("X")

        # retrieve the mutation variables
        Xp = self._do(problem, X, **kwargs)

        # the likelihood for a mutation on the individuals
        prob = get(self.prob, size=n_mut)
        mut = np.random.random(size=n_mut) <= prob
        # print(mut, '$$$$$$$$$$$')

        # store the mutated individual back to the population
        pop[mut].set("X", Xp[mut])
        
        ### sadman
        # for i in range(len(pop)):
        #     if pop[i].age == 0:
        #         print('x' * 500)
        
        for i in range(len(pop)):
            if mut[i] == True and pop[i].age != 1: 
                pop[i].age = pop[i].age - 1
            
            if int(pop[i].age) == 0:
                print('\n\n\n')
                print('o' * 100)
                print('\n\n\n')
                

        return pop

    def _do(self, problem, X, **kwargs):
        return X

    def get_prob_var(self, problem, **kwargs):
        prob_var = self.prob_var if self.prob_var is not None else min(0.5, 1 / problem.n_var)
        return get(prob_var, **kwargs)

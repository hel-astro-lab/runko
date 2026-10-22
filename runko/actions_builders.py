# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

# -*- coding: utf-8 -*-


from itertools import chain
from types import MethodType
import runko_cpp_bindings.actions as actions

quote = actions.quote
sequence = actions.sequence
parallel = actions.parallel


class meta_actions_builder:
    def __init__(self, car, method_name: str, prog_wrapper=lambda x: (quote, x)):
        self._car = car
        self._given_progs = []
        self._prog_wrapper = prog_wrapper

        def method(*vargs):
            if len(vargs) > 0:
                self._given_progs.append(vargs)

            return self


        setattr(self, method_name, method)


    def get(self):
        # Remove potentially nested builders.
        is_mab = lambda x: type(x[0]) == meta_actions_builder
        nested_removed = (x[0].get() if is_mab(x) else x for x in self._given_progs)
        return tuple(chain((self._car,), map(self._prog_wrapper, nested_removed)))




def sequence_builder():
    return meta_actions_builder(sequence, "then")


def parallel_builder():
    return meta_actions_builder(parallel, "add")

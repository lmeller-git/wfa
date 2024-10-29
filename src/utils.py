from time import perf_counter
from typing import Any


def timeit(func: object) -> object:

    def inner(*args, **kwargs) -> Any:
        start = perf_counter()
        returned_values = func(*args, **kwargs)
        end = perf_counter()
        print(f"time taken by function {func.__name__}: {end - start:.3f} s")
        return returned_values

    return inner

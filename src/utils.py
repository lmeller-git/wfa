from time import perf_counter


def timeit(func: object) -> None:

    def inner(*args, **kwargs) -> None:
        start = perf_counter()
        returned_values = func(*args, **kwargs)
        end = perf_counter()
        print(f"time taken by function {func.__name__}: {end - start:.3f}")
        return returned_values

    return inner

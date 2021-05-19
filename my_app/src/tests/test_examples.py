#!/usr/bin/env python
# Скрипт обрабатывает количество итераций и выполнение критериев остановки изменения
# Импортирование необходимые модули:
from __future__ import print_function
from pathlib import Path

import pytest

import topy


@pytest.mark.benchmark(
    group="param:filename", max_time=10, min_rounds=1,
)
@pytest.mark.parametrize(
    "filename", (str(filename) for filename in Path("examples").rglob("*.tpd"))
)
def test_optimise(filename, benchmark):
    # тип: (str) -> None
    # Оптимизация файла по адресу "filename"
    print("Optimizing file '%s'..." % filename)
    # Настройка топологии:
    t = topy.Topology()
    t.load_tpd_file(filename)
    t.set_top_params()
    benchmark.pedantic(topy.optimise, args=(t,), rounds=1, iterations=1)

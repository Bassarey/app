#!/usr/bin/env python

# Скрипт обрабатывает количество итераций и выполнение критериев остановки изменения
# Импортирование необходимых модулей:
from sys import argv
import topy


def optimise(fname):
    # Настройка
    t = topy.Topology()
    t.load_tpd_file(fname)
    t.set_top_params()
    topy.optimise(t)


if __name__ == '__main__':
    optimise(argv[1])


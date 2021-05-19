# Матрицы жесткости конечных элементов.
# Каталог 'data' для определения собственных конечных элементов
from os import path

from numpy import array, linspace, unique, sqrt, round, load
from numpy.linalg import eigvalsh

from .utils import get_logger
from .data.matlcons import _a, _nu, _E

logger = get_logger(__name__)

__all__ = ['Q4', 'Q5B',  'Q4a5B',  'Q4T',\
           'H8', 'H18B', 'H8T']

# Ошибки
MSG0 = 'Матрица жесткости из конечных элементов.'
MSG1 = 'Матрицы жесткости элементов не существуют. \ n Создано ... Пожалуйста, перезапустите \
твоя последняя попытка.'

# Путь к папке с данными:
pth = path.join(path.split(__file__)[0], 'data')

# 2D элементы
fname = path.join(pth, 'Q4bar.K')
try:
    Q4bar = load(fname)
except IOError:
    logger.info('It seems as though all or some of the element stiffness matrices')
    logger.info('do not exist. Creating them...')
    logger.info('This is usually only required once and may take a few minutes.')
    from topy.data import Q4bar_K
    Q4bar = load(fname)

# Матрица жесткости билинейного квадратного 4-узлового плоского напряженного элемента
fname = path.join(pth, 'Q4.K')
try:
    Q4 = load(fname)
except IOError:
    from topy.data import Q4_K
    Q4 = load(fname)

# Матрица жесткости квадратного 4-узлового плоского элемента напряжения '5-beta'
fname = path.join(pth, 'Q5B.K')
try:
    Q5B = load(fname)
except IOError:
    from topy.data import Q5B_K
    Q5B = load(fname)

# Матрица элемента, используемого в двухмерных тепловых задачах
fname = path.join(pth, 'Q4T.K')
try:
    Q4T = load(fname)
except IOError:
    from topy.data import Q4T_K
    Q4T = load(fname)

# Матрица жесткости квадратного элемента Q4a5B с четырьмя узлами 
# Этот элемент основан на предполагаемом элементе напряжения 5-бета для плоскости напряжение, но элементарные параметры вводятся и выбираются так, чтобы паразитные моды с нулевой энергией не вводятся, для которых расследование
# Необходимо количество характеристических уравнений матрицы элементарной жесткости.
# Толщина элемента установлена ​​= 1. Подробности см. В De Klerk and Groenwold.
# Символьное значение alpha_opt для гибки:
alpha2D = (2 * _a**2 * (1 - _nu) * (2 * _nu**2 - _nu + 1)) \
/ (3 * (_nu + 1) * _E**2)
Q4a5B = Q4 - alpha2D * _E * Q4bar  # матрица жесткости

# 3D элементы
# Матрица жесткости для трехлинейного трехмерного элемента шестигранник с 8 узлами
fname = path.join(pth, 'H8.K')
try:
    H8 = load(fname)
except IOError:
    from topy.data import H8_K
    H8 = load(fname)

# Матрица жесткости кубического 8-узлового элемента '18 -beta '
fname = path.join(pth, 'H18B.K')
try:
    H18B = load(fname)
except IOError:
    from topy.data import H18B_K
    H18B = load(fname)

# Матрица жесткости для трехлинейного трехмерного элемента шестигранник с 8 узлами для тепловые проблемы.
fname = path.join(pth, 'H8T.K')
try:
    H8T = load(fname)
except IOError:
    from topy.data import H8T_K
    H8T = load(fname)

# EOF elements.py

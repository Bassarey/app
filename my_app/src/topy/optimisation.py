from os import path, makedirs
from time import time

from numpy import array

from .utils import get_logger
from .visualisation import *
from .topology import *

logger = get_logger(__name__)


__all__ = ['optimise']

def optimise(topology, save=True, dir='./iterations'):
    # тип: (Topology, bool, str) -> None
    if not path.exists(dir):
        makedirs(dir)
    etas_avg = []

# Функция оптимизации
    def _optimise(t):
        t.fea()
        t.sens_analysis()
        t.filter_sens_sigmund()
        t.update_desvars_oc()
        # Ввод информации и создание изображения или геометрии
        if t.nelz:
            params = {
                'prefix': t.probname,
                'iternum': t.itercount,
                'time': 'none',
                'dir': dir
            }
            if save:
                create_3d_geom(t.desvars, **params)
        else:
            params = {
                'prefix': t.probname,
                'iternum': t.itercount,
                'time': 'none',
                'filetype': 'png',
                'dir': dir
            }
            if save:
                create_2d_imag(t.desvars, **params)

        
        str_ = '%4i  | %3.6e | %3.3f | %3.4e | %3.3f | %3.3f |  %1.3f  |  %3.3f '
        format_ = (t.itercount, t.objfval, t.desvars.mean(),\
            t.change, t.p, t.q, t.eta.mean(), t.svtfrac)
        logger.info(str_ % format_)
        # Составление списока средних этов
        etas_avg.append(t.eta.mean())


    # Создание начальной области проектирования
    logger.info('\n' + '='*80)
    # Запуск оптимизации
    str_ = '%5s | %11s | %5s | %10s | %5s | %5s | %7s | %5s '
    format_ = ('Iter', 'Obj. func.  ', 'Vol. ', 'Change    ', \
        'P_FAC', 'Q_FAC', 'Ave ETA', 'S-V frac.')
    logger.info(str_ % format_)
    logger.info('-'*80)
    ti = time()

    try:
        while topology.change > topology.chgstop:
            _optimise(topology)
    except AttributeError:
        for i in range(topology.numiter):
            _optimise(topology)
    te = time()

    # Ввод информации о соотношении сплошных и пустотных пространств:
    logger.info('\nSolid plus void to total elements fraction = %3.5f' %\
        (topology.svtfrac))
    # Отображение информации об итерациях

    logger.info('%d iterations took %3.3f minutes (%3.3f seconds/iteration)'\
        %(topology.itercount, (te - ti) / 60, (te - ti) / topology.itercount))
    logger.info('Average of all ETA\'s = %3.3f (average of all a\'s = %3.3f)' \
        % (array(etas_avg).mean(), 1/array(etas_avg).mean() - 1))




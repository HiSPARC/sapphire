"""Update or extend local data with latest data from public database

This package contains modules for updating local data:

:mod:`~sapphire.data.extend_local_data`
    add additional data

:mod:`~sapphire.data.update_local_data`
    bring already included data up to date

"""
from . import extend_local_data
from . import update_local_data


__all__ = ['extend_local_data',
           'update_local_data']

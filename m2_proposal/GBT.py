#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created initially on 8 July 2021 and modified in August 2022

Plots the visibilities at the GBT of either a list of objects on one date or a list of dates for one object.

@author: Emily Moravec
"""
from datetime import datetime
import pandas as pd
import pytz

from astropy.table import Table

import observing_prep as obs
# -------------------------------------------------------------------------------------------- #
# define dates of observation and relevant time arrays
# https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
tz = pytz.timezone("America/New_York")
# 23A
sem_23A = [datetime(2023, 2, 1, tzinfo=tz), 
        datetime(2023, 2, 15, tzinfo=tz), 
        datetime(2023, 3, 1, tzinfo=tz), 
        datetime(2023, 4, 1, tzinfo=tz)]

# 23B
sem_23B = [datetime(2023, 10, 15, tzinfo=tz), 
        datetime(2023, 11, 1, tzinfo=tz), 
        datetime(2023, 12, 1, tzinfo=tz), 
        datetime(2024, 1, 1, tzinfo=tz), 
        datetime(2024, 2, 1, tzinfo=tz)]

# 24A
sem_24A = [datetime(2024, 2, 1, tzinfo=tz), 
        datetime(2024, 2, 15, tzinfo=tz), 
        datetime(2024, 3, 1, tzinfo=tz), 
        datetime(2024, 4, 1, tzinfo=tz)]

# 23B + 24A
sem_23Bp24A = [datetime(2023, 10, 1, tzinfo=tz), 
        datetime(2023, 11, 1, tzinfo=tz), 
        datetime(2023, 12, 1, tzinfo=tz), 
        datetime(2024, 1, 1, tzinfo=tz), 
        datetime(2024, 2, 1, tzinfo=tz),
        datetime(2024, 3, 1, tzinfo=tz), 
        datetime(2024, 4, 1, tzinfo=tz)]

### Execute script
# for each day, plot the visibilities for a list of objects
obs.visibilities_loop_days_and_objects(table_path=location+'23B/23B_targets.txt',
                                        table=None,
                                        telescope='GBT',
                                        time_zone='ET',
                                        dates=sem_23B,
                                        semester='23B',
                                        plot_suffix='_test_targets_',
                                        file_save_location=location)

# for object, plot the visibilities for a list of dates
obs.visibilities_loop_days_for_one_object(table_path=location+'proposals/GBT/GBT22B/tables/GBT22B_targ.csv',
                                        table=None,
                                        telescope='GBT',
                                        time_zone='ET',
                                        dates=sem_23A,
                                        semester='23A',
                                        plot_suffix='_date_range_test',
                                        file_save_location=location+'proposals/GBT/')

""" Shorten worksheet name to <= 31 characters

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-10-10
:Copyright: 2019, Karr Lab
:License: MIT
"""

import obj_tables
import openpyxl


def transform(filename):
    # read
    wb = openpyxl.load_workbook(filename=filename)

    # transform
    ws = wb['!!Initial species concentrations']
    ws.title = '!!Init species concentrations'

    # save
    wb.save(filename)

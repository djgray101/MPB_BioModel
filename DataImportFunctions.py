import xlrd
import numpy as np
import pandas as pd
import random

def mpbdata(rowstart, rowend, colstart, colend, sheetnum, file_location):
    """
    :param rowstart: Start of the desired row in Excel data matrix
    :param rowend: End of the desired row in Excel data matrix
    :param sheetnum: Specifies the sheet number to pull the data from {1,2,3,4,5}-->{Volumes,Stems,Ratios,Greenattack,Infest ID}
    """
    #file_location = "C:/Users/djgra/Desktop/Thesis_MPB/MPB_Data/GridExperiments/Grid1_DataFiles.gdb/GridStems_DissolveBeetleJoin.xlsx"
    #file_location = "C:/Users/djgra/Desktop/Thesis_MPB/MPB_Data/GridExperiments/Grid1_DataFiles.gdb/GridStems_DissolveBeetleJoin_Experimenting.xlsx"
    workbook = xlrd.open_workbook(file_location)
    sheet = workbook.sheet_by_index(sheetnum)

    foo = [[sheet.cell_value(r, c) for r in range(rowstart, rowend + 1)] for c in range(colstart, colend + 1)]
    foo2 = np.array(foo).T.tolist()
    return foo2


def mpb_growthrates(dimension, periods, rowstart, rowend, column):
    # file_location = "C:/Users/djgra/Desktop/Thesis_MPB/MPB_Data/Hinton_Climate_Data/GrowthRates_All.xlsx"
    file_location = "C:/Users/djgra/Desktop/Thesis_MPB/MPB_Data/Hinton_Climate_Data/GrowthRates_All.xlsx"
    workbook = xlrd.open_workbook(file_location)
    sheet = workbook.sheet_by_index(0)

    value = [sheet.cell_value(r, column) for r in range(rowstart, rowend + 1)]
    foo2 = np.array(value).T.tolist()

    return foo2
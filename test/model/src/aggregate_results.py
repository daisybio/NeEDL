#/////////////////////////////////////////////////////////////////////////////#
#                                                                             #
#   Copyright (C) 2020 by David B. Blumenthal                                 #
#                                                                             #
#   This file is part of GenEpiSeeker.                                        #
#                                                                             #
#   GenEpiSeeker is free software: you can redistribute it and/or modify it   #
#   under the terms of the GNU General Public License as published by         #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   GenEpiSeeker is distributed in the hope that it will be useful,           #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with GenEpiSeeker. If not, see <http://www.gnu.org/licenses/>.      #
#                                                                             #
#/////////////////////////////////////////////////////////////////////////////#

##
# @file plot_p_values_cs_times.py
# @brief Plots mean p-values vs. mean evaluation times.
#
# @details 
# Usage: 
# ```sh
# $ python aggregate_results.py [--help] [-h] [--csv-directories <arg>] [--title <arg>] [--save-as <arg>] [--recursive]
# ```
#
# For more information, execute `$ python plot_p_values_cs_times.py --help`.

'''Plots mean p-values vs. mean evaluation times..'''

import os.path
import pandas as pd
import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams

print("\n**************************************************")
print("                   GenEpiSeeker                   ")
print("         aggregate results of model tests         ")
print("**************************************************")

parser = argparse.ArgumentParser(description="Aggregate results of model tests.")
parser.add_argument("--csv-directories", nargs="+", help="directories of the CSV files for which the plot should be generaed", required=True)
parser.add_argument("--save-as", help="name of CSV file where the aggregated results should be saved", required=True)
parser.add_argument("--recursive", action="store_true")
args = parser.parse_args()

res = []
corr_res = []
for csv_directory in args.csv_directories:
    if (args.recursive):
        res = res + [pd.read_csv(filename) for filename in glob.glob(csv_directory+"/**/*_AGG.csv", recursive=True)]
        corr_res = corr_res + [pd.read_csv(filename, index_col=0) for filename in glob.glob(csv_directory+"/**/*_CORR.csv", recursive=True)]
    else:
        res = res + [pd.read_csv(filename) for filename in glob.glob(csv_directory+"*_AGG.csv")]
        corr_res = corr_res + [pd.read_csv(filename, index_col=0) for filename in glob.glob(csv_directory+"*_CORR.csv")]
num_models = len(res[0].model)
models = [res[0].model[i] for i in range(num_models)]
corrs = pd.DataFrame(index=models,columns=models,data=0,dtype=float)
for cor in corr_res:
    corrs += cor
corrs /= float(len(corr_res))
fast_models = ["BAYESIAN_MODEL[]", "VARIANCE_MODEL[]", "PENETRANCE_MODEL[--score NLL]"]
shortcuts = {"BAYESIAN_MODEL[]" : "BM", "VARIANCE_MODEL[]" : "VM", "PENETRANCE_MODEL[--score NLL]" : "MLM", "PENETRANCE_MODEL[--score AIC]" : "PM_AIC",
                   "REGRESSION_MODEL[--score NLL]" : "QRMNLL", "REGRESSION_MODEL[--score AIC]" : "QRMAIC", "REGRESSION_MODEL[--score NLL-GAIN]" : "QRMNLLGAIN", "REGRESSION_MODEL[--score AIC-GAIN]" : "QRMAICGAIN"}
corrs_flat = pd.DataFrame(columns=["x","y","corr","x_name","y_name"],data=[[0,0,0,"",""] for i in range((num_models - 1) * (num_models - 1))])

index = 0
x = 0
for model_1 in models:
    if model_1 == "PENETRANCE_MODEL[--score AIC]":
        continue
    y = 0
    for model_2 in models:
        if model_2 == "PENETRANCE_MODEL[--score AIC]":
            continue
        corrs_flat.loc[index,"x"] = x
        corrs_flat.loc[index,"y"] = y
        corrs_flat.loc[index,"corr"] = corrs.loc[model_1,model_2]
        corrs_flat.loc[index,"x_name"] = shortcuts[model_1]
        corrs_flat.loc[index,"y_name"] = shortcuts[model_2]
        index += 1
        y += 1
    x += 1
corrs_flat.to_csv("{}_CORR.csv".format(args.save_as[:-4]), index=False)

all_p_values = [df.loc[i,"neg-log-p-value"] for df in res for i in range(num_models)]
max_real_p_values = max([val for val in all_p_values if val < 708.39])
p_values = {models[i] : [min(df.loc[i,"neg-log-p-value"], max_real_p_values + 2) for df in res] for i in range(num_models)}
best_p_values = [max([p_values[model][i] for model in models if model != "PENETRANCE_MODEL[--score AIC]"]) for i in range(len(res))]
best_fast_p_values = [max([p_values[model][i] for model in fast_models]) for i in range(len(res))]
num_best_p_values = {model : 0 for model in models}
num_best_fast_p_values = {model : 0 for model in models}
for i in range(len(res)):
    for model in models:
        if p_values[model][i] >= best_p_values[i] - 0.0001:
            num_best_p_values[model] += 1
    for model in fast_models:
        if p_values[model][i] >= best_fast_p_values[i] - 0.0001:
            num_best_fast_p_values[model] += 1
means_p_values = {model : np.mean(p_values[model]) for model in models}
std_p_values = {model : np.std(p_values[model]) for model in models}
evaluation_times = {models[i] : [df.loc[i,"evaluation-time"] for df in res] for i in range(num_models)}
means_evaluation_times = {model : np.mean(evaluation_times[model]) for model in models}
std_evaluation_times = {model : np.std(evaluation_times[model]) for model in models}
columns = []
dtypes = {}
for model in models:
    if model == "PENETRANCE_MODEL[--score AIC]":
        continue
    columns = columns + [shortcuts[model] + "_time_mean", shortcuts[model] + "_time_std", shortcuts[model] + "_p_value_mean", shortcuts[model] + "_p_value_std"]
dummy_data = [[0 for column in columns]]
df = pd.DataFrame(data=dummy_data, columns=columns, dtype=float)
df_counts = pd.DataFrame(columns=["best_p_value", "best_fast_p_value"],index=[shortcuts[model] for model in models if model != "PENETRANCE_MODEL[--score AIC]"])
for model in models:
    if model == "PENETRANCE_MODEL[--score AIC]":
        continue
    df[shortcuts[model] + "_time_mean"][0] = means_evaluation_times[model]
    df[shortcuts[model] + "_time_std"][0] = std_evaluation_times[model]
    df[shortcuts[model] + "_p_value_mean"][0] = means_p_values[model]
    df[shortcuts[model] + "_p_value_std"][0] = std_p_values[model]
    df_counts.loc[shortcuts[model],"best_p_value"] = num_best_p_values[model]
    df_counts.loc[shortcuts[model],"best_fast_p_value"] = num_best_fast_p_values[model]
df.to_csv(args.save_as, index=False, float_format="%g")
df_counts.to_csv("{}_COUNTS.csv".format(args.save_as[:-4]), index_label="model")
print("**************************************************\n")
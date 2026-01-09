# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: ACTHYDRO
#     language: python
#     name: python3
# ---

# %%
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objs as go

import baseflow
import baseflow.separation

import baseflow.utils


# %%
def plot_baseflow_separation(Q, baseflows, method_names=None, interactive=False):
    if method_names is None:
        method_names = list(baseflows.keys())

    if interactive:
        fig = make_subplots()
        fig.add_trace(go.Scatter(x=Q.index, y=Q.values, name='Streamflow', line=dict(width=2)))

        for method in method_names:
            if method in baseflows:
                fig.add_trace(go.Scatter(x=Q.index, y=baseflows[method], name=f'{method} Baseflow'))

        fig.update_layout(
            title='Baseflow Separation',
            xaxis_title='Date',
            yaxis_title='Flow',
            legend_title='Methods',
            hovermode="x unified",
            width=1100,
            height=600
        )

        return fig
    else:
        plt.figure(figsize=(12, 6))
        plt.plot(Q.index, Q.values, label='Observed Streamflow', alpha=0.7)

        for method in method_names:
            if method in baseflows:
                plt.plot(Q.index, baseflows[method], label=f'{method} Baseflow', alpha=0.7)

        plt.xlabel('Date')
        plt.ylabel('Flow')
        plt.title('Baseflow Separation')
        plt.legend()
        plt.grid(True)

        return plt.gcf()


# %%
# Read WOD data
df_wod = pd.read_csv('Data/SWIOID/WOD_410761.csv', index_col=0)
df_wod.index = pd.to_datetime(df_wod.index, format='mixed', dayfirst=True, errors='coerce')
df_wod = df_wod[~df_wod.index.isna()]
Q_wod = df_wod[df_wod.columns[0]]

# Read empirical flow data
df_empirical = pd.read_csv('Data/SWIOID/410761_empirical_flow_1992_2022.csv', index_col=0)
df_empirical.index = pd.to_datetime(df_empirical.index, format='mixed', dayfirst=True, errors='coerce')
df_empirical = df_empirical[~df_empirical.index.isna()]
Q_empirical = df_empirical[df_empirical.columns[0]]

Q_wod, Q_empirical

# %%
Q_wod.values, Q_empirical.values

# %%
BF_wod=baseflow.single(Q_wod)
BF_empirical=baseflow.single(Q_empirical)
BF_wod, BF_empirical

# %%
BF_wod_df = BF_wod[0]
BF_empirical_df = BF_empirical[0]

# %%
BF_wod_df

# %%
BF_empirical_df

# %%
merged_wod_df = pd.concat([Q_wod, BF_wod_df], axis=1)
merged_empirical_df = pd.concat([Q_empirical, BF_empirical_df], axis=1)

merged_wod_df.to_csv('temp/merged_baseflow_results_wod.csv')
merged_empirical_df.to_csv('temp/merged_baseflow_results_empirical.csv')


# %%
merged_wod_df

# %%
merged_empirical_df

# %%
import plotly.graph_objects as go

# Plot WOD
flow_wod = Q_wod
ukih_wod = BF_wod[0]['UKIH']

fig = go.Figure()

fig.add_trace(go.Scatter(
    x=flow_wod.index,
    y=flow_wod.values,
    mode='lines',
    name='WOD Flow',
    line=dict(color='black')
))

fig.add_trace(go.Scatter(
    x=ukih_wod.index,
    y=ukih_wod.values,
    mode='lines',
    name='WOD UKIH Baseflow',
    line=dict(color='royalblue')
))

# Plot Empirical
flow_empirical = Q_empirical
ukih_empirical = BF_empirical[0]['UKIH']

fig.add_trace(go.Scatter(
    x=flow_empirical.index,
    y=flow_empirical.values,
    mode='lines',
    name='Empirical Flow',
    line=dict(color='darkorange')
))

fig.add_trace(go.Scatter(
    x=ukih_empirical.index,
    y=ukih_empirical.values,
    mode='lines',
    name='Empirical UKIH Baseflow',
    line=dict(color='green')
))

fig.update_layout(
    title='WOD and Empirical Flow vs UKIH Baseflow',
    xaxis_title='Date',
    yaxis_title='Flow',
    legend=dict(x=0, y=1),
    height=800,
)

fig.show()


# %%
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Ensure indices are datetime and drop NaT
def ensure_dt(s):
    idx = s.index
    if not idx.inferred_type == 'datetime64':
        s.index = pd.to_datetime(idx, dayfirst=True, errors='coerce')
    s = s[~s.index.isna()]
    return s

ukih_series = ensure_dt(ukih_wod.copy())

# Set up WOD and empirical for superimposed KDEs
# ukih_series: baseflow (UKIH) for full WOD period
ukih_empirical = BF_empirical[0]['UKIH']
ukih_empirical = ensure_dt(ukih_empirical)
ukih_wod = ensure_dt(ukih_series)

months = [
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December"
]
ncols = 4
nrows = 3
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(22, 16))
axes = axes.flatten()

for m in range(1, 13):
    ax = axes[m-1]
    # Main WOD full record
    vals_wod = ukih_wod[ukih_wod.index.month == m]
    # Empirical (short period)
    vals_emp = ukih_empirical[ukih_empirical.index.month == m]

    # Plot WOD histogram
    if len(vals_wod) > 0:
        sns.histplot(vals_wod, bins=30, ax=ax, color='cornflowerblue', stat="density", edgecolor='white', alpha=0.4, label="WOD Hist")

    # KDE: WOD and empirical
    for vals, lab, col, lw in [
        (vals_wod, "WOD KDE", 'blue', 2),
        (vals_emp, "Empirical KDE", 'orange', 2),
    ]:
        x = vals.values
        x = x[~np.isnan(x)]
        if len(x) >= 2:
            sns.kdeplot(x, ax=ax, label=lab, color=col, linewidth=lw)

    # Median, P10, P90 for both WOD and Empirical
    stats_lines = []
    if len(vals_wod) > 0:
        median_val = vals_wod.median()
        p10 = vals_wod.quantile(0.10)
        p90 = vals_wod.quantile(0.90)
        stats_lines.append(
            f"WOD\n  Median: {median_val:.2f}\n  10th %: {p10:.2f}\n  90th %: {p90:.2f}"
        )
    if len(vals_emp) > 0:
        median_emp = vals_emp.median()
        p10_emp = vals_emp.quantile(0.10)
        p90_emp = vals_emp.quantile(0.90)
        stats_lines.append(
            f"Empirical\n  Median: {median_emp:.2f}\n  10th %: {p10_emp:.2f}\n  90th %: {p90_emp:.2f}"
        )
    if stats_lines:
        stats_text = "\n\n".join(stats_lines)
        bbox_props = dict(boxstyle="round,pad=0.4", facecolor="aliceblue", edgecolor="none")
        ax.text(
            0.98, 0.95,
            stats_text,
            ha='right', va='top',
            transform=ax.transAxes,
            fontsize=12,
            bbox=bbox_props,
            color='navy'
        )

    ax.set_title(months[m-1])
    ax.set_xlabel('UKIH Baseflow')
    ax.set_ylabel('Density')
    ax.legend(loc="upper left", fontsize=10)

plt.tight_layout()
plt.suptitle('Monthly Histograms and UKIH Baseflow KDEs', fontsize=22, y=1.02)
plt.show()


# %%
import plotly.graph_objs as go
import plotly.colors
import numpy as np
import calendar

# Ensure index is datetime
ukih_series = ukih_wod.copy()
if not ukih_series.index.inferred_type == 'datetime64':
    ukih_series.index = pd.to_datetime(ukih_series.index, dayfirst=True, errors='coerce')

# Drop rows where index couldn't be converted to datetime
ukih_series = ukih_series[~ukih_series.index.isna()]

months = [
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December"
]

# Use tab10 color palette from plotly
palette = plotly.colors.qualitative.Plotly if len(plotly.colors.qualitative.Plotly) >= 12 else plotly.colors.qualitative.Dark24

traces = []

for m in range(1, 13):
    month_vals = ukih_series[ukih_series.index.month == m]
    # Convert from ML/month to ML/day
    month_daily_vals = []
    for year in month_vals.index.year.unique():
        mask = (month_vals.index.month == m) & (month_vals.index.year == year)
        vals_this_year_month = month_vals[mask]
        if len(vals_this_year_month) > 0:
            days_in_month = calendar.monthrange(year, m)[1]
            daily_vals_this_month = vals_this_year_month / days_in_month
            month_daily_vals.append(daily_vals_this_month)
    if month_daily_vals:
        month_daily_vals = pd.concat(month_daily_vals)
    else:
        month_daily_vals = pd.Series([], dtype=float)
    # Add to traces if not empty
    if len(month_daily_vals) > 0:
        x = month_daily_vals.values
        x = x[~np.isnan(x)]  # Remove NaN, if any
        # KDE, use scipy.gaussian_kde for compatibility with plotly
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(x)
        # Set maximum x-axis extent to 200
        x_min = np.min(x)
        x_max = min(200, np.max(x))
        x_grid = np.linspace(x_min, x_max, 400)
        y_kde = kde(x_grid)
        traces.append(
            go.Scatter(
                x=x_grid, y=y_kde, 
                mode='lines',
                name=months[m-1],
                line=dict(color=palette[(m-1) % len(palette)], width=2),
                opacity=0.8
            )
        )

layout = go.Layout(
    title='Monthly KDEs of UKIH Baseflow WOD (ML/day)',
    xaxis=dict(title='UKIH Baseflow (ML/day)', range=[None, 200]),
    yaxis=dict(title='Density', range=[0, 0.15]),
    legend=dict(title='Month'),
    width=1000,
    height=600,
    margin=dict(l=80, r=250, t=80, b=80)
)

fig = go.Figure(data=traces, layout=layout)
fig.show()


# %%
import plotly.graph_objs as go
import plotly.colors
import numpy as np
import calendar

# Ensure index is datetime
ukih_series = ukih_empirical.copy()
if not ukih_series.index.inferred_type == 'datetime64':
    ukih_series.index = pd.to_datetime(ukih_series.index, dayfirst=True, errors='coerce')

# Drop rows where index couldn't be converted to datetime
ukih_series = ukih_series[~ukih_series.index.isna()]

months = [
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December"
]

# Use tab10 color palette from plotly
palette = plotly.colors.qualitative.Plotly if len(plotly.colors.qualitative.Plotly) >= 12 else plotly.colors.qualitative.Dark24

traces = []

for m in range(1, 13):
    month_vals = ukih_series[ukih_series.index.month == m]
    # Convert from ML/month to ML/day
    month_daily_vals = []
    for year in month_vals.index.year.unique():
        mask = (month_vals.index.month == m) & (month_vals.index.year == year)
        vals_this_year_month = month_vals[mask]
        if len(vals_this_year_month) > 0:
            days_in_month = calendar.monthrange(year, m)[1]
            daily_vals_this_month = vals_this_year_month / days_in_month
            month_daily_vals.append(daily_vals_this_month)
    if month_daily_vals:
        month_daily_vals = pd.concat(month_daily_vals)
    else:
        month_daily_vals = pd.Series([], dtype=float)
    # Add to traces if not empty
    if len(month_daily_vals) > 0:
        x = month_daily_vals.values
        x = x[~np.isnan(x)]  # Remove NaN, if any
        # KDE, use scipy.gaussian_kde for compatibility with plotly
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(x)
        # Set maximum x-axis extent to 200
        x_min = np.min(x)
        x_max = min(200, np.max(x))
        x_grid = np.linspace(x_min, x_max, 400)
        y_kde = kde(x_grid)
        traces.append(
            go.Scatter(
                x=x_grid, y=y_kde, 
                mode='lines',
                name=months[m-1],
                line=dict(color=palette[(m-1) % len(palette)], width=2),
                opacity=0.8
            )
        )

layout = go.Layout(
    title='Monthly KDEs of UKIH Baseflow EMPIRICAL (ML/day)',
    xaxis=dict(title='UKIH Baseflow (ML/day)', range=[None, 200]),
    yaxis=dict(title='Density', range=[0, 0.15]),
    legend=dict(title='Month'),
    width=1000,
    height=600,
    margin=dict(l=80, r=250, t=80, b=80)
)

fig = go.Figure(data=traces, layout=layout)
fig.show()


# %%
ukih_wod

# %%

# %%
merged_wod_df

# %%
# Create a new dataframe with 410761_flow and UKIH columns,
# and a new column for the proportion of UKIH (baseflow) to total flow
ukih_prop_wod_df = merged_wod_df[['410761_flow', 'UKIH']].copy()
ukih_prop_wod_df['UKIH_proportion'] = ukih_prop_wod_df['UKIH'] / ukih_prop_wod_df['410761_flow']
ukih_prop_wod_df


# %%
merged_empirical_df

# %%
# Create a new dataframe with 410761_flow and UKIH columns,
# and a new column for the proportion of UKIH (baseflow) to total flow
ukih_prop_empirical_df = merged_empirical_df[['Discharge_ML_day', 'UKIH']].copy()
ukih_prop_empirical_df['UKIH_proportion'] = ukih_prop_empirical_df['UKIH'] / ukih_prop_empirical_df['Discharge_ML_day']
ukih_prop_empirical_df

# %%
import matplotlib.pyplot as plt
import numpy as np
import calendar
from scipy.stats import gaussian_kde

# Add a Month column if not present
def ensure_month_column(df):
    if 'Month' not in df.columns:
        df = df.copy()
        df['Month'] = df.index.month
    return df

wod_df = ensure_month_column(ukih_prop_wod_df)
emp_df = ensure_month_column(ukih_prop_empirical_df)

months = [calendar.month_abbr[m] for m in range(1, 13)]

fig, axes = plt.subplots(3, 4, figsize=(18, 12), sharex=True, sharey=True)
axes = axes.flatten()

for m in range(1, 13):
    ax = axes[m-1]
    # Get subset for each source
    wod_subset = wod_df[wod_df['Month'] == m]['UKIH_proportion'].dropna()
    emp_subset = emp_df[emp_df['Month'] == m]['UKIH_proportion'].dropna()
    
    x_grid = np.linspace(0, 1.05, 300)
    has_plotted = False

    # Plot WOD KDE
    if not wod_subset.empty and len(wod_subset) > 1:
        wod_kde = gaussian_kde(wod_subset)
        wod_y = wod_kde(x_grid)
        ax.plot(x_grid, wod_y, color='tab:blue', lw=2, label='WOD')
        has_plotted = True

    # Plot Empirical KDE
    if not emp_subset.empty and len(emp_subset) > 1:
        emp_kde = gaussian_kde(emp_subset)
        emp_y = emp_kde(x_grid)
        ax.plot(x_grid, emp_y, color='tab:orange', lw=2, label='Empirical')
        has_plotted = True

    ax.set_title(months[m-1])

    # Include stats for both, with indication in the text
    box_props = dict(facecolor='white', alpha=0.8, edgecolor='gray')
    stats_ypos = 0.95

    if not wod_subset.empty:
        wod_mean = wod_subset.mean()
        wod_median = wod_subset.median()
        wod_p10 = wod_subset.quantile(0.10)
        wod_p90 = wod_subset.quantile(0.90)
        wod_stats = (
            f"WOD:\n"
            f"  Mean:   {wod_mean:.3f}\n"
            f"  Median: {wod_median:.3f}\n"
            f"  10th:   {wod_p10:.3f}\n"
            f"  90th:   {wod_p90:.3f}"
        )
        ax.text(0.05, stats_ypos, wod_stats, ha='left', va='top', fontsize=10,
                bbox=box_props, color='tab:blue', transform=ax.transAxes)
        stats_ypos -= 0.31

    if not emp_subset.empty:
        emp_mean = emp_subset.mean()
        emp_median = emp_subset.median()
        emp_p10 = emp_subset.quantile(0.10)
        emp_p90 = emp_subset.quantile(0.90)
        emp_stats = (
            f"Empirical:\n"
            f"  Mean:   {emp_mean:.3f}\n"
            f"  Median: {emp_median:.3f}\n"
            f"  10th:   {emp_p10:.3f}\n"
            f"  90th:   {emp_p90:.3f}"
        )
        ax.text(0.35, 0.95, emp_stats, ha='left', va='top', fontsize=10,
                bbox=box_props, color='tab:orange', transform=ax.transAxes)

    # Only add legend if something plotted
    if has_plotted:
        ax.legend(loc="upper right", fontsize=9)

for ax in axes:
    ax.set_xlim(0, 1.05)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('UKIH Baseflow Proportion')
    ax.set_ylabel('Density')

fig.suptitle('Distribution of UKIH Baseflow Proportion to Total Flow By Month\n(Comparison: WOD vs Empirical)', fontsize=18)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()


# %%
import matplotlib.pyplot as plt

# Group by month and calculate median for each DataFrame
wod_monthly_median = ukih_prop_wod_df.groupby(ukih_prop_wod_df.index.month)["UKIH_proportion"].median()
empirical_monthly_median = ukih_prop_empirical_df.groupby(ukih_prop_empirical_df.index.month)["UKIH_proportion"].median()

# Prepare x labels as month names
import calendar
months = range(1, 13)
month_names = [calendar.month_abbr[m] for m in months]

# Create comparison plot
plt.figure(figsize=(10, 6))
plt.plot(months, wod_monthly_median[months], marker='o', label="Without Development")
plt.plot(months, empirical_monthly_median[months], marker='o', label="Empirical Data")

plt.xticks(months, month_names)
plt.xlabel("Month")
plt.ylabel("Median UKIH Baseflow Proportion")
plt.title("Median UKIH Baseflow Proportion by Month")
plt.legend()
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.7, alpha=0.7)  # Add both horizontal and vertical grid lines
plt.tight_layout()
plt.show()

# %%
import matplotlib.pyplot as plt

# Group by month and calculate median of UKIH baseflow for each DataFrame
wod_monthly_median = ukih_prop_wod_df.groupby(ukih_prop_wod_df.index.month)["UKIH"].median()
empirical_monthly_median = ukih_prop_empirical_df.groupby(ukih_prop_empirical_df.index.month)["UKIH"].median()

# Prepare x labels as month names
import calendar
months = range(1, 13)
month_names = [calendar.month_abbr[m] for m in months]

# Create comparison plot
plt.figure(figsize=(10, 6))
plt.plot(months, wod_monthly_median[months], marker='o', label="Without Development")
plt.plot(months, empirical_monthly_median[months], marker='o', label="Empirical Data")

plt.xticks(months, month_names)
plt.xlabel("Month")
plt.ylabel("Median UKIH Baseflow")
plt.title("Median UKIH Baseflow by Month")
plt.legend()
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.7, alpha=0.7)
plt.tight_layout()
plt.show()

# %%
import matplotlib.pyplot as plt

# Group Q_wod and Q_empirical series by month and calculate median for each month
wod_monthly_median = Q_wod.groupby(Q_wod.index.month).median()
empirical_monthly_median = Q_empirical.groupby(Q_empirical.index.month).median()

# @file_context_0: Group by month and calculate median of UKIH baseflow for each DataFrame
wod_ukih_monthly_median = ukih_prop_wod_df.groupby(ukih_prop_wod_df.index.month)["UKIH"].median()
empirical_ukih_monthly_median = ukih_prop_empirical_df.groupby(ukih_prop_empirical_df.index.month)["UKIH"].median()

# Prepare x labels as month names
import calendar
months = range(1, 13)
month_names = [calendar.month_abbr[m] for m in months]

# Create comparison plot
plt.figure(figsize=(10, 6))
plt.plot(months, wod_monthly_median[months], marker='o', label="Without Development $Q_\mathrm{total}$")
plt.plot(months, empirical_monthly_median[months], marker='o', label="Empirical Data $Q_\mathrm{total}$")

# Add UKIH baseflow data from @file_context_0 as dashed lines
plt.plot(months, wod_ukih_monthly_median[months], marker='s', linestyle='--', label="Without Dev. UKIH Baseflow (dashed)")
plt.plot(months, empirical_ukih_monthly_median[months], marker='s', linestyle='--', label="Empirical UKIH Baseflow (dashed)")

# Add a horizontal line at 150 ML/d
plt.axhline(150, color='r', linestyle=':', linewidth=1.5, label='150 ML/d threshold')

plt.xticks(months, month_names)
plt.xlabel("Month")
plt.ylabel("Median Flow (ML/day)")
plt.title("Median Total Flows and UKIH Baseflow by Month @ 410761")
plt.legend()
plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.7, alpha=0.7)
plt.tight_layout()
plt.show()

# %%
Q_wod

# %%

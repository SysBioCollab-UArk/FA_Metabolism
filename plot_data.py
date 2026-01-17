import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
from pathlib import Path
import uuid


def make_expt_id(cell_line: str, n_chars: int = 6) -> str:
    token = uuid.uuid4().hex[:n_chars]
    return f"{cell_line}_{token}"


def _read_metadata(filepath, metadata_rows=9):
    """Read Seahorse metadata block (typically first 9 Excel rows, including time row)."""
    meta_raw = pd.read_excel(filepath, nrows=metadata_rows, header=None)
    meta = meta_raw.dropna(how="all").set_index(0)[1].to_dict()
    if 'Units' in meta.keys():
        meta['Data Units'] = meta['Units']
        meta.pop('Units')
    # add time units to metadata dict
    idx = next((i for i, s in enumerate(list(meta.keys())) if 'Time' in s), None)
    if idx is not None:
        s = list(meta.keys())[idx]
        time_units = s[s.find('(') + 1: s.find(')')]
        meta['Time Units'] = time_units
        meta.pop(s)

    return meta


def plot_atp_rate_assay_data(ax, atp_df, data_info, fs_labels=None):

    expt_id = atp_df['Expt_ID'].unique()
    if len(expt_id) > 1:
        raise ValueError('Multiple Expt_IDs detected.')

    labels = [data_info['mutation'] + '-', 'Corrected', 'Healthy']
    columns = ['Signal_Mut', 'Signal_Corr', 'Signal_Ctrl']
    glyco_mean = atp_df[atp_df['Signal_Type'] == 'glycoATP'][columns].to_numpy().flatten()
    mito_mean = atp_df[atp_df['Signal_Type'] == 'mitoATP'][columns].to_numpy().flatten()
    columns = ['SEM_Mut', 'SEM_Corr', 'SEM_Ctrl']
    glyco_sem = atp_df[atp_df['Signal_Type'] == 'glycoATP'][columns].to_numpy().flatten()
    mito_sem = atp_df[atp_df['Signal_Type'] == 'mitoATP'][columns].to_numpy().flatten()

    x = np.arange(len(labels))
    error_kw = {'elinewidth': 2, 'capthick': 2}
    ax.bar(x, glyco_mean, yerr=glyco_sem, capsize=6, label='glycoATP', color='r', error_kw=error_kw)
    ax.bar(x, mito_mean, bottom=glyco_mean, capsize=6, yerr=mito_sem, label='mitoATP', color='dodgerblue',
           error_kw=error_kw)

    ax.set_ylabel('ATP Production Rate\n(%s)' % data_info['sig_units'], fontsize=fs_labels)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)


def _flag_outliers_iqr(df, group_cols, value_col='Value', k=1.5):
    """
    Vectorized IQR outlier flagging within each group_cols group.
    Returns a boolean Series aligned to df.index.
    """
    q1 = df.groupby(group_cols, observed=True, sort=False)[value_col].transform('quantile', 0.25)
    q3 = df.groupby(group_cols, observed=True, sort=False)[value_col].transform('quantile', 0.75)
    iqr = q3 - q1

    lower = q1 - k * iqr
    upper = q3 + k * iqr

    # If a group has <2 points, quantiles may be NaN; treat as no outliers.
    is_outlier = (df[value_col] < lower) | (df[value_col] > upper)
    return is_outlier.fillna(False)


def add_outlier_flags(df, assay_kind, k=1.5, value_col='Value'):
    """
    Add an 'is_outlier' column to df, based on assay type.

    assay_kind:
      - 'kinetic' (OCR/ECAR/PER): group by ['Condition','Time']
      - 'atp_rate' (ATP Rate): group by ['Condition','Component']
    """
    df = df.copy()

    if assay_kind == 'kinetic':
        group_cols = ['Condition', 'Time']
    elif assay_kind == 'atp_rate':
        group_cols = ['Condition', 'Component']
    else:
        raise ValueError(f"Unknown assay_kind={assay_kind!r}")

    df['is_outlier'] = _flag_outliers_iqr(df, group_cols=group_cols, value_col=value_col, k=k)
    return df


def remove_outliers(df, assay_kind, k=1.5, value_col='Value'):
    """
    Return (df_clean, df_flagged, removed) where:
      - df_flagged contains an 'is_outlier' boolean column
      - df_clean excludes outliers
      - removed contains only outliers
    """
    df_flagged = add_outlier_flags(df, assay_kind=assay_kind, k=k, value_col=value_col)
    removed = df_flagged[df_flagged['is_outlier']].copy()
    df_clean = df_flagged[~df_flagged['is_outlier']].copy()
    return df_clean, df_flagged, removed


def _load_seahorse_kinetic_table(filepath, header_row=8, well_row=9):
    """
    Kinetic-style exports (OCR/ECAR/PER): multi-row header with group headers then well IDs.
    Excel row 9 = headers (Time (minutes), PD20, RV, Healthy...), Excel row 10 = well labels.
    """
    df_raw = pd.read_excel(filepath, header=[header_row, well_row])

    # Normalize the time column name to ('Time','') so downstream logic is consistent
    df_raw.columns = [
        ('Time', '') if c[0] in ('Time (min)', 'Time (minutes)', 'Time') else c
        for c in df_raw.columns
    ]

    # Melt using a scalar var_name (works across pandas versions)
    df_long = df_raw.melt(id_vars=[('Time', '')], var_name='variable',
                          value_name='Value').rename(columns={('Time', ''): 'Time'})

    # variable will be tuples like ('PD20','A1')
    df_long[['Condition', 'Well']] = pd.DataFrame(df_long['variable'].tolist(), index=df_long.index)
    df_long = df_long.drop(columns='variable')

    df_long['Time'] = pd.to_numeric(df_long['Time'], errors='coerce')
    df_long['Value'] = pd.to_numeric(df_long['Value'], errors='coerce')

    # Drop non-numeric / blank rows
    df_long = df_long.dropna(subset=['Time', 'Value'])

    return df_long


def _load_seahorse_atp_rate_table(filepath, header_row=8):
    """
    ATP Rate Assay export (basal rates table):
    columns: Group, Well, glycoATP Production Rate, mitoATP Production Rate
    Excel row 9 = header row (pandas header_row=8).
    Returns long format: Condition, Well, Component (glycoATP/mitoATP), Value.
    """
    df = pd.read_excel(filepath, header=header_row)

    # Robust column finding (handles small naming variations)
    col_map = {}
    for c in df.columns:
        c_str = str(c).strip()
        if c_str.lower() == 'group':
            col_map[c] = 'Condition'
        elif c_str.lower() == 'well':
            col_map[c] = 'Well'
        elif 'glycoatp' in c_str.lower():
            col_map[c] = 'glycoATP'
        elif 'mitoatp' in c_str.lower():
            col_map[c] = 'mitoATP'

    df = df.rename(columns=col_map)

    required = {'Condition', 'Well', 'glycoATP', 'mitoATP'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"ATP Rate parser: missing expected columns: {missing}. Found: {list(df.columns)}")

    df = df[['Condition', 'Well', 'glycoATP', 'mitoATP']].copy()
    df['glycoATP'] = pd.to_numeric(df['glycoATP'], errors='coerce')
    df['mitoATP'] = pd.to_numeric(df['mitoATP'], errors='coerce')

    df_long = df.melt(id_vars=['Condition', 'Well'], value_vars=['glycoATP', 'mitoATP'], var_name='Component',
                      value_name='Value').dropna(subset=['Value'])

    return df_long


def load_seahorse_excel(filepath, metadata_rows=9, kinetic_header_row=8, kinetic_well_row=9, atp_header_row=8,
                        remove_outliers_flag=False, outlier_k=1.5):
    """
    Unified Seahorse loader that auto-detects assay type from metadata and routes to the
    appropriate parser. Optionally flags/removes outliers.

    Returns
    -------
    df : pandas.DataFrame
        If remove_outliers_flag=False:
          - kinetic: columns ['Time','Condition','Well','Value']
          - atp_rate: columns ['Condition','Well','Component','Value']
        If remove_outliers_flag=True:
          - df is cleaned (outliers removed)
    meta : dict
        Metadata from the top of the file.
    extras : dict
        Only returned if remove_outliers_flag=True:
          {
            'df_flagged': DataFrame with is_outlier column,
            'removed': DataFrame containing removed points
          }
    """
    meta = _read_metadata(filepath, metadata_rows=metadata_rows)
    data_type = str(meta.get('Data Type', '')).strip().lower()

    if 'atp production rate' in data_type:
        assay_kind = 'atp_rate'
        df = _load_seahorse_atp_rate_table(filepath, header_row=atp_header_row)
    else:
        assay_kind = 'kinetic'
        df = _load_seahorse_kinetic_table(filepath, header_row=kinetic_header_row, well_row=kinetic_well_row)

    if not remove_outliers_flag:
        return df, meta, dict()

    df_clean, df_flagged, removed = remove_outliers(df, assay_kind=assay_kind, k=outlier_k)
    df_clean['Condition'] = df_clean['Condition'].astype(str).str.strip()
    df_clean['Well'] = df_clean['Well'].astype(str).str.strip()

    extras = {
        'assay_kind': assay_kind,
        'df_flagged': df_flagged,
        'removed': removed
    }

    return df_clean, meta, extras


def plot_seahorse_data(datafile, expt_ids='all', datatype='all', scale_to_first_point=False, save_plot=True,
                       show_plot=False, squeeze_plots=False, **kwargs):

    # process kwargs
    figsize_adjust = kwargs.get('figsize_adjust', (1, 1))
    xlim = kwargs.get('xlim', (None, None))
    ylim = kwargs.get('ylim', (None, None))
    title_type = kwargs.get('title', 'line_mut_type')
    title_kw = kwargs.get('title_kw', {})
    suffix = kwargs.get('suffix', None)
    suptitle = kwargs.get('suptitle', None)
    sharey = kwargs.get('sharey', False)

    fontsizes = kwargs.get('fontsizes', {})
    fs_labels = fontsizes.get('labels', None)
    fs_ticks = fontsizes.get('ticks', None)
    fs_legend = fontsizes.get('legend', None)
    fs_title = fontsizes.get('title', None)
    fs_suptitle = fontsizes.get('suptitle', None)

    preferred_fig_order = ['OCR', 'ECAR', 'PER', 'ATP']

    all_data = pd.read_csv(datafile)

    if expt_ids == 'all':
        expt_ids = all_data['Expt_ID'].unique()
    elif isinstance(expt_ids, str):
        expt_ids = [expt_ids]
    print(expt_ids)

    # only keep data for specified expt IDs
    all_data = all_data[all_data['Expt_ID'].isin(expt_ids)]

    if datatype != 'all':
        if isinstance(datatype, str):
            datatype = [datatype]
        all_data = all_data[all_data['Signal_Type'].isin(datatype)]
    print(all_data.columns)

    # get signal types and numbers of plots per signal type
    sig_types = all_data['Signal_Type'].unique()
    sig_types_dict = dict(
        zip(
            sig_types,
            [all_data[all_data['Signal_Type'] == sig_type]['Expt_ID'].nunique() for sig_type in sig_types]
        )
    )
    sig_types_dict['ATP'] = (sig_types_dict.get('mitoATP', 0) + sig_types_dict.get('glycoATP', 0)) // 2

    # replace 'mitoATP' and 'glycoATP' with 'ATP'
    for i in range(len(sig_types)):
        if re.search(r'ATP$', sig_types[i]):
            sig_types[i] = 'ATP'
    sig_types = list(set(sig_types))  # remove duplicates

    # one column per signal type
    ncols = len(sig_types)

    # make sure signal types match known types
    if any(sig_type not in preferred_fig_order for sig_type in sig_types):
        raise ValueError('Unrecognized signal type detected in %s. Supported types are %s.' %
                         (str(sig_types), str(preferred_fig_order)))
    # order 'fig_order' based on 'preferred_order'
    fig_order = sorted(sig_types, key=lambda item: preferred_fig_order.index(item))
    print(fig_order)

    # get number of rows
    if not squeeze_plots:
        nrows = len(expt_ids)  # default is to put plots for each expt on a separate row
    else:
        nrows = max([sig_types_dict[sig_type] for sig_type in sig_types])
        for sig_type in sig_types:
            sig_types_dict[sig_type] = nrows

    fig = plt.figure(constrained_layout=True,
                     figsize=(6.4 * 0.8 * ncols * figsize_adjust[0], 4.8 * 0.65 * nrows * figsize_adjust[1]))

    row = 0
    ax_col_dict = {col: [] for col in range(ncols)}
    for i, expt_id in enumerate(expt_ids):
        print(expt_id)
        data_expt = all_data[all_data['Expt_ID'] == expt_id]
        sig_types = data_expt['Signal_Type'].unique()
        atp_assays = [sig_type for sig_type in sig_types if 'ATP' in sig_type]
        sig_types = [[sig_type] for sig_type in sig_types if 'ATP' not in sig_type]
        if len(atp_assays) > 1:
             sig_types += [atp_assays]

        for sig_type in sig_types:
            print('   ', sig_type)
            col = fig_order.index(sig_type[0] if len(sig_type) == 1 else 'ATP')

            if squeeze_plots:
                sig_types_dict[sig_type[0] if len(sig_type) == 1 else 'ATP'] -= 1
                row = nrows - sig_types_dict[sig_type[0] if len(sig_type) == 1 else 'ATP'] - 1

            ax = fig.add_subplot(nrows, ncols, row * ncols + col + 1)
            ax_col_dict[col].append(ax)

            data_expt_sig = data_expt[data_expt['Signal_Type'].isin(sig_type)]

            data_info = {
                'mutation': data_expt_sig['Mutation'].unique(),
                'cell_line': data_expt_sig['Cell_Line'].unique(),
                'cell_type': data_expt_sig['Cell_Type'].unique(),
                'time_units': data_expt_sig['Time_Units'].unique(),
                'sig_units': data_expt_sig['Signal_Units'].unique()
            }

            # make sure all metadata is unique
            for key in data_info.keys():
                if len(data_info[key]) > 1:
                    raise ValueError('Multiple %s values detected for expt %s: %s' %
                                     (key, expt_id, str(data_info[key])))
                else:
                    data_info[key] = data_info[key][0]

            # Plot ATP Rate Assay data
            if len(sig_type) > 1:
                plot_atp_rate_assay_data(ax, data_expt_sig, data_info, fs_labels=fs_labels)
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
                if fs_ticks is not None:
                    ax.yaxis.get_offset_text().set_fontsize(0.8 * fs_ticks)

            # Plot Mito and Glyco Rate Assay data
            else:
                sig_type = sig_type[0]
                scale_factors = [data_expt_sig[n].iloc[0] if scale_to_first_point else 1
                                 for n in ['Signal_Mut', 'Signal_Corr', 'Signal_Ctrl']]
                sig_mutant = data_expt_sig['Signal_Mut'] / scale_factors[0]
                sig_corrected = data_expt_sig['Signal_Corr'] / scale_factors[1]
                sig_healthy = data_expt_sig['Signal_Ctrl'] / scale_factors[2]
                sem_mutant = data_expt_sig['SEM_Mut'] / scale_factors[0]
                sem_corrected = data_expt_sig['SEM_Corr'] / scale_factors[1]
                sem_healthy = data_expt_sig['SEM_Ctrl'] / scale_factors[2]

                ax.errorbar(data_expt_sig['Time'], sig_mutant, yerr=sem_mutant, marker='o', ls='-', ms=10,
                            capsize=6, label=data_info['mutation'] + '-')
                ax.errorbar(data_expt_sig['Time'], sig_corrected, yerr=sem_corrected, marker='s', ls='-', ms=10,
                            capsize=6, label='Corrected')
                ax.errorbar(data_expt_sig['Time'], sig_healthy, yerr=sem_healthy, marker='^', ls='-', ms=10,
                            capsize=6, label='Healthy')

                ax.set_xlabel('Time (%s)' % data_info['time_units'], fontsize=fs_labels)
                ax.set_ylabel('%s (%s)' % (sig_type, data_info['sig_units']) \
                                  if not scale_to_first_point else 'Relative %s' % sig_type, fontsize=fs_labels)

                ax.set_xlim(xlim)

            # Set the plot title. Currently supported plot title types are: 'none', 'line_mut_type', and 'expt_it'.
            # The default is 'line_mut_type'.
            extra_space = '\n' if row > 0 else ''
            if title_type != 'none':
                title = extra_space
                if title_type == 'line_mut_type':
                    title += '%s, %s-deficient, %s' % \
                            (data_info['cell_line'], data_info['mutation'], data_info['cell_type'])
                elif title_type == 'expt_id':
                    title += expt_id
                else:
                    raise ValueError('Unrecognized title type: %s' % title_type)
                ax.set_title(title, fontsize=fs_title, **title_kw)

            ax.tick_params(axis='both', which='major', labelsize=fs_ticks)
            ax.legend(loc='best', fontsize=fs_legend, frameon=False, borderpad=0.1)

        row += 1

    # share y-axes, if requested
    if sharey:
        for col in range(ncols):
            bottom = min([ax.get_ylim()[0] for ax in ax_col_dict[col]])
            top = max([ax.get_ylim()[1] for ax in ax_col_dict[col]])
            for ax in ax_col_dict[col]:
                ax.set_ylim((bottom, top))

    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=fs_suptitle, fontweight='bold')

    if save_plot is not False:
        outdir = '.' if save_plot is True else save_plot
        filename = 'Seahorse_plots%s.png' % ('_' + suffix if suffix is not None else '')
        fig.savefig(os.path.join(outdir, filename), dpi=300)

    if show_plot:
        plt.show()


def export_seahorse_kinetic_to_csv(df: pd.DataFrame, meta: dict, out_csv: str | Path, *, expt_id: str,
                                   condition_map: dict, mutation: str, cell_line: str, cell_type: str,
                                   write_mode: str = "append", float_format: str | None = None
                                   ) -> pd.DataFrame:
    """
    Export OCR/ECAR/PER kinetic Seahorse data (tidy format) to a standardized CSV schema.

    Parameters
    ----------
    df:
        Tidy kinetic dataframe with columns ['Time','Condition','Well','Value'].
    meta:
        Metadata dict from loader. Should contain 'Data Type' and 'Units' when available.
    out_csv:
        CSV path to write/append.
    condition_map:
        dict mapping roles -> condition labels in df:
            {'mut': 'PD20', 'corr': 'RV', 'ctrl': 'Healthy'}
        Values must match df['Condition'] values.
    mutation, cell_line, cell_type:
        Manual fields you want in every exported row.
    expt_id:
        If None, will try meta['Expt_ID'] then meta['Assay Name'] (fallback).
        In practice, you likely want to pass this explicitly (e.g., derived from filename).
    time_units:
        If None, defaults to 'min'.
    signal_type:
        If None, uses meta['Data Type'] (e.g., OCR/ECAR/PER).
    signal_units:
        If None, uses meta['Units'] (e.g., pmol/min/1000 Cells) and normalizes a bit.

    write_mode:
        "append" -> create if missing else append without header
        "overwrite" -> overwrite existing file and write header
    float_format:
        Optional, e.g., "%.6f"

    Returns
    -------
    exported_df:
        The exact dataframe written/appended (one row per timepoint).
    """
    REQUIRED_COLS = ["Time", "Time_Units", "Mutation", "Cell_Line", "Cell_Type", "Signal_Type", "Signal_Units",
                     "Signal_Mut", "Signal_Corr", "Signal_Ctrl", "SEM_Mut", "SEM_Corr", "SEM_Ctrl", "Expt_ID",]

    out_csv = Path(out_csv)

    # --- validate input dataframe ---
    needed = {"Time", "Condition", "Well", "Value"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"df is missing required columns: {missing}. Found: {list(df.columns)}")

    # --- info from metadata ---
    time_units = str(meta['Time Units']).strip()
    signal_type = str(meta["Data Type"]).strip()  # OCR/ECAR/PER expected
    signal_units = str(meta["Data Units"]).strip()

    # --- compute per-timepoint mean and standard error for each Condition ---
    # summary = df.groupby(["Condition", "Time"], sort=False, observed=True).agg(mean=("Value", "mean")).reset_index()
    summary = df.groupby(["Condition", "Time"], sort=False, observed=True).\
        agg(mean=("Value", "mean"), sem=("Value", "sem"),  n=("Value", "count")).reset_index()

    # --- pivot into columns by role (mut/corr/ctrl) ---
    # pivot = summary.pivot(index="Time", columns="Condition", values="mean")
    pivot_mean = summary.pivot(index="Time", columns="Condition", values="mean")
    pivot_sem = summary.pivot(index="Time", columns="Condition", values="sem")

    # Pull out the three signals using the mapping
    def _get_cond(role: str) -> str:
        if role not in condition_map:
            raise ValueError(f"condition_map must contain key {role!r}. Got keys: {list(condition_map)}")
        return condition_map[role]

    mut_cond = _get_cond("mut")
    corr_cond = _get_cond("corr")
    ctrl_cond = _get_cond("ctrl")

    for cond in [mut_cond, corr_cond, ctrl_cond]:
        if cond not in pivot_mean.columns:
            raise ValueError(
                f"Condition {cond!r} not found in data. "
                f"Available conditions: {list(pivot_mean.columns)}"
            )
        if cond not in pivot_sem.columns:
            raise ValueError(
                f"Condition {cond!r} not found in data. "
                f"Available conditions: {list(pivot_sem.columns)}"
            )

    exported = pd.DataFrame({
        "Time": pivot_mean.index.astype(float),
        "Time_Units": time_units,
        "Mutation": mutation,
        "Cell_Line": cell_line,
        "Cell_Type": cell_type,
        "Signal_Type": signal_type,
        "Signal_Units": signal_units,
        "Signal_Mut": pivot_mean[mut_cond].to_numpy(),
        "Signal_Corr": pivot_mean[corr_cond].to_numpy(),
        "Signal_Ctrl": pivot_mean[ctrl_cond].to_numpy(),
        "SEM_Mut": pivot_sem[mut_cond].to_numpy(),
        "SEM_Corr": pivot_sem[corr_cond].to_numpy(),
        "SEM_Ctrl": pivot_sem[ctrl_cond].to_numpy(),
        "Expt_ID": expt_id,
    }).reset_index(drop=True)

    # Ensure column order
    exported = exported[REQUIRED_COLS]

    # --- write/append logic ---
    write_mode = write_mode.lower().strip()
    if write_mode not in {"append", "overwrite"}:
        raise ValueError("write_mode must be 'append' or 'overwrite'.")

    if write_mode == "overwrite":
        exported.to_csv(out_csv, index=False, mode="w", header=True, float_format=float_format)
    else:
        # append
        file_exists = out_csv.exists()
        exported.to_csv(out_csv, index=False, mode="a", header=not file_exists, float_format=float_format)

    return exported


def export_seahorse_atp_rate_to_csv(df: pd.DataFrame, meta: dict, out_csv: str | Path, *, expt_id: str,
                                    condition_map: dict, mutation: str, cell_line: str, cell_type: str,
                                    write_mode: str = "append", float_format: str | None = None
                                    ) -> pd.DataFrame:
    """
    Export Seahorse ATP Rate Assay data (tidy format) to the same standardized CSV schema used
    for OCR/ECAR/PER exports.

    Expected input df columns:
        ['Condition', 'Well', 'Component', 'Value']
    where Component is typically 'mitoATP' and 'glycoATP'.

    Output:
        Two rows (one per Component), with Time and Time_Units left blank.

    Notes
    -----
    - Signal_Type will be set to the Component value (e.g., mitoATP, glycoATP).
    - Signal_Units will be taken from meta["Data Units"].
    """

    REQUIRED_COLS = ["Time", "Time_Units", "Mutation", "Cell_Line", "Cell_Type", "Signal_Type", "Signal_Units",
                     "Signal_Mut", "Signal_Corr", "Signal_Ctrl", "SEM_Mut", "SEM_Corr", "SEM_Ctrl", "Expt_ID",]

    out_csv = Path(out_csv)

    # --- validate input dataframe ---
    needed = {"Condition", "Well", "Component", "Value"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"df is missing required columns: {missing}. Found: {list(df.columns)}")

    # --- normalize labels (protect against Excel whitespace) ---
    df = df.copy()
    df["Condition"] = df["Condition"].astype(str).str.strip()
    df["Component"] = df["Component"].astype(str).str.strip()

    # --- info from metadata ---
    # ATP rate rows have no time axis in the export, but we keep the same schema with blanks.
    signal_units = str(meta["Data Units"]).strip()

    # --- compute per-component mean and standard error for each Condition ---
    summary = df.groupby(["Component", "Condition"], sort=False, observed=True). \
        agg(mean=("Value", "mean"), sem=("Value", "sem"), n=("Value", "count")).reset_index()

    # --- pivot into columns by role (mut/corr/ctrl) for mean and sem ---
    pivot_mean = summary.pivot(index="Component", columns="Condition", values="mean")
    pivot_sem = summary.pivot(index="Component", columns="Condition", values="sem")

    # Pull out the three signals using the mapping
    def _get_cond(role: str) -> str:
        if role not in condition_map:
            raise ValueError(f"condition_map must contain key {role!r}. Got keys: {list(condition_map)}")
        return condition_map[role]

    mut_cond = _get_cond("mut")
    corr_cond = _get_cond("corr")
    ctrl_cond = _get_cond("ctrl")

    for cond in [mut_cond, corr_cond, ctrl_cond]:
        if cond not in pivot_mean.columns:
            raise ValueError(
                f"Condition {cond!r} not found in data. "
                f"Available conditions: {list(pivot_mean.columns)}"
            )
        if cond not in pivot_sem.columns:
            raise ValueError(
                f"Condition {cond!r} not found in data. "
                f"Available conditions: {list(pivot_sem.columns)}"
            )

    # Ensure we output exactly these two components (order matches your example)
    components_out = ["mitoATP", "glycoATP"]
    missing_components = [c for c in components_out if c not in pivot_mean.index]
    if missing_components:
        raise ValueError(
            f"Component(s) {missing_components} not found in data. "
            f"Available components: {list(pivot_mean.index)}"
        )

    exported = pd.DataFrame({
        "Time": ["", ""],
        "Time_Units": ["", ""],
        "Mutation": [mutation, mutation],
        "Cell_Line": [cell_line, cell_line],
        "Cell_Type": [cell_type, cell_type],
        "Signal_Type": components_out,
        "Signal_Units": [signal_units, signal_units],
        "Signal_Mut": [pivot_mean.loc[c, mut_cond] for c in components_out],
        "Signal_Corr": [pivot_mean.loc[c, corr_cond] for c in components_out],
        "Signal_Ctrl": [pivot_mean.loc[c, ctrl_cond] for c in components_out],
        "SEM_Mut": [pivot_sem.loc[c, mut_cond] for c in components_out],
        "SEM_Corr": [pivot_sem.loc[c, corr_cond] for c in components_out],
        "SEM_Ctrl": [pivot_sem.loc[c, ctrl_cond] for c in components_out],
        "Expt_ID": [expt_id, expt_id],
    }).reset_index(drop=True)

    # Ensure column order
    exported = exported[REQUIRED_COLS]

    # --- write/append logic ---
    write_mode = write_mode.lower().strip()
    if write_mode not in {"append", "overwrite"}:
        raise ValueError("write_mode must be 'append' or 'overwrite'.")

    if write_mode == "overwrite":
        exported.to_csv(out_csv, index=False, mode="w", header=True, float_format=float_format)
    else:
        # append
        file_exists = out_csv.exists()
        exported.to_csv(out_csv, index=False, mode="a", header=not file_exists, float_format=float_format)

    return exported


if __name__ == '__main__':
    import os

    data_path = 'Excel_Data'  # '.'
    datafile = os.path.join(data_path, 'fa_seahorse_export.csv')  # 'seahorse_data.csv')

    cell_lines = ['PD20', 'PD220']
    suptitles = ['FANCD2-deficient (PD20)', 'FANCA-deficient (PD220)']

    # expt_ids = ['PD20_121625', 'PD20_121825', 'PD20_122625', 'PD20_122825']
    # expt_ids = ['PD220_112225', 'PD220_112625_r1', 'PD220_112625_r2', 'PD220_120525_r1', 'PD220_120525_r2',
    #             'PD220_121325', 'PD220_121425']
    # suptitle = 'FANCA-deficient (PD220)'
    # 'all'  # ['PD20_050125', 'PD20_122124', 'PD20_031325', 'PD20_050125']

    kwargs = {
        'figsize_adjust': (1.8, 2.2),
        'title': 'none',  # 'expt_id' 'line_mut_type' 'none'
        'title_kw': {'fontweight': 'normal', 'color': 'k'},
        'xlim': (-2, 102),
        'fontsizes': {'labels': 32, 'ticks': 32, 'legend': 20, 'title': 20, 'suptitle': 40},
        'sharey': True
    }

    df = pd.read_csv(datafile)
    for cell_line, suptitle in zip(cell_lines, suptitles):
        expt_ids = df.loc[df['Expt_ID'].str.startswith(cell_line), 'Expt_ID'].unique().tolist()
        kwargs['suptitle'] = suptitle
        plot_seahorse_data(datafile, expt_ids=expt_ids, datatype='all', scale_to_first_point=False, show_plot=True,
                           squeeze_plots=True, **kwargs)

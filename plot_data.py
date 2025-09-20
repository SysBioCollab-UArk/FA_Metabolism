import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re


def plot_atp_rate_assay_data(ax, atp_df, data_info):

    expt_id = atp_df['Expt_ID'].unique()
    if len(expt_id) > 1:
        raise ValueError('Multiple Expt_IDs detected.')

    labels = [data_info['mutation'] + '-', 'Corrected', 'Healthy']
    columns = ['Signal_Mut', 'Signal_Corr', 'Signal_Ctrl']
    glyco_data = atp_df[atp_df['Signal_Type'] == 'glycoATP'][columns].to_numpy().flatten()
    mito_data = atp_df[atp_df['Signal_Type'] == 'mitoATP'][columns].to_numpy().flatten()

    if len(glyco_data) > 3:
        raise ValueError('More than 3 glycoATP values detected: %s' % str(glyco_data))
    if len(mito_data) > 3:
        raise ValueError('More than 3 mitoATP values detected: %s' % str(mito_data))

    x = np.arange(len(labels))
    ax.bar(x, glyco_data, label='glycoATP', color='r')
    ax.bar(x, mito_data, bottom=glyco_data, label='mitoATP', color='dodgerblue')

    ax.set_ylabel('ATP Production Rate\n(%s)' % data_info['sig_units'], fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)


def plot_seahorse_data(datafile, expt_ids='all', datatype='all', scale_to_first_point=False, save_plot=True,
                       show_plot=False, **kwargs):

    preferred_fig_order = ['OCR', 'ECAR', 'PER', 'ATP']

    all_data = pd.read_csv(datafile)

    if expt_ids == 'all':
        expt_ids = all_data['Expt_ID'].unique()
    elif isinstance(expt_ids, str):
        expt_ids = [expt_ids]
    print(expt_ids)

    if datatype != 'all':
        if isinstance(datatype, str):
            datatype = [datatype]
        all_data = all_data[all_data['Signal_Type'].isin(datatype)]
    print(all_data.columns)

    # get number of rows and columns for figure
    nrows = len(expt_ids)  # put plots for each expt on a separate row
    sig_types = all_data['Signal_Type'].unique()
    # replace 'mitoATP' and 'glycoATP' with 'ATP'
    for i in range(len(sig_types)):
        if re.search(r'ATP$', sig_types[i]):
            sig_types[i] = 'ATP'
    sig_types = list(set(sig_types))  # remove duplicates
    ncols = len(sig_types)
    # make sure signal types match known types
    if any(sig_type not in preferred_fig_order for sig_type in sig_types):
        raise ValueError('Unrecognized signal type detected in %s. Supported types are %s.' %
                         (str(sig_types), str(preferred_fig_order)))
    # order 'fig_order' based on 'preferred_order'
    fig_order = sorted(sig_types, key=lambda item: preferred_fig_order.index(item))
    print(fig_order)

    fig = plt.figure(constrained_layout=True, figsize=(6.4 * 0.8 * ncols, 4.8 * 0.65 * nrows))

    row = 0
    for i, expt_id in enumerate(expt_ids):
        print(expt_id)
        data_expt = all_data[all_data['Expt_ID'] == expt_id]
        sig_types = data_expt['Signal_Type'].unique()
        atp_assays = [sig_type for sig_type in sig_types if 'ATP' in sig_type]
        sig_types = [[sig_type] for sig_type in sig_types if 'ATP' not in sig_type]
        if len(atp_assays) > 1:
             sig_types += [atp_assays]

        # col = 0
        for sig_type in sig_types:
            print('   ', sig_type)
            col = fig_order.index(sig_type[0] if len(sig_type) == 1 else 'ATP')

            ax = fig.add_subplot(nrows, ncols, row * ncols + col + 1)

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
                plot_atp_rate_assay_data(ax, data_expt_sig, data_info)
                ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            # Plot Mito and Glyco Rate Assay data
            else:
                sig_type = sig_type[0]
                scale_factors = [data_expt_sig[n].iloc[0] if scale_to_first_point else 1
                                 for n in ['Signal_Mut', 'Signal_Corr', 'Signal_Ctrl']]
                sig_mutant = data_expt_sig['Signal_Mut'] / scale_factors[0]
                sig_corrected = data_expt_sig['Signal_Corr'] / scale_factors[1]
                sig_healthy = data_expt_sig['Signal_Ctrl'] / scale_factors[2]

                ax.plot(data_expt_sig['Time'], sig_mutant, 'o-', ms=10, label=data_info['mutation'] + '-')
                ax.plot(data_expt_sig['Time'], sig_corrected, 's-', ms=10, label='Corrected')
                ax.plot(data_expt_sig['Time'], sig_healthy, '^-', ms=10, label='Healthy')

                ax.set_xlabel('Time (%s)' % data_info['time_units'], fontsize=16)
                ax.set_ylabel('%s (%s)' % (sig_type, data_info['sig_units']) \
                                  if not scale_to_first_point else 'Relative %s' % sig_type, fontsize=16)

                ax.set_xlim(kwargs.get('xlim', (None, None)))
                ax.set_ylim(kwargs.get('ylim', (None, None)))

            # Set the plot title. Currently supported plot title types are: 'none', 'line_mut_type', and 'expt_it'.
            # The default is 'line_mut_type'.
            extra_space = '\n' if row > 0 else ''
            title_type = kwargs.get('title', 'line_mut_type')
            if title_type != 'none':
                title = extra_space
                if title_type == 'line_mut_type':
                    title += '%s, %s-deficient, %s' % \
                            (data_info['cell_line'], data_info['mutation'], data_info['cell_type'])
                elif title_type == 'expt_id':
                    title += expt_id
                else:
                    raise ValueError('Unrecognized title type: %s' % title_type)
                title_kw = kwargs.get('title_kw', {})
                ax.set_title(title, fontsize=16, **title_kw)

            ax.tick_params(axis='both', which='major', labelsize=16)
            ax.legend(loc='best', fontsize=12)

            # col += 1
        row += 1

    if save_plot is not False:
        outdir = '.' if save_plot is True else save_plot
        suffix = kwargs.get('suffix', None)
        filename = 'Seahorse_plots%s.png' % ('_' + suffix if suffix is not None else '')
        fig.savefig(os.path.join(outdir, filename), dpi=300)

    if show_plot:
        plt.show()

if __name__ == '__main__':
    import os

    data_path = '.'
    datafile = os.path.join(data_path, 'seahorse_data.csv')
    expt_ids = 'all'  # ['PD20_050125', 'PD20_122124', 'PD20_031325', 'PD20_050125']

    plot_seahorse_data(datafile, expt_ids=expt_ids, datatype='all', scale_to_first_point=True, show_plot=True,
                       title='expt_id', title_kw={'fontweight': 'normal', 'color': 'k'}, xlim=(-2, 102))
    # title='expt_id' 'line_mut_type'
    # ['mitoATP', 'glycoATP']

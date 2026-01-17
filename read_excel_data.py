from plot_data import load_seahorse_excel, export_seahorse_kinetic_to_csv, export_seahorse_atp_rate_to_csv, make_expt_id
import pandas as pd
import os
import glob
import re

base_path = 'Excel_Data'
directories = ['PD20', 'PD220']

OVERWRITE = True

cell_lines = ['PD20', 'PD220']
mutations = ['FANCD2', 'FANCA']
cell_types = ['Fibroblast', 'Fibroblast']
corrected = ['RV', 'LV']

labels_dict = {
    'PD220': 'FANCA-',
    'PD20': 'FANCD2-',
    'LV': 'Corrected',
    'RV': 'Corrected'
}

out_csv = os.path.join(base_path, "fa_seahorse_export.csv")
if OVERWRITE and os.path.exists(out_csv):
    os.remove(out_csv)

for directory, cell_line, mutation, cell_type, corr in zip(directories, cell_lines, mutations, cell_types, corrected):

    path = os.path.join(base_path, directory)

    files = glob.glob(os.path.join(path, '*.xlsx'))
    for file in files:
        print('filename:', file)

        df_clean, meta, extras = load_seahorse_excel(file, remove_outliers_flag=True, outlier_k=1.5) # 1.5) # 3

        # for key in meta.keys():
        #     print('%s: %s' % (key, meta[key]))
        # print(df_clean.columns)
        # quit()

        # get date of the experiment, if it's in the assay name
        expt_date = re.search(r'(\d+\.\d+\.\d+)', meta['Assay Name'])
        suffix = '' if expt_date is None else str('_%s' % expt_date.group(1)).replace('.', '')

        condition_order = list(dict.fromkeys(df_clean['Condition']))
        # Convert to an ordered categorical
        df_clean['Condition'] = pd.Categorical(df_clean['Condition'], categories=condition_order, ordered=True)

        # Map three conditions to the mutant/corrected/control roles
        condition_map = {"mut": cell_line, "corr": corr, "ctrl": "Healthy"}

        print('Outliers:')
        if 'Time' in df_clean.columns:
            if extras.get('removed', None) is not None:
                print(extras['removed'][['Condition', 'Time', 'Well', 'Value']])

            # Export data to a CSV file
            exported = export_seahorse_kinetic_to_csv(
                df=df_clean,
                meta=meta,
                out_csv=out_csv,
                condition_map=condition_map,
                mutation=mutation,
                cell_line=cell_line,
                cell_type=cell_type,
                expt_id=make_expt_id(cell_line + suffix),
                write_mode='append',  # "append" or "overwrite"
                float_format="%.8f"
            )
        else:
            if extras.get('removed', None) is not None:
                print(extras['removed'][['Condition', 'Component', 'Well', 'Value']])

            # Export data to a CSV file
            exported = export_seahorse_atp_rate_to_csv(
                df=df_clean,
                meta=meta,
                out_csv=out_csv,
                condition_map=condition_map,
                mutation=mutation,
                cell_line=cell_line,
                cell_type=cell_type,
                expt_id=make_expt_id(cell_line + suffix),
                write_mode='append',  # "append" or "overwrite"
                float_format="%.8f"
            )

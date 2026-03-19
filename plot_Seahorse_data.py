import pandas as pd
from read_excel_data import read_excel_data
from plot_utils import plot_seahorse_data

# USER-DEFINED VARIABLES
data_path = 'Excel_Data'  # directory where the Excel spreadsheets are located
cell_lines = ['PD20', 'PD220']
mutations = ['FANCD2', 'FANCA']
cell_types = ['Fibroblast', 'Fibroblast']
corrected = ['RV', 'LV']

OVERWRITE = True  # overwrite existing `seahorse_export.csv` file?
suptitles = ['FANCD2-deficient (PD20)', 'FANCA-deficient (PD220)']  # title for each figure

scale_to_first_point = True  # set to False for RAW data, set to True for RELATIVE data

kwargs = {
    'figsize_adjust': (1.8, 2.2),  # increase or decrease figure size
    'title': 'none',  # 'expt_id' 'line_mut_type' 'none'
    'title_kw': {'fontweight': 'normal', 'color': 'k'},  # for individual plot titles
    'xlim': (-2, 102),  # x-axis limits
    'fontsizes': {'labels': 30, 'ticks': 30, 'legend': 20, 'title': 30, 'suptitle': 40},
    'sharey': True  # share y-axes among all plots?
}

# CODE TO PLOT DATA (DON'T MODIFY)
datafile = read_excel_data(data_path, cell_lines, cell_lines, mutations, cell_types, corrected, overwrite=OVERWRITE)
df = pd.read_csv(datafile)
for cell_line, suptitle in zip(cell_lines, suptitles):
    expt_ids = df.loc[df['Expt_ID'].str.startswith(cell_line), 'Expt_ID'].unique().tolist()
    kwargs['suptitle'] = suptitle
    kwargs['suffix'] = '%s_%s' % (cell_line, 'RAW' if not scale_to_first_point else 'REL')
    plot_seahorse_data(datafile, expt_ids=expt_ids, datatype='all', scale_to_first_point=scale_to_first_point,
                       show_plot=True, squeeze_plots=True, **kwargs)
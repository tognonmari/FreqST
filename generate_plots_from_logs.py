import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd
import parser_for_logs as pfl
import numpy as np
import pickle
import matplotlib.transforms as mtrans
import matplotlib.patches as mpatches
RESULTING_STATS_PKL = './resulting_statistics.pkl'
AGGREGATED_STATS_PKL = './aggregated_runs.pkl'
DATASETS = ['athens_small', 'berlin_10', 'synth_100_100_5', 'chicago_4','chicago']#, 'chicago_4', 'chicago', 'berlin_10']

MODES = ['aided_means_random','aided_means','means']
#MODES = ['means', 'centers', 'aided_means', 'aided_centers', 'aided_means_random]

DATASETS_ROOT_FOLDERS = {'athens_small' : './athens_small', 'chicago': "./chicago", 'berlin_10' : './berlin_10', 'chicago_4' : './chicago_4', 'berlin' : './berlin', 'synth_100_100_5' : './synth_100_100_5'}
COSTS_CLASSIFICATION = ['LOW', 'MEDIUM']
COSTS = {'athens_small': {'means' : [[1,0.00003,128], [1,0.0003,128]], 'aided_means' : [[1,0.00003,128], [1,0.0003,128]], 'centers': [[1, 0.003, 128], [1,0.03, 128]], 'aided_centers': [[1, 0.003, 128], [1,0.03, 128]], 'aided_means_random': [[1, 0.00003, 128], [1, 0.0003, 128]]},
         'chicago' : {'means' : [[1,0.00003,888], [1,0.0003,888]], 'aided_means' : [[1,0.00003,888], [1,0.0003,888]], 'centers': [[1, 0.003, 888], [1,0.03, 888]], 'aided_centers': [[1, 0.003, 888], [1,0.03, 888]],'aided_means_random' : [[1,0.00003,888], [1,0.0003,888]]},
         'chicago_4' :  {'means' : [[1,0.00003,222], [1,0.0003,222]], 'aided_means' : [[1,0.00003,222], [1,0.0003,222]], 'centers': [[1, 0.1, 222], [1,1, 222]], 'aided_centers': [[1, 0.1, 222], [1,1, 222]], 'aided_means_random' : [[1,0.00003,222], [1,0.0003,222]]},
         'berlin_10' : {'means' : [[1,0.00003,2717], [1,0.0003,2717]], 'aided_means' : [[1,0.00003,2717], [1,0.0003,2717]], 'centers': [[1, 0.1,2717], [1,1, 2717]], 'aided_centers': [[1, 0.1, 2717], [1,1, 2717]], 'aided_means_random':  [[1,0.00003,2717], [1,0.0003,2717]]},
         'berlin' : {'means':[[1,0.00003,27188]], 'aided_means':  [[1,0.00003,27188]]},
         'synth_100_100_5' :{'aided_means_random' : [[1.0, 0.0003, 100],[1.0, 0.0001, 100]], 'means': [[1.0, 0.0003, 100],[1.0, 0.0001, 100]], 'aided_means': [[1.0, 0.0003, 100],[1.0, 0.0001, 100]]}
           
        }
MODES_CODES = {'means': 0, 'centers' : 1, 'aided_means':2, 'aided_centers' : 3, 'aided_means_random': 4}
#SIMPLIFICATION_FACTORS = {'means': [0.2,0.0], 'centers': [0.2,0.0], 'aided_means' : [0.0], 'aided_centers' :[0.0]}

SIMPLIFICATION_FACTORS = {'means' : [0.2], 'aided_means': [0.2,0.0],  'aided_means_random' :[0.2, 0.0]}

SEEDS = [0,1,2,3,4]

FREQUENCY_THRESHOLDS = [0.1]

K_RANDOM = [3,5]

def generate_tracker_logs_file_path(dataset_name, mode, costs, simplification, threshold, k, seed):

    folder = f"{DATASETS_ROOT_FOLDERS[dataset_name]}/logs"

    if not os.path.exists(folder):
        os.makedirs(folder)
    if(mode == 'aided_means_random'):
        return f"{folder}/tracker_{mode}_{'_'.join(str(cost) for cost in costs)}_{simplification}_{threshold}_{k}_{seed}.txt"
    elif (mode=='aided_means'):
        return f"{folder}/tracker_{mode}_{'_'.join(str(cost) for cost in costs)}_{simplification}_{threshold}.txt"
    # means as is 
    return f'{folder}/tracker_{mode}_{'_'.join(str(cost) for cost in costs)}_{simplification}.txt'

def retrieve_mode_from_list(arg_list):
    if arg_list[0]== 'means':
        return 'means', 1
    if arg_list[2] == 'random':
        return 'aided_means_random', 3
    return 'aided_means',2

def fill_dict_with_info_from_name(filename):
    statistics_dict = {'dataset_size': 0,'max_rss_gb': 0.0, 'dataset': '', 'run_id': '', 'mode' : '', 'frequency_threshold': -1.0, 'distances': [], 'fp_time': [], 'uncovered_pts_progression' :[], 'distances_progression' :[],  'efficacy': 0.0, 'total_time' : 0.0, 'costs' : [], 'k': -1, 'seed': -1, 'simplification_factor' : 0.2}
    statistics_dict['dataset'] = os.path.normpath(filename).split(os.sep)[-3]
    filename = os.path.basename(filename)

    statistics_dict['run_id']=filename
    arg_list = filename.replace(".txt", "").split('_')
    if('tracker' in arg_list):
        arg_list.remove('tracker')
    
    mode, first_cost_index = retrieve_mode_from_list(arg_list)
    
    statistics_dict['mode'] = mode
    statistics_dict['costs'] = tuple([float(a) for a in arg_list[first_cost_index:first_cost_index+3]])
    statistics_dict['disambiguation'] = 'LOW' if statistics_dict['costs'][1]<0.00029 else 'MEDIUM'

    statistics_dict['simplification_factor'] = float(arg_list[first_cost_index+3])
    if mode=='means':
        return statistics_dict
    statistics_dict['frequency_threshold'] = float(arg_list[first_cost_index+4])
    if mode=='aided_means':
        return statistics_dict
    statistics_dict['k']  = int(arg_list[first_cost_index+5])
    statistics_dict['seed']= int(arg_list[first_cost_index+6])
    
    return statistics_dict


def extract_statistics_from_file(filename)-> dict:

    print(f'Analyzing {filename}')
    statistics_dict = fill_dict_with_info_from_name(filename)

    statistics_dict['fp_selected_progression'] = []
    with open(filename, 'r') as logs_file:
        distance_to_fill = -1 
        head = [next(logs_file) for _ in range(2)]
        statistics_dict['dataset_size'] = int(head[1].rstrip())
        for line in logs_file:
            if line.startswith('Initializing'):
                statistics_dict['distances'].append(float(line.rstrip().split(" ")[-1]))
            elif line.strip().startswith('**************'):
                
                statistics_dict['fp_selected_progression'].append(True)
            elif line.startswith("TIME  FOR FP AT SQ DIST : "):
                if len(statistics_dict['fp_time'])==0:
                    statistics_dict['fp_time'] = len(statistics_dict['distances']) * [0.0]
                statistics_dict['fp_time'][distance_to_fill] = float(line.rstrip().split(" ")[-1])
                distance_to_fill = distance_to_fill-1    
            elif line.startswith("UNCOVERED POINTS"):
                statistics_dict['uncovered_pts_progression'].append(int(line.rstrip().split(" ")[-1]))
                if(len(statistics_dict['uncovered_pts_progression'])> len(statistics_dict['fp_selected_progression'])):
                    statistics_dict['fp_selected_progression'].append(False)
            elif line.startswith(' Efficacy:'):
                statistics_dict['efficacy'] = float(line.rstrip().split(" ")[-1])
            elif line.startswith('Total time:'):
                statistics_dict['total_time'] = float(line.rstrip().split(" ")[-1])
            elif line.startswith('	Maximum resident set size (kbytes)'):
                statistics_dict['max_rss_gb'] = float(line.rstrip().split(" ")[-1]) * (10 ** (-6))
            elif line.startswith('Best cluster has distance:'): #selected distances progressions
                statistics_dict['distances_progression'].append(float(line.rstrip().replace(',','').split(" ")[4]))
            else:
                continue
    assert(len(statistics_dict['uncovered_pts_progression'])== len(statistics_dict['fp_selected_progression']))
    return statistics_dict

def display_uncovered_pts_progression():
    # For each dataset and for each s.f. and for each cost vector, for fixed k and seed 
    # STEP 0: LOAD THE DF AND KEEP THE NEEDED ROWS (ONE K AND ONE SEED FOR THE RANDOM + AIDED MEANS)
    pickle_in = open(RESULTING_STATS_PKL, 'rb')
    aided_statistics_df = pickle.load(pickle_in)
    pickle_in.close()
    aided_statistics_df = aided_statistics_df[aided_statistics_df['mode']!= 'means']
    condition = ((aided_statistics_df['mode'] != 'aided_means_random' ) | (
    (aided_statistics_df['mode'] == 'aided_means_random') & (aided_statistics_df['seed'] == 0) & (aided_statistics_df['k'] == 5))) & aided_statistics_df['simplification_factor']>0
    aided_statistics_df = aided_statistics_df[condition]
    
    # STEP 1: LOAD THE TRACKER FILES FOR MEANS 
    tracker_files_lists = []
    for ds in DATASETS:
        for costs in COSTS[ds]['means']:
            for sf in SIMPLIFICATION_FACTORS['means']:
                file_path = generate_tracker_logs_file_path(ds,'means', costs,sf,0,1,0)
                tracker_files_lists.append(file_path)

    tracker_dict_list = []
    for file_path in tracker_files_lists:
        tracker_dict_list.append(extract_statistics_from_file(file_path))

    tracker_statistics_df = pd.DataFrame(tracker_dict_list)
    
    # STEP 1.1: MERGE THE TWO DATAFRAMES
    df =  pd.concat([tracker_statistics_df, aided_statistics_df], ignore_index=True)
    print(df[df['dataset']=='synth_100_100_05'])
    # STEP 2: PLOT ALL STUFF RELATIVE TO ONE DATASET AND ONE COST VECTOR ON THE SAME SUBPLOT.
    
    
    fig, axes = plt.subplots(len(DATASETS), 2, figsize=(12, 4 * len(DATASETS)), sharex=False, sharey=False)
    for i, dataset in enumerate(DATASETS):
        
        for j, cost_class in enumerate(COSTS_CLASSIFICATION):
            ax = axes[i, j]
            subset = df[(df['dataset'] == dataset) & ((df['disambiguation']==cost_class))]
            print(f"NOW LOOKING AT {dataset} for CLASS {cost_class}")
            print(f"THE SUBSET HAS {subset['mode']}")
            assert(len(subset.index)==3)
            dataset_size =  df[df['dataset']==dataset]['dataset_size'].iloc[0]
            
            for mode in subset['mode'].unique():
                mode_data = subset[subset['mode'] == mode]
                if not mode_data.empty:
                    # Assume one row per mode-cost_vector-dataset combo
                    uncovered = mode_data['uncovered_pts_progression'].iloc[0]
                    
                    uncovered = [int(u) for u in uncovered]
                    uncovered.insert(0,dataset_size)
                    covered = dataset_size - uncovered
                    
                    steps = np.arange(len(covered))
                    ax.plot(steps, covered, label=mode,marker='d', markersize=(3 if len(steps)<100 else 0))
            ax.axhline(y=dataset_size, color='k', linestyle='--')
            ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')
            # Titles and labels
            ax.set_title(f"{dataset} ")
            ax.set_xlabel("Step")
            ax.set_ylabel("Covered points")
            ax.legend()

    plt.tight_layout()
    fig.savefig('./images/coverage_per_iteration.png')

    # STEP 3: SAVE THE FIGURE 
    pass

def barcharts_with_costs_and_times():
    plt.rcParams['figure.dpi'] = 300
    pickle_in = open(AGGREGATED_STATS_PKL, 'rb')
    df = pickle.load(pickle_in)
    pickle_in.close()
    #fig, axes = plt.subplots(4,len(DATASETS), figsize=(5 * 4, 3 * len(DATASETS)), squeeze=False)
    fig, axes = plt.subplots(nrows=4,ncols=(len(DATASETS)+1), figsize=(5 * 4,( 3 * len(DATASETS))), squeeze=False)
    fig.subplots_adjust(left=0.1,  top=0.9, bottom=0, right=1)
    for i, identifier in enumerate(['LOW','MEDIUM']):
        for j, dataset in enumerate(DATASETS):
            # Filter data for this identifier and dataset
            
            subset = df[(df['disambiguation'] == identifier) & (df['dataset'] == dataset) & (df['simplification_factor']>0.1)]
            
            # Get the two axes: quality (even row), time (odd row)
            ax_quality = axes[2 * i][j]
            ax_time = axes[2 * i + 1][j]

            # Sort by mode for consistency
            subset = subset.sort_values(by='mode')
            yerr = []
            # Plot quality
            color='skyblue'
            sns.barplot(subset, hue='mode', x='mode', y='mean_efficacy',palette='Set2',ax=ax_quality,edgecolor='k') 
            ax_quality.grid(True, linestyle='--', linewidth=0.5, color='gray')
            ax_quality.set_xticks(range(len(subset['mode'])))
            ax_quality.set_xticklabels(['AM','AMR','M'])
            #ax_quality.set_title(f"Cost - {identifier} - {dataset}")
            ax_quality.set_ylabel("Cost")
            ax_quality.set_xlabel("Method")
            color_map = dict(zip(subset['mode'],sns.color_palette('Set2') ))
            hatch_map= { 'Frequent Pathlets\n Computation Time' : 'xx', 'Greedy Algorithm \nComputation Time':''}
            abbreviations= dict(zip(subset['mode'],['AM','AMR','M'] ))
            mode_patches = [
                mpatches.Patch(color=color_map[mode], label=abbreviations[mode])
                    for mode in subset['mode']
            ]
            
            '''ax_quality.bar(subset['mode'],subset['mean_efficacy'], yerr= subset['stddev_efficacy'],capsize=0, color=color)
            ax_quality.set_frame_on(True)'''
           
            # Plot timefor idx, row in subset.iterrows():
            #subset['mean_total_time'] = subset['mean_total_time'].apply(np.log10)
            
            #subset['mean_fp_time'] = subset['mean_fp_time'].apply(np.log10)
            bar2 = sns.barplot(subset,x='mode',hue='mode', y='mean_total_time', palette='Set2', ax=ax_time, edgecolor='k')
            bar1 = sns.barplot(subset, x='mode', hue='mode', y='mean_fp_time',palette='Set2',ax=ax_time, hatch='xx', edgecolor='k')
            bar2.set_yscale('log')
            bar1.set_yscale('log')
            #ax_time.set_title(f"Time - {identifier} - {dataset}")
            ax_time.set_ylabel("Time (s)")
            ax_time.set_xlabel("Method")
            ax_time.set_xticks(range(len(subset['mode'])))
            ax_time.set_xticklabels(['AM','AMR','M'])
            ax_time.grid(True, linestyle='--', linewidth=0.5, color='gray')
            abbreviations= dict(zip(subset['mode'],['AM','AMR','M'] ))
    mode_patches = [
        mpatches.Patch(color=color_map[mode], label=abbreviations[mode])
            for mode in subset['mode']
    ]
    
    hatch_patches=[]
    hatch_patches = [
        mpatches.Patch(facecolor='white', edgecolor='black', hatch=hatch_map[ft], label=ft)
        for ft in hatch_map.keys()
    ]
    
    legend_items = tuple(mode_patches + hatch_patches)
             

    for i in np.arange(0,4):
        axes[i][-1].set_frame_on(False)
        axes[i][-1].axis('off')
    fig.legend(handles=tuple(legend_items), title="Legend", loc='center right',fontsize='x-large',title_fontsize='large')
    plt.tight_layout()
    
    fig.subplots_adjust(left=0.1, top=0.9, right=1)
    column_labels = [' '.join(ds.split('_')).capitalize() for ds in DATASETS]
    row_labels = [r'Cost ($c_2=3\cdot 10^{-5}$)', r'Running Time($c_2=3\cdot 10^{-5}$)',r'Cost ($c_2=3\cdot 10^{-4}$)', r'Running Time ($c_2=3\cdot 10^{-4}$)']

    # Add column labels on top
    for col in range(0,len(DATASETS)):
        ax = axes[0, col]  # Top row, column `col`
        pos = ax.get_position()
        x = pos.x0 + pos.width / 2
        y = pos.y1 + 0.02  # Slightly above the top of the Axes
        fig.text(x, y, column_labels[col], ha='center', va='bottom', fontsize=14, )

    # Add row labels on the left
    for row in range(0,4):
        ax = axes[row, 0]  # Leftmost column, row `row`
        pos = ax.get_position()
        x = pos.x0 - 0.04  # Slightly to the left
        y = pos.y0 + pos.height / 2
        fig.text(x, y, row_labels[row], ha='center', va='center', fontsize=14, rotation=90)

    third_col_ax = axes[0, 2]  
    pos = third_col_ax.get_position()
    x = pos.x0 + pos.width / 2
    y = pos.y1 + 0.045 
    fig.text(x, y, 'Dataset', ha='center', va='bottom', fontsize=16, weight='medium')
    fig.text(x,y+ 0.025, 'Comparison with Means Clustering Baseline',ha='center', va='bottom', fontsize=16, weight='semibold')
    top_row_pos = axes[1, 0].get_position()
    second_row_pos = axes[2, 0].get_position()

    # Horizontal center of all columns:
    left = axes[0, 0].get_position().x0
    
    x = left - 0.06 

    # Vertical position: halfway between bottom of row 0 and top of row 1
    y = second_row_pos.y1 + (top_row_pos.y0 - second_row_pos.y1) / 2

    fig.text(x, y, 'Experimental Results', ha='center', va='center', fontsize=16, weight='medium',rotation=90)
    fig.savefig('./images/quality_time_modes_simpl.png')
    plt.clf()
  
    fig1, axes1 = plt.subplots(nrows=4,ncols=(len(DATASETS)+1), figsize=(5 * 4,( 3 * len(DATASETS))), squeeze=False)
    fig1.subplots_adjust(left=0.1,  top=0.9, bottom=0, right=1)
    for i, identifier in enumerate(['LOW','MEDIUM']):
        for j, dataset in enumerate(DATASETS):
            # Filter data for this identifier and dataset
            
            subset = df[(df['disambiguation'] == identifier) & (df['dataset'] == dataset) & (df['simplification_factor']!=0.2)]
            
            # Get the two axes: quality (even row), time (odd row)
            ax_quality = axes1[2 * i][j]
            ax_time = axes1[2 * i + 1][j]

            # Sort by mode for consistency
            subset = subset.sort_values(by='mode')
            yerr = []
            # Plot quality
            color='skyblue'
            sns.barplot(subset, hue='mode', x='mode', y='mean_efficacy',palette='Set2',ax=ax_quality,edgecolor='k') 
            ax_quality.grid(True, linestyle='--', linewidth=0.5, color='gray')
            ax_quality.set_xticks(range(len(subset['mode'])))
            ax_quality.set_xticklabels(['AM','AMR','M'])
            #ax_quality.set_title(f"Cost - {identifier} - {dataset}")
            ax_quality.set_ylabel("Cost")
            ax_quality.set_xlabel("Method")
            color_map = dict(zip(subset['mode'],sns.color_palette('Set2') ))
            hatch_map= { 'Frequent Pathlets\n Computation Time' : 'xx', 'Greedy Algorithm \nComputation Time':''}
            abbreviations= dict(zip(subset['mode'],['AM','AMR','M'] ))
            mode_patches = [
                mpatches.Patch(color=color_map[mode], label=abbreviations[mode])
                    for mode in subset['mode']
            ]
            
            '''ax_quality.bar(subset['mode'],subset['mean_efficacy'], yerr= subset['stddev_efficacy'],capsize=0, color=color)
            ax_quality.set_frame_on(True)'''
           
            # Plot timefor idx, row in subset.iterrows():
            #subset['mean_total_time'] = subset['mean_total_time'].apply(np.log10)
            
            #subset['mean_fp_time'] = subset['mean_fp_time'].apply(np.log10)
            bar2 = sns.barplot(subset,x='mode',hue='mode', y='mean_total_time', palette='Set2', ax=ax_time, edgecolor='k')
            bar1 = sns.barplot(subset, x='mode', hue='mode', y='mean_fp_time',palette='Set2',ax=ax_time, hatch='xx', edgecolor='k')
            bar2.set_yscale('log')
            bar1.set_yscale('log')
            #ax_time.set_title(f"Time - {identifier} - {dataset}")
            ax_time.set_ylabel("Time (s)")
            ax_time.set_xlabel("Method")
            ax_time.set_xticks(range(len(subset['mode'])))
            ax_time.set_xticklabels(['AM','AMR','M'])
            ax_time.grid(True, linestyle='--', linewidth=0.5, color='gray')
            abbreviations= dict(zip(subset['mode'],['AM','AMR','M'] ))
    mode_patches = [
        mpatches.Patch(color=color_map[mode], label=abbreviations[mode])
            for mode in subset['mode']
    ]
    
    hatch_patches=[]
    hatch_patches = [
        mpatches.Patch(facecolor='white', edgecolor='black', hatch=hatch_map[ft], label=ft)
        for ft in hatch_map.keys()
    ]
    
    legend_items = tuple(mode_patches + hatch_patches)
             

    
    
    for i in np.arange(0,4):
        axes1[i][-1].set_frame_on(False)
        axes1[i][-1].axis('off')
    fig1.legend(handles=tuple(legend_items), title="Legend", loc='center right',fontsize='x-large',title_fontsize='large')
    plt.tight_layout()
    
    fig1.subplots_adjust(left=0.1, top=0.9, right=1)
    column_labels = [' '.join(ds.split('_')).capitalize() for ds in DATASETS]
    row_labels = [r'Cost ($c_2=3\cdot 10^{-5}$)', r'Running Time($c_2=3\cdot 10^{-5}$)',r'Cost ($c_2=3\cdot 10^{-4}$)', r'Running Time ($c_2=3\cdot 10^{-4}$)']

    # Add column labels on top
    for col in range(0,len(DATASETS)):
        ax = axes1[0, col]  # Top row, column `col`
        pos = ax.get_position()
        x = pos.x0 + pos.width / 2
        y = pos.y1 + 0.02  # Slightly above the top of the Axes
        fig1.text(x, y, column_labels[col], ha='center', va='bottom', fontsize=14, )

    # Add row labels on the left
    for row in range(0,4):
        ax = axes1[row, 0]  # Leftmost column, row `row`
        pos = ax.get_position()
        x = pos.x0 - 0.04  # Slightly to the left
        y = pos.y0 + pos.height / 2
        fig1.text(x, y, row_labels[row], ha='center', va='center', fontsize=14, rotation=90)

    third_col_ax = axes1[0, 2]  
    pos = third_col_ax.get_position()
    x = pos.x0 + pos.width / 2
    y = pos.y1 + 0.045 
    fig1.text(x, y, 'Dataset', ha='center', va='bottom', fontsize=16, weight='medium')
    fig1.text(x,y+ 0.025, 'Comparison with Means Clustering Baseline',ha='center', va='bottom', fontsize=16, weight='semibold')
    top_row_pos = axes1[1, 0].get_position()
    second_row_pos = axes1[2, 0].get_position()

    # Horizontal center of all columns:
    left = axes1[0, 0].get_position().x0
    
    x = left - 0.06 

    # Vertical position: halfway between bottom of row 0 and top of row 1
    y = second_row_pos.y1 + (top_row_pos.y0 - second_row_pos.y1) / 2

    fig1.text(x, y, 'Experimental Results', ha='center', va='center', fontsize=16, weight='medium',rotation=90)

    fig1.savefig('./images/quality_time_modes.png')
    
    
def frequent_pathlet_selection_per_iter():
    # Gather Tracker info
    tracker_files = []
    for dataset in DATASETS:
        for simplification in [0.2]:
            for costs in COSTS[dataset]['means']:
                tracker_files.append(generate_tracker_logs_file_path(dataset,'means', costs,simplification, 0,0,0))

    for dataset in DATASETS:
        for simplification in [0.2]:
            for costs in COSTS[dataset]['aided_means_random']:
                tracker_files.append(generate_tracker_logs_file_path(dataset,'aided_means_random', costs,simplification, 0.1,5,0))
    tracker_stats_dict = []
    for file in tracker_files:
        tracker_stats_dict.append(extract_statistics_from_file(file))

    #repeat the same plots of the iteration vs coverage, but no marker and scatter markers separately
    df = pd.DataFrame(tracker_stats_dict)
    fig, axes = plt.subplots(len(DATASETS), 2, figsize=(12, 4 * len(DATASETS)), sharex=False, sharey=False)

    for i, dataset in enumerate(DATASETS):
        for j, cost_class in enumerate(COSTS_CLASSIFICATION):
            ax = axes[i, j]
            subset = df[(df['dataset'] == dataset) & ((df['disambiguation']==cost_class)) & df['simplification_factor']>0]
            #print(f"NOW LOOKING AT {dataset} for CLASS {cost_class}")
            #print(f"THE SUBSET HAS {subset['mode']}")
            assert(len(subset.index)==2)
            dataset_size =  df[df['dataset']==dataset]['dataset_size'].iloc[0]
            
            for mode in subset['mode'].unique():
                mode_data = subset[subset['mode'] == mode]
                if not mode_data.empty:
                    # Assume one row per mode-cost_vector-dataset combo
                    uncovered = mode_data['uncovered_pts_progression'].iloc[0]
                    
                    uncovered = [int(u) for u in uncovered]
                    #uncovered.insert(0,dataset_size)
                    covered = dataset_size - uncovered
                    
                    steps = np.arange(len(covered))
                    ax.plot(steps, covered, label=mode , markersize=(3 if len(steps)<100 else 0))
                    condition_highlight = np.array(mode_data['fp_selected_progression'].iloc[0])
                    print((condition_highlight))
                    ax.scatter(x=steps[condition_highlight], y=covered[condition_highlight], marker='o')
            ax.axhline(y=dataset_size, color='k', linestyle='--')
            ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')
            # Titles and labels
            ax.set_title(f"{dataset} ")
            ax.set_xlabel("Step")
            ax.set_ylabel("Covered points")
            ax.legend()

    plt.tight_layout()
    fig.savefig('./images/fp_selection_per_iteration.png')
    pass

def pile_up_of_frequent_pathlets():
    
    pass

def main():
    plt.rcParams['figure.dpi'] = 300
    #display_uncovered_pts_progression()
    barcharts_with_costs_and_times()
    #frequent_pathlet_selection_per_iter()


if __name__=='__main__':
    main()
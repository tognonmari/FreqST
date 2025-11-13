import pandas as pd
import os
import tqdm
import pickle
import numpy
DATASETS = ['athens_small', 'berlin_10', 'synth_100_100_5', 'chicago_4','chicago']#, 'chicago_4', 'chicago', 'berlin_10']

MODES = ['aided_means_random','aided_means','means']
#MODES = ['means', 'centers', 'aided_means', 'aided_centers', 'aided_means_random]

DATASETS_ROOT_FOLDERS = {'athens_small' : './athens_small', 'chicago': "./chicago", 'berlin_10' : './berlin_10', 'chicago_4' : './chicago_4', 'berlin' : './berlin', 'synth_100_100_5' : './synth_100_100_5'}

COSTS = {'athens_small': {'means' : [[1,0.00003,128], [1,0.0003,128]], 'aided_means' : [[1,0.00003,128], [1,0.0003,128]], 'centers': [[1, 0.003, 128], [1,0.03, 128]], 'aided_centers': [[1, 0.003, 128], [1,0.03, 128]], 'aided_means_random': [[1, 0.00003, 128], [1, 0.0003, 128]]},
         'chicago' : {'means' : [[1,0.00003,888], [1,0.0003,888]], 'aided_means' : [[1,0.00003,888], [1,0.0003,888]], 'centers': [[1, 0.003, 888], [1,0.03, 888]], 'aided_centers': [[1, 0.003, 888], [1,0.03, 888]],'aided_means_random' : [[1,0.00003,888], [1,0.0003,888]]},
         'chicago_4' :  {'means' : [[1,0.00003,222], [1,0.0003,222]], 'aided_means' : [[1,0.00003,222], [1,0.0003,222]], 'centers': [[1, 0.1, 222], [1,1, 222]], 'aided_centers': [[1, 0.1, 222], [1,1, 222]], 'aided_means_random' : [[1,0.00003,222], [1,0.0003,222]]},
         'berlin_10' : {'means' : [[1,0.00003,2717], [1,0.0003,2717]], 'aided_means' : [[1,0.00003,2717], [1,0.0003,2717]], 'centers': [[1, 0.1,2717], [1,1, 2717]], 'aided_centers': [[1, 0.1, 2717], [1,1, 2717]], 'aided_means_random':  [[1,0.00003,2717], [1,0.0003,2717]]},
         'berlin' : {'means':[[1,0.00003,27188]], 'aided_means':  [[1,0.00003,27188]]},
         'synth_100_100_5' :{'aided_means_random' : [[1.0, 0.0003, 100],[1.0, 0.0001, 100]], 'means': [[1.0, 0.0003, 100],[1.0, 0.0001, 100]], 'aided_means': [[1.0, 0.0003, 100],[1.0, 0.0001, 100]]}
           
        }
MODES_CODES = {'means': 0, 'centers' : 1, 'aided_means':2, 'aided_centers' : 3, 'aided_means_random': 4}
#SIMPLIFICATION_FACTORS = {'means': [0.2,0.0], 'centers': [0.2,0.0], 'aided_means' : [0.0], 'aided_centers' :[0.0]}

SIMPLIFICATION_FACTORS = {'means' : [0.0,0.2], 'aided_means': [0.2,0.0],  'aided_means_random' :[0.2, 0.0]}

SEEDS = [0,1,2,3,4]

FREQUENCY_THRESHOLDS = [0.1]

K_RANDOM = [3,5]


def generate_logs_file_path(dataset_name, mode, costs, simplification, threshold, k, seed):

    folder = f"{DATASETS_ROOT_FOLDERS[dataset_name]}/logs"

    if not os.path.exists(folder):
        os.makedirs(folder)
    if(mode == 'aided_means_random'):
        return f"{folder}/{mode}_{'_'.join(str(cost) for cost in costs)}_{simplification}_{threshold}_{k}_{seed}.txt"
    elif (mode=='aided_means'):
        return f"{folder}/{mode}_{'_'.join(str(cost) for cost in costs)}_{simplification}_{threshold}.txt"
    # means as is 
    return f'{folder}/{mode}_{'_'.join(str(cost) for cost in costs)}_{simplification}.txt'

def generate_out_file_path(dataset_name, mode, costs, simplification, threshold, k, seed):

    folder = f"{DATASETS_ROOT_FOLDERS[dataset_name]}/out"

    if not os.path.exists(folder):
        os.makedirs(folder)
    if(mode == 'aided_means_random'):
        return f"{folder}/{mode}_{'_'.join(str(cost) for cost in costs)}_{simplification}_{threshold}_{k}_{seed}.txt"
    elif (mode=='aided_means'):
        return f"{folder}/{mode}_{'_'.join(str(cost) for cost in costs)}_{simplification}_{threshold}.txt"
    # means as is 
    return f'{folder}/{mode}_{'_'.join(str(cost) for cost in costs)}_{simplification}.txt'

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
    
    mode, first_cost_index = retrieve_mode_from_list(arg_list)
    statistics_dict['mode'] = mode
    statistics_dict['costs'] = tuple([float(a) for a in arg_list[first_cost_index:first_cost_index+3]])
    statistics_dict['simplification_factor'] = float(arg_list[first_cost_index+3])
    statistics_dict['disambiguation'] = 'LOW' if statistics_dict['costs'][1]<0.00025 else 'MEDIUM'
    
    if mode=='means':
        return statistics_dict
    statistics_dict['frequency_threshold'] = float(arg_list[first_cost_index+4])
    if mode=='aided_means':
        return statistics_dict
    statistics_dict['k']  = int(arg_list[first_cost_index+5])
    statistics_dict['seed']= int(arg_list[first_cost_index+6])
    
    return statistics_dict


def extract_statistics_from_file(filename)-> dict:

    
    statistics_dict = fill_dict_with_info_from_name(filename)

    with open(filename, 'r') as logs_file:
        distance_to_fill = -1 
        head = [next(logs_file) for _ in range(2)]
        statistics_dict['dataset_size'] = int(head[1].rstrip())
        for line in logs_file:
            if line.startswith('Initializing'):
                statistics_dict['distances'].append(float(line.rstrip().split(" ")[-1]))
                
            elif line.startswith("TIME  FOR FP AT SQ DIST : "):
                if len(statistics_dict['fp_time'])==0:
                    statistics_dict['fp_time'] = len(statistics_dict['distances']) * [0.0]
                statistics_dict['fp_time'][distance_to_fill] = float(line.rstrip().split(" ")[-1])
                distance_to_fill = distance_to_fill-1
                
            elif line.startswith("UNCOVERED POINTS"):
                statistics_dict['uncovered_pts_progression'].append(int(line.rstrip().split(" ")[-1]))
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
    statistics_dict['total_fp_time'] = sum(statistics_dict['fp_time'])/1000
    return statistics_dict




def main():
    files_to_examine = []
    for dataset in DATASETS:
        for mode in MODES:
            for costs in COSTS[dataset][mode]:
                for simplification in SIMPLIFICATION_FACTORS[mode]:
                    if mode !='means':
                        for frequency_threshold in FREQUENCY_THRESHOLDS:

                            if mode !='aided_means':
                                for k in K_RANDOM:
                                    for seed in SEEDS:
                                        files_to_examine.append(generate_logs_file_path(dataset, mode, costs,simplification,frequency_threshold,k,seed))
                            else:
                                files_to_examine.append(generate_logs_file_path(dataset, mode, costs,simplification,frequency_threshold,0,0))
                    else:
                        files_to_examine.append(generate_logs_file_path(dataset, mode, costs,simplification,0,0,0))
    
    stat_dict_list = []

    for i in tqdm.tqdm(range(0, len(files_to_examine)), desc="Analyzing files : "):
        try:
            
            stat_dict_list.append(extract_statistics_from_file(files_to_examine[i]))
        except FileNotFoundError:
            print(f'something went wrong with {files_to_examine[i]}\n')
    
    
    stats_df = pd.DataFrame(stat_dict_list)
    stats_df.to_csv("./resulting_statistics.csv",float_format="%.5f")
    pickle_out = open("./resulting_statistics.pkl", "wb")
        
    pickle.dump(stats_df, pickle_out)
    pickle_out.close()
    
    #Compute means and stddev of times and costs of runs with mode aidedmean random that share the same dataset, cost vector, ft, sf
    grouped_df = stats_df[stats_df['mode']=='aided_means_random']
    grouped_df = grouped_df.groupby(by=['dataset', 'simplification_factor','frequency_threshold', 'costs','mode','disambiguation'])
    
    
    mean_score = grouped_df['efficacy'].mean()
    stddev_score = grouped_df['efficacy'].sem()

    mean_time = grouped_df['total_time'].mean()
    stddev_time = grouped_df['total_time'].sem()
    mean_fp_time=grouped_df['total_fp_time'].mean()
    
    summary_df = pd.DataFrame({'mean_efficacy': mean_score, 'stddev_efficacy': stddev_score, 'mean_total_time': mean_time, 'stddev_total_time': stddev_time, 'mean_fp_time': mean_fp_time}).reset_index()

    remaining_runs = stats_df[stats_df['mode']!='aided_means_random'].groupby(by=['dataset', 'simplification_factor','frequency_threshold', 'costs','mode','disambiguation'])
    
    remaining_summary_df = pd.DataFrame()
    remaining_summary_df['mean_efficacy'] = remaining_runs['efficacy'].mean()
    remaining_summary_df['mean_total_time'] = remaining_runs['total_time'].mean()
    remaining_summary_df['stddev_efficacy'] = 0.0
    remaining_summary_df['stddev_total_time'] = 0.0
    remaining_summary_df['mean_fp_time'] = remaining_runs['total_fp_time'].mean()
    remaining_summary_df = remaining_summary_df.reset_index()
    
    merged_df = pd.concat([remaining_summary_df, summary_df], ignore_index=True).sort_values(by='dataset')
    
    merged_df.to_csv('./aggregated_runs.csv', float_format="%.2f")
    merged_df.drop(columns=['disambiguation','frequency_threshold']).to_latex('./aggregated_runs.tex', float_format="%.2f", index=False)
    
    for dataset in DATASETS:
            filtered_df= merged_df[(merged_df['dataset']==dataset)&(merged_df['simplification_factor']<0.1)]
            filtered_df.sort_values(by='mode').to_latex(f'./tables/{dataset}_{simplification}.tex', index=False)

    for dataset in DATASETS:
            filtered_df= merged_df[(merged_df['dataset']==dataset)&(merged_df['simplification_factor']>0.1)]
            filtered_df.sort_values(by='mode').to_latex(f'./tables/{dataset}_{simplification}.tex', index=False)
            
    filtered_df= merged_df[(merged_df['simplification_factor']<0.1)& (merged_df['dataset']=='synth_100_100_5') ]
    filtered_df.sort_values(by='mode').to_latex(f'./tables/synth_100_100_5_0.0.tex', index=False)
    pickle_out = open("./aggregated_runs.pkl", "wb")
        
    pickle.dump(merged_df, pickle_out)
    pickle_out.close()
    
    
if __name__ == '__main__':
    main()
    

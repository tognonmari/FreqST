import subprocess
import sys
from pathlib import Path
import os
import logging


SAMPLER_EXECUTABLE = "./chernoff_sampler"

PM_EXECUTABLE = "./helloworld"

EPSILON_LIST = [0.05, 0.1, 0.15, 0.2]

DELTA_LIST = [0.025, 0.05, 0.1]

DATASET_FOLDERS = ["./athens_small", "./chicago", "./berlin"]

RADIUS_LIST = [50, 100, 200]

FREQUENCY_LIST = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]

LOG_FILE = "./experiment_runner.log"

OUTPUT_REDIRECT_FOLDER = "stdout"

SAMPLE_FOLDER_NAME = "samples"

RESULT_FOLDER_NAME = "results"

def main():

    logger = logging.getLogger('ALGO ST EXPERIMENTS')
    logger.setLevel(logging.INFO)
    #Generate the necessary samples
    for epsilon in EPSILON_LIST:
        for delta in DELTA_LIST:
            for data_folder in DATASET_FOLDERS:

                if not os.path(f"{data_folder}/{SAMPLE_FOLDER_NAME}").is_dir():
                    os.makedirs(f"{data_folder}/{SAMPLE_FOLDER_NAME}")
            # Generate one sample for now (seed 0)
                try:
                    str_exec = f"{SAMPLER_EXECUTABLE} -e {epsilon} -d {delta} {data_folder}/merged.txt {data_folder}/{SAMPLE_FOLDER_NAME}"
                    child_process = subprocess.run(str_exec, capture_output = True, text = True)

                except subprocess.CalledProcessError as e:
                    logger.error(f"Error executing {str_exec} : {e}")

                    continue
    
    # Not all the samples guarantee to be created, because chernoff smaple is loose

    # check utility dirs exist

    # Execute the pattern matching on the samples

    for dataset in DATASET_FOLDERS:

        p = f"{dataset}/{SAMPLE_FOLDER_NAME}"

        for sample in Path(p).iterdir():

            if(sample.is_dir()):

                logger.info("there's an unexpected directory in the samples.")
                continue

            for freq in FREQUENCY_LIST:

                for rad in RADIUS_LIST:

                    try:
                        #ensure out dir and file exist
                        

                        res_file_name = f"{sample.replace(".txt", "")}_{rad}_{freq}.txt"
                        outfilepath = Path(f"{dataset}/{OUTPUT_REDIRECT_FOLDER}/{res_file_name}")
                        
                        if not outfilepath.exists():
                            handler = open(outfilepath, 'w')
                            handler.close()

                        str_exec = f"{PM_EXECUTABLE} -r {rad} -f {freq} {p}/{sample} {dataset}/merged.txt {dataset}/{RESULT_FOLDER_NAME}/{res_file_name} > {dataset}/{OUTPUT_REDIRECT_FOLDER}/{res_file_name}"
                        child_process = subprocess.run(str_exec, capture_output=True, text=True)
                            

                    except subprocess.CalledProcessError as e:
                        logger.error(f"Error executing {str_exec} : {e}")
                        logger.error(child_process.stderr)



if __name__ == "__main__":

    main()
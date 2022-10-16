import os
import multiprocessing
import argparse
import shutil
import time
import subprocess
import logging
import string
import random

'''
This module contains miscellaneous functions

Most of the codes are inspired by Unicycler according to the terms of the GNU General Public License
https://github.com/rrwick/Unicycler/blob/master/unicycler/misc.py
'''

def random_string(n=6):
    """
    Generate random string with acsii lowercase letter

    Args:
        n (int, optional): length of the random string. Defaults to 6.
    """    
    return ''.join(random.choices(string.ascii_lowercase, k=n))

def current_time():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    
def default_threads_number():
    set_maximum_threads_number = 8
    return min(int(multiprocessing.cpu_count()/2), set_maximum_threads_number)

def write_log(screen_log, log_file, level=logging.INFO):
    """
    Print log while write into log_file

    Args:
        screen_log (str): Information to log
        log_file (str): log file name
        level (logging, optional): _description_. Defaults to logging.INFO.
    """    
    
    print(screen_log)
    logging.basicConfig(filename=log_file, filemode="a", level=logging.INFO)
    logging.log(level=level, msg=screen_log)

    return

def cmd(command, description, log_file):

    screen_log = "Start {description}\n".format(time=current_time(), description=description)
    write_log(screen_log, log_file)
    write_log(command, log_file)

    # Run command
    ret = subprocess.run(command,shell=True)
    if ret.returncode == 0:
        # Print end
        screen_log = "\nFinish {description}\n".format(description=description)
        write_log(screen_log, log_file)
    
    else:
        screen_log = "ERROR: {description} failed, for more details, please check the log {log}".format(description=description, log=log_file)
        
        logging.basicConfig(filename=log_file, filemode="a", level=logging.ERROR)
        logging.info(screen_log)
        exit(screen_log)
   
    return

# pipeline_utils.py
"""
Utility functions for pipeline management and configuration handling.
This module provides functions to log status messages, run shell commands,
validate project names, and load configuration files in various formats.
It is designed to be used in a pipeline context where these utilities can help
manage the execution of steps and maintain a clear log of operations.
"""
import subprocess
from datetime import datetime
import os
import re
import json

try:
    import yaml
except ImportError:
    yaml = None
try:
    import tomllib  # Python 3.11+ for TOML support
except ImportError:
    tomllib = None

def log_status(log_file, message):
    """
    Append a timestamped status message to the log file.
    Opens the log file in append mode for each write to avoid conflicts in parallel runs.
    """
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log_file, 'a') as log:
        log.write(f"[{timestamp}] {message}\n")
        log.flush()

def run_command(command, step_name, log_file, critical=False):
    """
    Run a shell command and log its status (SUCCESS/FAILURE).
    If the command fails and `critical` is True, exit the program.
    """
    try:
        subprocess.run(command, shell=True, check=True)
        log_status(log_file, f"{step_name}: SUCCESS")
    except subprocess.CalledProcessError:
        log_status(log_file, f"{step_name}: FAILURE")
        print(f"Error: {step_name} failed. Check {log_file} for details.")
        if critical:
            exit(1)

def is_valid_project_name(project_name):
    """
    Validate project name: must consist of letters, numbers, underscores, 
    and not contain the substring 'NODE' (reserved for internal use).
    """
    return bool(re.match(r'^[A-Za-z0-9_]+$', project_name)) and "NODE" not in project_name

def load_config(config_path):
    """
    Load configuration parameters from a YAML, JSON, or TOML file.
    Returns a dictionary of config values. Requires PyYAML for .yaml files.
    """
    with open(config_path, 'r') as f:
        if config_path.endswith(('.yaml', '.yml')):
            if yaml is None:
                raise ImportError("PyYAML is not installed. Please install it to use YAML config files.")
            config = yaml.safe_load(f)
        elif config_path.endswith('.toml'):
            if tomllib is None:
                raise ImportError("TOML support is not available. Use Python 3.11+ or install a toml library.")
            config = tomllib.load(f)
        else:
            # Default to JSON
            config = json.load(f)
    return config

# pipeline_utils.py
"""
Utility functions for pipeline management, configuration handling, and
result packaging.

This module provides functions to log status messages, run shell commands,
validate project names, load configuration files in various formats, and bundle
pipeline outputs into a portable artifact with an integrity manifest.
"""
import subprocess
from datetime import datetime
import os
import re
import json
import hashlib
import shutil
import tempfile
import zipfile
from pathlib import Path

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


def _hash_file(file_path, chunk_size=65536):
    """Return the SHA-256 checksum for the provided file path."""
    digest = hashlib.sha256()
    with open(file_path, 'rb') as handle:
        for block in iter(lambda: handle.read(chunk_size), b''):
            digest.update(block)
    return digest.hexdigest()


def _collect_manifest_entries(root_dir, exclude_paths=None):
    """Return manifest entries for all files under ``root_dir`` except excluded ones."""
    entries = []
    exclude_paths = set(exclude_paths or [])
    for path in sorted(Path(root_dir).rglob('*')):
        if path.is_file():
            rel_path = path.relative_to(root_dir)
            rel_str = rel_path.as_posix()
            if rel_str in exclude_paths:
                continue
            entries.append(
                {
                    'path': rel_str,
                    'size': path.stat().st_size,
                    'sha256': _hash_file(path),
                }
            )
    return entries


def package_outputs(items, bundle_path, bundle_format='directory', metadata=None, manifest_name='manifest.json', overwrite=False):
    """Collect pipeline outputs into a directory or zip archive with a manifest.

    Parameters
    ----------
    items : sequence of tuple[str, str]
        Iterable of ``(relative_path, source_path)`` pairs describing which
        files or directories to include and where to place them within the
        bundle.
    bundle_path : str or Path
        Destination directory (for ``bundle_format='directory'``) or archive
        file. Parent directories are created automatically.
    bundle_format : {'directory', 'zip'}
        Output format. ``'directory'`` preserves the bundled layout on disk,
        while ``'zip'`` creates a compressed archive.
    metadata : dict, optional
        Arbitrary metadata to include in the manifest file.
    manifest_name : str, optional
        Name of the manifest file to generate inside the bundle.
    overwrite : bool, optional
        When ``True`` existing bundle directories or archives are removed
        before writing new contents. Otherwise a ``FileExistsError`` is raised.

    Returns
    -------
    Path
        Path to the resulting directory or archive.
    """

    if bundle_format not in {'directory', 'zip'}:
        raise ValueError("bundle_format must be either 'directory' or 'zip'")

    metadata = metadata or {}
    bundle_path = Path(bundle_path)

    if bundle_format == 'directory':
        if bundle_path.exists():
            if overwrite and bundle_path.is_dir():
                shutil.rmtree(bundle_path)
            elif overwrite and bundle_path.is_file():
                bundle_path.unlink()
            elif bundle_path.exists():
                raise FileExistsError(f"Bundle destination '{bundle_path}' already exists")
        bundle_root = bundle_path
        bundle_root.mkdir(parents=True, exist_ok=True)
        cleanup = None
    else:
        # Prepare temporary directory and ensure parent directory exists.
        if bundle_path.exists():
            if overwrite:
                if bundle_path.is_dir():
                    shutil.rmtree(bundle_path)
                else:
                    bundle_path.unlink()
            else:
                raise FileExistsError(f"Bundle destination '{bundle_path}' already exists")
        bundle_path.parent.mkdir(parents=True, exist_ok=True)
        temp_dir = tempfile.mkdtemp(prefix='sprout_bundle_')
        bundle_root = Path(temp_dir)
        cleanup = bundle_root

    for relative_path, source_path in items:
        src = Path(source_path)
        if not src.exists():
            raise FileNotFoundError(f"Bundle source '{source_path}' does not exist")
        dest = bundle_root / relative_path
        dest.parent.mkdir(parents=True, exist_ok=True)
        if src.is_dir():
            shutil.copytree(src, dest, dirs_exist_ok=True)
        else:
            shutil.copy2(src, dest)

    manifest_entries = _collect_manifest_entries(bundle_root)
    manifest = {
        'generated_at': datetime.utcnow().isoformat() + 'Z',
        'bundle_format': bundle_format,
        'metadata': metadata,
        'files': manifest_entries,
    }

    manifest_path = bundle_root / manifest_name
    with open(manifest_path, 'w', encoding='utf-8') as handle:
        json.dump(manifest, handle, indent=2)

    if bundle_format == 'zip':
        with zipfile.ZipFile(bundle_path, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(bundle_root.rglob('*')):
                archive.write(path, arcname=path.relative_to(bundle_root).as_posix())
        if cleanup and cleanup.exists():
            shutil.rmtree(cleanup)
        return bundle_path

    return bundle_root

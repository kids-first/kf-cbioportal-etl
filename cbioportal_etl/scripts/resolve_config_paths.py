import os

def resolve_config_paths(config, tool_dir):
    """
    Resolve paths dynamically based on assumptions:
    - Paths starting with 'scripts/' or 'REFS/' are relative to the tool directory.
    """
    for key, value in config.items():
        if isinstance(value, dict):
            resolve_config_paths(value, tool_dir)
        elif isinstance(value, str) and value.startswith(("REFS/", "scripts/", "external_scripts/")):
            config[key] = os.path.abspath(os.path.join(tool_dir, value))

    return config
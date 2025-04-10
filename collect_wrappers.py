import json
import os


def find_wrappers_recursive(current_dir="."):
    dirs = []
    for d in os.listdir(current_dir):
        subdir = os.path.join(current_dir, d)
        if os.path.isdir(subdir) and not d.startswith("."):
            if os.path.exists(os.path.join(subdir, "wrapper.py")):
                dirs += [subdir]
            else:
                dirs += find_wrappers_recursive(current_dir=subdir)
    return dirs


with open("workflows.json", "w") as f:
    json.dump({"workflow": find_wrappers_recursive()}, f)

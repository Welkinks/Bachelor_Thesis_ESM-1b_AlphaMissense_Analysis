import sys
from pathlib import Path

# sys.path - list of directories that Python searches for modules
def setup_project_root():
    project_root = Path.cwd()
    if "src" not in [p.name for p in project_root.iterdir()]:
        project_root = project_root.parent
    if str(project_root) not in sys.path:
        sys.path.append(str(project_root))


# In every notebook, run this setup first
"""

# --- Project Setup ---
from setup_notebook import setup_project_root
setup_project_root()

"""

### --- Explanations for Josef --- ###
# sys.path is a list of directories that Python searches for modules.
# notebook_path = os.getcwd() is the Current Working Directory.
# use if __name__ == "__main__": to ensure the script runs only when executed directly, not when imported.
#         setup_project_root()
#         print("✓ Project root added to sys.path")

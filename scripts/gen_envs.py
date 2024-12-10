import os
from pathlib import Path
import shutil

# Prompt the user for replacement choice
replace_choice = input(
    "Do you want to replace existing .env files with .env.template? (y - replace/n - copy nonexisting): ").strip().lower()

# Iterate through all .env.template files
for template_path in Path(".").rglob(".env.template"):
    dir_path = template_path.parent
    env_file = dir_path / ".env"

    if env_file.exists():
        if replace_choice in ["y", "yes"]:
            shutil.copy(template_path, env_file)
            print(f"Replaced existing .env with .env.template in {dir_path}")
        else:
            print(f".env already exists in {dir_path}, skipping...")
    else:
        shutil.copy(template_path, env_file)
        print(f"Copied .env.template to .env in {dir_path}")

# Get the home directory path
home = str(Path.home())
variable_name = os.path.join(home, ".nolabs")

# Update the variable_name value in .env files
for env_file in Path(".").rglob(".env"):
    with env_file.open("r") as file:
        lines = file.readlines()

    with env_file.open("w") as file:
        for line in lines:
            if line.startswith(f"NOLABS_HOME="):
                file.write(f"NOLABS_HOME={variable_name}\n")
                print(f"Updated NOLABS_HOME in {env_file} to {variable_name}")
            else:
                file.write(line)

# Ensure the directory exists
if not os.path.exists(variable_name):
    os.makedirs(variable_name)
    print(f"Created directory at {variable_name}")
else:
    print(f"Directory {variable_name} already exists, skipping creation.")

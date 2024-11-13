#!/bin/bash

# Create the directory if it does not already exist
if [ ! -d "$NOLABS_HOME" ]; then
    mkdir -p "$NOLABS_HOME"
    echo "Created directory at $NOLABS_HOME"
else
    echo "Directory $NOLABS_HOME already exists, skipping creation."
fi

# Prompt to ask if the user wants to replace all .env files with .env.template if they already exist
read -p "Do you want to replace existing .env files with .env.template? (y - replace/n - copy nonexisting): " replace_choice

# Find and process .env.template files
find . -type f -name ".env.template" | while read -r template_path; do
    dir=$(dirname "$template_path")
    env_file="$dir/.env"

    if [ -f "$env_file" ]; then
        # If the .env file exists and user chose to replace, copy the template
        if [[ "$replace_choice" == "y" || "$replace_choice" == "Y" ]]; then
            cp "$template_path" "$env_file"
            echo "Replaced existing .env with .env.template in $dir"
        else
            echo ".env already exists in $dir, skipping..."
        fi
    else
        # If .env does not exist, just copy the template
        cp "$template_path" "$env_file"
        echo "Copied .env.template to .env in $dir"
    fi
done

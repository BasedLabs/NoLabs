#!/bin/bash

read -p "Do you want to replace existing .env files with .env.template? (y - replace/n - copy nonexisting): " replace_choice

find . -type f -name ".env.template" | while read -r template_path; do
    dir=$(dirname "$template_path")
    env_file="$dir/.env"

    if [ -f "$env_file" ]; then
        if [[ "$replace_choice" == "y" || "$replace_choice" == "Y" ]]; then
            cp "$template_path" "$env_file"
            echo "Replaced existing .env with .env.template in $dir"
        else
            echo ".env already exists in $dir, skipping..."
        fi
    else
        cp "$template_path" "$env_file"
        echo "Copied .env.template to .env in $dir"
    fi
done

home=$(realpath "$HOME")
variable_name=$HOME/.nolabs
sed -i "s|^${variable_name}=.*|${variable_name}=${home}|" '.env'
    echo "Updated ${variable_name} in '.env' to ${home}"

if [ ! -d "$NOLABS_HOME" ]; then
    mkdir -p "$NOLABS_HOME"
    echo "Created directory at $NOLABS_HOME"
else
    echo "Directory $NOLABS_HOME already exists, skipping creation."
fi

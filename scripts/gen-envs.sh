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

if [ ! -d "${variable_name}" ]; then
    mkdir -p "${variable_name}"
    echo "Created directory at ${variable_name}"
else
    echo "Directory ${variable_name} already exists, skipping creation."
fi

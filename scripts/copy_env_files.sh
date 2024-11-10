#!/bin/bash

find . -type f -name ".env.template" | while read -r template_path; do
    dir=$(dirname "$template_path")
    env_file="$dir/.env"
    if [ ! -f "$env_file" ]; then
        cp "$template_path" "$env_file"
        echo "Copied .env.template to .env in $dir"
    else
        echo ".env already exists in $dir, skipping..."
    fi
done

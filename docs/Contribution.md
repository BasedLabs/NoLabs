# NoLabs Contributor Guide

## Forking the Repository
1. Visit the [NoLabs GitHub repository](https://github.com/BasedLabs/NoLabs).
2. Click the "Fork" button in the top-right corner of the page.
3. This will create a copy of the repository in your GitHub account.

## Making Changes
### Adding New Models
1. To add new models, create a folder within the microservices/ directory.
2. The folder should contain:
   - A build directory with a Dockerfile inside it.
   - A directory for the source code of the microservice.
   - A `client/` directory, generated using a command like this:
   ```shell
    npx @openapitools/openapi-generator-cli generate \
       -i https://127.0.0.1:5737/openapi.json \
       -g python \
       -o ./client \
       --additional-properties=packageName=diffdock_microservice
   ```
   
3. Ensure the microservice includes FastAPI's surface, and Uvicorn should start upon Docker run, similar to the example files:
   Dockerfile: https://github.com/BasedLabs/NoLabs/blob/master/microservices/esmfold/build/Dockerfile
   API file: https://github.com/BasedLabs/NoLabs/blob/master/microservices/esmfold/esmfold/api.py
4. Once everything is working, add the microservice to the docker-compose.yaml file in the main directory.
5. If changes are made to the microservice, update the version of the image in the docker-compose file for the respective microservice:
   ```yaml
   esmfold:
     image: 'ghcr.io/basedlabs/esmfold:1.0.0' # Update to higher version
   ```
## Changing Backend/Frontend in NoLabs
   Pull requests (PRs) should be either bug fixes, meaningful refactoring, or new features. 
   Reach out to [jaktenstid](https://github.com/jaktenstid) or [timurishmuratov7](https://github.com/timurishmuratov7) for significant new features or refactoring.
   Update the image version of NoLabs in the docker-compose file:
   ```yaml
   nolabs:
     image: 'ghcr.io/basedlabs/nolabs:1.2.1' # Update to higher version
   ```
### Opening a Pull Request
   1. Commit your changes using 
   ```
   git commit -m "Your commit message here"
   ```
   
   2. Push your changes to your forked repository with 
   ```
   git push origin <branch_name>
   ```
   
   3. Visit your forked repository on GitHub.
   4. Click on the "Pull Request" button.
   5. Provide a clear title and description for your pull request, explaining the changes you've made.
   6. Submit the pull request.

## Best Practices
- Use Conventional Commits: write descriptive commit messages using Conventional Commits. 
- Branch Naming: Name your branches as "feature/<feature_name>", "fix/<fix_name>", or "refactor/<refactor_name>".
- Python Formatting: Follow PEP8 formatting guidelines for Python code.
- Recommended Tools: We strongly recommend using PyCharm and Webflow when developing NoLabs.
- GitHub Workflows: For new microservices, add them to GitHub workflows (master template and PR template).

## Collaborating
   Feel free to discuss your changes with the community in the pull request comments.
   Be open to feedback and iterate on your changes based on suggestions.
   Once your pull request is approved, it will be merged into the main repository.

#!/usr/bin/env fish
# Initialise git and push to a PRIVATE GitHub repository.
# Requires git; the push step uses the GitHub CLI (gh). Install gh with:
#   sudo apt install gh    then    gh auth login
cd (dirname (status --current-filename))

git init -b main
git add .
git commit -m "Initial analysis: DESeq2 results, Wnt target and skeletal gene set tests"

if type -q gh
    gh repo create craniofacial-wnt-rnaseq --private --source . --push
else
    echo "gh not found. Create an empty private repo on github.com, then run:"
    echo "  git remote add origin git@github.com:YOUR_USER/craniofacial-wnt-rnaseq.git"
    echo "  git push -u origin main"
end

https://techwritingmatters.com/how-to-update-your-forked-repository-on-github


In the command line, move to your repository folder: cd repositoryname.
Run git fetch upstream. This fetches all the changes from the original repo.
Run git checkout develop. This assumes the main branch of your repository is called develop. It may be called master, so type the appropriate command. This switches you to that branch in case you’re not there already.
Run git merge upstream/develop. Again, substitute develop for whatever your main branch is. This merges all the changes from the original repo to your local clone.
Run git push. This pushes the changes from your local clone to your forked repo in GitHub. If for some reason this doesn’t work, try git push origin develop (or master).
Refresh your Github page and your fork should now be even (in sync) with the original repo.

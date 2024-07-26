# Globnut env
Repo contains code for the article "Eurasian plant diversity will benefit from lowering N deposition and preserving N/P variation" (in prep).

## Reproducibility
Users are encouraged to use the renv package to make sure they all use the same package versions.
renv will run checks and let you know if the packages you are using are the same versions as your collaborators'. If you install additional packages, let renv save that package name and version by using renv::snapshot() and collaborators will receive this information through git. renv will then recommend collaborators to run renv::restore() to install missing packages.

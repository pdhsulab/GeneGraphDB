# Developer Environment 

## Docker

Use docker to build a consistent development environment.  Docker is a lightweight virtual machine ("container"),
that lets your run code in a performant sandbox with well defined software dependencies.
Please skim some online intros (medium, youtube, official site...) until you're satisfied you get the gist.

Docker ports cleanly between a developer laptop and a cloud production server.
Minimizing the gap between these environments reduces the surface area for deployment bugs.

Docker makes sharing code between developers easier (e.g. if your new commit requires SQLAlchemy to work,
instead of requesting everyone install SQLAlchemy over Slack, you can simply update the docker image and assume
everyone will pull the new image developers when your code is merged.)
Building a consistent local development environment allows for faster iteration and reduces developer friction.

Whenever possible, "pin" the version of dependencies  (e.g. python==3.9.9, or numpy==1.22.0).
This reduces the chance of bugs caused by dependency updates (e.g. breaking API change from numpy 2.0.0).
Although it's a generally a good policy to keep dependencies up-to-date (mostly for security reasons),
leaving dependencies "floating" (to always grab latest version) can discourage developers from updating the docker
image for fear of tripping over dependency-update bugs.


## setup.sh

In this directory, setup.sh includes aliases for
running a local docker container (`genegraphdb_docker_run`),
building a docker image locally (`genegraphdb_docker_build`),
or running a local jupyter notebook (`genegraphdb_docker_jupyter`).

To access these aliases, please run `source ./developer_env/setup.sh` 
and consider adding this line to your .bashrc/.zshrc or similar.

Sourcing setup.sh will also change your local git settings to automatically strip output from ipython notebooks.
This decreases the size of files committed to github, decreases the risk that sensitive information is committed,
and decreases the noise of commits (ipython notebook outputs can change without any code changing).


## CloudBuild

CloudBuild is a tool offered from GCP for CI/CD ("continuous integration continuous deploy").
For the typical engineering team, this means building/compiling code, running tests, and deploying the built code to
e.g. a web server.  For now, we're just using it to build docker images.

### Building a new base image

```
git tag 'build_base' -f
git push origin build_base -f
```

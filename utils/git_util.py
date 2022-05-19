import os
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


class RunGitCommandError(Exception):
    pass


def get_branch_name():
    return _run_git_command("rev-parse --abbrev-ref HEAD")


def get_default_commit_sha():
    """A wrapper for `get_commit_sha` that catches errors

    Returns the results of `git_commit_sha` unless the git command raises an error,
    then the function returns a default string.
    """
    try:
        sha_string = get_commit_sha()
    except RunGitCommandError:
        sha_string = ">no git"

    return sha_string


def get_commit_sha():
    return _run_git_command("rev-parse --short HEAD")


def get_full_commit_sha():
    return _run_git_command("rev-parse HEAD")


def does_commit_exist_on_remote(commit_sha):
    """Checks if a commit exists in the remote repository

    `git branch -r --contains {commit_sha}` returns a string of all branches
    that contains the commit corresponding to the sha.

    If the git command returns an empty string (meaning: the commit doesn't
    exist on any branch in the remote repository) then this function returns
    False.

    Git will throw an error if the commit does not exist in the remote nor
    the local repository. Git will throw an error if the commit sha is malformed.
    """
    s = _run_git_command("branch -r --contains {}".format(commit_sha)).strip()
    return bool(s)


def _run_git_command(command):
    git_command = "git {}".format(command)
    cmd = "cd {} && {}".format(SCRIPT_DIR, git_command)

    try:
        output = subprocess.getoutput(cmd)
    except subprocess.CalledProcessError:
        raise RunGitCommandError('Git command "{}" threw an error'.format(git_command))

    return output.rstrip()

import os


def create_dir(dir_name):
    os.system("mkdir " + dir_name)


def download_dir(bucket_dir):
    bucket_path = "gs://jluo_bucket/" + bucket_dir
    command = "gsutil cp -m " + bucket_path + " ."
    os.system(command)

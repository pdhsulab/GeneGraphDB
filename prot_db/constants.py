import datetime

GCP_PROJECT_ID = "durrant"
GCS_BUCKET_NAME = "gs://durrant_prot_pred_db"
LOCAL_DATA_DIR = "/GeneGraphDB/data/durrant_prot_pred_db"

TIMESTAMP_FORMAT = "%Y%m%dT%H:%M:%S"


def get_timestamp(dt: datetime.datetime = None):
    if dt is None:
        dt = datetime.datetime.utcnow()
    return dt.utcnow().strftime(TIMESTAMP_FORMAT)

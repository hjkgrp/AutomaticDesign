from molSimplifyAD.utils.pymongo_tools import push2db
from .job_converter import loop_convert_jobs


def push_job_manager_jobs(basedir, database, collection, tag, subtag,
                          publication=False, user=False, pwd=False,
                          host="localhost", port=27017, auth=False,
                          all_runs_pickle=False, update_fields=False,
                          outpath='/home/db_backup'):
    print("converting jobs to dftrun objects...")
    runs = loop_convert_jobs(basedir)
    print("pushing results...")
    push2db(database, collection, tag, subtag,
            publication=publication, user=user,
            pwd=pwd, host=host,
            port=port, auth=auth,
            all_runs_pickle=all_runs_pickle,
            update_fields=update_fields,
            all_runs_list=runs,
            outpath=outpath)

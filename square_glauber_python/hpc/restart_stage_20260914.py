"""Targeted five-stream checkpoint extension; preserves earlier stages."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from square_glauber.restart_campaign import frozen_case, run_sample_worker, consolidate_case
from square_glauber.perfect_campaign import _atomic_json
from square_glauber.perfect import MonotoneCFTPSampler, CFTPNotCoalesced
from square_glauber.restart import sample_with_checkpoint

ROOT=Path('/nobackup/fvkl37/square_glauber_perfect_covariance_20260910')
MANIFEST=ROOT/'campaign_manifest.csv'
STAGE=ROOT/'restart_streams_1073741824_20260914'
CASES=[(670,1250),(682,1314),(855,812),(889,1015)]
OLD_CAP=536870912
NEW_CAP=1073741824

def idle():
    if subprocess.check_output(['squeue','-h','-u',os.environ['USER']],text=True).strip():
        raise RuntimeError('require empty queue before this stage')

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    action=sys.argv[1]
    if action=='prepare':
        idle()
        summary=json.loads((ROOT/'campaign_summary.json').read_text())
        assert summary['complete_environments']==23796 and summary['censored_environments']==4
        assert not STAGE.exists()
        plan=[]
        for task,index in CASES:
            case,env,metadata,seeds,backend=frozen_case(MANIFEST,task,ROOT,index)
            assert not (case/'pair.json').exists()
            for sample,seed in [(1,seeds.chain1_seed),(2,seeds.chain2_seed)]:
                if (case/f'sample{sample}.json').exists(): continue
                checkpoint=case/f'sample{sample}_restart.json'
                before=digest(checkpoint)
                sampler=MonotoneCFTPSampler(env,backend=backend,chunk_attempts=4000000)
                try:
                    sample_with_checkpoint(sampler,seed,checkpoint_path=checkpoint,max_horizon_sweeps=OLD_CAP)
                except CFTPNotCoalesced as error:
                    assert error.horizon_sweeps==OLD_CAP
                else: raise AssertionError('expected completed failed horizon')
                assert digest(checkpoint)==before
                plan.append(dict(task_id=task,environment_index=index,sample_index=sample,
                    environment_id=metadata['environment_id'],cftp_seed=seed,
                    environment_sha256=metadata['environment_sha256'],
                    checkpoint_path=str(checkpoint),checkpoint_sha256=before))
        assert len(plan)==5
        STAGE.mkdir()
        for name in ['campaign_summary.json','campaign_status.csv','center_measurements.csv',
                     'spatial_measurements.csv','chain_diagnostics.csv','censored_environment_ids.txt',
                     'incomplete_task_ids.txt','incomplete_slurm_array_spec.txt']:
            shutil.copy2(ROOT/name,STAGE/name)
        for row in plan:
            shutil.copy2(row['checkpoint_path'],STAGE/f"{row['environment_id']}-sample{row['sample_index']}-checkpoint.json")
        _atomic_json(STAGE/'plan.json',plan)
        print(json.dumps(plan,indent=2))
    elif action=='submit':
        idle()
        plan=json.loads((STAGE/'plan.json').read_text()); assert len(plan)==5
        assert not (STAGE/'jobs.json').exists()
        for row in plan:
            assert digest(Path(row['checkpoint_path']))==row['checkpoint_sha256']
        analysis=ROOT/'analysis_frozen_restart_20260914'; assert not analysis.exists()
        os.environ.update(OUTPUT_DIR=str(ROOT),MANIFEST=str(MANIFEST),ANALYSIS_DIR=str(analysis))
        job=subprocess.check_output(['sbatch','--parsable','--array=0-4%5',
            '--export=ALL','hpc/restart_stream_task_20260914.slurm'],text=True).strip()
        _atomic_json(STAGE/'jobs.json',dict(worker_job=job))
        finish=subprocess.check_output(['sbatch','--parsable','--dependency=afterany:'+job,
            '--export=ALL','hpc/restart_stream_finish_20260914.slurm'],text=True).strip()
        _atomic_json(STAGE/'jobs.json',dict(worker_job=job,finish_job=finish))
        print(job,finish)
    elif action=='worker':
        row=json.loads((STAGE/'plan.json').read_text())[int(os.environ['SLURM_ARRAY_TASK_ID'])]
        assert digest(Path(row['checkpoint_path']))==row['checkpoint_sha256']
        result=run_sample_worker(MANIFEST,row['task_id'],ROOT,row['environment_index'],row['sample_index'],
            max_horizon_sweeps=NEW_CAP,chunk_attempts=4000000)
        print(json.dumps(result,indent=2))
    elif action=='consolidate':
        for task,index in CASES:
            print(consolidate_case(MANIFEST,task,ROOT,index))
    else: raise ValueError(action)

if __name__=='__main__': main()

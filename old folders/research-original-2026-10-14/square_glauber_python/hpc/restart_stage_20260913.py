"""One-off frozen 11-stream extension. No old worker may overlap this stage."""
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
STAGE=ROOT/'restart_streams_536870912_20260913'
CASES=[(98,1908),(149,727),(516,484),(670,1250),(682,1314),
       (736,99),(831,669),(855,812),(889,1015)]

def main():
    action=sys.argv[1]
    if action=='prepare':
        if subprocess.check_output(['squeue','-h','-u',os.environ['USER']],text=True).strip():
            raise RuntimeError('require empty job queue before legacy-to-restart handoff')
        summary=json.loads((ROOT/'campaign_summary.json').read_text())
        assert summary['complete_environments']==23791 and summary['censored_environments']==9
        assert not STAGE.exists()
        plan=[]
        for task,index in CASES:
            case,env,metadata,seeds,backend=frozen_case(MANIFEST,task,ROOT,index)
            assert not (case/'pair.json').exists()
            for sample,seed in [(1,seeds.chain1_seed),(2,seeds.chain2_seed)]:
                if (case/f'sample{sample}.json').exists(): continue
                legacy=case/'censoring'/f'sample{sample}_cap-268435456.json'
                # Validate legacy history without running a single update.
                assert not (case/f'sample{sample}_restart.json').exists()
                sampler=MonotoneCFTPSampler(env,backend=backend,chunk_attempts=4000000)
                try:
                    sample_with_checkpoint(sampler,seed,checkpoint_path=case/f'sample{sample}_restart.json',
                        max_horizon_sweeps=268435456,legacy_censor_path=legacy)
                except CFTPNotCoalesced as error:
                    assert error.horizon_sweeps==268435456
                else: raise AssertionError('expected validated censor')
                plan.append(dict(task_id=task,environment_index=index,sample_index=sample,
                    environment_id=metadata['environment_id'],cftp_seed=seed,
                    environment_sha256=metadata['environment_sha256'],legacy_censor_path=str(legacy)))
        assert len(plan)==11
        STAGE.mkdir()
        for name in ['campaign_summary.json','campaign_status.csv','center_measurements.csv',
                     'spatial_measurements.csv','chain_diagnostics.csv','censored_environment_ids.txt',
                     'incomplete_task_ids.txt','incomplete_slurm_array_spec.txt']:
            shutil.copy2(ROOT/name,STAGE/name)
        _atomic_json(STAGE/'plan.json',plan)
        print(json.dumps(plan,indent=2))
    elif action=='submit':
        plan=json.loads((STAGE/'plan.json').read_text()); assert len(plan)==11
        assert not (STAGE/'jobs.json').exists()
        if subprocess.check_output(['squeue','-h','-u',os.environ['USER']],text=True).strip():
            raise RuntimeError('unexpected jobs at submission')
        analysis=ROOT/'analysis_frozen_restart_20260913'; assert not analysis.exists()
        os.environ.update(OUTPUT_DIR=str(ROOT),MANIFEST=str(MANIFEST),ANALYSIS_DIR=str(analysis))
        job=subprocess.check_output(['sbatch','--parsable','--array=0-10%11',
            '--export=ALL','hpc/restart_stream_task.slurm'],text=True).strip()
        _atomic_json(STAGE/'jobs.json',dict(worker_job=job))
        finish=subprocess.check_output(['sbatch','--parsable','--dependency=afterany:'+job,
            '--export=ALL','hpc/restart_stream_finish.slurm'],text=True).strip()
        _atomic_json(STAGE/'jobs.json',dict(worker_job=job,finish_job=finish))
        print(job,finish)
    elif action=='worker':
        row=json.loads((STAGE/'plan.json').read_text())[int(os.environ['SLURM_ARRAY_TASK_ID'])]
        result=run_sample_worker(MANIFEST,row['task_id'],ROOT,row['environment_index'],row['sample_index'],
            max_horizon_sweeps=536870912,chunk_attempts=4000000,
            legacy_censor_path=row['legacy_censor_path'])
        print(json.dumps(result,indent=2))
    elif action=='consolidate':
        for task,index in CASES:
            print(consolidate_case(MANIFEST,task,ROOT,index))
    else: raise ValueError(action)

if __name__=='__main__': main()

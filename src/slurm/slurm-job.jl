using DrWatson

const DEFAULT_HOST_NAME = "vacc-user2.uvm.edu"
const DEFAULT_USER = "bcase"
const DEFAULT_EMAIL = "bcase@uvm.edu"
const DEFAULT_EMAIL_TYPE = "ALL"
const DEFAULT_TIME = "00:15:00"
const DEFAULT_JULIA_EXEC = "/users/b/c/bcase/julia-1.7.1/bin/julia"

last_job = nothing

get_last_job() = last_job === nothing ? ErrorException("No last script") : last_job

"""
Write and submit a given julia `script` with some #SBATCH commands\n
this script should use `addprocs` to use `srun` with the specified resources
"""
function submit_slurm_job(
    script=get_last_job(), nodes=1, tasks_per_node=2;
    host_name=DEFAULT_HOST_NAME,
    user=DEFAULT_USER,
    email=DEFAULT_EMAIL,
    email_type=DEFAULT_EMAIL_TYPE,
    time=DEFAULT_TIME,
    julia_exec=DEFAULT_JULIA_EXEC
)
    # job_name = split(script, ".jl")[1]
    f = """
    #!/bin/bash                                                                    
    #SBATCH --nodes=$nodes                                
    #SBATCH --ntasks=$(nodes*tasks_per_node)                        
    #SBATCH --ntasks-per-node=$tasks_per_node 
    #SBATCH --cpus-per-task=1                                                                 
    #SBATCH --mem-per-cpu=4G                                                       
    #SBATCH --time=$time
    #SBATCH --job-name=$job_name                                            
    #SBATCH --output=$job_name/$job_name.out                                           
    #SBATCH --mail-user=$email                                     
    #SBATCH --mail-type=$email_type

    $julia_exec -t 4 $(scriptsdir())/$script \$SLURM_NTASKS
    """

    last_job = job_name
    dest = "$user@$host_name"

    run(`git add slurm-job scripts/$script`)
    run(`git commit -m "sending $script to cluster; git push`)
    run(`ssh $dest "git pull; sbatch $job_name/slurm-job; exit"`)
    # run(`ssh $dest "mkdir $job_name; exit"`)
    # run(`scp $script $(dest):$job_name/$script`) # send julia script over to new directory ~/$job_name
    # write("slurm-job", f) # local slurm script
    # run(`scp slurm-job $(dest):$job_name/slurm-job`) # send slurm script over
    # run(`ssh $dest "sbatch $job_name/slurm-job; exit"`) # submit job
    # rm("slurm-job")
end

# TODO up to user to make sure fout is unique atm
function get_slurm_output(
    fout, job=get_last_job(), outdir="results"; cleanup=false, host_name=DEFAULT_HOST_NAME, user=DEFAULT_USER
)
    dest = "$user@$host_name"
    run(`scp $(dest):$fout $(outdir)/$fout`)
    if cleanup 
        run(`ssh $dest "rm -r julia-*.out $job/; exit"`)
    end
end
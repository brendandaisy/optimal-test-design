using Distributed
                                                                                                             
hosts = []                                                                     
pids = []                                                                      
for i in workers()                                                             
    host, pid, t = fetch(                                                      
          @spawnat i (gethostname(), getpid())
)                                                                              
    println("$host: $pid with t=$t")
    push!(hosts, host)
    push!(pids, pid)
end
                                                                               
# The Slurm resource allocation is released when all the workers have          
# exited                                                                       
for i in workers()                                                             
    rmprocs(i)                                                                 
end
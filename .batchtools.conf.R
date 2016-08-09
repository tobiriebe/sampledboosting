cluster.functions = batchtools::makeClusterFunctionsSlurm("/home/hpc/ua341/di25koz/sampledboosting/slurm_lmulrz.tmpl", 
clusters = "serial")
default.resources = list(walltime = 300L, memory = 512L, ntasks = 1L)

max.concurrent.jobs = 1000L

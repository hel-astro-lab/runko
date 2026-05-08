# for dry-run add -n as an option
rsync -ahmr \
    --verbose \
    --info=progress2 \
    --exclude='*/restart/' \
    --exclude='*/ram-usage/' \
    --include='*/' \
    --include-from=include.rsync  \
    --exclude='*' \
    --exclude=jobs \
    lumi:/users/jnattila/scratch_26/runko-gpu/projects/pic-shock/ .

    #turso:/home/jnattila/runko-gpu/projects/pic-shock/ .

    #hile:/wrk-kappa/users/jnattila/runko/projects/pic-shocks/ .
    #jnattila@turso.cs.helsinki.fi:/wrk-kappa/users/jnattila/runko/projects/pic-shocks/ .
    #jnattila@gateway.flatironinstitute.org:ceph/runko/projects/shocks3d/ .
    #-e 'ssh -p 61022' \
    #jnattila@gateway.flatironinstitute.org:ceph/runko/projects/shocks3d/ .
    #jnattila@rusty:ceph/runko/projects/shocks3d/ .
    #jnattila@beskow.pdc.kth.se:/cfs/klemming/nobackup/j/jnattila/runko/projects/turbulence/ .



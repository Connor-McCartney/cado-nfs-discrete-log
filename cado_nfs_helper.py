import subprocess
cado_path = "/home/connor/Documents/cado-nfs/cado-nfs.py"

def cado_dlog_unknown_base(p, ell, target):
    out = subprocess.run(["python", f"{cado_path}", "-dlp", 
                          f"{ell=}", f"{target=}", f"{p}"], 
                         capture_output=True)
    for line in out.stderr.decode().splitlines():
        logbase_ = line.split("logbase = ")
        logtarget_ = line.split("log(target) = ")
        if len(logbase_) > 1:
            logbase = int(logbase_[-1])
        if len(logtarget_) > 1:
            logtarget = int(logtarget_[-1].split()[0])

    assert 1 == pow(pow(logbase, logtarget, p) * pow(target, -1, p), (p-1)//ell, p)
    return logbase, logtarget

def cado_dlog(p, ell, a, base):
    lb1, log_a = cado_dlog_unknown_base(p, ell, a)
    lb2, log_b = cado_dlog_unknown_base(p, ell, base)
    assert lb1 == lb2
    dlog = (log_a * pow(log_b, -1, ell)) % ell
    return dlog

def pohlig_hellman():
    #todo

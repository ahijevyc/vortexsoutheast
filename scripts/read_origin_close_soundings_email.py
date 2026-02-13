import datetime
import re
import sys

"""
Input format example
Earl98; upshr BMX 9/3/00z, dnshr JAX 9/3/12z
Frances98; upshr BRO 9/12/00z, dnshr LCH 9/12/00z

Output format example

Earl,19980903T00,upshr,BMX
Earl,19980903T12,dnshr,JAX
Frances,19980912T00,upshr,BRO
Frances,19980912T00,dnshr,LCH

"""
ifile = sys.argv[1]
with open(ifile, encoding="utf-8") as f:
    for line in f:
        line = line.strip()  # no newline
        s = re.search(r"^\w+;\s+", line)
        stormyearsemicolon = s.group()
        line = line.replace(stormyearsemicolon, "")
        stormyear = stormyearsemicolon.rstrip("; ")
        storm = stormyear[:-2]
        year = int(stormyear[-2:])
        if year < 70:
            year += 2000
        else:
            year += 1900

        parts = line.split(",")
        for p in parts:
            p = p.strip()
            s = re.search(r"([a-z]+) ([A-Z][A-Z][A-Z]) (\S+)", p)
            shr, stid, t = s.groups()
            assert shr == "dnshr" or shr == "upshr", f"unexpected shear direction {shr}"
            month, day, hour = t.split("/")
            assert hour.endswith("z")
            hour = hour.rstrip("z")
            t = datetime.datetime(year, int(month), int(day), int(hour))
            tstr = t.strftime("%Y%m%dT%H")
            print(f"{storm},{tstr},{shr},{stid}")

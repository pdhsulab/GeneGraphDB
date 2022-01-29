def revcomp(s):
    o = ""
    for c in s[::-1]:
        if c == "A":
            o += "T"
        elif c == "T":
            o += "A"
        elif c == "C":
            o += "G"
        elif c == "G":
            o += "C"
        else:
            o += c
    return o

import os
BASE = 'F:\\nano-everything\\mini-automation-theory\\2. mini-time-domain-analysis\\mini-steady-state-error'
SRC = os.path.join(BASE, "src")
def w(fn, content):
    with open(os.path.join(SRC, fn), "w") as f: f.write(content)
    return content.count(chr(10)) + 1

# Auto-generated source files
total = 0
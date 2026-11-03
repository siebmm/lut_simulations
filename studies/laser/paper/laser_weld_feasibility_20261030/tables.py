"""Generate paper tables from the portable, unrounded evidence."""
from figures import *
def num(v,n=2):return '--' if v is None else f'{v:.{n}f}'
def texpath(p):return '$'+label(p).replace(' | ',r'\mid ')+'$'
def write(name,rows):
    (HERE/name).write_text('\n'.join(' & '.join(r)+r' \\' for r in rows)+'\n',encoding='utf-8')
def build():
    rows=[]
    for p in CORE:
        f,pa=[metric(k,p) for k in ['full','partial']]
        rows.append([texpath(p)]+[num(m[k],1 if k=='near_support_pct' else 2) for k in ['cmax_db','c99_db','near_support_pct'] for m in [f,pa]])
    write('table_core.tex',rows)
    rows=[]
    for name,title in [('baseline','Equal / 0.5--5'),('band1','Equal / 1--5'),('weighted','Response-weighted / 0.5--5'),('band1_weighted','Response-weighted / 1--5')]:
        m=metric('partial','1-2_SSS',name);rows.append([title,num(m['near']['um_s']),num(m['cmax_db']),num(m['c99_db']),num(m['peak']['distance_mm'],3)])
    write('table_improvement.tex',rows)
    rows=[]
    for p in CORE:
        m=metric('full',p,'shoulders_plate');rows.append(['Full '+texpath(p),'Shoulders; plates',num(m['near_support_pct'],1),num(m['cmax_db']),num(m['c99_db'])])
    for name,title in [('weighted_plate','All five; plates'),('weighted_central_plate','Floor three; plates'),('weighted_shoulders','Shoulders; all')]:
        m=metric('partial','1-2_SSS',name);rows.append(['Partial '+texpath('1-2_SSS'),title,num(m['near_support_pct'],1),num(m['cmax_db']),num(m['c99_db'])])
    write('table_joint.tex',rows)
    # CSV remains useful outside a LaTeX installation.
    import csv
    rows=[]
    for r in INDEX:
        for case,m in r['metrics'].items():
            rows.append({k:r[k] for k in ['uid','stage','kind','path','name']}|{'case':case}|{k:m[k] for k in ['support_pct','near_support_pct','cmax_db','c99_db']}|{'near_um_s':m['near']['um_s'] if m['near'] else None,'peak_distance_mm':m['peak']['distance_mm'] if m['peak'] else None})
    with (HERE/'all_results.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
if __name__=='__main__':build()

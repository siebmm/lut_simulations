"""Typeset measured values directly from the frozen-study analysis."""
from study import *
AN={k:json.loads((OUT/k/'analysis.json').read_text()) for k in DATA}
def run(k,p,n='baseline'):return next(r for r in AN[k]['runs'] if r['path']==p and r['name']==n)
def m(k,p,n='baseline',c='lof'):return run(k,p,n)['signals'][c]
def f(v,n=2):return '--' if v is None else f'{v:.{n}f}'
def val(a,key):return a[key]['um_s'] if a[key] else None
def pl(p):return '$'+label(p).replace(' | ',r'\mid ')+'$'
def output(name,rows):
    (HERE/('table_'+name+'.tex')).write_text('\n'.join(' & '.join(str(v) for v in row)+r' \\' for row in rows)+'\n',encoding='utf-8')

for kind in DATA:
    rows=[];controls=[];fine=[]
    for p in PATHS:
        a=m(kind,p)
        rows.append([pl(p),f(a['support_pct'],1),f(a['near_support_pct'],1),f(val(a,'near'),3),f(a['peak']['distance_mm'] if a['peak'] else None,3),f(a['cmax_db']),f(a['c99_db'])])
        b=m(kind,p,c='reference');c=m(kind,p,c='difference')
        controls.append([pl(p),f(val(a,'near'),3),f(val(b,'near'),3),f(val(c,'near'),3),f(c['near_support_pct'],1),f(c['cmax_db']),f(c['c99_db'])])
        if kind=='partial' or p in CORE:
            b=m(kind,p,'fine');fine.append([pl(p),f(val(a,'peak'),3),f(val(b,'peak'),3),f(b['peak']['distance_mm'] if b['peak'] else None,3),f(a['c99_db']),f(b['c99_db'])])
    output(kind+'_baseline',rows);output(kind+'_controls',controls);output(kind+'_fine',fine)
rows=[]
for p in PATHS:
    a,b=[m(k,p) for k in DATA]
    rows.append([pl(p),f(a['near_support_pct'],1),f(b['near_support_pct'],1),f(a['cmax_db']),f(b['cmax_db']),f(a['c99_db']),f(b['c99_db'])])
output('combined',rows)
rows=[]
for r in json.loads((OUT/'matched_comparison.json').read_text()):
    a,b=r['full'],r['partial'];rows.append([pl(r['path']),f(a['support_pct'],1),f(a['near_support_pct'],1),f(a['cmax_db']),f(b['cmax_db']),f(a['c99_db']),f(b['c99_db'])])
output('matched',rows)
choices=[('coarse','Baseline'),('no_coherence','Coherence off'),('screen_off','Direct screening off'),('screen_ssss','Add SSSS screening'),('screen_reflected','Add all reflected screening'),('gate10',r'Angular gate 10\%'),('gate30',r'Angular gate 30\%'),('offset1','Exclusion 1 mm'),('offset2','Exclusion 2 mm'),('receivers4','Minimum 4 receivers'),('receivers16','Minimum 16 receivers'),('one_source','Minimum 1 source'),('three_sources','Minimum 3 sources'),('band15','Band 1--5 MHz'),('tail20','Pulse tail 20 dB'),('tail40','Pulse tail 40 dB'),('shoulders','Shoulder sources only'),('floor_sources','Fill-floor sources only')]
rows=[]
for n,title in choices:
    a=m('partial','1-2_SSS',n);rows.append([title,f(a['support_pct'],1),f(a['near_support_pct'],1),f(a['cmax_db']),f(a['c99_db'])])
output('sensitivity',rows)
rows=[]
for p in PATHS:
    a=m('full',p,'coarse');b=m('full',p,'exclude_repairs');rows.append([pl(p),f(a['support_pct'],1),f(b['support_pct'],1),f(a['c99_db']),f(b['c99_db']),f(b['peak']['distance_mm'] if b['peak'] else None,3)])
output('repairs',rows)
for kind in DATA:
    rows=[]
    for p in WIDE:
        a=m(kind,p,'overview');q=a['peak'];rows.append([pl(p),f(a['support_pct'],1),f(q['x_mm'] if q else None),f(q['depth_mm'] if q else None),f(val(a,'peak'),3),f(a['cmax_db']),f(a['c99_db'])])
    output(kind+'_wide',rows)
print('Tables generated')

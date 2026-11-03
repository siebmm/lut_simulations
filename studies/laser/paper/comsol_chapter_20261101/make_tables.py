"""Generate every reported score from archived, untruncated image arrays."""
from pathlib import Path
import json
import numpy as np
import make_artifacts as f
HERE=Path(__file__).resolve().parent;DATA=HERE/'data'

def rows(name):return json.loads((DATA/name).read_text())['rows']
def metric(rr,id,case='lof'):return next(r for r in rr if r['id']==id and r['case']==case)
def num(x,d=2):return '--' if x is None else f'{x:.{d}f}'
def pathlabel(p):return '$'+f.label(p).replace(' | ',r'\mid ')+'$'
def textfile(name,lines):
    path=HERE/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main():
    full=rows('full_confirmation_metrics.json');partial=rows('partial_confirmation_metrics.json');sel=rows('full_selected_metrics.json')
    primary=[('Full',metric(sel,'all'),'0--5; 257'),('Partial',metric(partial,'baseline__1-2_SSS__fine'),'0.5--5; 657')]
    macros=[]
    for kind,r,_ in primary:
        for name,v,dec in [('Contrast',r['cmax_db'],2),('Percentile',r['c99_db'],2),('Coverage',r['near_support_pct'],1),('Peak',r['peak']['um_s'],2),('Distance',r['peak']['distance_mm'],3)]:
            macros.append(r'\newcommand{\comsol'+kind+name+'}{'+num(v,dec)+'}')
    selected=metric(partial,'baseline__1-2_SSS__fine')
    for case,prefix in [('lof',''),('reference','Pristine'),('difference','Difference')]:
        r=metric(partial,'baseline__1-2_SSS__fine',case)
        for name,v,dec in [('PeakExact',r['peak']['um_s'],3),('X',r['peak']['x_mm'],2),
            ('Z',r['peak']['depth_mm'],2),('DistanceExact',r['peak']['distance_mm'],4),
            ('ContrastExact',r['cmax_db'],2)]:
            macros.append(r'\newcommand{\comsolPartial'+prefix+name+'}{'+num(v,dec)+'}')
    for prefix,id in [('Finer','baseline__1-2_SSS__finer'),('Wide','baseline__1-2_SSS__wide'),
        ('FloorSource','baseline__1-2_SSS__central'),('Plate','baseline__1-2_SSS__plate'),
        ('FloorReceiver','baseline__1-2_SSS__floor'),('BandOneFive','band1_5__1-2_SSS__fine'),
        ('ShortBand','bp_short_matched_guard__1-2_SSS__fine'),('ShortLow','lp_short_matched_guard__1-2_SSS__fine')]:
        r=metric(partial,id)
        for name,v,dec in [('Contrast',r['cmax_db'],2),('Coverage',r['near_support_pct'],1)]:
            macros.append(r'\newcommand{\comsolPartial'+prefix+name+'}{'+num(v,dec)+'}')
    macros.append(r'\newcommand{\comsolPartialDisplayReference}{'+f'{f.PARTIAL_REFERENCE:g}'+'}')
    textfile('comsol_values.tex',macros)
    table=[]
    for kind,r,band in primary:
        table.append(f'{kind} & '+('$S\\mid S$' if kind=='Full' else '$S\\mid SS$')+' & '+band+' & '+num(r['cmax_db'])+' & '+num(r['c99_db'])+' & '+num(r['near_support_pct'],1)+' & '+num(r['peak']['distance_mm'],3)+r' \\')
    textfile('tables/comsol/primary.tex',table)
    table=[]
    for p in f.CORE:
        r=metric(full,'baseline__'+p+'__fine')
        table.append(pathlabel(p)+' & '+num(r['cmax_db'])+' & '+num(r['c99_db'])+' & '+num(r['near_support_pct'],1)+' & '+num(r['peak']['distance_mm'],3)+r' \\')
    textfile('tables/comsol/full_common.tex',table)
    table=[]
    for title,source,id in [('Full $S\\mid S$',sel,''),('Full $L\\mid L$',full,'baseline__1-1_LL__'),('Partial $S\\mid SS$',partial,'baseline__1-2_SSS__')]:
        for group,disp in [('all','All sources and receivers'),('shoulders','Shoulder sources only'),('central','Weld-surface sources only'),('plate','Side-surface receivers only'),('floor','Weld-surface receivers only'),('shoulders_plate','Shoulders and side-surface receivers')]:
            query=group if source is sel else id+('fine' if group=='all' else group)
            r=metric(source,query);table.append(title+' & '+disp+' & '+num(r['cmax_db'])+' & '+num(r['near_support_pct'],1)+r' \\')
        table.append(r'\addlinespace')
    textfile('tables/comsol/acquisition.tex',table)
    survey=rows('survey_metrics.json');protocol=json.loads((DATA/'survey_protocol.json').read_text())
    for kind in ['full','partial']:
        table=[]
        for name in [v['name'] for v in protocol['variants']]:
            for p in protocol['paths'][kind]:
                r=next(r for r in survey if r['kind']==kind and r['variant']==name and r['path']==p and r['case']=='lof')
                # Machine-readable all-case metrics are in the evidence bundle;
                # a compact per-filter best-case summary avoids 308 printed rows.
        for name in [v['name'] for v in protocol['variants']]:
            subset=[r for r in survey if r['kind']==kind and r['variant']==name and r['case']=='lof']
            good=[r for r in subset if r['cmax_db'] is not None]
            r=max(good,key=lambda r:r['cmax_db'])
            npeaks=sum(bool(v['peak'] and v['peak']['distance_mm']<=.5) for v in subset)
            table.append(f.NAMES[name].replace('%',r'\%').replace('0.5-','0.5--').replace('0-5','0--5').replace('1-5','1--5')+' & '+pathlabel(r['path'])+' & '+num(r['cmax_db'])+' & '+num(r['near_support_pct'],1)+' & '+str(npeaks)+r' \\')
        textfile('tables/comsol/'+kind+'_filters.tex',table)
    r=metric(partial,'baseline__1-2_SSS__fine');r1=metric(partial,'band1_5__1-2_SSS__fine')
    table=[]
    for name in ['baseline','band1_5','bp_short_matched_guard','lp_short_matched_guard','lp_short_initial_rest']:
        for p in ['1-2_SSS','1-1_SS','1-1_LL','3-2_SSSSS']:
            r=metric(partial,name+'__'+p+'__fine')
            table.append(f.NAMES[name]+' & '+pathlabel(p)+' & '+num(r['cmax_db'])+' & '+num(r['near_support_pct'],1)+r' \\')
    textfile('tables/comsol/partial_confirmation.tex',table)
    checks={'full_selected':metric(sel,'all'),'partial_selected':metric(partial,'baseline__1-2_SSS__fine'),
        'full_finer':metric(full,'lp_short_matched_guard__1-1_SS__finer'),
        'partial_finer':metric(partial,'baseline__1-2_SSS__finer'),
        'full_wide':metric(sel,'wide'),'partial_wide':metric(partial,'baseline__1-2_SSS__wide'),
        'full_gate15':metric(sel,'gate15'),'full_gate25':metric(sel,'gate25')}
    f.write(DATA/'reported_values.json',checks)
    print('Generated tables and values from numerical metrics')

if __name__=='__main__':main()

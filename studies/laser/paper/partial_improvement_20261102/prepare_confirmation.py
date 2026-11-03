"""Record follow-up checks after viewing the two exploratory stages."""
import copy
import json
from run import HERE, OUT, s, base

choices = {
    'baseline': base(),
    'weighted_1_5': dict(base(), processing='weighted', low_MHz=1.),
    'equal_2_5': dict(base(), low_MHz=2.),
}
jobs=[]
for name,selection in choices.items():
    for field,changes in [
        ('fine',dict(spacing_mm=.05)),
        ('wide',dict(spacing_mm=.25,x_min_mm=-45,x_max_mm=45,depth_min_mm=.25,depth_max_mm=25.25)),
    ]:
        jobs.append(dict(id=f'{name}__{field}',selection=dict(selection,**changes)))
candidate=choices['weighted_1_5']
for parity in ('even','odd'):
    jobs.append(dict(id=f'weighted_1_5__receivers_{parity}',selection=dict(candidate,aperture=parity)))
for omitted in range(5):
    jobs.append(dict(id=f'weighted_1_5__omit_source{omitted}',
                     selection=dict(candidate,sources=[i for i in range(5) if i!=omitted])))
p=dict(reason='Response-weighted 1-5 MHz retains two-source admission and 95.8% near support at 8.59 dB contrast; equal-mean 2-5 MHz is a simpler alternative at 7.99 dB. Confirm localization, field dependence, absolute amplitude and source/receiver dependence.',
       interpretation='These checks reuse the same acquisition and known flaw. They are consistency/sensitivity checks, not independent specimen validation. No stronger coherence exponent or nonlinear display is used.',
       previous_metrics_sha256={stage:s.sha(OUT/stage/'metrics.json') for stage in ('survey','frequency')},
       requests=jobs)
path=HERE/'confirmation_protocol.json'
if path.exists():assert json.loads(path.read_text())==p
else:s.write(path,p)
print('Recorded',len(jobs),'confirmation configurations',flush=True)

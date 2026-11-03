"""Independent waveform focusing through the new branch dispatch."""
from verify import *
from functools import lru_cache

class RoutedReference(PartialReference):
    @lru_cache(maxsize=None)
    def cached_ray(self,a,b,modes,outgoing,via):
        if not via:return super().ray(a,b,modes,outgoing)
        value=oracle(self,a,b,via,modes[0],outgoing)
        if value is None:return None
        d=np.diff(value['points'],axis=0)
        value.update(first=np.arctan2(d[0,0],abs(d[0,1])),last=np.arctan2(d[-1,0],abs(d[-1,1])))
        return value
    def ray(self,a,b,modes,outgoing=False):
        via=self.route['after_via' if outgoing else 'before_via']
        return self.cached_ray(tuple(a),tuple(b),modes,outgoing,via)

def main():
    catalogue=json.loads((OUT/'verification/native/path_catalogue.json').read_text())['paths']
    ids=['exp_corner_3_S_tx_id_2','exp_corner_4_S_rx_id_2','exp_corner_4_L_rx_id_2',
         'exp_corner_4_S_tx_direct_2','exp_wall_2_S_tx_pixel-id_1',
         'exp_wall_4_S_rx_pixel-id_2','exp_wall_2_L_tx_pixel-id_1','exp_wall_2_S_tx_direct_1']
    points=np.array([[-5,10],[-5,20],[0,21],[5,10],[5,20],[2.7,19.5]])/1000
    jobs=[]
    for path in ids:
        c=copy.deepcopy(BASE);c.update(path=path,cases=['lof','reference','difference'],processing='weighted',coherence=False,pulse_screening='off',min_offset_mm=2)
        jobs.append(dict(selection=c,points_m=points.tolist()))
    request=OUT/'verification/focus_request.json';write(request,dict(requests=jobs))
    run_request(request,OUT/'verification/focus')
    ref=RoutedReference(DATA);rows=[]
    for i,path in enumerate(ids):
        ref.route=next(p for p in catalogue if p['id']==path)
        before,after=ref.route['before'],ref.route['after']
        expected=ref.focus(f'{len(before)}-{len(after)}_{before+after}',points)
        errors=[]
        with np.load(OUT/'verification/focus'/str(i)/'images.npz') as a:
            for case in expected['maps']:
                np.testing.assert_allclose(a[case+'_maps'][2],expected['maps'][case],rtol=2e-6,atol=2e-13,equal_nan=True)
                np.testing.assert_allclose(a[case+'_means'][2],expected['means'][case],rtol=2e-6,atol=2e-13)
                np.testing.assert_array_equal(a[case+'_counts'][2],expected['source_count'])
                np.testing.assert_array_equal(a[case+'_receiver_usage'][2],expected['usage'])
                errors.append(float(np.max(abs(a[case+'_means'][2]-expected['means'][case]))))
        rows.append(dict(path=path,max_complex_mean_error_m_s=max(errors),supported_points=int(np.isfinite(expected['maps']['lof']).sum())))
        print(rows[-1],flush=True)
    write(OUT/'verification/focus_verification.json',dict(passed=True,points_mm=(points*1000).tolist(),paths=rows))
if __name__=='__main__':main()

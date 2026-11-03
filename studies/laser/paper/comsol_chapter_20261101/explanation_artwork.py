"""Explanatory figures using recorded data and explicit illustrative ray geometry."""
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import make_artifacts as f
from figure_layout import legend_strip


def arrow_segment(ax, a, b, color, ls='-'):
    ax.plot([a[0],b[0]],[a[1],b[1]],color=color,ls=ls,lw=1.2)
    p=np.asarray(a)*.6+np.asarray(b)*.4
    q=np.asarray(a)*.4+np.asarray(b)*.6
    ax.annotate('',xy=q,xytext=p,arrowprops=dict(arrowstyle='->',color=color,lw=1.0))


def walkthrough():
    v=np.load(f.DATA/'recordings/full.npz')
    d=f.DATASETS['full']; cs=d['material']['cs_m_s']/1000
    source=np.array([-12.,0.]); receiver=np.array([float(v['receiver_x_mm'][94]),0.])
    point=np.array([0.,20.])
    lengths=np.array([np.linalg.norm(point-source),np.linalg.norm(receiver-point)])
    time=.03+lengths.sum()/cs
    trace=v['lof_trace_um_s']; t=v['time_us']
    sample=np.interp(time,t,trace.real)+1j*np.interp(time,t,trace.imag)
    fig=plt.figure(figsize=(7.15,3.35),layout='constrained')
    gs=fig.add_gridspec(2,3,width_ratios=[1.15,1.35,.72],height_ratios=[1,.22])
    ax=fig.add_subplot(gs[0,0])
    ax.add_patch(Rectangle((-26,0),34,25.4,fc='#f4f4f4',ec='#777777',lw=.6))
    f.outline(ax,'full',marker='none')
    arrow_segment(ax,source,point,'#d55e00');arrow_segment(ax,point,receiver,'#d55e00')
    ax.scatter(*source,marker='v',c='black',s=30,zorder=5)
    ax.scatter(*receiver,marker='s',c='black',s=23,zorder=5)
    ax.scatter(*point,c='black',s=18,zorder=5)
    ax.text(source[0],-1.1,'S0',ha='center',va='bottom',fontsize=7)
    ax.text(receiver[0],-1.1,'R94',ha='center',va='bottom',fontsize=7)
    ax.annotate('Trial pixel',xy=point,xytext=(-20,22.9),fontsize=7,
                arrowprops=dict(arrowstyle='-',lw=.5))
    ax.set(xlim=(-26,8),ylim=(27,-4),xlabel='$x$ (mm)',ylabel='Depth (mm)',
           title='(a) Predict the travel time',xticks=[-20,-10,0],yticks=[0,10,20,25.4])
    ax.set_aspect('equal')
    ax=fig.add_subplot(gs[0,1])
    ax.plot(t,trace.real,c='#444444',lw=.8,label='Signed waveform')
    ax.plot(t,abs(trace),c='#8c564b',ls='--',lw=.8,label='Envelope')
    ax.axvline(time,c='black',ls=':',lw=1,label=f'Predicted time: {time:.2f} '+r'$\mu$s')
    ax.scatter(time,sample.real,c='black',s=19,zorder=5)
    ids=(t>=14)&(t<=23); bound=np.max(abs(trace[ids]))*1.15
    ax.set(xlim=(14,23),ylim=(-bound,bound),xlabel=r'Time ($\mu$s)',
           ylabel=r'Filtered $v_y$ ($\mu$m/s)',title='(b) Sample the recording')
    # All waveform keys and timing text occupy a separate row outside the data.
    legend_strip(fig,gs[1,:2],*ax.get_legend_handles_labels(),ncol=3,fontsize=7,columnspacing=1.3)
    ax=fig.add_subplot(gs[0,2])
    row=f.row('full','baseline__1-1_SS__fine')
    im=f.heat(ax,row,reference=10,title='(c) Repeat and combine',marker='none')
    ax.set_title('(c) Form the image',loc='left',pad=10)
    ax.plot(*point,marker='o',mfc='none',mec='white',ms=4,lw=0)
    cax=fig.add_subplot(gs[1,2])
    fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-40,-20,0],label='Image (dB)')
    cax.set_box_aspect(.12)
    f.save(fig,'saft_walkthrough')
    result=dict(status='passed',illustrative_pixel_mm=point.tolist(),source_mm=source.tolist(),
                receiver_mm=receiver.tolist(),segment_lengths_mm=lengths.tolist(),
                shear_speed_mm_per_us=cs,arrival_time_us=float(time),
                sampled_complex_velocity_um_s=[float(sample.real),float(sample.imag)],
                recording='S0/R94; existing full-weld LOF recording, 0.5-5 MHz',
                image_uid=row['uid'],image_uses_all_five_sources=True,
                note='The pixel is an explanatory example, not a prescribed flaw coordinate or a fixed reconstruction aperture.')
    (f.HERE/'verification/walkthrough_check.json').write_text(json.dumps(result,indent=2)+'\n')


def paths():
    fig=plt.figure(figsize=(7.15,2.75),layout='constrained')
    gs=fig.add_gridspec(2,4,height_ratios=[1,.19])
    source=np.array([-12.,0.]);receiver=np.array([12.,0.]);point=np.array([0.,20.]);height=25.4
    bounce_left=np.array([-12+(0+12)*height/(2*height-point[1]),height])
    bounce_right=np.array([12*(height-point[1])/(2*height-point[1]),height])
    definitions=[('L | L',[source,point],[point,receiver],'L'),
                 ('S | S',[source,point],[point,receiver],'S'),
                 ('S | SS',[source,point],[point,bounce_right,receiver],'S'),
                 ('SS | S',[source,bounce_left,point],[point,receiver],'S')]
    for i,(name,inc,rec,mode) in enumerate(definitions):
        ax=fig.add_subplot(gs[0,i]);ax.set_facecolor('#f6f6f6')
        ax.axhline(0,c='#777777',lw=.7);ax.axhline(height,c='#777777',lw=.7)
        color='#0072b2' if mode=='L' else '#d55e00'
        for branch in [inc,rec]:
            for a,b in zip(branch[:-1],branch[1:]):arrow_segment(ax,a,b,color)
        ax.scatter(*source,marker='v',c='black',s=27,zorder=5)
        ax.scatter(*receiver,marker='s',c='black',s=18,zorder=5)
        ax.scatter(*point,marker='o',c='black',s=15,zorder=5)
        ax.set(xlim=(-16,16),ylim=(28,-3),xlabel='$x$ (mm)',title=f'({chr(97+i)}) {name}',
               xticks=[-12,0,12],yticks=[0,10,20,25.4])
        ax.set_aspect('equal')
        if i==0:ax.set_ylabel('Depth (mm)')
        else:ax.set_yticklabels([])
    handles=[Line2D([],[],color='#0072b2',label='L segment'),Line2D([],[],color='#d55e00',label='S segment')]
    handles += [Line2D([],[],ls='',marker=m,c='black',markersize=4,label=l)
                for m,l in [('v','Source'),('s','Receiver'),('o','Trial pixel')]]
    legend_strip(fig,gs[1,:],handles,ncol=5,fontsize=7,columnspacing=1.5)
    f.save(fig,'path_guide')


if __name__=='__main__':
    f.INDEX=json.loads((f.DATA/'images.json').read_text())
    f.DATASETS={k:json.loads((f.DATA/(k+'_dataset.json')).read_text()) for k in ['full','partial']}
    walkthrough();paths()

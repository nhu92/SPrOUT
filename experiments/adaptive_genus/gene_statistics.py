"""Shared-gene resampling; no uncalibrated biological presence p-values."""
from __future__ import annotations
from dataclasses import asdict,dataclass
import numpy as np
from scipy.stats import norm


@dataclass(frozen=True)
class Settings:
    z_threshold: float=.9
    min_margin: float=.1
    min_genes: int=3
    min_comparable: int=10
    lower_threshold: float=.02
    lower_quantile: float=.1
    bootstrap_replicates: int=1000
    outer_replicates: int=100
    inner_replicates: int=200
    seed: int=20261015

    def validate(self):
        if self.min_margin!=.1:raise ValueError('Feature cache fixes margin at 0.1; re-extract with a new schema for another margin')
        if not np.isfinite(self.z_threshold):raise ValueError('Z cutoff must be finite')
        if not 0<self.lower_quantile<.5 or not 0<=self.lower_threshold<=1:raise ValueError('Invalid stability settings')
        if min(self.min_genes,self.min_comparable,self.bootstrap_replicates,self.outer_replicates,self.inner_replicates)<1:raise ValueError('Counts must be positive')


def zscores(scores):
    sd=np.std(scores,axis=-1,keepdims=True)
    return np.divide(scores-np.mean(scores,axis=-1,keepdims=True),sd,out=np.zeros_like(scores,dtype=float),where=sd>0)


def bootstrap_counts(rng,reps,n):
    indices=rng.integers(0,n,size=(reps,n))
    counts=np.zeros((reps,n),dtype=int)
    np.add.at(counts,(np.repeat(np.arange(reps),n),indices.ravel()),1)
    return counts


def wilson(success,n,confidence=.95):
    """Monte Carlo binomial interval, not biological confidence of membership."""
    if n<1:return None,None
    p=success/n;z=norm.ppf((1+confidence)/2);den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return (0. if success==0 else float(max(0,center-half)),
            1. if success==n else float(min(1,center+half)))


def states(z,k,n,raw,lower,settings):
    native=z>settings.z_threshold
    rescue=(k>=settings.min_genes)&(n>=settings.min_comparable)&(lower>=settings.lower_threshold)
    hybrid=native|rescue
    unconfirmed=hybrid&~rescue&(n>=2)&(raw<2)
    return native,rescue,hybrid,hybrid&~unconfirmed,unconfirmed


def analyze(data,settings=Settings()):
    settings.validate()
    a,w,e,k,v=[np.asarray(data[name],dtype=float) for name in ['acs','weights','eligible','informative','votes']]
    if not(a.shape==w.shape==e.shape==k.shape==v.shape) or a.ndim!=2:raise ValueError('Feature matrices must align')
    g,f=a.shape
    if g!=len(data['genes']) or f!=len(data['families']):raise ValueError('Gene/family labels do not align')
    if g<2 or f<1:raise ValueError('At least two genes and one family are required')
    if not all(np.isfinite(x).all() for x in [a,w,e,k,v]):raise ValueError('Nonfinite feature')
    if ((w<0)|(w>1)).any() or (a<0).any():raise ValueError('Invalid score bounds')
    if any(not np.isin(x,[0,1]).all() for x in [e,k,v]) or (w[e==0]!=0).any():raise ValueError('Invalid coverage or support masks')
    rng=np.random.default_rng(settings.seed)
    counts=bootstrap_counts(rng,settings.bootstrap_replicates,g)
    den=counts@e;boot=np.divide(counts@w,den,out=np.zeros_like(den),where=den>0)
    lower=np.quantile(boot,settings.lower_quantile,axis=0)
    z=zscores(a.sum(axis=0));n=e.sum(axis=0);support=k.sum(axis=0);raw=v.sum(axis=0)
    native,rescue,hybrid,selected,unconfirmed=states(z,support,n,raw,lower,settings)
    acs_boot=counts@a;acs_z_boot=zscores(acs_boot)
    # Full outer rule, including a fresh inner gene bootstrap on each outer draw.
    # Repeated draws are bootstrap observations; they do not change observed
    # counts in the exported panel. This is a conditional stability diagnostic.
    outer=bootstrap_counts(np.random.default_rng(settings.seed+1),settings.outer_replicates,g)
    pass_count=np.zeros(f,dtype=int)
    for multiplicity in outer:
        indices=np.repeat(np.arange(g),multiplicity)
        inner=bootstrap_counts(rng,settings.inner_replicates,g)
        ed=e[indices];wd=w[indices]
        iden=inner@ed
        ilower=np.quantile(np.divide(inner@wd,iden,out=np.zeros_like(iden),where=iden>0),settings.lower_quantile,axis=0)
        present=(multiplicity>0).astype(float)
        result=states(zscores(multiplicity@a),present@k,present@e,present@v,ilower,settings)
        pass_count+=result[3]
    # Delete-one influence uses score means; it does not pretend to be the
    # complete nested selection rule or an independent biological replicate.
    sums=w.sum(axis=0);means=np.divide(sums,n,out=np.zeros_like(sums),where=n>0)
    leave_den=n[None,:]-e
    leave=np.divide(sums[None,:]-w,leave_den,out=np.zeros_like(w),where=leave_den>0)
    max_shift=np.max(np.abs(leave-means),axis=0)
    share=np.divide(w.max(axis=0),sums,out=np.zeros_like(sums),where=sums>0)
    rows=[]
    for j,family in enumerate(data['families']):
        low,high=wilson(int(pass_count[j]),settings.outer_replicates)
        reason=('retained_insufficient_coverage_to_judge' if selected[j] and not rescue[j] and n[j]<2 else
                'selected' if selected[j] else 'unconfirmed_no_repeated_gene_support' if unconfirmed[j] else 'not_admitted')
        if not hybrid[j]:
            reason=('insufficient_comparable_coverage' if n[j]<settings.min_comparable else
                'too_few_discriminating_genes' if support[j]<settings.min_genes else 'weak_or_unstable_margin')
        rows.append(dict(family=family,parent_order=data['family_parent'].get(family,''),
            acs=float(a[:,j].sum()),z=float(z[j]),native_selected=bool(native[j]),rescue_selected=bool(rescue[j]),
            hybrid_selected=bool(hybrid[j]),selected=bool(selected[j]),unconfirmed=bool(unconfirmed[j]),reason=reason,
            admission_route='both' if native[j] and rescue[j] else 'acs' if native[j] else 'gene' if rescue[j] else 'neither',
            raw_gene_votes=int(raw[j]),discriminating_genes=int(support[j]),comparable_genes=int(n[j]),
            support_fraction=float(support[j]/n[j]) if n[j] else None,
            reference_species=data['reference_species'].get(family,0),
            missing_sibling_genes=data['missing_sibling_genes'].get(family,0),
            mean_margin=float(means[j]) if n[j] else None,bootstrap_lower=float(lower[j]),
            margin_interval=[float(x) for x in np.quantile(boot[:,j],[.025,.975])],
            acs_interval=[float(x) for x in np.quantile(acs_boot[:,j],[.025,.975])],
            z_interval=[float(x) for x in np.quantile(acs_z_boot[:,j],[.025,.975])],
            acs_bootstrap_frequency=float(np.mean(acs_z_boot[:,j]>settings.z_threshold)),
            selection_frequency=float(pass_count[j]/settings.outer_replicates),
            selection_frequency_mc_interval=[low,high],max_gene_share=float(share[j]),
            max_delete_one_mean_shift=float(max_shift[j]),presence_p=None,adjusted_p=None,
            p_status='not_calibrated_no_valid_absence_null'))
    rows.sort(key=lambda r:(-r['selection_frequency'],-r['bootstrap_lower'],-r['z'],r['family']))
    return dict(schema=1,project=data['project'],settings=asdict(settings),rows=rows,genes=g,
        truth_used=False,input_sha256=data['input_sha256'],
        z_origin='User reports historical NPV approximately 90% at Z > 0.9; not a probability calibration or universal guarantee',
        stability_scope='Conditional on observed recovered genes and reference panel; nested C rule with fixed settings; bootstrap multiplicities weight scores, while support/coverage gates count distinct observed gene identities',
        interval_scope='Margin/ACS/Z: gene percentile bootstrap; frequency interval: Monte Carlo Wilson precision only',
        p_scope='Biological absence null unavailable; no biological presence p-values or FDR guarantees are inferred')

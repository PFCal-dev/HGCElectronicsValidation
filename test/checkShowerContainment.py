import ROOT
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use([hep.style.CMS, hep.style.firamath])

url='/eos/cms/store/cmst3/group/hgcal/CMG_studies/psilva/DigiTester/SingleK0LEGun_eta2p5_13_0_0_pre4_D99.root'
if len(sys.argv)>1:
    url=sys.argv[1]

rdf=ROOT.RDataFrame('ana/hits',url)
rdf=rdf.Filter('crossCalo && gdradius<100 && isSci==0 && qsim>12')
df=pd.DataFrame(rdf.AsNumpy(columns=['gdradius','z','gpt','layer','qsim','cce']))

fig,ax=plt.subplots(1,2,figsize=(12,6))
for reg,mask in [
        ('CE-E',(np.abs(df['z'])<362.18)),
        ('CE-Hfine',(np.abs(df['z'])>362.18) & (np.abs(df['z'])<437.87)),
        ('CE-Hcoarse',(np.abs(df['z'])>437.87))
]:
    if mask.sum()<100: continue
    h,e=np.histogram(np.clip(df[mask]['gdradius'],0,20),np.linspace(0,20,100),weights=df[mask]['qsim']*df[mask]['cce'])
    h=h/h.sum()
    hep.histplot(h,e,label=reg,ax=ax[0])
    hep.histplot(np.cumsum(h),e,ax=ax[1])
ax[0].legend()
for i in range(2):
    ax[i].grid()
    ax[i].set_xlabel('$\Delta R$ [cm]')
ax[0].set_ylabel('PDF')
ax[1].set_ylabel('CDF')
ax[1].plot([0,20],[0.68,0.68],ls='--',color='gray')
fig.tight_layout()
plt.savefig(os.path.basename(url).replace('.root','.png'))
plt.show()

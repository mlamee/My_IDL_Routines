; Monte Carlo Simulation of Intergalactic Medium absorption along several lines of sight. we follow Inoue and Iwata (2008), with the updated parameters from Inoue et al. (2011). 
; We use the provided probability density functions and generate IGM absorbers along each sightline with variety of neutral hydrogen column density NHi, doppler broadening parameter b and redshift z. 
; The results are saved in separate files. The routine "tra23.pro" reads these files and create absorption spectra for each sightline. 
;Another routine "av3.35.pro" takes into account all the simulated sightlines and calculates the mean and the median IGM absorption spectrum of a source at z=zs. 

; z is the redshift of the IGM absorbers
Function  f,z
   z1=1.2
   z2=4.0
   y1=0.2
   y2=2.5
   y3=4.0
   if z ge 0 and z le z1 then return,500.*((1+z)/(1+z1))^y1
   if z gt z1 and z le z2 then return,500.*((1+z)/(1+z1))^y2
   if z gt z2 then return, 500.*((1+z2)/(1+z1))^y2+((1+z)/(1+z2))^y3
end

; b is the doppler broadening parameter
function h,b
   b0=23. ;km/s
   return, 4.*b0^4*exp(-(b0^4/b^4))/b^5.
end

;NHi is the column density of neutral hydrogen
function g,NHi
   B1=1.7
   B2=1.2
   Nc=1.6e17 ;1/cm^2
   N1=1.e12 ;1/cm^2
   Nu=1.e22 ;1/cm^2
   if NHi ge N1 and NHi lt Nc then return, (NHi/Nc)^(-B1)
   if NHi ge Nc and NHi le Nu then return, (NHi/Nc)^(-B2)
end
;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;
pro igm
;spawn,'date'
; Number of simulated sightlines
num=5000
z1=1.2
z2=4.0
y1=0.2
y2=2.5
y3=4.0

count=0L
; Redshift of the emitting source
zs=3.33
bd=lindgen(1200)/10.+0.1
h=bd*0.
for i=0, n_elements(bd)-1 do h[i]=h(bd[i])
hh=h*0.
ff=(1.e14-1.e12)/20000.
nh11=lindgen(20000L)*ff+1.e12
bb=(1.6e17-1.e14)/100000. 
nh12=lindgen(100001L)*bb+1.e14  
remove,[0,100000],nh12
nh1=[nh11,nh12] 
cc=(1.e22-1.6e17)/100000.
nh2=lindgen(100001L)*cc+1.6e17
nh=[nh1,nh2] 
g=nh*0.
for i=0,n_elements(nh)-1 do g[i]=g(nh[i])
gg=lindgen(n_elements(g)-1)*0.
norm=tsum(nh,g)
for i=1, n_elements(h)-1 do hh[i]=int_tabulated(bd[0:i],h[0:i],/double) 
for i=1,n_elements(nh1)-1 do gg[i]=tsum(nh1[0:i],g[0:i])/norm
for i=n_elements(nh1)+1,n_elements(nh)-1 do gg[i-1]=gg[n_elements(nh1)-1]+tsum(nh[n_elements(nh1):i],g[n_elements(nh1):i])/norm
spawn,'date'

; Generating "num" absorption spectra of a source at redshift zs for different sightlines and saving each spectrum
for w=0,num-1 do begin
   name=strcompress('~/IGM/z3.33/los3.33.'+string(w+1)+'.dat',/remove_all)
   openw,lun,name,/get_lun
   z=0.0
   za=0.
   nhi=0.
   bdop=0.
   while z le zs and z ge 0 do begin
      f1=f(z)
      f1=f1[0]
      dz=lindgen(zs*10000L+1-fix(10000.*z))/10000.
      phi=fltarr(n_elements(dz))
      p1=f1*exp(-f1*dz)   
      a=randomu(s); for phi
      a2=randomu(s) ;for hh
      a3=randomu(s) ;for nh
      ii=0
      while phi[ii] le a and ii+1 lt n_elements(dz) do begin
         ii=ii+1 
         phi[ii]=int_tabulated(z+dz[0:ii],p1[0:ii],/double)
      endwhile

      if a lt phi[ii]  then begin      
         j=ii-1
         b=(phi[j+1]-phi[j])/(dz[j+1]-dz[j]) ; Linear interpolation
         c=phi[j]-b*(z+dz[j])
         z=(a-c)/b
      
         k2=where(hh lt a2)
         j2=k2[n_elements(k2)-1]
         if j2 eq n_elements(hh)-1 then j2=j2-1
         b2=(hh[j2+1]-hh[j2])/(bd[j2+1]-bd[j2])
         c2=hh[j2]-b2*bd[j2]
         b_dop=(a2-c2)/b2
      
         k3=where(gg lt a3)
         j3=k3[n_elements(k3)-1]
         if j3 eq n_elements(gg)-1 then j3=j3-1
         b3=(gg[j3+1]-gg[j3])/(nh[j3+1]-nh[j3])
         c3=gg[j3]-b3*nh[j3]
         n_hi=(a3-c3)/b3
   
         printf,lun,z,b_dop,n_hi
         if z lt zs then za=[za,z]           ;redshift of the first absorber 
         bdop=[bdop,b_dop]
         nhi=[nhi,n_hi]
      endif else begin
         print,'No absorber'
         print, 'a',a, '   phi',phi[n_elements(phi)-1]
         z=10*zs
      endelse
   endwhile
free_lun,lun
endfor

spawn,'date'
;stop
end

;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;
;This reads the files created by inoue.pro and creates absorption spectra for each sightline by taking into account absorption in the continuum (below 912 ̊A)
;and in the Lyman series.
; Atomic parameters fi (oscilator strengh of the ith transition)and Aij (Einestein coefficients) are used to calculate the Lyman series cross sections "sigi" assuming a Voigt line profile.
pro tra23
spawn,'date'

; Lyman continuum cross section
sigll=6.30e-18 ;1/cm^2
h0=6.63e-27 ;erg/s
c0=3.e10 ;cm/s
; frequency of photons with 13.6 ev energy
vll=c0/(912.324e-8)
; mass of electron
me=9.11e-28 ;g
; electron charge
e=4.80e-10
;redshift of the source
zs=3.33
zlim=700.*(1+zs)/1300.-1.
; Lyman series wavelength in cm 
lami=1.e-8*[1215.67, 1025.72, 992.537, 949.743, 937.803, 930.748, 926.226, 923.150, 920.963, 919.352, 918.129, 917.181, 916.429, 915.824, 915.329, 914.919, 914.576, 914.286, 914.039, 913.826, 913.641, 913.480, 913.339, 913.215, 913.104, 913.006, 912.918, 912.839, 912.768, 912.703, 912.645, 912.592, 912.543, 912.499, 912.458, 912.420, 912.385, 912.353, 912.324]
; oscilator strengths
fi=[0.4162, 7.910e-2, 2.899e-2, 1.394e-2, 7.799e-3, 4.814e-3, 3.183e-3, 2.216e-3, 1.605e-3, 1.201e-3, 9.214e-4, 7.227e-4, 5.774e-4, 4.686e-4, 3.856e-4, 3.211e-4, 2.702e-4, 2.296e-4, 1.967e-4, 1.698e-4, 1.476e-4, 1.291e-4, 1.136e-4, 1.005e-4, 8.928e-5, 7.970e-5, 7.144e-5, 6.429e-5, 5.806e-5, 5.261e-5, 4.782e-5, 4.360e-5, 3.986e-5, 3.653e-5, 3.357e-5, 3.092e-5, 2.854e-5, 2.640e-5, 2.446e-5]
;fi=double(fi)
; Einestein coefficients
Aij=[4.699e8, 5.575e7, 1.278e7, 4.125e6, 1.644e6, 7.568e5, 3.869e5, 2.143e5, 1.263e5, 7.834e4, 5.066e4, 3.393e4, 2.341e4, 1.657e4, 1.200e4, 8858.,6654., 5077., 3928., 3077., 2438., 1952., 1578., 1286., 1057., 875.3, 729.7, 612.2, 516.7, 438.6, 374.2, 320.8, 276.3, 239., 207., 181., 158.4, 139.1, 122.6]
; frequency of Lyman series
vi=c0/lami
; Generating  a continuum spectrum with fixed flux density of 1  
lamb=1e-8*(lindgen(6000L)/10.+700.)
vs=c0/lamb
cc0=16*c0*!pi^2
pii=sqrt(!pi)
set_plot,'ps'

; w represents the sightlines. we can run it for everything or for a small patch of sightlines we have produced.
for w=901,1000 do begin 
   name=strcompress('~/IGM/z3.33/los3.33.'+string(w+1)+'.dat',/remove_all)
   name2=strcompress('~/IGM/z3.33/Trans3.33.'+string(w+1)+'.dat',/remove_all)
   name3=strcompress('Trans3.33.'+string(w+1)+'.ps',/remove_all)
   readcol,name,z,b,nhi
   ind=where(z lt zlim)
   remove,ind,z
   remove,ind,b
   remove,ind,nhi
   b=b*1.e5 ;cm/s
   vd=fltarr(n_elements(vi),n_elements(b))
   for i=0,n_elements(vi)-1 do for j=0,n_elements(b)-1 do vd[i,j]=vi[i]*b[j]/c0  

   vz=fltarr(n_elements(z),n_elements(vs))
   siglc=vz*0.
   for m=0, n_elements(z)-1 do begin
      vz[m,*]=vs*(1+z[m])/(1.+zs)  ;redshifting the frequency
      q=where(vz[m,*] ge vll)
      if q[0] ne -1 then siglc[m,q]=sigll*(vll/vz[m,q])^3.
   endfor

   phi=fltarr(n_elements(vs),n_elements(vi),n_elements(b))
   ; The Voigt profile
   for i=0,n_elements(vi)-1 do for k=0, n_elements(b)-1 do phi[*,i,k]=voigt(total(Aij[0:i])*(lami[i]^2/(cc0*vd[i,k])) ,(vz[k,*]-vi[i])/vd[i,k])/(pii*vd[i,k])
   sigi=fltarr(n_elements(b),n_elements(vs))
   for k=0, n_elements(b)-1 do for j=0,n_elements(vs)-1 do for i=0,n_elements(vi)-1 do sigi[k,j]=sigi[k,j]+fi[i]*phi[j,i,k]
   ; The Lyman series cross section
   sigi=(!pi*e^2)/(me*c0)*sigi
   openw,lun,name2,/get_lun

   ; Calculating the effective opacity at different wavelength for each sightline
   Teff=fltarr(n_elements(vs))
   for j=0, n_elements(vs)-1 do begin
      for k=0, n_elements(b)-1 do Teff[j]=Teff[j]+nhi[k]*(siglc[k,j]+sigi[k,j])
      printf,lun,lamb[j]*1.e8, Teff[j]
   endfor

   device,file=name3,/color,/encapsul,ysize=10,xsize=20
   plot,lamb*1.e8,exp(-Teff),yr=[0,1.1],xr=[700,1300],xtitle='Wavelength (' +cgsymbol('Angstrom') + ')',ytitle='IGM Transmission',title='IGM Transmission for zs=3.33 rest frame',charsize=1.2,color=get_colour_by_name('black'), thick=1
   device,/close
   free_lun,lun
endfor
set_plot,'x' 
spawn,'date'
;stop
end
;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;
; Reads the output of "tra23.pro" and calculate the average mean, median and 16-84 percentiles of the IGM absorption for a source at redshift z=zs 
pro av2

tav=fltarr(6000)
tam=fltarr(6000)
tap=fltarr(2,6000)

readcol,'~/IGM/z3.35/tran3.33.dat',name,format='a'
nn=n_elements(name)
t=fltarr(nn,6000)
openw,lun,'IGM-z3.33.dat',/get_lun
for w=0,nn-1 do begin
   readcol,name[w],lam,tt
   t[w,*]=tt
endfor

for i=0,5999 do begin 
   tav[i]=mean(exp(-t[0:nn-1,i]))
   tam[i]=median(exp(-t[0:nn-1,i]))
   tap[*,i]=percentiles(reform(exp(-t[0:nn-1,i])),value=[0.16,0.84])  
   printf,lun,lam[i],tav[i], tam[i]
endfor
free_lun,lun

lam1=[indgen(700),lam]
tav1=[fltarr(700),tav]

set_plot,'ps'
device,file='IGM-z3.33.ps',/color,/encapsul,ysize=10,xsize=20
plot,lam,tav,yr=[0,1.1],xtitle='Wavelength (' +cgsymbol('Angstrom') + ')',ytitle='IGM Transmission',title='IGM Transmission for zs=3.33 rest frame',charsize=1.2,color=get_colour_by_name('black'), thick=1,xr=[600,1100],xstyle=1,charthick=3,xthick=3,ythick=3;,/nodata
oplot,lam,tam,thick=1,color=get_colour_by_name('blue')
for j=0,5900,100 do oplot,[lam[j],lam[j]],[tap[0,j],tap[1,j]];
device,/close

set_plot,'x'
;stop 
end

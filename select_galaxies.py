"""
Creating a mock catalogue from the SHARK runs on pawsey.
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM

from read import read_lightcone, read_photometry_data_hdf5, read_filter_names
from write import CatalogueDetails, stamp_preamble


cosmo = FlatLambdaCDM(Om0=0.3, H0=67.51) # Change this to freeze in our cosmology.
cat_details = CatalogueDetails(area = 107.889, mag_filter='', mag_cut=23, redshift_cut=3, version='b0.0.1')


galaxy_fields = {'galaxies': ('dec', 'ra', 'zobs',
                           'id_galaxy_sky','sfr_burst','sfr_disk','mstars_bulge','mstars_disk','rstar_bulge_apparent',
                           'rstar_disk_apparent','id_group_sky','dc', 'mvir_hosthalo', 'type')}

sed_fields = {'SED/ap_dust': ('total', 'bulge_t')}
lightcone_dir = '/scratch/pawsey0119/clagos/Stingray/output/medi-SURFS/Shark-TreeFixed-ReincPSO-kappa0p002/deep-optical-final/'
subdir = 'split/'
sed_file = "Sting-SED-VST-eagle-rr14"
sub_volumes = np.arange(2)

galaxy_data = read_lightcone(lightcone_dir, subdir, galaxy_fields, sub_volumes, 'mock')
sed_ids, sed_data = read_photometry_data_hdf5(lightcone_dir, subdir, sed_fields, sub_volumes, sed_file)





def prepare_data(phot_data, ids_sed, hdf5_data, hdf5_data_mvir, subvols, lightcone_dir,  nbands):

    (dec, ra, zobs, idgal, sfrb, sfrd, mstarb, mstard, rsb, rsd, id_group, dc, mvir_host, typeg) = hdf5_data
    (mvirz0, idgal, snap, subv) = hdf5_data_mvir

    #components of apparent magnitudes:
    #(len(my_data), 2, 2, 5, nbands)
    #0: disk instability bulge
    #1: galaxy merger bulge
    #2: total bulge
    #3: disk
    #4: total
    SEDs_dust   = phot_data[0]
    SEDs_dust_bulge = phot_data[1]

    mstartot = np.log10((mstarb+mstard)/cosmo.h)
    sfrtot = np.log10((sfrb+sfrd)/cosmo.h/1e9)
    re = (rsb*mstarb + mstard*rsd) / (mstarb+mstard)
    BT = mstarb / (mstarb+mstard)
    mvir = np.log10(mvir_host/cosmo.h)
    mvirz0 = np.log10(mvirz0/cosmo.h)

    fluxSEDs = 10.0**(SEDs_dust/(-2.5)) * 3631.0 * 1e3 #in mJy

    dgal = 4.0 * np.pi * pow((1.0+zobs) * dc/cosmo.h, 2.0)

    #S850 microns selected to have a flux >0.05mJy
    ind = np.where((SEDs_dust[6,:] > 0) & (SEDs_dust[6,:] <= 23) & (zobs <= 3))
    SEDs_dustin = SEDs_dust[:,ind]
    SEDs_dustin = SEDs_dustin[:,0,:]
   

    cat_details = CatalogueDetails(area = 107.889, mag_filter='app_Z_VISTA', mag_cut=23, redshift_cut=3, version='0.0.1')

    writeon = True
    if(writeon == True):
       with open('/scratch/pawsey0119/clagos/Stingray/output/medi-SURFS/Shark-TreeFixed-ReincPSO-kappa0p002/deep-optical/Shark-deep-opticalLightcone-WAVES.txt', 'w', encoding='utf-8') as file:
            stamp_preamble(file)
            cat_details.stamp_details(file)

            file.write("#units\n")
            file.write("#mstar and mvir in [Msun]\n")
            file.write("#sfr[Msun/yr]\n")
            file.write("#re[arcsec]\n")
            file.write("#magnitudes AB\n")
            file.write("#type_galaxy: =0 for central, and >0 for satellites\n")
            file.write("#id_group_sky = -1 if a galaxy is the only one in its host halo\n")
            file.write("#\n")
            file.write("#dec ra redshift log10(mstar) log10(sfr) re B/T app_u_VST app_g_VST app_r_VST app_i_VST app_Z_VISTA app_Y_VISTA app_J_VISTA app_H_VISTA app_K_VISTA id_group_sky log10(mvir) log10(mvir_z0) type_galaxy\n")
            for a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,r,s,t,u in zip(dec[ind], ra[ind], zobs[ind], mstartot[ind], sfrtot[ind], re[ind], BT[ind], SEDs_dustin[2,:], SEDs_dustin[3,:], SEDs_dustin[4,:], SEDs_dustin[5,:], SEDs_dustin[6,:], SEDs_dustin[7,:], SEDs_dustin[8,:], SEDs_dustin[9,:], SEDs_dustin[10,:], id_group[ind], mvir[ind], mvirz0[ind], typeg[ind]):
                file.write("%5.10f %5.10f %5.10f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f  %5.2f %5.2f %5.2f %5.2f %5.2f %10.0f %5.2f %5.2f %5.2f\n" % (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,r,s,t,u))

def main():

    lightcone_dir = '/scratch/pawsey0119/clagos/Stingray/output/medi-SURFS/Shark-TreeFixed-ReincPSO-kappa0p002/deep-optical-final/'
    outdir= '/sratch/pawsey0119/clagos/Stingray/output/medi-SURFS/Shark-TreeFixed-ReincPSO-kappa0p002/deep-optical-final/Plots/'
    obsdir= '/software/projects/pawsey0119/clagos/shark/data/'


    subvols = range(64)
    sed_file = "Sting-SED-VST-eagle-rr14"
    filter_names = read_filter_names(lightcone_dir, 'split/', sed_file)
    fields_sed = {'SED/ap_dust': ('total', 'bulge_t')}

    ids_sed, seds = read_photometry_data_hdf5(lightcone_dir, 'split/', fields_sed, subvols, sed_file)

    fields = {'galaxies': ('dec', 'ra', 'zobs',
                           'id_galaxy_sky','sfr_burst','sfr_disk','mstars_bulge','mstars_disk','rstar_bulge_apparent',
                           'rstar_disk_apparent','id_group_sky','dc', 'mvir_hosthalo', 'type')}
    fields_mvir = {'galaxies': ('mvir_z0','id_galaxy_sam','snapshot','subvolume')}

    hdf5_data = read_lightcone(lightcone_dir, 'split/', fields, subvols, "mock")

    hdf5_data_mvir = read_lightcone(lightcone_dir, 'split/', fields_mvir, subvols, "final_mvir")

    nbands = len(seds[0])
    prepare_data(seds, ids_sed, hdf5_data, hdf5_data_mvir, subvols, lightcone_dir, nbands)


if __name__ == '__main__':
    main()

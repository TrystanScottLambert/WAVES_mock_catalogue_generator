# Properties 
## Galaxy Properties 
| Property | Description | 
| --- | --- | 
| `coline_flux_int_1` | Integrated CO(1-0) line flux [W/m^2] | 
| `coline_flux_int_10` | Integrated CO(10-9) line flux [W/m^2] | 
| `coline_flux_int_2` | Integrated CO(2-1) line flux [W/m^2] | 
| `coline_flux_int_3` | Integrated CO(3-2) line flux [W/m^2] | 
| `coline_flux_int_4` | Integrated CO(4-3) line flux [W/m^2] | 
| `coline_flux_int_5` | Integrated CO(5-4) line flux [W/m^2] | 
| `coline_flux_int_6` | Integrated CO(6-5) line flux [W/m^2] | 
| `coline_flux_int_7` | Integrated CO(7-6) line flux [W/m^2] | 
| `coline_flux_int_8` | Integrated CO(8-7) line flux [W/m^2] | 
| `coline_flux_int_9` | Integrated CO(9-8) line flux [W/m^2] | 
| `coline_flux_int_vel_1` | Velocity-integrated CO(1-0) line flux [Jy km/s] | 
| `coline_flux_int_vel_10` | Velocity-integrated CO(10-9) line flux [Jy km/s] | 
| `coline_flux_int_vel_2` | Velocity-integrated CO(2-1) line flux [Jy km/s] | 
| `coline_flux_int_vel_3` | Velocity-integrated CO(3-2) line flux [Jy km/s] | 
| `coline_flux_int_vel_4` | Velocity-integrated CO(4-3) line flux [Jy km/s] | 
| `coline_flux_int_vel_5` | Velocity-integrated CO(5-4) line flux [Jy km/s] | 
| `coline_flux_int_vel_6` | Velocity-integrated CO(6-5) line flux [Jy km/s] | 
| `coline_flux_int_vel_7` | Velocity-integrated CO(7-6) line flux [Jy km/s] | 
| `coline_flux_int_vel_8` | Velocity-integrated CO(8-7) line flux [Jy km/s] | 
| `coline_flux_int_vel_9` | Velocity-integrated CO(9-8) line flux [Jy km/s] | 
| `hiline_flux_peak` | [s/km] normalised peak HI line flux density of inclined galaxy (multiply by hiline_flux_int_vel to get Jy values) | 
| `hiline_flux_central` | [s/km] normalised central HI line flux density of inclined galaxy (multiply by hiline_flux_int_vel to get Jy values) | 
| `hiline_width_peak` | [km/s] HI line-width between flux peaks of inclined galaxy in rest-frame velocity units | 
| `hiline_width_50` | [km/s] HI line-width at 50% of the peak flux of inclined galaxy in rest-frame velocity units | 
| `hiline_width_20` | [km/s] HI line-width at 20% of the peak flux of inclined galaxy in rest-frame velocity units | 
| `hiline_flux_peak_eo` | [s/km] normalised peak HI line flux density of edge-on galaxy (multiply by hiline_flux_int_vel to get Jy values) | 
| `hiline_flux_central_eo` | [s/km] normalised central HI line flux density of edge-on galaxy (multiply by hiline_flux_int_vel to get Jy values) | 
| `hiline_width_peak_eo` | [km/s] HI line-width between flux peaks of edge-on galaxy in rest-frame velocity units | 
| `hiline_width_50_eo` | [km/s] HI line-width at 50% of the peak flux of edge-on galaxy in rest-frame velocity units | 
| `hiline_width_20_eo` | [km/s] HI line-width at 20% of the peak flux of edge-on galaxy in rest-frame velocity units | 
| `h2line_flux_peak` | [s/km] normalised peak molecular line flux density of inclined galaxy (multiply by velocity-integrated flux to get Jy) | 
| `h2line_flux_central` | [s/km] normalised central molecular line flux density of inclined galaxy (multiply by velocity-integrated flux to get Jy) | 
| `h2line_width_peak` | [km/s] molecular line-width between flux peaks of inclined galaxy in rest-frame velocity units | 
| `h2line_width_50` | [km/s] molecular line-width at 50% of the peak flux of inclined galaxy in rest-frame velocity units | 
| `h2line_width_20` | [km/s] molecular line-width at 20% of the peak flux of inclined galaxy in rest-frame velocity units | 
| `h2line_flux_peak_eo` | [s/km] normalised peak molecular line flux density of edge-on galaxy (multiply by velocity-integrated flux to get Jy) | 
| `h2line_flux_central_eo` | [s/km] normalised central molecular line flux density of edge-on galaxy (multiply by velocity-integrated flux to get Jy) | 
| `h2line_width_peak_eo` | [km/s] molecular line-width between flux peaks of edge-on galaxy in rest-frame velocity units | 
| `h2line_width_50_eo` | [km/s] molecular line-width at 50% of the peak flux of edge-on galaxy in rest-frame velocity units | 
| `h2line_width_20_eo` | [km/s] molecular line-width at 20% of the peak flux of edge-on galaxy in rest-frame velocity units | 
| `snapshot` | snapshot index | 
| `subvolume` | subvolume index | 
| `tile` | tile index in tiling array | 
| `zobs` | redshift in observer-frame | 
| `zcmb` | redshift in CMB frame | 
| `zcos` | cosmological redshift without peculiar motions | 
| `dc` | [Mpc/h] comoving distance | 
| `ra` | [deg] right ascension | 
| `dec` | [deg] declination | 
| `id_galaxy_sky` | unique galaxy ID in mock sky, from 1 to n | 
| `id_galaxy_sky_smart` | unique galaxy ID in mock sky, smart format | 
| `id_group_sky` | unique group ID if galaxy is in a group, -1 otherwise | 
| `id_galaxy_sam` | galaxy ID in SAM | 
| `id_halo_sam` | host halo ID in SAM | 
| `type` | galaxy type (0=central, 1=satellite in halo, 2=orphan) | 
| `inclination` | [deg] inclination = angle between line-of-sight and spin axis | 
| `pa` | [deg] position angle from north to east | 
| `mag` | apparent magnitude (generic: M/L ratio of 1, no k-correction) | 
| `vpec_x` | [proper km/s] x-component of peculiar velocity | 
| `vpec_y` | [proper km/s] y-component of peculiar velocity | 
| `vpec_z` | [proper km/s] z-component of peculiar velocity | 
| `vpec_r` | [proper km/s] line-of-sight peculiar velocity | 
| `mstars_disk` | [Msun/h] stellar mass in the disk | 
| `mstars_bulge` | [Msun/h] stellar mass in the bulge | 
| `mgas_disk` | [Msun/h] gas mass in the disk | 
| `mgas_bulge` | [Msun/h] gas mass in the bulge | 
| `matom_disk` | [Msun/h] atomic gas mass in the disk | 
| `matom_bulge` | [Msun/h] atomic mass in the bulge | 
| `mmol_disk` | [Msun/h] molecular gas mass in the disk | 
| `mmol_bulge` | [Msun/h] molecular gas mass in the bulge | 
| `l_x` | [Msun pMpc km/s] x-component of total angular momentum | 
| `l_y` | [Msun pMpc km/s] y-component of total angular momentum | 
| `l_z` | [Msun pMpc km/s] z-component of total angular momentum | 
| `jdisk` | [cMpc/h km/s] specific angular momentum of the disk | 
| `jbulge` | [cMpc/h km/s] specific angular momentum of the bulge | 
| `mvir_hosthalo` | [Msun/h] host halo mass | 
| `mvir_subhalo` | [Msun/h] subhalo mass | 
| `zgas_disk` | metallicity of the gas in the disk | 
| `zgas_bulge` | metallicity of the gas in the bulge | 
| `sfr_disk` | [Msun/h/Gyr] star formation rate in the disk | 
| `sfr_burst` | [Msun/h/Gyr] star formation rate in the bulge | 
| `mbh` | [Msun/h] black hole mass | 
| `mbh_acc_hh` | [Msun/h/Gyr] black hole accretion rate in the hot halo mode | 
| `mbh_acc_sb` | [Msun/h/Gyr] black hole accretion rate in the starburst mode | 
| `vvir_hosthalo` | [km/s] host halo virial velocity | 
| `vvir_subhalo` | [km/s] subhalo virial velocity | 
| `cnfw_subhalo` | NFW concentration parameter of subhalo | 
| `rstar_disk_apparent` | [arcsec] apparent semi-major axis of half-mass ellipse of stellar disk | 
| `rstar_bulge_apparent` | [arcsec] apparent semi-major axis of half-mass ellipse of stellar bulge | 
| `rgas_disk_apparent` | [arcsec] apparent semi-major axis of half-mass ellipse of gas disk | 
| `rgas_bulge_apparent` | [arcsec] apparent semi-major axis of half-mass ellipse of gas bulge | 
| `rstar_disk_intrinsic` | [cMpc/h] intrinsic half-mass radius of stellar disk | 
| `rstar_bulge_intrinsic` | [cMpc/h] intrinsic half-mass radius of stellar bulge | 
| `rgas_disk_intrinsic` | [cMpc/h] intrinsic half-mass radius of gas disk | 
| `rgas_bulge_intrinsic` | [cMpc/h] intrinsic half-mass radius of gas bulge | 
| `hiline_flux_int` | [W/m^2] integrated HI line flux | 
| `hiline_flux_int_vel` | [Jy km/s] velocity-integrated HI line flux | 
## Group Properties 
| Property | Description | 
| --- | --- | 
| `snapshot` | snapshot index | 
| `subvolume` | subvolume index | 
| `tile` | tile index in tiling array | 
| `zobs` | redshift in observer-frame | 
| `zcmb` | redshift in CMB frame | 
| `zcos` | cosmological redshift without peculiar motions | 
| `dc` | [Mpc/h] comoving distance | 
| `ra` | [deg] right ascension | 
| `dec` | [deg] declination | 
| `vpec_x` | [proper km/s] x-component of peculiar velocity | 
| `vpec_y` | [proper km/s] y-component of peculiar velocity | 
| `vpec_z` | [proper km/s] z-component of peculiar velocity | 
| `vpec_r` | [proper km/s] line-of-sight peculiar velocity | 
| `sigma_3d_all` | [proper km/s] 3D peculiar velocity dispersion of ALL group members, including non-detections | 
| `sigma_los_detected` | [proper km/s] line-of-sight peculiar velocity dispersion of detected group members | 
| `id_group_sky` | unique parent halo ID in mock sky | 
| `id_halo_sam` | parent halo ID in SAM | 
| `mvir` | [Msun/h] virial mass | 
| `n_galaxies_total` | total number of galaxies that live in the same group (host halo) | 
| `n_galaxies_selected` | number of galaxies that live in the same group (host halo) and are present in the mock sky | 
| `flag` | group flag: 0 if group complete, >0 if truncated by survey edge (+1), snapshot limit (+2), tile edge (+4), shell edge (+8) | 

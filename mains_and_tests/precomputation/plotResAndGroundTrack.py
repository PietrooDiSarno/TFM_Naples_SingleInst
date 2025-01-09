from pySPICElib.kernelFetch import kernelFetch
from pySPICElib.SPICEtools import *

import PSOA
import numpy as np
import spiceypy as spice
from FuturePackage import Instrument
from FuturePackage import ROIDataBase


target_body = "GANYMEDE"  # Can be a list of strings or a single string

if target_body == "GANYMEDE":
    METAKR = ['https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_crema_5_1_150lb_23_1_default_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_crema_5_1_150lb_23_1_comms_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_crema_5_1_150lb_23_1_conjctn_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_crema_5_1_150lb_23_1_flybys_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_crema_5_1_150lb_23_1_baseline_v03.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_ptr_soc_pcw2_s01p00_v02.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_ptr_soc_pcw2_s02p00_v02.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_ptr_soc_lega_s07p00_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_ptr_soc_s007_01_s06p00_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_attc_000060_230414_240531_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_attm_000059_240817_240827_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_lpbooms_f160326_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_magboom_f160326_v04.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_majis_scan_zero_v02.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_swi_scan_zero_v02.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sa_crema_5_1_150lb_23_1_default_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sa_crema_5_1_150lb_23_1_baseline_v04.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sa_ptr_soc_s007_01_s02p00_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_mga_crema_5_1_150lb_23_1_default_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_mga_crema_5_1_150lb_23_1_baseline_v04.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_mga_ptr_soc_s007_01_s02p00_v01.bc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_v40.tf',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_sci_v17.tf',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_ops_v10.tf',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_dsk_surfaces_v11.tf',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_roi_v02.tf',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_events_crema_5_1_150lb_23_1_v02.tf',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_stations_topo_v01.tf',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/rssd0002.tf',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/earth_topo_050714.tf',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/earthfixediau.tf',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/estrack_v04.tf',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_europa_plasma_torus_v03.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_io_plasma_torus_v05.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_jup_ama_gos_ring_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_jup_halo_ring_v04.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_jup_main_ring_v04.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_jup_the_ring_ext_v01.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_jup_the_gos_ring_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_fixed_v01.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_bus_v07.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_gala_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_janus_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_jmc1_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_jmc2_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_navcam1_v01.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_navcam2_v01.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_lpb1_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_lpb2_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_lpb3_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_lpb4_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_rwi_v04.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_scm_v03.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_mag_v06.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_majis_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_mga_apm_v04.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_mga_dish_v04.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_pep_jdc_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_pep_jei_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_pep_jeni_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_pep_jna_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_pep_nim_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_rimemx_v03.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_rimepx_v03.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_sapy_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_samy_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_str1_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_str2_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_str3_v02.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_swi_v03.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/dsk/juice_sc_uvs_v01.bds',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_gala_v05.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_janus_v08.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_jmc_v02.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_jmag_v02.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_majis_v08.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_navcam_v01.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_pep_v14.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_radem_v03.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_rime_v04.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_rpwi_v03.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_str_v01.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_swi_v07.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_uvs_v06.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ik/juice_aux_v02.ti',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/lsk/naif0012.tls',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/pck/pck00011.tpc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/pck/de-403-masses.tpc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/pck/gm_de431.tpc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/pck/inpop19a_moon_pa_v01.bpc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/pck/earth_070425_370426_predict.bpc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/pck/juice_jup011.tpc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/pck/juice_roi_v01.tpc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/sclk/juice_fict_160326_v02.tsc',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_sci_v04.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_struct_v21.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_struct_internal_v01.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_cog_v00.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_cog_000060_230416_240516_v01.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_roi_v02.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/mar085_20200101_20400101.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/earthstns_fx_050714.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/estrack_v04.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_earthstns_v01.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/jup365_19900101_20500101.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/jup343_19900101_20500101.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/jup344-s2003_j24_19900101_20500101.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/jup346_19900101_20500101.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/de432s.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/inpop19a_19900101_20500101.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/noe-5-2017-gal-a-reduced_20200101_20380902.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_crema_5_1_150lb_23_1_plan_v01.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_orbc_000060_230414_310721_v01.bsp',
          'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_orbm_000059_240817_240827_v01.bsp']

kf = kernelFetch()
kf.ffList(urlKernelL=METAKR, forceDownload=False)


ROIs_filename = "../../data/roi_info/ganymede_roi_info.txt"  # Can be a list of strings or a single string

target_radii = spice.bodvrd(target_body, "RADII", 3)[1][1] # [km]

ifov = 15e-6
npix = 1735
imageRate = 10.
fs = 20.
instrument = Instrument(ifov, npix, imageRate, fs)

DB = ROIDataBase(ROIs_filename, target_body)
roi = DB.getROIs()

observer = 'JUICE'
flybys_ganymede =[["2033 NOV 26 18:22:11", "2033 NOV 27 18:22:11"],
                  ["2034 JAN 14 06:38:51", "2034 JAN 15 06:38:51"],
                  ["2034 JUN 05 18:53:51", "2034 JUN 06 18:53:51"],
                  ["2034 JUL 11 19:50:31", "2034 JUL 12 19:50:31"],
                  ["2034 SEP 07 06:03:51", "2034 SEP 08 06:03:51"],
                  ["2034 SEP 28 18:48:51", "2034 SEP 29 18:48:51"],
                  ["2034 NOV 18 09:58:51", "2034 NOV 19 09:58:51"]]

flybys_callisto = [["2033 JAN 13 13:32:12", "2033 JAN 14 13:32:12"],
                   ["2032 SEP 27 04:26:36", "2032 SEP 28 04:26:36"],
                   ["2032 NOV 16 05:53:16", "2032 NOV 17 05:53:16"],
                   ["2032 DEC 27 21:06:34", "2032 DEC 28 21:06:34"]]
                   #["2032 AUG 14 10:09:17", "2032 AUG 15 10:09:17"],
                   #["2032 SEP 10 11:52:29", "2032 SEP 11 11:52:29"],
                   #["2032 SEP 27 04:26:36", "2032 SEP 28 04:26:36"],
                   #["2032 OCT 13 20:53:36", "2032 OCT 14 20:53:36"],
                   #["2032 OCT 30 14:21:22", "2032 OCT 31 14:21:22"],
                   #["2032 NOV 16 05:53:16", "2032 NOV 17 05:53:16"],
                   #["2032 DEC 27 21:06:34", "2032 DEC 28 21:06:34"],
                   #["2033 JAN 13 13:32:12", "2033 JAN 14 13:32:12"],
                   #["2033 MAR 04 14:44:09", "2033 MAR 05 14:44:09"],
                   #["2033 APR 06 23:29:30", "2033 APR 07 23:29:30"],
                   #["2033 MAY 10 08:21:32", "2033 MAY 11 08:21:32"],
                   #["2033 JUL 16 02:11:44", "2033 JUL 17 02:11:44"]]


if target_body == "CALLISTO":
    flybys = flybys_callisto
else:
    flybys = flybys_ganymede

import cv2

filename = "../../data/mosaics/Ganymede_mosaic.jpg"
img = cv2.imread(filename, cv2.IMREAD_COLOR)
img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)  # Convert from BGR to RGB format

# Plot the image
fig, ax = plt.subplots(figsize=(10, 6))
ax.imshow(img_rgb, extent=[-180, 180, -90, 90])
scatter_lons = []
scatter_lats = []
scatter_alts = []
scatter_res = []

for i in range(len(roi)):
    vert_x = [0]*(len(roi[i].vertices) +1)
    vert_y = [0]*(len(roi[i].vertices) +1)
    for k in range(len(roi[i].vertices)):
        if roi[i].vertices[k][0] < 0:
            vert_x[k] = roi[i].vertices[k][0] + 180
        elif roi[i].vertices[k][0] > 0:
             vert_x[k] = roi[i].vertices[k][0] - 180
        vert_y[k] =  roi[i].vertices[k][1]
    vert_x[-1] = vert_x[0]
    vert_y[-1] = vert_y[0]
    if np.abs(vert_x[0] - vert_x[1]) > 100:
        continue
    else:
        ax.plot(np.array(vert_x), np.array(vert_y), c = 'white', linewidth = 1, linestyle = '--', alpha = 0.35)

for num, fb in enumerate(flybys):
    et0 = spice.utc2et(fb[0])
    et1 = spice.utc2et(fb[1])
    tstep = 3
    et = np.arange(et0, et1 + tstep, tstep)

    for t in et:

        alt = altitude(observer, target_body, t)
        gtlon, gtlat = PSOA.groundtrack(observer, t, target_body)
        resolution = pointres(ifov, [gtlon, gtlat], t, target_body, observer)
        #if alt < 30000:
        scatter_lats.append(gtlat)
        scatter_lons.append(gtlon)
        scatter_alts.append(alt)
        scatter_res.append(resolution)


# Create the scatter plot using a colormap
scatter = ax.scatter(scatter_lons, scatter_lats, marker= 's', c=scatter_alts, cmap='turbo_r', s=5)
scatter_2 = ax.scatter(scatter_lons, scatter_lats, marker= 's', c=scatter_res, cmap='turbo_r', s=5)

# Color bar for altitude
#cbar = plt.colorbar(scatter, ax=ax, orientation='horizontal', extend='max', pad=0.15, shrink=0.6)
#cbar.set_label('JUICE Altitude (km)', fontsize=15)
cbar_2 = plt.colorbar(scatter_2, ax=ax, orientation='horizontal', extend='max', pad=0.10, shrink=0.6)
cbar_2.set_label('JANUS Pixel Scale (km/px)', fontsize=15)
cbar_2.ax.xaxis.set_ticks_position('bottom')
cbar_2.ax.xaxis.set_label_position('bottom')
# Set axis titles and plot details
ax.set_title(f"JANUS Resolution across {target_body}", fontsize=20)
ax.set_xlim([-180, 180])
ax.set_ylim([-90, 90])
ax.set_xlabel('Longitude (°)', fontsize=18)
ax.set_ylabel('Latitude (°)', fontsize=18)
ax.set_aspect('equal', 'box')  # Use 'datalim' to avoid aspect ratio conflicts with colorbar

plt.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.07)

plt.show()


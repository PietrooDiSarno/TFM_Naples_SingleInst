import os
import random

import PSOA
import matplotlib
from pySPICElib import etToAxisStrings
from pySPICElib.kernelFetch import kernelFetch
import matplotlib.pyplot as plt
import numpy as np
import spiceypy as spice

from FuturePackage import ROIDataBase
from mosaic_algorithms.auxiliar_functions.plot.trgobsvec import trgobsvec

target_body = "GANYMEDE"  # Can be a list of strings or a single string

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

flyby = ["2034 JUN 05 18:53:51", "2034 JUN 06 18:53:51"]
start = spice.str2et(flyby[0])
end = spice.str2et(flyby[-1])

etvec = np.linspace(start, end, num = 1000)
distances = []

for t in etvec:
    pos, _ = spice.spkpos('GANYMEDE', t, 'ECLIPJ2000', 'NONE', 'JUICE')
    distances.append(np.linalg.norm(pos))

fig_, ax = plt.subplots()
plt.plot(etvec, distances, color = 'r')
plt.xlabel('Time instant')
plt.ylabel('Distance [km]')
etv, ets = etToAxisStrings(etvec, ndates=8, accurate = True)

plt.xticks(etv, ets)
plt.xticks(rotation=15)
plt.title('Distance of JUICE from Ganymede.')
plt.show()
plt.close(fig_)

ROIs_filename = "./data/roi_info/ganymede_roi_info.txt"  # Can be a list of strings or a single string
DB = ROIDataBase(ROIs_filename, target_body)
rois = DB.getROIs()
roinames = DB.getnames()


for roi in rois:
     angles = []
     distances = []
     for t in etvec:
        angles.append(PSOA.emissionang(roi.centroid, t, roi.body, 'JUICE'))
        distances.append(trgobsvec(roi.centroid, t, 'GANYMEDE', 'JUICE')[1])
     matplotlib.use('TkAgg')
     fig, ax1 = plt.subplots()
     plt.xlabel('Date')
     plt.title("Emission angle and distance from the centroid of " + roi.name)
     ax1.plot(etvec, angles, 'r', label='Emission angle [°]')
     ax2 = ax1.twinx()
     ax2.plot(etvec, distances, 'b', label='Distance [km]')

     etv, ets = etToAxisStrings(etvec, 8, accurate = True)
     ax1.set_xticks(etv)
     ax1.set_xticklabels(ets, rotation=15)

     lines1, labels1 = ax1.get_legend_handles_labels()
     lines2, labels2 = ax2.get_legend_handles_labels()
     lines = lines1 + lines2
     labels = labels1 + labels2
     ax1.legend(lines, labels)
     #plt.show()

     path = 'C:\\Users\\pietr\\ROI_makespans\\trial'
     fname = os.path.join(path, roi.name + ".tiff")
     fig.tight_layout()
     plt.subplots_adjust(left=0.1, right=0.9, bottom=0.20, top=0.92)
     plt.savefig(fname, format="tiff", bbox_inches="tight")


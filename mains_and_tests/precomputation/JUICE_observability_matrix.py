import seaborn as sns

import pandas as pd

import os
from pySPICElib.kernelFetch import kernelFetch

from pySPICElib.SPICEtools import *
from FuturePackage import ROIDataBase
from plotonMap import plotonMap

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

#################################################################################################################
ROIs_filename = "../../data/roi_info/ganymede_roi_info.txt"
target_body = "GANYMEDE"  # Can be a list of strings or a single string
DB = ROIDataBase(ROIs_filename, target_body)
rois = DB.getROIs()
roinames = DB.getnames()


flybys = [["2033 NOV 26 18:22:11", "2033 NOV 27 18:22:11"],
          ["2034 JAN 14 06:38:51", "2034 JAN 15 06:38:51"],
          ["2034 JUN 05 18:53:51", "2034 JUN 06 18:53:51"],
          ["2034 JUL 11 19:50:31", "2034 JUL 12 19:50:31"],
          ["2034 SEP 07 06:03:51", "2034 SEP 08 06:03:51"],
          ["2034 SEP 28 18:48:51", "2034 SEP 29 18:48:51"],
          ["2034 NOV 18 09:58:51", "2034 NOV 19 09:58:51"]]
matrix1 = []
a = []
names = []
anot = []
rfb = []
ismin = []
globalMat = []
previous = 'None'
for i, name in enumerate(roinames):
    names.append(name)
    roinum = roinames[i][14]
    a = [np.nan]*7
    for j, fb in enumerate(flybys):
        patron = f"output_{roinames[i]}_fb{str(j)}.txt"
        for file in os.listdir("../../data/roi_files"):
            if file == patron:
                rfb.append(name)
                with open("../../data/roi_files/"+file, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if 'minRes = ' in line:
                            minres = line.split('=')
                            #print(float(minres[1].strip()))
                            a[j] = float(minres[1].strip())
                            #print(matrix[i][j])
    matrix1.append(a)
    globalMat.append(a)


df = pd.DataFrame(globalMat)
df = df.dropna(how='all')

valid_rows = df.index
filtered_names = [roinames[i] for i in valid_rows]

rows = len(df)
cols = len(df.columns)

mid_index = rows // 2
df1 = df.iloc[:mid_index]
df2 = df.iloc[mid_index:]

names1 = filtered_names[:mid_index]
names1 = [n[-6:] for n in names1]
names2 = filtered_names[mid_index:]
names2 = [n[-6:] for n in names2]

cmap = sns.color_palette("coolwarm_r", as_cmap=True)

plt.figure(figsize=(cols, mid_index))
sns.heatmap(df1, cmap=cmap, annot=True, fmt=".2f", linewidth=.5, vmin=0)
plt.gca().set_facecolor('grey')
plt.grid(False)
plt.gca().set_yticklabels(names1, rotation=0)
plt.xlabel('Flyby number')
plt.show()

plt.figure(figsize=(cols, rows - mid_index))
sns.heatmap(df2, cmap=cmap, annot=True, fmt=".2f", linewidth=.5, vmin=0)
plt.gca().set_facecolor('grey')
plt.grid(False)
plt.gca().set_yticklabels(names2, rotation=0)
plt.xlabel('Flyby number')
plt.ylabel('ROI')
plt.show()


plotonMap(globalMat, rois)

import copy
import random
import sys
from pySPICElib.kernelFetch import kernelFetch
from FuturePackage.roiDataBase import ROIDataBase
import subprocess
from pySPICElib.SPICEtools import *

if len(sys.argv) < 2:
    nAgents = input('Number of agents: ')

#################################################################################################################
np.random.seed(123)
random.seed(123)

target_body = ["GANYMEDE"] # can be a list of strings or a single list

if target_body == ["GANYMEDE"]:
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

else:
      if target_body == ["CALLISTO"]:
          METAKR = ['https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_crema_5_1_150lb_23_1_b2_default_v01.bc',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_crema_5_1_150lb_23_1_b2_comms_v01.bc',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_crema_5_1_150lb_23_1_b2_conjctn_v01.bc',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_crema_5_1_150lb_23_1_b2_flybys_v01.bc',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sc_crema_5_1_150lb_23_1_b2_baseline_v01.bc',
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
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sa_crema_5_1_150lb_23_1_b2_default_v01.bc',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sa_crema_5_1_150lb_23_1_b2_baseline_v01.bc',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_sa_ptr_soc_s007_01_s02p00_v01.bc',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_mga_crema_5_1_150lb_23_1_b2_default_v01.bc',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_mga_crema_5_1_150lb_23_1_b2_baseline_v01.bc',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/juice_mga_ptr_soc_s007_01_s02p00_v01.bc',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_v40.tf',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_sci_v17.tf',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_ops_v11.tf',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_dsk_surfaces_v11.tf',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_roi_v02.tf',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/juice_events_crema_5_1_150lb_23_1_b2_v02.tf',
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
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_crema_5_1_150lb_23_1_b2_v01.bsp',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_orbc_000060_230414_310721_v01.bsp',
                'https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_orbm_000059_240817_240827_v01.bsp']

kf = kernelFetch()
kf.ffList(urlKernelL=METAKR, forceDownload=False)

# INPUTS

#   a) ROI INFO
#   a.1) Raw data info
if target_body == ['GANYMEDE']:
    ROIs_filename = "../../../data/roi_info/ganymede_roi_info.txt"  # Can be a list of strings or a single string
    ROIs_antijovian = "../../../data/roi_info/ganymede_anti_jovian.txt"
else:
      if target_body ==["CALLISTO"]:
          ROIs_filename = "../../../data/roi_info/callisto_roi_info.txt"  # Can be a list of strings or a single string
          ROIs_antijovian = "../../../data/roi_info/callisto_anti_jovian.txt"

DB = ROIDataBase(ROIs_filename, target_body)
rois = DB.getROIs()
roinames = DB.getnames()


n_ROIs = int(nAgents)
k = len(roinames) // n_ROIs
i_start = np.array(range(k)) * n_ROIs
i_end = np.array(range(1, k + 1)) * n_ROIs
if i_end[-1] < len(roinames):
      i_start = np.append(i_start, i_start[-1] + n_ROIs)
      i_end = np.append(i_end, len(roinames))

### Start of subprocesses
for index in range(len(i_start)):
      c = []
      for i, roiname in enumerate(roinames[i_start[index]:i_end[index]]):
          #if i == 104 or i == 12: continue #Smth wrong with these ones (maybe end-of-the-world Ganymede rois)
          c.append([sys.executable, 'checkOneROI.py'] + [roiname])
      #roiname = roinames[0]
      #c.append([sys.executable, 'checkOneROI.py'] + [roiname])

      proc = []  # list of p
      finished = []  # finished procs

      # spawn all procs
      for cmd in c:
          proc.append(subprocess.Popen(cmd))
          finished.append(0)
      # wait for completion
      nf = 0  # number of finished processes
      while True:
          for a in range(0, len(proc)):
              ret = proc[a].poll()
              if ret is not None and finished[a] == 0:
                  finished[a] = 1
                  nf = nf + 1
                  print('** command', c[a], 'finished ', ret, 'nf=', nf, 'of ', len(proc))
          if nf == len(proc):
              print('** all finished')
              break



import PMOT as pm
import sys
import spiceypy as spice
from pySPICElib.kernelFetch import kernelFetch

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

if len(sys.argv) > 1:
    print(sys.argv[1])
    pathToData = sys.argv[1]
    caseName = sys.argv[2]

    GApop = input('Population: ')
    GAgen = input('Generations: ')
    nm = input('Mutants: ')
    nd = input('Descendants: ')
    ne = input('Elites: ')
    ncm = input('Can mutate: ')
    savingIte = input('Start saving at generation: ')
    loadingIte = input('Start loading at generation: ')
    info = input('Info: ')
    nAgents = input('Number of agents: ')
else:
    GApop = input('Population: ')
    GAgen = input('Generations: ')
    nm = input('Mutants: ')
    nd = input('Descendants: ')
    ne = input('Elites: ')
    ncm = input('Can mutate: ')
    savingIte = input('Start saving at generation: ')
    loadingIte = input('Start loading at generation: ')
    info = input('Info: ')
    nAgents = input('Number of agents: ')
    pathToData = 'JUICE_multiAgent'
    caseName = 'JUICE'

ads = pm.agentDataSharing(True, pathToData, caseName)
lastAgentFound = ads.getLastKnownAgent()
#print('lastAgentFound=', lastAgentFound)
#arguments = [GApop, GAgen, nm, nd, ne, ncm, savingIte, loadingIte, info]
#firstAgent = lastAgentFound + 1
#if True:
#    sp = pm.spawnAgents('JUICE_schedule_spawn', range(firstAgent, firstAgent + int(nAgents)), None, pathToData, caseName,
#                        arguments=arguments)
#    sp.spawn()

agent, generation, fit, ind = ads.findBest()
print('Best result obtained: agent =', agent, 'generation =', generation, 'fit =', fit, 'ind =', ind)
ads.plotAllHistory()

for i in range(len(ind.stol)):
    start = spice.et2utc(ind.stol[i], 'C', 0)
    end = spice.et2utc(ind.stol[i] + ind.obsLength[i], 'C', 0)
    print('ROI number: ', i,'from: ', start, 'to', end, 'res:', ind.qroi[i],'cov:', ind.croi[i])


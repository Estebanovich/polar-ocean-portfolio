# Verbatim copies of the MSc thesis files (Exp. 4)

Copied unmodified (`cp -p`) from the read-only thesis tree. Original absolute path and SHA-256 of each file; the copy has the same hash. Thesis binary inputs (> 1 MB) are not copied: Task 3 regenerates them or they are read from the original path. Note: the thesis runs used MNC output; Exp. 4 changes only the I/O to MDS.

Helper modules of the Output_IWBay scripts: `visualizar_espectros_metodologia.py` (imported by `visualizar_espectros_bay_vs_nobay.py`), `mapas_psd_dominio.py` (imported by `mapas_psd_dominio_bay_vs_nobay.py`), `coords_iwbay.py` and `MIT_DatProcesW.py` (imported by `mapas_psd_dominio.py`). `MIT_DatProcesW.py` also needs `MITgcmutils`, `cmocean`, `pandas`, `netCDF4` at import time.

| Original path | SHA-256 |
|---|---|
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/code/DIAGNOSTICS_SIZE.h` | `48a35a2137a47ed37828e5dba9c18e0115c83dba4e7eec7c1b8c60f30828d385` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/code/MNC_SIZE.h` | `6b3aada224ad55da2402e87ecdbeb287313d60687e3b208d7b3ee67e61f86580` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/code/OBCS_OPTIONS.h` | `e8ef561a1975d34e1d6c789606c8728170a9008d6e361cfbe50c4b9c59952641` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/code/packages.conf` | `fa04854368e5f62a245f603f5eb00a93465b3cdd3207652f994e568499833d70` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/code/SIZE_ORIGINAL.h` | `4bd99b4c58ab490c243c510a5e054afe8e2146f78c3ef6a98507ceb0658c5ede` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/code/SIZE.h` | `47e3cccf181d0d4ad2d80d4259919b3037e9d9cefb6e850eaa4d874de9c07ed0` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/code/SIZE.h_mpi` | `dd11947e9de3a4c13e3fbcd02fb1572dd1e89cd6ee1020b0de74a41085a9a4c4` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/input/bay_domain.py` | `7e6e9d7a679a9c08234cfaa457e5fefee7c5a2cbe6a2e8cb91baa12d6f39c8b2` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/input/make_stretched_grid.py` | `2dd8efbf0ebb5fd4e8bfcaca1f052c4b14d4d458cbd1911119b4b74bdffeaf5b` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/input/make_T_S_bin_exp.py` | `984f76663e88ac466ae7037947259204b43224cd65e9575f1ab364afa9ebcccc` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/input/wind_patch.py` | `3241d5cbc17df32ad83b1bd9cb0b82bdbfe120c875f5170095a10eb62adc41b5` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand_nobay/data.diagnostics` | `c615b7656452c73e2528316ed56716ce838f8a47689ef516adee75a1c9b360d3` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand_nobay/data.obcs` | `f93e9ae5e45590b703c83d6da01626393030780cde3a7010b4833fc976883670` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand_nobay/data.pkg` | `1c86295be40a64cc957e5fe1f1dd93006ab339fa192c1650fed2898205f7eaf4` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand_nobay/data.stage1` | `238597f3cfd12ca66d0ca4093f1eee66302a396cc1e2813ec51e621b40574e6b` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand_nobay/data.stage2` | `205ac82e05d49a5e92dadae5c79d0c9f10c1564c95b35732daa5a04ca5801439` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand_nobay/eedata` | `c37f927b330c60c784f7fea12273be3b0df665497b0f56f673f6b0012d55b66d` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand/data.diagnostics` | `c615b7656452c73e2528316ed56716ce838f8a47689ef516adee75a1c9b360d3` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand/data.obcs` | `f93e9ae5e45590b703c83d6da01626393030780cde3a7010b4833fc976883670` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand/data.pkg` | `1c86295be40a64cc957e5fe1f1dd93006ab339fa192c1650fed2898205f7eaf4` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand/data.stage1` | `847ab7ad02069994bb0a735ead335360d96000e6221968b41cfaa3b0fa35c6c5` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand/data.stage2` | `5b881aa71ca2b1ed754b2c83a7564a1f8d65a64adba6e6535605b31b4e79ea0d` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand/eedata` | `c37f927b330c60c784f7fea12273be3b0df665497b0f56f673f6b0012d55b66d` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay/coords_iwbay.py` | `557dfd11fc01f6e85caaea26b384bcb42f96e253d10025efcb6d3aa21025c2d5` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay/energy_flux_IWBay.py` | `ed39dc901781a6bc84688156bb2cb2eae94675bb532c4512a6b5c3ef32ec5068` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay/mapas_psd_dominio_bay_vs_nobay.py` | `91910cdad918d9e6b1e1d68d1b6a18e29a2154fc06f424b49586f15736381039` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay/mapas_psd_dominio.py` | `0616e8dbd1b78357a16d22b8dbb5bf5147be3bc0a4c2cfb4e1bb22bda71f884d` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay/MIT_DatProcesW.py` | `99d1747184ae67213779f13c968fecc1799a0a0f84443d1953b38c7fcffbc4a8` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay/plot_flujo_boca_bt_bc.py` | `83f15198cd61dd4b9e5386bb07a32b584e732535bd657b196e3dc5c18a00fba3` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay/visualizar_espectros_bay_vs_nobay.py` | `90799e98591294042cb337ef92913c2ff454793e174ad66884a9f0755ecee078` |
| `/Volumes/ESTEBAN_MIT/MITgcm/verification/Output_IWBay/visualizar_espectros_metodologia.py` | `07418338dfca4d7712da9515ea8aeec0ced69d8ac288a4bf41caa0bc1602a46e` |

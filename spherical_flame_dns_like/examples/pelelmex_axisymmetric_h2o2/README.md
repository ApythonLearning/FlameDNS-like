# PeleLMeX Axisymmetric H2/O2 Case Template

This directory contains a PeleLMeX `Exec` case template for a 2D RZ premixed
H2/O2 spherical hot-kernel simulation. It is a template because PeleLMeX is not
installed in the current machine environment detected by
`scripts/check_dependencies.py`.

Files:

- `GNUmakefile`: 2D GNUmake build configuration.
- `inputs.axisym_h2o2`: run-time inputs with all requested parameters.
- `pelelmex_prob.H`: problem parameters and spherical high-temperature kernel initialization.
- `pelelmex_prob.cpp`: `ParmParse` reader for `prob.*` inputs.

Generate the 1D Cantera profile first:

```bash
cd spherical_flame_dns_like
python scripts/compute_1d_h2o2_flame.py --h2-volume-fraction 0.10 --output-dir profiles/h2o2_phi_auto
```

Then copy this directory into `PeleLMeX/Exec/RegTests/AxisymmetricH2O2`, update
`PELE_HOME` in `GNUmakefile` if needed, and build:

```bash
make -j4 TPL
make -j4
./PeleLMeX2d.gnu.MPI.ex inputs.axisym_h2o2
```

If your PeleLMeX checkout uses different state-index names or a different H2/O2
mechanism, adjust `URHO`, `UMX`, `UMY`, `UEDEN`, `UTEMP`, `UFS`, `H2_ID`, and
`O2_ID` in `pelelmex_prob.H` to match the generated headers in that checkout.
The template intentionally does not generate or postprocess synthetic DNS data.

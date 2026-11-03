"""Print the EFFECTIVE numba codegen and environment as JSON (no JIT compilation happens)."""
import hashlib
import json
import os
import platform
import sys

import llvmlite
import llvmlite.binding as ll
import numba
import numpy
import scipy
from numba.core import config
from numba.core.codegen import get_host_cpu_features

cpu = config.CPU_NAME if config.CPU_NAME is not None else ll.get_host_cpu_name()
feat = config.CPU_FEATURES if config.CPU_FEATURES is not None else get_host_cpu_features()
try:
    import pyamg
    pyamg_v = pyamg.__version__
except Exception as e:  # noqa: BLE001
    pyamg_v = f"unavailable: {e}"
try:
    from threadpoolctl import threadpool_info
    blas = threadpool_info()
except Exception as e:  # noqa: BLE001
    blas = f"unavailable: {e}"
try:
    svml = bool(config.USING_SVML)
except Exception:
    svml = None
out = {
    "python": sys.version, "platform": platform.platform(), "machine": platform.machine(),
    "numba": numba.__version__, "llvmlite": llvmlite.__version__,
    "llvm": ".".join(map(str, ll.llvm_version_info)),
    "triple": ll.get_process_triple(),
    "config_CPU_NAME": config.CPU_NAME, "config_CPU_FEATURES": config.CPU_FEATURES,
    "effective_cpu_name": cpu,
    "effective_features_sha256": hashlib.sha256(feat.encode()).hexdigest(),
    "effective_features_n": len([f for f in feat.split(",") if f]),
    "ENABLE_AVX": config.ENABLE_AVX, "OPT": str(config.OPT), "USING_SVML": svml,
    "numpy": numpy.__version__, "scipy": scipy.__version__, "pyamg": pyamg_v,
    "env": {k: v for k, v in os.environ.items()
            if k.startswith(("NUMBA_", "OMP_", "OPENBLAS_", "MKL_", "PYTHON", "PROV_"))},
    "blas": blas,
}
print(json.dumps(out, indent=1, sort_keys=True, default=str))

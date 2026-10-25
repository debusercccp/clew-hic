"""Fixture con il .mcool sintetico di ``clew.synthetic`` (TAD e compartimenti noti)."""
from __future__ import annotations

import pytest

from clew.synthetic import (  # noqa: F401  (riesportati per i test)
    CHROMSIZES, COMP_BP, RESOLUTIONS, TAD_BP, make_synthetic_mcool, true_compartment,
)


@pytest.fixture(scope="session")
def synthetic_mcool(tmp_path_factory) -> str:
    """Percorso di un .mcool bilanciato a 10, 20 e 40 kb."""
    return make_synthetic_mcool(tmp_path_factory.mktemp("hic"))


@pytest.fixture(scope="session")
def synthetic_cool(synthetic_mcool) -> str:
    """Percorso del .cool a risoluzione singola (non bilanciato)."""
    return synthetic_mcool.replace(".mcool", ".cool")

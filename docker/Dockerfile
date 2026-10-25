# Immagine unica per JupyterLab, job batch e test. Le dipendenze vengono da
# pyproject.toml + uv.lock, le stesse di `uv sync` in locale.

# --- build: uv risolve dal lockfile e compila cooltools (solo sorgente su PyPI) ---
FROM python:3.11.16-slim-bookworm AS build
COPY --from=ghcr.io/astral-sh/uv:0.12.11 /uv /bin/uv
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential \
 && rm -rf /var/lib/apt/lists/*
ENV UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_PYTHON_DOWNLOADS=never \
    UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy
WORKDIR /build
COPY pyproject.toml uv.lock ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --all-groups

# --- runtime: solo il venv, niente compilatori ---
FROM python:3.11.16-slim-bookworm

# Utente con lo stesso UID/GID dell'host: i file scritti in notebooks/ e
# results/ restano di proprietà tua.
ARG HOST_UID=1000
ARG HOST_GID=1000
RUN groupadd -o -g "${HOST_GID}" clew \
 && useradd -o -m -u "${HOST_UID}" -g "${HOST_GID}" -s /bin/bash clew

COPY --from=build /opt/venv /opt/venv
ENV PATH=/opt/venv/bin:$PATH \
    PYTHONPATH=/work \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib
USER clew
WORKDIR /work
EXPOSE 8888

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--ServerApp.root_dir=/work"]
